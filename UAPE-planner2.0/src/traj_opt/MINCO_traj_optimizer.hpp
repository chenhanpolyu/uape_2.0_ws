/*
    MIT License

    Copyright (c) 2021 Zhepei Wang (wangzhepei@live.com)

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

#ifndef MINCOOPT_HPP
#define MINCOOPT_HPP

#include "trajectory.hpp"
#include <call_states/ros_communicate.h>
#include "mlmap.h"
// @article{WANG2021GCOPTER,
//     title={Geometrically Constrained Trajectory Optimization for Multicopters},
//     author={Wang, Zhepei and Zhou, Xin and Xu, Chao and Gao, Fei},
//     journal={arXiv preprint arXiv:2103.00190},
//     year={2021}
// }

// class ConstraintPoints
// {
// public:
//     int cp_size; // deformation points
//     Eigen::MatrixXd points;
//     std::vector<std::vector<Eigen::Vector3d>> base_point; // The point at the statrt of the direction vector (collision point)
//     std::vector<std::vector<Eigen::Vector3d>> direction;  // Direction vector, must be normalized.
//     std::vector<bool> flag_temp;                          // A flag that used in many places. Initialize it everytime before using it.

//     void resize_cp(const int size_set)
//     {
//         cp_size = size_set;

//         base_point.clear();
//         direction.clear();
//         flag_temp.clear();

//         points.resize(3, size_set);
//         base_point.resize(cp_size);
//         direction.resize(cp_size);
//         flag_temp.resize(cp_size);
//     }

//     void segment(ConstraintPoints &buf, const int start, const int end)
//     {
//         if (start < 0 || end >= cp_size || points.rows() != 3)
//         {
//             ROS_ERROR("Wrong segment index! start=%d, end=%d", start, end);
//             return;
//         }

//         buf.resize_cp(end - start + 1);
//         buf.points = points.block(0, start, 3, end - start + 1);
//         buf.cp_size = end - start + 1;
//         for (int i = start; i <= end; i++)
//         {
//             buf.base_point[i - start] = base_point[i];
//             buf.direction[i - start] = direction[i];
//         }
//     }

//     static inline int two_thirds_id(Eigen::MatrixXd &points, const bool touch_goal)
//     {
//         return touch_goal ? points.cols() - 1 : points.cols() - 1 - (points.cols() - 2) / 3;
//     }

//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
// };

class MINCO_S3_G // G FOR General
{
public:
    MINCO_S3_G() = default;
    ~MINCO_S3_G() { A.destroy(); }
    double pena = 0.0;
    double pena_ball = 0.0;
    double pena_sta = 0.0;
    double pena_acc = 0.0;
    double cost_total;
    double dyn_coe = 15.0;
    size_t opt_iter_num = 0;
    mlmap::Ptr MLmap;
    Eigen::MatrixXd global_bound_;
    // my
    // double compute_time = 0;
private:
    int N;

    Eigen::Matrix3d headPVA;
    Eigen::Matrix3d tailPVA;
    Eigen::VectorXd T1;
    BandedSystem A;
    Eigen::MatrixXd b;

    // Temp variables
    Eigen::VectorXd T2;
    Eigen::VectorXd T3;
    Eigen::VectorXd T4;
    Eigen::VectorXd T5;
    Eigen::MatrixXd gdC;
    double *cost_block;

    std::vector<std::pair<size_t, Eigen::Vector3d>> contorl_pairs;

private:
    template <typename EIGENVEC>
    inline void addGradJbyT(EIGENVEC &gdT) const
    {
        for (int i = 0; i < N; i++)
        {
            gdT(i) += 36.0 * b.row(6 * i + 3).squaredNorm() +
                      288.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T1(i) +
                      576.0 * b.row(6 * i + 4).squaredNorm() * T2(i) +
                      720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T2(i) +
                      2880.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T3(i) +
                      3600.0 * b.row(6 * i + 5).squaredNorm() * T4(i);
        }
        return;
    }

    template <typename EIGENMAT>
    inline void addGradJbyC(EIGENMAT &gdC) const
    {
        for (int i = 0; i < N; i++)
        {
            gdC.row(6 * i + 5) += 240.0 * b.row(6 * i + 3) * T3(i) +
                                  720.0 * b.row(6 * i + 4) * T4(i) +
                                  1440.0 * b.row(6 * i + 5) * T5(i);
            gdC.row(6 * i + 4) += 144.0 * b.row(6 * i + 3) * T2(i) +
                                  384.0 * b.row(6 * i + 4) * T3(i) +
                                  720.0 * b.row(6 * i + 5) * T4(i);
            gdC.row(6 * i + 3) += 72.0 * b.row(6 * i + 3) * T1(i) +
                                  144.0 * b.row(6 * i + 4) * T2(i) +
                                  240.0 * b.row(6 * i + 5) * T3(i);
        }
        return;
    }

    inline void solveAdjGradC(Eigen::MatrixXd &gdC) const
    {
        A.solveAdj(gdC);
        return;
    }

    template <typename EIGENVEC>
    inline void addPropCtoT(const Eigen::MatrixXd &adjGdC, EIGENVEC &gdT) const
    {
        Eigen::MatrixXd B1(6, 3), B2(3, 3);

        Eigen::RowVector3d negVel, negAcc, negJer, negSnp, negCrk;

        for (int i = 0; i < N - 1; i++)
        {
            negVel = -(b.row(i * 6 + 1) +
                       2.0 * T1(i) * b.row(i * 6 + 2) +
                       3.0 * T2(i) * b.row(i * 6 + 3) +
                       4.0 * T3(i) * b.row(i * 6 + 4) +
                       5.0 * T4(i) * b.row(i * 6 + 5));
            negAcc = -(2.0 * b.row(i * 6 + 2) +
                       6.0 * T1(i) * b.row(i * 6 + 3) +
                       12.0 * T2(i) * b.row(i * 6 + 4) +
                       20.0 * T3(i) * b.row(i * 6 + 5));
            negJer = -(6.0 * b.row(i * 6 + 3) +
                       24.0 * T1(i) * b.row(i * 6 + 4) +
                       60.0 * T2(i) * b.row(i * 6 + 5));
            negSnp = -(24.0 * b.row(i * 6 + 4) +
                       120.0 * T1(i) * b.row(i * 6 + 5));
            negCrk = -120.0 * b.row(i * 6 + 5);

            B1 << negSnp, negCrk, negVel, negVel, negAcc, negJer;

            gdT(i) += B1.cwiseProduct(adjGdC.block<6, 3>(6 * i + 3, 0)).sum();
        }

        negVel = -(b.row(6 * N - 5) +
                   2.0 * T1(N - 1) * b.row(6 * N - 4) +
                   3.0 * T2(N - 1) * b.row(6 * N - 3) +
                   4.0 * T3(N - 1) * b.row(6 * N - 2) +
                   5.0 * T4(N - 1) * b.row(6 * N - 1));
        negAcc = -(2.0 * b.row(6 * N - 4) +
                   6.0 * T1(N - 1) * b.row(6 * N - 3) +
                   12.0 * T2(N - 1) * b.row(6 * N - 2) +
                   20.0 * T3(N - 1) * b.row(6 * N - 1));
        negJer = -(6.0 * b.row(6 * N - 3) +
                   24.0 * T1(N - 1) * b.row(6 * N - 2) +
                   60.0 * T2(N - 1) * b.row(6 * N - 1));

        B2 << negVel, negAcc, negJer;

        gdT(N - 1) += B2.cwiseProduct(adjGdC.block<3, 3>(6 * N - 3, 0)).sum();

        return;
    }

    template <typename EIGENMAT>
    inline void addPropCtoP(const Eigen::MatrixXd &adjGdC, EIGENMAT &gdInP) const
    {
        for (int i = 0; i < N - 1; i++)
        {
            gdInP.col(i) += adjGdC.row(6 * i + 5).transpose(); // the last row for this segment
        }
        return;
    }

    inline void normalizeFDF(const Eigen::Vector3d &x,
                             Eigen::Vector3d &xNor,
                             Eigen::Matrix3d &G) const
    {
        const double a = x(0), b = x(1), c = x(2);
        const double aSqr = a * a, bSqr = b * b, cSqr = c * c;
        const double ab = a * b, bc = b * c, ca = c * a;
        const double xSqrNorm = aSqr + bSqr + cSqr;
        const double xNorm = sqrt(xSqrNorm);
        const double den = xSqrNorm * xNorm;
        xNor = x / xNorm;
        G(0, 0) = bSqr + cSqr;
        G(0, 1) = -ab;
        G(0, 2) = -ca;
        G(1, 0) = -ab;
        G(1, 1) = aSqr + cSqr;
        G(1, 2) = -bc;
        G(2, 0) = -ca;
        G(2, 1) = -bc;
        G(2, 2) = aSqr + bSqr;
        G /= den;
        return;
    }

    inline void positiveSmoothedL1(const double &x, double &f, double &df) const
    {
        const double pe = 1.0e-4;
        const double half = 0.5 * pe;
        const double f3c = 1.0 / (pe * pe);
        const double f4c = -0.5 * f3c / pe;
        const double d2c = 3.0 * f3c;
        const double d3c = 4.0 * f4c;

        if (x < pe)
        {
            f = (f4c * x + f3c) * x * x * x;
            df = (d3c * x + d2c) * x * x;
        }
        else
        {
            f = x - half;
            df = 1.0;
        }

        return;
    }
    inline void getFoot(int seg, const Eigen::Vector3d &pos, Eigen::Vector3d &foot, const std::vector<Eigen::Vector3d> &path)
    {
        int choose_id;
        double min_dist = INFINITY, len, proj_on_line, dot_res, choose_len, dist2line;
        Eigen::Vector3d base;
        for (size_t i = 0; i < path.size() - 1; i++)
        {
            base = path[i + 1] - path[i];
            len = (base).norm();
            proj_on_line = (pos - path[i]).dot(base);         //|a| * |b| * cos(theta)
            if (proj_on_line > 0 && proj_on_line < len * len) // foot on this line
            {
                dist2line = ((pos - path[i]).cross(base)).norm() / len;
                if (dist2line < min_dist)
                {
                    min_dist = dist2line;
                    choose_id = i;
                    dot_res = proj_on_line;
                    choose_len = len;
                }
            }
        }
        if (!std::isinf(min_dist))
        {
            foot = path[choose_id] + dot_res / choose_len * (path[choose_id + 1] - path[choose_id]) / choose_len;
        }
        else
        {
            base = path[seg + 1] - path[seg];
            foot = path[seg] + ((pos - path[seg]).dot(base)) / (base).squaredNorm() * base;
        }
    }

    inline bool in_global_bound(const Eigen::Vector3d &pos, float safe_margin)
    {
        return (pos(0) > this->global_bound_(0, 0) &&
                pos(0) < this->global_bound_(0, 1) &&
                pos(1) > this->global_bound_(1, 0) &&
                pos(1) < this->global_bound_(1, 1) &&
                pos(2) > this->global_bound_(2, 0) &&
                pos(2) < this->global_bound_(2, 1));
    }
    template <typename EIGENVEC>
    inline void addTimeIntPenalty(const Eigen::VectorXi cons,
                                  const Eigen::VectorXi &idxHs,
                                  const std::vector<Eigen::MatrixXd> &cfgHs,
                                  const Eigen::Vector3d &ellipsoid,
                                  const double safeMargin,
                                  const double vMax,
                                  const double thrAccMin,
                                  const double thrAccMax,
                                  const double bdrMax,
                                  const double gAcc,
                                  const Eigen::Vector4d ci,
                                  double &cost,
                                  EIGENVEC &gdT,
                                  Eigen::MatrixXd &gdC,
                                  const double plan_t,
                                  const dynobs_tmp *dynobs_pointer,
                                  const std::vector<Eigen::Vector3d> &path)
    {
        pena = 0.0;
        pena_ball = 0.0;
        pena_sta = 0.0;
        pena_acc = 0.0;
        double vMaxSqr;
        if (opt_iter_num % 10 == 0)
        {
            contorl_pairs.clear();
        }

        if (dynobs_pointer->ball_number > 0)
        {
            if (vMax <= 2)
            {
                vMaxSqr = 32.0;
            }
            else
            {
                vMaxSqr = vMax * vMax * 8;
            }
        }
        // else if (dynobs_pointer->dyn_number >0)
        //  {vMaxSqr = vMax * vMax*1.9;}
        else
        {
            vMaxSqr = vMax * vMax;
        }

        // const double thrAccMinSqr = thrAccMin * thrAccMin;
        const double thrAccMaxSqr = thrAccMax * thrAccMax;
        // const double bdrMaxSqr = bdrMax * bdrMax;
        std::vector<bool> mask_dyn_idx(dynobs_pointer->dyn_number, false);
        Eigen::Vector3d pos, vel, acc, jer, sna;
        double step, alpha;
        double s1, s2, s3, s4, s5;
        Eigen::Matrix<double, 6, 1> beta0, beta1, beta2, beta3, beta4;
        int K;
        // Eigen::Matrix3d rotM;
        // double signedDist, signedDistSqr, signedDistCub;
        // double gradSignedDt;
        // Eigen::Vector3d h, zB, czB, xC, yB, xB;
        // Eigen::Matrix3d dnczB, dzB, cdzB, dyB, dxB;
        Eigen::Vector3d outerNormal, point;
        // double eNorm;
        // Eigen::Vector3d eNormGd;
        // Eigen::Matrix3d gradSdTxyz;
        // Eigen::Vector3d gradSdT;
        // double outerNormaldVel;
        // Eigen::Matrix<double, 6, 3> gradSdCx, gradSdCy, gradSdCz, gradSdC;
        // Eigen::Matrix<double, 6, 3> beta2dOuterNormalTp, beta0dOuterNormalTp;

        double violaVel, violaThrl, violaThrh;
        double violaVelPenaD, violaThrlPenaD, violaThrhPenaD;
        double violaVelPena, violaThrlPena, violaThrhPena;
        Eigen::Matrix<double, 6, 3> gradViolaVc, gradViolaThrlc, gradViolaThrhc;
        double gradViolaVt, gradViolaThrlt, gradViolaThrht;
        double fThr, sqrMagThr;
        Eigen::Vector3d dfThr, dSqrMagThr;
        Eigen::Vector3d h;
        Eigen::Vector3d ct_center, foot;
        double omg, violaPos, violaPosPenaD, violaPosPena;
        double t_gap, t_now;
        double wei_dyn = ci(0) * dyn_coe;
        double wei_dyn_acc = ci(0) * 0.0025; // 0.0005; //
        double wei_ball = ci(0) * 24;
        double current_pena;
        int innerLoop, idx;
        constexpr double inv_a2 = 1 / 2.0 / 2.0, inv_b2 = 1.0;
        double inv_x, inv_y, inv_z;
        double dist2line, gradViolaPt;
        size_t cur_id_check = 0;
        Eigen::Vector3d gradP, direction;
        size_t total_samp_45 = static_cast<size_t>((cons.sum() + N) * 4 / 5);
        size_t current_samp = 0;
        for (int i = 0; i < N; i++)
        {
            const auto &c = b.block<6, 3>(i * 6, 0);
            s1 = 0.0;
            step = T1(i) / cons(i);
            innerLoop = cons(i) + 1;

            for (int j = 0; j < innerLoop; j++)
            {
                current_samp = innerLoop * i + j;
                if (i > 0)
                {
                    t_now = plan_t + T1.head(i).sum() + s1;
                }
                else
                {
                    t_now = plan_t + s1;
                }
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0 << 1.0, s1, s2, s3, s4, s5;
                // cout << "mk3200" <<endl;
                beta1 << 0.0, 1.0, 2.0 * s1, 3.0 * s2, 4.0 * s3, 5.0 * s4;
                beta2 << 0.0, 0.0, 2.0, 6.0 * s1, 12.0 * s2, 20.0 * s3;
                beta3 << 0.0, 0.0, 0.0, 6.0, 24.0 * s1, 60.0 * s2;
                beta4 << 0.0, 0.0, 0.0, 0.0, 24.0, 120.0 * s1;
                // cout << "mk3201" <<endl;
                alpha = 1.0 / cons(i) * j;
                pos = c.transpose() * beta0; // c_0->c_5
                vel = c.transpose() * beta1;
                acc = c.transpose() * beta2;
                jer = c.transpose() * beta3;
                sna = c.transpose() * beta4;

                h = acc;
                h(2) += gAcc;
                fThr = h.norm();
                dfThr = h / fThr;
                sqrMagThr = fThr * fThr;
                dSqrMagThr = 2 * h;

                violaThrh = sqrMagThr - thrAccMaxSqr;

                omg = (j == 0 || j == innerLoop - 1) ? 0.5 : 1.0;

                // idx = idxHs(i);
                // K = cfgHs[idx].cols();
                // for (int k = 0; k < K; k++)
                // {
                //     outerNormal = cfgHs[idx].col(k).head<3>();
                //     violaPos = outerNormal.dot(pos - cfgHs[idx].col(k).tail<3>()) + 1.5 * safeMargin; // projection on the direction of normal vector
                //     if (violaPos > 0.0 && T1(i) - step > s1 && step < s1)
                //     {
                //         positiveSmoothedL1(violaPos, violaPosPena, violaPosPenaD);
                //         gdC.block<6, 3>(i * 6, 0) += omg * step * ci(0) * violaPosPenaD * beta0 * outerNormal.transpose();
                //         gdT(i) += omg * (ci(0) * violaPosPenaD * alpha * outerNormal.dot(vel) * step +
                //                          ci(0) * violaPosPena / cons(i));
                //         current_pena = omg * step * violaPosPena;
                //         pena += current_pena * ci(0);
                //         pena_sta += current_pena;
                //     }
                // }
                if (current_samp <= total_samp_45)
                {
                    if (opt_iter_num % 10 == 0)
                    {
                        if (MLmap->getOccupancy(pos, safeMargin) == MLmap->OCCUPIED)
                        {
                            getFoot(i, pos, foot, path);
                            contorl_pairs.emplace_back(current_samp, foot);
                            // path[i] path[i + 1]
                        }
                    }

                    else if (!in_global_bound(pos, safeMargin))
                    {
                        getFoot(i, pos, foot, path);
                        contorl_pairs.emplace_back(current_samp, foot);
                        // cout<<"out of the global bound: "<<pos.transpose()<<endl;
                        // path[i] path[i + 1]
                    }
                }
                if (contorl_pairs.size() > cur_id_check && (current_samp) == contorl_pairs[cur_id_check].first)
                {
                    direction = pos - contorl_pairs[cur_id_check].second;
                    dist2line = direction.norm();
                    // double dist2line2 = dist2line*dist2line;
                    double costp = ci(0) * pow(dist2line, 2);
                    pena_sta += costp * omg * step / ci(0);
                    pena += costp * omg * step;
                    gradP = ci(0) * 2 * direction; // cost increase along gradient direction
                    gdC.block<6, 3>(i * 6, 0) += beta0 * gradP.transpose() * omg * step;
                    gradViolaPt = alpha * gradP.transpose() * vel;
                    gdT(i) += omg * (costp / innerLoop + step * gradViolaPt);
                    cur_id_check++;
                }
                if (violaThrh > 0.0)
                {
                    positiveSmoothedL1(violaThrh, violaThrhPena, violaThrhPenaD);
                    gradViolaThrhc = beta2 * dSqrMagThr.transpose();
                    gradViolaThrht = alpha * dSqrMagThr.transpose() * jer;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(2) * violaThrhPenaD * gradViolaThrhc;
                    gdT(i) += omg * (ci(2) * violaThrhPenaD * gradViolaThrht * step +
                                     ci(2) * violaThrhPena / cons(i));
                    pena += omg * step * ci(2) * violaThrhPena;
                }

                for (int m = 0; m < dynobs_pointer->dyn_number && current_samp <= total_samp_45; m++)
                {
                    if (m == 0)
                    {
                        t_gap = t_now - dynobs_pointer->time_stamp;
                    }
                    // t_gap = t_now - dynobs_pointer->time_stamp;
                    // cout << "mk321-1" <<endl<<t_gap<<endl;
                    double obj_prop_conv = pow(dynobs_pointer->max_accs[m](1) + dynobs_pointer->max_accs[m](2) * t_gap * t_gap, 0.5);
                    obj_prop_conv = obj_prop_conv > 0.5 ? 0.5 : obj_prop_conv;
                    Eigen::Vector3d conv_vec = {obj_prop_conv, obj_prop_conv, 0.0};
                    // conv_vec *= 0.5;
                    Eigen::Vector3d acc_vec = {dynobs_pointer->max_accs[m](0) * t_gap * t_gap / 2, dynobs_pointer->max_accs[m](0) * t_gap * t_gap / 2, 0};
                    ct_center = dynobs_pointer->centers[m] + t_gap * dynobs_pointer->vels[m];
                    Eigen::Vector3d check_vec = ((ct_center - pos).cwiseAbs() - dynobs_pointer->obs_sizes[m] * 0.5 - conv_vec);
                    Eigen::Vector3d check_vec_acc = check_vec - acc_vec;
                    double grad_prev_t = 0.0;
                    if ((check_vec_acc.array() < 0.0).all())
                    {
                        // Eigen::Vector3d dist_vec = pos - ct_center;
                        //   cout << "max_accs" <<endl<<dynobs_pointer->max_accs[m]<<endl;
                        Eigen::Vector3d dist_vec = pos - ct_center;
                        Eigen::Vector3d dJ_dP;
                        if ((check_vec.array() < safeMargin).all() && !mask_dyn_idx[m]) //
                        {
                            mask_dyn_idx[m] = true;
                            Eigen::Vector3d half_len = (dynobs_pointer->obs_sizes[m] * 0.5 + conv_vec).array() + safeMargin + 0.05;
                            Eigen::Vector3d axis = {dist_vec(0) > 0 ? ct_center(0) + half_len(0) : ct_center(0) - half_len(0),
                                                    dist_vec(1) > 0 ? ct_center(1) + half_len(1) : ct_center(1) - half_len(1),
                                                    dist_vec(2) > 0 ? ct_center(2) + half_len(2) : ct_center(2) - half_len(2)};
                            Eigen::Vector3d dist = pos - axis;
                            double dist_err2 = dist.squaredNorm();
                            current_pena = dist_err2 * omg * step;
                            pena += wei_dyn * current_pena;
                            pena_ball += current_pena;
                            dJ_dP = omg * step * wei_dyn * 2 * dist; // gradient!
                            gdC.block<6, 3>(i * 6, 0) += beta0 * dJ_dP.transpose();
                            gdT(i) += omg * wei_dyn * dist_err2 / innerLoop + dJ_dP.dot(vel) * j / innerLoop;

                            // original cost and gradient:
                            //  double ellip_dist2 = dist_vec(2) * dist_vec(2) * inv_z + dist_vec(0) * dist_vec(0) * inv_x + dist_vec(1) * dist_vec(1) * inv_y;
                            //  double dist2_err = 2.5 - ellip_dist2;
                            //  double dist2_err2 = dist2_err * dist2_err;
                            //  double dist2_err3 = dist2_err2 * dist2_err;
                            //  current_pena = wei_dyn * dist2_err3 * omg * step;
                            //  pena += current_pena;
                            //  pena_ball += current_pena;
                            //  dJ_dP = omg * step * wei_dyn * 3 * dist2_err2 * (-2) * Eigen::Vector3d(inv_x * dist_vec(0), inv_y * dist_vec(1), inv_z * dist_vec(2)); // gradient!
                            //  gdC.block<6, 3>(i * 6, 0) += beta0 * dJ_dP.transpose();
                            //  // gdT(i) += omg * wei_dyn * dist2_err3 / innerLoop + dJ_dP.dot(vel - dynobs_pointer->vels[m]) * j / innerLoop;
                            //  gdT(i) += omg * wei_dyn * dist2_err3 / innerLoop + dJ_dP.dot(vel ) * j / innerLoop;
                            grad_prev_t += dJ_dP.dot(-dynobs_pointer->vels[m]);
                        }
                        else if (t_gap < 3.0)
                        {
                            Eigen::Vector3d half_len = dynobs_pointer->obs_sizes[m] * 0.5 + conv_vec + acc_vec;
                            inv_z = 1 / (half_len(2) * half_len(2));
                            inv_x = 1 / (half_len(0) * half_len(0));
                            inv_y = 1 / (half_len(1) * half_len(1));
                            double ellip_dist2 = dist_vec(2) * dist_vec(2) * inv_z + dist_vec(0) * dist_vec(0) * inv_x + dist_vec(1) * dist_vec(1) * inv_y;
                            double dist2_err = 2.5 - ellip_dist2;
                            double dist2_err2 = dist2_err * dist2_err;
                            double dist2_err3 = dist2_err2 * dist2_err;
                            pena += wei_dyn_acc * dist2_err3 * omg * step;
                            pena_acc += wei_dyn_acc * dist2_err3 * omg * step;
                            dJ_dP = wei_dyn_acc * 3 * dist2_err2 * (-2) * Eigen::Vector3d(inv_x * dist_vec(0), inv_y * dist_vec(1), inv_z * dist_vec(2)); // gradient!
                            gdC.block<6, 3>(i * 6, 0) += beta0 * dJ_dP.transpose() * omg * step;
                            // gdT(i) += omg * (wei_dyn_acc * dist2_err3 / innerLoop + step * dJ_dP.dot(vel - dynobs_pointer->vels[m]) * j / innerLoop);
                            gdT(i) += omg * (wei_dyn_acc * dist2_err3 / innerLoop + step * dJ_dP.dot(vel) * j / innerLoop);
                            grad_prev_t += dJ_dP.dot(-dynobs_pointer->vels[m]);
                        }

                        if (i > 0)
                        {
                            gdT.head(i).array() += omg * step * grad_prev_t;
                        }
                    }
                    // cout << "mk322" <<endl;
                }

                for (int m = 0; m < dynobs_pointer->ball_number; m++)
                {
                    if (m == 0)
                    {
                        t_gap = t_now - dynobs_pointer->ball_time_stamp;
                    }
                    // t_gap = t_now - dynobs_pointer->time_stamp;
                    // cout << "mk321" <<endl<<t_gap<<endl<<m<<"  j: "<<j<<"  innerloop:  "<<innerLoop<<endl;
                    ct_center = dynobs_pointer->ballpos[m] + t_gap * dynobs_pointer->ballvel[m] + 0.5 * t_gap * t_gap * dynobs_pointer->ballacc[m];
                    Eigen::Vector3d check_vec = ((ct_center - pos).cwiseAbs() - dynobs_pointer->ball_sizes[m] * 0.5);
                    if ((ct_center - pos).norm() < dynobs_pointer->ball_sizes[m](0) * 0.5) //((check_vec.array()<0.0).all())
                    {
                        double sa = dynobs_pointer->ball_sizes[m].squaredNorm() / 4;
                        // Eigen::Vector3d dist_vec = pos - ct_center;
                        //  cout<<"check_vec,ct_center,pos:\n"<<check_vec<<",   \n"<<ct_center<<",   \n"<<pos<<endl<<t_gap<<endl<<dynobs_pointer->ball_sizes[m](0)*0.5<<endl<<(ct_center - pos).norm()<<endl;
                        Eigen::Vector3d dist_vec = pos - ct_center;
                        double ellip_dist2 = dist_vec(2) * dist_vec(2) * inv_a2 + (dist_vec(0) * dist_vec(0) + dist_vec(1) * dist_vec(1)) * inv_b2;
                        double dist2_err = sa - ellip_dist2;
                        double dist2_err2 = dist2_err * dist2_err;
                        double dist2_err3 = dist2_err2 * dist2_err;
                        current_pena = wei_ball * dist2_err3 * omg * step;
                        pena += current_pena;
                        pena_ball += current_pena;
                        Eigen::Vector3d dJ_dP = wei_ball * 3 * dist2_err2 * (-2) * Eigen::Vector3d(inv_b2 * dist_vec(0), inv_b2 * dist_vec(1), inv_a2 * dist_vec(2)); // gradient!
                        gdC.block<6, 3>(i * 6, 0) += beta0 * dJ_dP.transpose() * omg * step;
                        gdT(i) += omg * (wei_ball * dist2_err3 / innerLoop + step * dJ_dP.dot(vel - dynobs_pointer->ballvel[m]) * j / innerLoop);
                        // cout<<"gradient: "<<beta0 * dJ_dP.transpose()* omg * step<<"\n"<<omg * (wei_ball * dist2_err3 / innerLoop + step * dJ_dP.dot(vel - dynobs_pointer->ballvel[m])*j/innerLoop)<<endl;
                        double grad_prev_t = dJ_dP.dot(-dynobs_pointer->ballvel[m]);
                        if (i > 0)
                        {
                            gdT.head(i).array() += omg * step * grad_prev_t;
                        }
                    }
                    // cout << "mk322" <<endl;
                }
                violaVel = vel.squaredNorm() - vMaxSqr;
                if (violaVel > 0.0)
                {
                    
                    positiveSmoothedL1(violaVel, violaVelPena, violaVelPenaD);
                    // cout<<"Add vel pena: "<< violaVel<<"  "<<violaVelPena<<endl;
                    gradViolaVc = 2.0 * beta1 * vel.transpose();
                    gradViolaVt = 2.0 * alpha * vel.transpose() * acc;
                    gdC.block<6, 3>(i * 6, 0) += omg * step * ci(1) * violaVelPenaD * gradViolaVc;
                    gdT(i) += omg * ci(1) * (violaVelPenaD * gradViolaVt * step + violaVelPena / cons(i));
                    pena += omg * step * ci(1) * violaVelPena;
                }
                s1 += step;
            }
        }

        cost += pena;
        opt_iter_num++;
        return;
    }

public:
    inline void reset(const Eigen::Matrix3d &headState,
                      const Eigen::Matrix3d &tailState,
                      const int &pieceNum)
    {
        N = pieceNum;
        headPVA = headState;
        tailPVA = tailState;
        T1.resize(N);
        A.create(6 * N, 6, 6);
        b.resize(6 * N, 3);
        gdC.resize(6 * N, 3);
        return;
    }

    inline Eigen::MatrixXd getInitConstraintPoints(const int K) const
    {
        Eigen::MatrixXd pts(3, N * K + 1);
        Eigen::Vector3d pos;
        Eigen::Matrix<double, 6, 1> beta0;
        double s1, s2, s3, s4, s5;
        double step;
        int i_dp = 0;

        for (int i = 0; i < N; ++i)
        {
            const auto &c = b.block<6, 3>(i * 6, 0);
            step = T1(i) / K;
            s1 = 0.0;
            double t = 0;
            // innerLoop = K;

            for (int j = 0; j <= K; ++j)
            {
                s2 = s1 * s1;
                s3 = s2 * s1;
                s4 = s2 * s2;
                s5 = s4 * s1;
                beta0 << 1.0, s1, s2, s3, s4, s5;
                pos = c.transpose() * beta0;
                pts.col(i_dp) = pos;

                s1 += step;
                if (j != K || (j == K && i == N - 1))
                {
                    ++i_dp;
                }
            }
        }

        return pts;
    }

    inline void generate(const Eigen::MatrixXd &inPs,
                         const Eigen::VectorXd &ts)
    {
        T1 = ts;
        T2 = T1.cwiseProduct(T1);
        T3 = T2.cwiseProduct(T1);
        T4 = T2.cwiseProduct(T2);
        T5 = T4.cwiseProduct(T1);

        A.reset();
        b.setZero();

        A(0, 0) = 1.0;
        A(1, 1) = 1.0;
        A(2, 2) = 2.0;
        b.row(0) = headPVA.col(0).transpose();
        b.row(1) = headPVA.col(1).transpose();
        b.row(2) = headPVA.col(2).transpose();

        for (int i = 0; i < N - 1; i++)
        {
            A(6 * i + 3, 6 * i + 3) = 6.0;
            A(6 * i + 3, 6 * i + 4) = 24.0 * T1(i);
            A(6 * i + 3, 6 * i + 5) = 60.0 * T2(i);
            A(6 * i + 3, 6 * i + 9) = -6.0;
            A(6 * i + 4, 6 * i + 4) = 24.0;
            A(6 * i + 4, 6 * i + 5) = 120.0 * T1(i);
            A(6 * i + 4, 6 * i + 10) = -24.0;
            A(6 * i + 5, 6 * i) = 1.0;
            A(6 * i + 5, 6 * i + 1) = T1(i);
            A(6 * i + 5, 6 * i + 2) = T2(i);
            A(6 * i + 5, 6 * i + 3) = T3(i);
            A(6 * i + 5, 6 * i + 4) = T4(i);
            A(6 * i + 5, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i) = 1.0;
            A(6 * i + 6, 6 * i + 1) = T1(i);
            A(6 * i + 6, 6 * i + 2) = T2(i);
            A(6 * i + 6, 6 * i + 3) = T3(i);
            A(6 * i + 6, 6 * i + 4) = T4(i);
            A(6 * i + 6, 6 * i + 5) = T5(i);
            A(6 * i + 6, 6 * i + 6) = -1.0;
            A(6 * i + 7, 6 * i + 1) = 1.0;
            A(6 * i + 7, 6 * i + 2) = 2 * T1(i);
            A(6 * i + 7, 6 * i + 3) = 3 * T2(i);
            A(6 * i + 7, 6 * i + 4) = 4 * T3(i);
            A(6 * i + 7, 6 * i + 5) = 5 * T4(i);
            A(6 * i + 7, 6 * i + 7) = -1.0;
            A(6 * i + 8, 6 * i + 2) = 2.0;
            A(6 * i + 8, 6 * i + 3) = 6 * T1(i);
            A(6 * i + 8, 6 * i + 4) = 12 * T2(i);
            A(6 * i + 8, 6 * i + 5) = 20 * T3(i);
            A(6 * i + 8, 6 * i + 8) = -2.0;

            b.row(6 * i + 5) = inPs.col(i).transpose();
        }

        A(6 * N - 3, 6 * N - 6) = 1.0;
        A(6 * N - 3, 6 * N - 5) = T1(N - 1);
        A(6 * N - 3, 6 * N - 4) = T2(N - 1);
        A(6 * N - 3, 6 * N - 3) = T3(N - 1);
        A(6 * N - 3, 6 * N - 2) = T4(N - 1);
        A(6 * N - 3, 6 * N - 1) = T5(N - 1);
        A(6 * N - 2, 6 * N - 5) = 1.0;
        A(6 * N - 2, 6 * N - 4) = 2 * T1(N - 1);
        A(6 * N - 2, 6 * N - 3) = 3 * T2(N - 1);
        A(6 * N - 2, 6 * N - 2) = 4 * T3(N - 1);
        A(6 * N - 2, 6 * N - 1) = 5 * T4(N - 1);
        A(6 * N - 1, 6 * N - 4) = 2;
        A(6 * N - 1, 6 * N - 3) = 6 * T1(N - 1);
        A(6 * N - 1, 6 * N - 2) = 12 * T2(N - 1);
        A(6 * N - 1, 6 * N - 1) = 20 * T3(N - 1);

        b.row(6 * N - 3) = tailPVA.col(0).transpose();
        b.row(6 * N - 2) = tailPVA.col(1).transpose();
        b.row(6 * N - 1) = tailPVA.col(2).transpose();

        A.factorizeLU();
        A.solve(b);

        return;
    }

    inline double getTrajJerkCost() const
    {
        double objective = 0.0;
        for (int i = 0; i < N; i++)
        {
            objective += 36.0 * b.row(6 * i + 3).squaredNorm() * T1(i) +
                         144.0 * b.row(6 * i + 4).dot(b.row(6 * i + 3)) * T2(i) +
                         192.0 * b.row(6 * i + 4).squaredNorm() * T3(i) +
                         240.0 * b.row(6 * i + 5).dot(b.row(6 * i + 3)) * T3(i) +
                         720.0 * b.row(6 * i + 5).dot(b.row(6 * i + 4)) * T4(i) +
                         720.0 * b.row(6 * i + 5).squaredNorm() * T5(i);
        }
        return objective;
    }

    template <typename EIGENVEC, typename EIGENMAT>
    inline void evalTrajCostGrad(const Eigen::VectorXi &cons,
                                 const Eigen::VectorXi &idxHs,
                                 const std::vector<Eigen::MatrixXd> &cfgHs,
                                 const Eigen::Vector3d &ellipsoid,
                                 const double &safeMargin,
                                 const double &vMax,
                                 const double &thrAccMin,
                                 const double &thrAccMax,
                                 const double &bdrMax,
                                 const double &gAcc,
                                 const Eigen::Vector4d &ci,
                                 double &cost,
                                 EIGENVEC &gdT,
                                 EIGENMAT &gdInPs,
                                 const double plan_t,
                                 const dynobs_tmp *dynobs_pointer,
                                 const std::vector<Eigen::Vector3d> &path)
    {
        gdT.setZero();
        gdInPs.setZero();
        gdC.setZero();

        cost = getTrajJerkCost();
        // cout << "cost for jerk: " << cost <<endl;  //<<"b:"<<b<<endl<<"T1: "<<T1<<endl;
        addGradJbyT(gdT);
        addGradJbyC(gdC);
        // cout << "mk32,cost: " <<cost<<endl;
        addTimeIntPenalty(cons, idxHs, cfgHs, ellipsoid, safeMargin,
                          vMax, thrAccMin, thrAccMax, bdrMax,
                          gAcc, ci, cost, gdT, gdC, plan_t, dynobs_pointer,
                          path);

        // calculate the cost and gradient w.r.t the time allocation T and polynomial param matrix C

        solveAdjGradC(gdC);
        addPropCtoT(gdC, gdT);
        addPropCtoP(gdC, gdInPs);
    }

    inline Trajectory getTraj(void) const
    {
        Trajectory traj;
        traj.reserve(N);
        for (int i = 0; i < N; i++)
        {
            traj.emplace_back(T1(i), b.block<6, 3>(6 * i, 0).transpose().rowwise().reverse());
        }
        return traj;
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

class POLY_OPT
{
private:
    // Use C2 or Cinf diffeo
    bool c2dfm;

    // Use soft time or not
    bool softT;

    // Weight for time regularization term
    double rho;

    // Fixed total time
    double sumT;

    // Minimum Jerk Optimizer
    MINCO_S3_G jerkOpt;
    mlmap::Ptr MLmap;
    // Temp variables for problem solving
    Eigen::MatrixXd iState;
    Eigen::MatrixXd fState;
    Eigen::MatrixXd global_bound_;
    // Each col of cfgHs denotes a facet (outter_normal^T,point^T)^T
    std::vector<Eigen::Vector3d> path;
    std::vector<Eigen::Vector3d> cld_points;
    std::vector<Eigen::MatrixXd> cfgVs;
    std::vector<Eigen::MatrixXd> cfgHs;
    Eigen::MatrixXd gdInPs;

    // Piece num for each polytope
    Eigen::VectorXi intervals;
    // Assignment vector for point in V-polytope
    Eigen::VectorXi idxVs;
    // Assignment vector for piece in H-polytope
    Eigen::VectorXi idxHs;

    int coarseN;
    int fineN;
    int dimFreeT;
    int dimFreeP;
    Eigen::VectorXd coarseT;
    Eigen::VectorXd fineT;
    Eigen::MatrixXd innerP;

    // Params for constraints
    Eigen::VectorXi cons;
    Eigen::Vector4d chi;

    Eigen::Vector3d ellipsoid;
    double safeMargin;
    double vMax;
    double thrAccMin;
    double thrAccMax;
    double bdrMax;
    double gAcc;

    // L-BFGS Solver Parameters
    lbfgs::lbfgs_parameter_t lbfgs_params;

    // for dynamic obstacles
    double plan_t;
    const dynobs_tmp *dynobs_pointer;

private:
    template <typename EIGENVEC>
    static inline void forwardT(const EIGENVEC &t,
                                Eigen::VectorXd &vecT,
                                bool soft,
                                const double &sT,
                                bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = vecT.size();
                for (int i = 0; i < M; i++)
                {
                    vecT(i) = t(i) > 0.0
                                  ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                                  : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
            }
            else
            {
                vecT = t.array().exp();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                for (int i = 0; i < Ms1; i++)
                {
                    vecT(i) = t(i) > 0.0
                                  ? ((0.5 * t(i) + 1.0) * t(i) + 1.0)
                                  : 1.0 / ((0.5 * t(i) - 1.0) * t(i) + 1.0);
                }
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            }
            else
            {
                int Ms1 = t.size();
                vecT.head(Ms1) = t.array().exp();
                vecT(Ms1) = 0.0;
                vecT /= 1.0 + vecT.sum();
                vecT(Ms1) = 1.0 - vecT.sum();
                vecT *= sT;
            }
        }
        return;
    }

    template <typename EIGENVEC>
    static inline void backwardT(const Eigen::VectorXd &vecT,
                                 EIGENVEC &t,
                                 bool soft,
                                 bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = vecT.size();
                for (int i = 0; i < M; i++)
                {
                    t(i) = vecT(i) > 1.0
                               ? (sqrt(2.0 * vecT(i) - 1.0) - 1.0)
                               : (1.0 - sqrt(2.0 / vecT(i) - 1.0));
                    /*

                    t  = 2/[(tau-1)^2+1] tau<0
                       = 1/2[(tau+1)^2+1] tau>=0
                    */
                }
            }
            else
            {
                t = vecT.array().log();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                t = vecT.head(Ms1) / vecT(Ms1);
                for (int i = 0; i < Ms1; i++)
                {
                    t(i) = t(i) > 1.0
                               ? (sqrt(2.0 * t(i) - 1.0) - 1.0)
                               : (1.0 - sqrt(2.0 / t(i) - 1.0));
                }
            }
            else
            {
                int Ms1 = t.size();
                t = (vecT.head(Ms1) / vecT(Ms1)).array().log();
            }
        }
        return;
    }

    template <typename EIGENVEC>
    static inline void forwardP(const EIGENVEC &p,
                                const Eigen::VectorXi &idVs,
                                const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                Eigen::MatrixXd &inP)
    {
        int M = inP.cols();
        Eigen::VectorXd q;
        int j = 0, k, idx;
        for (int i = 0; i < M; i++)
        {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;
            q = 2.0 / (1.0 + p.segment(j, k).squaredNorm()) * p.segment(j, k);
            inP.col(i) = cfgPolyVs[idx].rightCols(k) * q.cwiseProduct(q) +
                         cfgPolyVs[idx].col(0); // q.cwiseProduct(q) is all positive, so the recovered joint point inP must be inside the intersection part
            j += k;
        }
        return;
    }

    static inline double objectiveNLS(void *ptrPOBs,
                                      const double *x,
                                      double *grad,
                                      const int n)
    {
        const Eigen::MatrixXd &pobs = *(Eigen::MatrixXd *)ptrPOBs;
        Eigen::Map<const Eigen::VectorXd> p(x, n);
        Eigen::Map<Eigen::VectorXd> gradp(grad, n);

        double qnsqr = p.squaredNorm();
        double qnsqrp1 = qnsqr + 1.0;
        double qnsqrp1sqr = qnsqrp1 * qnsqrp1;
        Eigen::VectorXd r = 2.0 / qnsqrp1 * p;

        Eigen::Vector3d delta = pobs.rightCols(n) * r.cwiseProduct(r) +
                                pobs.col(1) - pobs.col(0);
        double cost = delta.squaredNorm();
        Eigen::Vector3d gradR3 = 2 * delta;

        Eigen::VectorXd gdr = pobs.rightCols(n).transpose() * gradR3;
        gdr = gdr.array() * r.array() * 2.0;
        gradp = gdr * 2.0 / qnsqrp1 -
                p * 4.0 * gdr.dot(p) / qnsqrp1sqr;

        return cost;
    }

    template <typename EIGENVEC>
    static inline void backwardP(const Eigen::MatrixXd &inP,
                                 const Eigen::VectorXi &idVs,
                                 const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                 EIGENVEC &p)
    {
        int M = inP.cols();
        int j = 0, k, idx;

        // Parameters for tiny nonlinear least squares
        double minSqrD;
        lbfgs::lbfgs_parameter_t nls_params;
        lbfgs::lbfgs_load_default_parameters(&nls_params);
        nls_params.g_epsilon = FLT_EPSILON;
        nls_params.max_iterations = 128;

        Eigen::MatrixXd pobs;
        for (int i = 0; i < M; i++)
        {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;
            p.segment(j, k).setConstant(1.0 / (sqrt(k + 1.0) + 1.0));
            pobs.resize(3, k + 2);
            pobs << inP.col(i), cfgPolyVs[idx];
            lbfgs::lbfgs_optimize(k,
                                  p.data() + j,
                                  &minSqrD,
                                  &POLY_OPT::objectiveNLS,
                                  nullptr,
                                  nullptr,
                                  &pobs,
                                  &nls_params);

            j += k;
        }

        return;
    }

    template <typename EIGENVEC>
    static inline void addLayerTGrad(const Eigen::VectorXd &t,
                                     EIGENVEC &gradT,
                                     bool soft,
                                     const double &sT,
                                     bool c2)
    {
        if (soft)
        {
            if (c2)
            {
                int M = t.size();
                double denSqrt;
                for (int i = 0; i < M; i++)
                {
                    if (t(i) > 0)
                    {
                        gradT(i) *= t(i) + 1.0;
                    }
                    else
                    {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        gradT(i) *= (1.0 - t(i)) / (denSqrt * denSqrt);
                    }
                }
            }
            else
            {
                int M = t.size();
                gradT.head(M).array() *= t.array().exp();
            }
        }
        else
        {
            if (c2)
            {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd dExpTau(Ms1);
                double expTauSum = 0.0, gFreeDotExpTau = 0.0;
                double denSqrt, expTau;
                for (int i = 0; i < Ms1; i++)
                {
                    if (t(i) > 0)
                    {
                        expTau = (0.5 * t(i) + 1.0) * t(i) + 1.0;
                        dExpTau(i) = t(i) + 1.0;
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    }
                    else
                    {
                        denSqrt = (0.5 * t(i) - 1.0) * t(i) + 1.0;
                        expTau = 1.0 / denSqrt;
                        dExpTau(i) = (1.0 - t(i)) / (denSqrt * denSqrt);
                        expTauSum += expTau;
                        gFreeDotExpTau += expTau * gFree(i);
                    }
                }
                denSqrt = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * dExpTau.array() / denSqrt -
                                  (gFreeDotExpTau - gTail * expTauSum) * dExpTau.array() / (denSqrt * denSqrt);
                gradT(Ms1) = 0.0;
            }
            else
            {
                int Ms1 = t.size();
                Eigen::VectorXd gFree = sT * gradT.head(Ms1);
                double gTail = sT * gradT(Ms1);
                Eigen::VectorXd expTau = t.array().exp();
                double expTauSum = expTau.sum();
                double denom = expTauSum + 1.0;
                gradT.head(Ms1) = (gFree.array() - gTail) * expTau.array() / denom -
                                  (gFree.dot(expTau) - gTail * expTauSum) * expTau.array() / (denom * denom);
                gradT(Ms1) = 0.0;
            }
        }
        return;
    }

    template <typename EIGENVEC_0, typename EIGENVEC_1>
    static inline void addLayerPGrad(EIGENVEC_0 &p,
                                     const Eigen::VectorXi &idVs,
                                     const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                                     const Eigen::MatrixXd &gradInPs,
                                     EIGENVEC_1 &grad)
    {
        int M = gradInPs.cols();

        int j = 0, k, idx;
        double qnsqr, qnsqrp1, qnsqrp1sqr;
        Eigen::VectorXd q, r, gdr;
        for (int i = 0; i < M; i++)
        {
            idx = idVs(i);
            k = cfgPolyVs[idx].cols() - 1;

            q = p.segment(j, k);
            qnsqr = q.squaredNorm();
            qnsqrp1 = qnsqr + 1.0;
            qnsqrp1sqr = qnsqrp1 * qnsqrp1; // the same in the backwardP()
            r = 2.0 / qnsqrp1 * q;
            gdr = cfgPolyVs[idx].rightCols(k).transpose() * gradInPs.col(i);
            gdr = gdr.array() * r.array() * 2.0;

            grad.segment(j, k) = gdr * 2.0 / qnsqrp1 -
                                 q * 4.0 * gdr.dot(q) / qnsqrp1sqr;

            j += k;
        }

        return;
    }

    static inline void splitToFineT(const Eigen::VectorXd &cT,
                                    const Eigen::VectorXi &intervs,
                                    Eigen::VectorXd &fT)
    {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++)
        {
            inverv = intervs(i);
            fT.segment(offset, inverv).setConstant(cT(i) / inverv);
            offset += inverv;
        }
        return;
    }

    static inline void mergeToCoarseGradT(const Eigen::VectorXi &intervs,
                                          Eigen::VectorXd &fineGdT)
    {
        int M = intervs.size();
        int offset = 0;
        int inverv;
        for (int i = 0; i < M; i++)
        {
            inverv = intervs(i);
            fineGdT(i) = fineGdT.segment(offset, inverv).mean();
            offset += inverv;
        }
        return;
    }

    static inline double objectiveFunc(void *ptrObj,
                                       const double *x,
                                       double *grad,
                                       const int n)
    {
        POLY_OPT &obj = *(POLY_OPT *)ptrObj;
        const int dimT = obj.dimFreeT;
        const int dimP = obj.dimFreeP;
        const double rh = obj.rho;
        Eigen::Map<const Eigen::VectorXd> t(x, dimT);
        Eigen::Map<const Eigen::MatrixXd> p(x + dimT, 3, obj.coarseN - 1); // add dimT to pointer means to shift the map start pos.
        Eigen::Map<Eigen::VectorXd> gradt(grad, dimT);
        Eigen::VectorXd proxyGradT(obj.fineN);
        Eigen::Map<Eigen::VectorXd> gradp(grad + dimT, dimP);
        // cout<<"gradt: "<<gradt<<endl; //<<"gradp:"<<gradp<<endl
        forwardT(t, obj.coarseT, obj.softT, obj.sumT, obj.c2dfm);
        obj.fineT = obj.coarseT;
        // splitToFineT(obj.coarseT, obj.intervals, obj.fineT);  //evenly divide the coarse T into fine T
        // p.resize(3, obj.coarseN - 1);
        obj.innerP = p;
        // forwardP(p, obj.idxVs, obj.cfgVs, obj.innerP);
        // use vertexes of the intersection Polyh to represent the traj joint point, so it is constrained inside the intersection part strictly.
        // p is the weights vector of each vertex
        obj.jerkOpt.cost_total = 0;
        obj.jerkOpt.generate(obj.innerP, obj.fineT);
        // cout << "mk31" <<endl;
        // get_collide_points();
        obj.jerkOpt.evalTrajCostGrad(obj.cons, obj.idxHs, obj.cfgHs, obj.ellipsoid,
                                     obj.safeMargin, obj.vMax, obj.thrAccMin,
                                     obj.thrAccMax, obj.bdrMax, obj.gAcc, obj.chi,
                                     obj.jerkOpt.cost_total, proxyGradT, obj.gdInPs,
                                     obj.plan_t, obj.dynobs_pointer, obj.path);

        obj.jerkOpt.cost_total += rh * obj.coarseT.sum();
        proxyGradT.array() += rh;

        // mergeToCoarseGradT(obj.intervals, proxyGradT);
        // grad of T
        // T,P->tau kesi
        addLayerTGrad(t, proxyGradT, obj.softT, obj.sumT, obj.c2dfm);
        // addLayerPGrad(p, obj.idxVs, obj.cfgVs, obj.gdInPs, gradp);  //convert the gradient w.r.t of the joint points to the actual optimized variables p.
        Eigen::MatrixXd gdInPs = obj.gdInPs;
        gdInPs.resize(dimP, 1);
        gradp = gdInPs;
        gradt = proxyGradT.head(dimT);
        // cout << "final cost: " << cost <<endl;
        return obj.jerkOpt.cost_total;
    }

public:
    inline void gridMesh(const Eigen::Matrix3d &iState,
                         const Eigen::Matrix3d &fState,
                         const std::vector<Eigen::MatrixXd> &cfgPolyVs,
                         const double &gridResolution,
                         Eigen::VectorXi &intervalsVec) const
    {
        int M = intervalsVec.size();

        int curInterval, k;
        Eigen::Vector3d lastP, curP;
        curP = iState.col(0);
        for (int i = 0; i < M - 1; i++)
        {
            lastP = curP;
            k = cfgPolyVs[2 * i + 1].cols() - 1;
            curP = cfgPolyVs[2 * i + 1].rightCols(k).rowwise().sum() / (1.0 + k) +
                   cfgPolyVs[2 * i + 1].col(0);
            curInterval = ceil((curP - lastP).norm() / gridResolution);
            intervalsVec(i) = curInterval > 0 ? curInterval : 1;
        }
        lastP = curP;
        curP = fState.col(0);
        curInterval = ceil((curP - lastP).norm() / gridResolution);
        intervalsVec(M - 1) = curInterval > 0 ? curInterval : 1;

        return;
    }

    inline bool setup(const double &rh,
                      const double &st,
                      const Eigen::MatrixXd &iniState,
                      const Eigen::MatrixXd &finState,
                      const mlmap::Ptr map,
                      std::vector<Eigen::Vector3d> FE_path,
                      const int &itgSpaces,
                      const double &margin, // 0
                      const double &vm,
                      const double &minThrAcc,
                      const double &maxThrAcc,
                      const double &bodyRateMax,
                      const double &g,
                      const Eigen::Vector4d &w,
                      bool c2diffeo,
                      const double plan_start_t,
                      const dynobs_tmp *dyn_pointer,
                      const Eigen::MatrixXd &global_bound)
    {
        // Setup for optimization parameters
        cout << "To setup" << endl;
        if (FE_path.size() < 3)
        {
            cout << "front-end waypoints number < 3!" << endl;
            if (FE_path.size() == 2)
                FE_path.insert(FE_path.begin() + 1, (FE_path[0] + FE_path[1]) / 2);
            else
                return false;
        }
        cout << "Initialize MINCO traj with waypoints: " << endl;
        for (auto pt : FE_path)
        {
            cout << pt.transpose() << endl;
        }
        MLmap = map;
        iState = iniState;
        fState = finState;
        plan_t = plan_start_t;
        dynobs_pointer = dyn_pointer;
        path = FE_path;
        coarseN = FE_path.size() - 1;

        intervals.resize(coarseN);

        c2dfm = c2diffeo;

        softT = rh > 0;
        if (softT)
        {
            rho = rh;
            sumT = 1.0; // optimize total time
        }
        else
        {
            rho = 0.0;
            sumT = st;
        }

        // cout << dynobs_pointer << endl;

        fineN = coarseN;
        cons.resize(fineN);
        cons.setConstant(itgSpaces); // how many sample points per segment to evaluate the cost

        dimFreeT = softT ? coarseN : coarseN - 1; // softT is true tau size is n
        dimFreeP = 3 * (coarseN - 1);

        chi = w;
        safeMargin = margin;
        vMax = vm;
        thrAccMin = minThrAcc;
        thrAccMax = maxThrAcc;
        bdrMax = bodyRateMax;
        gAcc = g;

        // Make a legal initial speed
        double tempNorm;
        tempNorm = iState.col(1).norm();
        iState.col(1) *= tempNorm > vMax ? (vMax / tempNorm - 0.01) : 1.0;
        tempNorm = fState.col(1).norm();
        fState.col(1) *= tempNorm > vMax ? (vMax / tempNorm - 0.01) : 1.0;

        // Setup for L-BFGS solver
        lbfgs::lbfgs_load_default_parameters(&lbfgs_params);

        // Allocate temp variables
        coarseT.resize(coarseN);
        fineT.resize(fineN);
        innerP.resize(3, fineN - 1);
        gdInPs.resize(3, fineN - 1);
        jerkOpt.reset(iniState, finState, fineN);
        jerkOpt.MLmap = this->MLmap;
        jerkOpt.global_bound_ = global_bound;
        //
        // hzcmy load param

        return true;
    }

    template <typename EIGENVEC>
    inline void setInitial(
        Eigen::VectorXd &vecT,
        EIGENVEC &vecInP) const
    {
        constexpr double maxSpeedForAllocatiion = 10.0;

        int M = vecT.size();
        Eigen::Vector3d lastP, curP, delta;
        curP = iState.col(0);
        for (int i = 0; i < M - 1; i++)
        {
            vecInP.col(i) = path[i + 1];
        }
        for (int i = 0; i < M - 1; i++)
        {
            lastP = curP;
            curP = vecInP.col(i);
            delta = curP - lastP;
            vecT(i) = delta.norm() / std::min(vMax, maxSpeedForAllocatiion);
        }
        lastP = curP;
        curP = fState.col(0);
        delta = curP - lastP;
        vecT(M - 1) = delta.norm() / std::min(vMax, maxSpeedForAllocatiion);
        return;
    }

    inline void get_collide_points(void)
    {
        ;
        // TODO;
    }
    inline double optimize(Trajectory &traj,
                           const double &relCostTol)
    {
        double *x = new double[dimFreeT + dimFreeP];
        Eigen::Map<Eigen::VectorXd> t(x, dimFreeT);
        Eigen::Map<Eigen::MatrixXd> innerP(x + dimFreeT, 3, (coarseN - 1));

        setInitial(coarseT, innerP); // initialize the waypoint innerP and time vector  coarseT.

        backwardT(coarseT, t, softT, c2dfm); // get virtual T
        // backwardP(innerP, idxVs, cfgVs, p);  //cannot find why use this function to get p

        double minObjectivePenalty;
        lbfgs_params.mem_size = 8; // default
        lbfgs_params.past = 3;
        // lbfgs_params.delta = 1.0e-4;
        lbfgs_params.g_epsilon = 1.0e-5; // default
        lbfgs_params.min_step = 1.0e-20; // default
        lbfgs_params.max_iterations = 100;
        lbfgs_params.abs_curv_cond = 0;
        lbfgs_params.delta = relCostTol;
        size_t re_opt_num = 0;
        bool vel_fail = false;
        // cout << "mk2" <<endl;
        // cout << dynobs_pointer << endl;
        cout << "gcopter dynamic obs number:" << dynobs_pointer->dyn_number << "///" << dynobs_pointer->ball_number << endl;
        while ((re_opt_num == 0 | jerkOpt.pena_sta > 0.4 | vel_fail) && re_opt_num < 5) // if optimized traj is not safe, re-optimize it with hot start
        {
            lbfgs::lbfgs_optimize(dimFreeT + dimFreeP,
                                  x,
                                  &minObjectivePenalty,
                                  &POLY_OPT::objectiveFunc,
                                  nullptr,
                                  nullptr,
                                  this,
                                  &lbfgs_params);
            re_opt_num++;
            if (jerkOpt.pena_sta > 0.4)
                chi(0) *= 1.5;                        // increase the weight for collision cost
            forwardT(t, coarseT, softT, sumT, c2dfm); // update final coarseT
            fineT = coarseT;
            jerkOpt.generate(innerP, fineT); // generate final trajectory
            traj = jerkOpt.getTraj();
            if (traj.getMaxVelRate() > vMax + 0.2)
            {
                chi(1) *= 2.0;
                vel_fail = true;
            }
            else
            {
                vel_fail = false;
            }
            cout << "Optimize times: " << re_opt_num << " obs penalty: " << jerkOpt.pena_sta << endl;
            cout << "Max traj speed: " << traj.getMaxVelRate() << " vMax: " << vMax << endl;
            
        }
        cout << "initial state: \n"<<iState<<"\n final state: \n"<<fState<<endl;
        // for (auto j = 0; j < dimFreeT; j++)
        // {
        //     cout<<"t and x"<<j<<" : "<<t(j)<<" "<<x[j]<<endl;
        // }
        // cout << "mk3" <<endl;

        // splitToFineT(coarseT, intervals, fineT);
        // forwardP(p, idxVs, cfgVs, innerP); // update the innerP

        cout << "Total cost (include time): " << jerkOpt.cost_total << "\npena: " << jerkOpt.pena << "\npena_staticObs: " << jerkOpt.pena_sta << "\npena_dynobjects: " << jerkOpt.pena_ball << "\npena_dyn_accdistrib: " << jerkOpt.pena_acc << endl;
        // std::cout<<"-------------------------------------\n";
        // std::cout<<"total time is "<<jerkOpt.compute_time<<" s";

        if (jerkOpt.pena_sta + jerkOpt.pena_ball > 0.5)
        {
            cout << "Optimized trajectory is not safe?" << endl;
        }
        delete[] x;

        return jerkOpt.getTrajJerkCost();
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

#endif
