#pragma once

// #include "am_traj.hpp"
// #include <planners/wayp_reader.h>
#include <traj_opt/traj_utils.h> //decomp_ros is also included in it
#include <traj_opt/se3gcopter_cpu.hpp>
#include <traj_opt/MINCO_traj_optimizer.hpp>
#include <call_states/ros_communicate.h>
#include <iostream>
#include <stdio.h>
#include <math.h>
#include <Eigen/Eigen>
#include <stdlib.h>
#include <Tools.h>
#include "mlmap.h"
using namespace std;
using namespace Eigen;
// using namespace min_jerk;
struct Config
{
    string odomFrame;

    // Params
    double scaleSI;
    // double mapHeight;
    Vector3d polyhedronBox, global_min, global_size;
    MatrixXd global_bound;
    double rho;
    double totalT;
    int qdIntervals;
    double horizHalfLen;
    double vertHalfLen;
    double safeMargin;
    double velMax;
    double thrustAccMin;
    double thrustAccMax;
    double horizontalAccMax;
    double bodyRateMax;
    double gravAcc;
    Vector4d penaltyPVTB;
    bool useC2Diffeo;
    bool if_debug;
    double optRelTol;
    double yaw_replan_t;
    double trajVizWidth;
    Vector3d trajVizRGB;
    string routeStoragePath;
    string ellipsoidPath;
    Vector4d ellipsoidVizRGBA;
    int max_yaw_range_n = 10;
    double yaw_w = 0.4;
    double yaw_gap_max = 0.55;
    double horizon = 8.0;
    double yaw_reso = 0.175;
    double delta_t_yaw = 0.5;
    double yaw_active_distance = 4.0;
    bool applyUnknownAwareYaw = true;
    bool yawplan = false;
    inline void loadParameters(const ros::NodeHandle &nh_priv)
    // {       cout << "mk3" << endl;
    {
        vector<double> vecPolyhedronBox, vecPenaltyPVTB, vecTrajVizRGB, vecEllipsoidVizRGBA, vecGlobal_min, vecGlobal_size;
        ;

        nh_priv.getParam("Yawplan", yawplan);
        nh_priv.getParam("MaxYawRange", max_yaw_range_n);
        nh_priv.getParam("YawGapWeight", yaw_w);
        nh_priv.getParam("YawGapMax", yaw_gap_max);
        nh_priv.getParam("YawResolution", yaw_reso);
        nh_priv.getParam("ApplyUnknownAwareYaw", applyUnknownAwareYaw);
        nh_priv.getParam("YawDeltaT", delta_t_yaw);
        nh_priv.getParam("YawReplanTimeInterval", yaw_replan_t);
        nh_priv.getParam("YawActiveDistance", yaw_active_distance);
        
        nh_priv.getParam("ScaleSI", scaleSI);
        nh_priv.getParam("PolyhedronBox", vecPolyhedronBox);
        // cout << "mk4" << endl<< scaleSI << endl << vecPolyhedronBox[2] << endl << vecPolyhedronBox[0] <<endl;
        polyhedronBox << vecPolyhedronBox[0], vecPolyhedronBox[1], vecPolyhedronBox[2];
        nh_priv.getParam("GlobalBox_min", vecGlobal_min);
        nh_priv.getParam("GlobalBox_size", vecGlobal_size);
        global_min << vecGlobal_min[0], vecGlobal_min[1], vecGlobal_min[2];
        global_size << vecGlobal_size[0], vecGlobal_size[1], vecGlobal_size[2];

        nh_priv.getParam("Rho", rho);
        nh_priv.getParam("TotalT", totalT);

        nh_priv.getParam("QdIntervals", qdIntervals);
        nh_priv.getParam("HorizHalfLen", horizHalfLen);

        nh_priv.getParam("VertHalfLen", vertHalfLen);
        nh_priv.getParam("SafeMargin", safeMargin);
        nh_priv.getParam("VelMax", velMax);
        
        nh_priv.getParam("ThrustAccMin", thrustAccMin);
        nh_priv.getParam("ThrustAccMax", thrustAccMax);
        nh_priv.getParam("BodyRateMax", bodyRateMax);
        nh_priv.getParam("GravAcc", gravAcc);
        nh_priv.getParam("PenaltyPVTB", vecPenaltyPVTB);
        penaltyPVTB << vecPenaltyPVTB[0], vecPenaltyPVTB[1], vecPenaltyPVTB[2], vecPenaltyPVTB[3];
        nh_priv.getParam("UseC2Diffeo", useC2Diffeo);
        nh_priv.getParam("OptRelTol", optRelTol);
        nh_priv.getParam("TrajVizWidth", trajVizWidth);
        nh_priv.getParam("TrajVizRGB", vecTrajVizRGB);
        nh_priv.getParam("if_debug", if_debug);
        nh_priv.getParam("search/horizon", horizon);
        trajVizRGB << vecTrajVizRGB[0], vecTrajVizRGB[1], vecTrajVizRGB[2];

        nh_priv.getParam("EllipsoidVizRGBA", vecEllipsoidVizRGBA);
        
        
        ellipsoidVizRGBA << vecEllipsoidVizRGBA[0], vecEllipsoidVizRGBA[1], vecEllipsoidVizRGBA[2], vecEllipsoidVizRGBA[3];

        global_bound.resize(3, 2);
        global_bound.col(0) = global_min.array() + safeMargin;
        global_bound.col(1) = (global_min + global_size).array() - safeMargin;
        horizontalAccMax = sqrt(thrustAccMax*thrustAccMax - gravAcc*gravAcc);

    }
};

/*
// Generate Trajectory by waypoints.txt
*/
class TrajectoryGenerator_fast
{
private:
    Matrix<double, 3, 5> camera_vertex_b;
    ros::NodeHandle nh_;
    double MaxVel, MaxVelCal;
    double adapt_Margin = 0.35;
    double Acc = 5.0;
    bool if_config = false; // if the waypoints and polyhedrons are updated
    Matrix3d iS, fS;        // xyz * pva        // intial & finial state
    bool last_traj_polyH_check;
    bool last_dyn_check_fail = false;
    bool yaw_timeout = true;
    double yaw_plan_tm = 0.0, plan_tm = 0.0;
    SE3GCOPTER GCOpt;
    POLY_OPT MINCOOpt;

    VectorXd allocateTime(const MatrixXd &wayPs, double vel, double acc);

    bool dyn_safe_check(Vector3d pt, double check_t);
    void sort_vec(const VectorXd &vec, VectorXd &sorted_vec, VectorXi &ind);
    void Traj_opt(const MatrixXd &iniState, const MatrixXd &finState, double plan_t);
    void Traj_opt_noPolyH(mlmap::Ptr map, const MatrixXd &iniState, const MatrixXd &finState, double plan_t);
    void Yaw_plan(double plan_t);
    bool inFOV(Matrix<double, 3, 5> camera_vertex, Vector3d ct_center);
    inline double cal_yaw_with_2dVec(Vector2d p);
    inline bool HaveToCheckUnkownPoint(double traj_t_start, Vector2d& first_unkown_point_traj, mlmap::Ptr map);
    
public:
    //         MatrixXd waypoints;    // pva * time
    TrajectoryGenerator_fast(Matrix<double, 3, 5> camera_vertex_bb)
    {
        camera_vertex_b = camera_vertex_bb;
    }
    States drone_state;
    Config config;
    vec_Vec3f *obs_pointer;
    dynobs_tmp *dynobs_pointer;
    vector<MatrixXd> hPolys;
    // vec_E<Polyhedron3D> decompPolys;
    Trajectory traj;
    VectorXd ts; // time for pieces
    VectorXf bound;
    vector<Vector3d> wPs;
    double total_t;
    double Max_traj_vel;
    Vec3f pt, fail_pt;
    uint check_sfc_ind;
    vector<double> yaw_plan;
    vector<double> yaw_plan_t;
    Vector3d last_check_pos;
    ros::Time last_check_time;
    // bool check_traj_safe(const MatrixXd &cd_c, const VectorXd &cd_r, const double start_t);
    bool inGlobalBound(const Vector3d &pos);
    void replan_traj_noPolyH(mlmap::Ptr map, Vector3d &start, Vector3d &vi, Vector3d &ai, MatrixXd &waypoints, vector<Vector3d> &start_end_divs, dynobs_tmp *dynobs, double plan_t, bool full_trip = true, bool if_reach = true);
    bool check_traj_safe(const double plan_t, Vector3d &start, mlmap::Ptr map, dynobs_tmp *dynobs, double start_t, bool &if_time_out, bool pcl_update);
    void read_param(const ros::NodeHandle *nh_priv);

    int N;
    void get_desire(double timee, Vector3d &p_d, Vector3d &v_d, Vector3d &a_d, Vector3d &p_d_yaw);
    void get_traj_samples(MatrixXd &sp_pos, MatrixXd &sp_vel, MatrixXd &sp_acc, double start_t);
    Vector2d getYaw(double t);
    Vector2d getYaw(double t,  mlmap::Ptr map);

};
