#include <planner_fast.h>

void TrajectoryGenerator_fast::read_param(const ros::NodeHandle *nh_priv)
{
  config.loadParameters(*nh_priv);
}

void TrajectoryGenerator_fast::replan_traj_noPolyH(mlmap::Ptr map, Vector3d &start,
                                                   Vector3d &vi, Vector3d &ai,
                                                   MatrixXd &waypoints,
                                                   vector<Vector3d> &start_end_divs,
                                                   dynobs_tmp *dynobs,
                                                   double plan_t,
                                                   bool full_trip,
                                                   bool if_reach)
{
  dynobs_pointer = dynobs;

  if (full_trip)
  {
    wPs.clear();
    if (waypoints.cols() < 3)
    {
      wPs.emplace_back(start);
      wPs.emplace_back(((start + waypoints.col(1)) / 2));
      wPs.emplace_back(waypoints.col(1));
    }
    else
    {
      wPs.emplace_back(start);
      for (uint i = 1; i < waypoints.cols(); i++)
      {
        if (((waypoints.col(i) - wPs.back()).norm() > 0.8) | i == waypoints.cols() - 1)
          wPs.emplace_back(waypoints.col(i));
      }
    }

    //  if_config = false;
  }
  Matrix3d iniState, finState;

  iniState.col(0) = start;
  iniState.col(1) = vi;
  iniState.col(2) = ai;
  finState.col(0) = wPs.back();
  finState.col(1) = Vector3d::Zero(3);
  finState.col(2) = Vector3d::Zero(3);

  Traj_opt_noPolyH(map, iniState, finState, plan_t);
  if (config.yawplan && dynobs_pointer->dyn_number > 0)
  {
    
    try
    {
      Yaw_plan(yaw_plan_tm);
      yaw_plan_tm = ros::Time::now().toSec();
    }
    catch (const std::exception &e)
    {
      std::cerr << "caught error: " << e.what() << '\n';
      yaw_plan.clear();
    }

    yaw_timeout = false;
  }
}


void TrajectoryGenerator_fast::Yaw_plan(double plan_t)
{
  // total_t = traj.getTotalDuration();
  yaw_plan.clear();
  if ((plan_t - plan_tm) > (total_t - config.delta_t_yaw - 1e-3))
  return;
  Matrix<double, 3, 5> camera_vertex;
  Vector3d sp_pos, sp_acc, sp_vel;
  Matrix3d Rota;
  Vector3d ct_center;
  double thrust, t_base, score;
  Vector2d v1, v2;
  double plan_duration = 2.0;
  int rows = (min(total_t, plan_t - plan_tm + plan_duration) - (plan_t - plan_tm)) / config.delta_t_yaw + 1;
  int cols = config.max_yaw_range / config.yaw_reso + 1;
  double v_psi;
  double M_yaw[rows][cols], M_vis_score[rows][cols], M_score[rows][cols];

  // vector<vector<Matrix<double, 3, 5>>> M_camera_vertex;
  // vector<Matrix<double, 3, 5>> fov_plan;
  int M_parent[rows][cols];
  double sp_yaw1, sp_theta, sp_phi;
  
  // Matrix<double, rows, cols> M_yaw, M_yaw1, M_vis_score, M_score;
  // Matrix<int, rows, cols> M_parent;
  v1 << 1.0, 0.0;
  int row = 0, col = 0;
  cout << "yaw plan begin" << endl;

  for (double ti = plan_t - plan_tm; ti < min(total_t, plan_t - plan_tm + plan_duration); ti += config.delta_t_yaw)
  {

    t_base = plan_tm - dynobs_pointer->time_stamp + ti;
    // yaw_plan_t.push_back(ti);
    sp_pos = traj.getPos(ti);
    sp_acc = traj.getAcc(ti);
    sp_vel = traj.getPos(total_t) - traj.getPos(ti); // traj.getVel(ti+0.1); //
    sp_acc(2) += G;
    thrust = sp_acc.norm();
    v2 = sp_vel.head(2);
    v_psi = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
    if (v2(1) < 0)
    {
      v_psi = -v_psi;
    }
    if (dynobs_pointer->dyn_number == 0)
    {
      yaw_plan.push_back(v_psi);
      // if (config.if_debug) std::cout<<"row,col:"<<row<<" "<<col<<" "<<rows<<" "<<cols<<std::endl;
      if (row == rows - 1)
      {
        // std::cout<<"No dynamic obs, use velocity yaw "<<std::endl;
        return;
      }
      row++;
      continue;
    }
    for (double sp_yaw = v_psi - config.max_yaw_range / 2; sp_yaw <= v_psi + config.max_yaw_range / 2; sp_yaw += config.yaw_reso)
    {
      if (row == 0)
      {
        M_vis_score[row][col] = 1;
        M_yaw[row][col] = drone_state.Euler(2); //
        M_score[row][col] = 3;
        continue;
      }

      if (sp_yaw > M_PI)
        sp_yaw1 = sp_yaw - 2 * M_PI;
      else if (sp_yaw < -M_PI)
        sp_yaw1 = sp_yaw + 2 * M_PI;
      else
        sp_yaw1 = sp_yaw;
      sp_theta = atan((sp_acc(0) + sp_acc(1) * tan(sp_yaw1)) / (sp_acc(2) * (cos(sp_yaw1) + sin(sp_yaw1) * tan(sp_yaw1))));
      sp_phi = acos(sp_acc(2) / thrust / cos(sp_theta));
      // AngleAxisd rollAngle(AngleAxisd(sp_phi, Vector3d::UnitX()));
      // AngleAxisd pitchAngle(AngleAxisd(sp_theta, Vector3d::UnitY()));
      // AngleAxisd yawAngle(AngleAxisd(sp_yaw1, Vector3d::UnitZ()));
      // Quaterniond quaternion = yawAngle * pitchAngle * rollAngle;
      // Rota = Quaternion2Rota(quaternion.normalized());
      Rota = Quaternion2Rota(Euler2Quaternion(Vector3d(sp_phi, sp_theta, sp_yaw1)));
      camera_vertex = (Rota * camera_vertex_b).array().colwise() + sp_pos.array();
      double score_vis = 0;
      // bool traj_vel_vis = inFOV(camera_vertex, (traj.getVel(ti).normalized() * 0.5 + sp_pos));
      //   score_traj_vel_vis = 0.0;
      // else
      //   score_traj_vel_vis = -5.0;
      // std::cout<<"Got the sample Rota: "<<Rota<<"row,col:"<<row<<" "<<col<<" "<<rows<<" "<<cols<<std::endl;
      for (int i = 0; i < dynobs_pointer->dyn_number; i++)
      {
        ct_center = dynobs_pointer->centers[i] + t_base * dynobs_pointer->vels[i];
        if (inFOV(camera_vertex, ct_center))
        {
          score_vis += 3 * max(0.0, double(4 - (ct_center - sp_pos).norm()));
          // cout << "dynobs in fov, score_vis:  " <<score_vis<< " sp_yaw: "<<sp_yaw<<endl;
        }
      }
      // score_vis += score_traj_vel_vis;
      M_vis_score[row][col] = score_vis;
      M_yaw[row][col] = sp_yaw1; //
      // M_camera_vertex[row][col] = camera_vertex;
      col += 1;
    }
    col = 0;
    row += 1;
  }
  // cout << "yaw plan begin" << endl;
  // std::cout<<"Got the vis score matrix: \n"<<M_vis_score<<std::endl;
  double max_total_score = 0;
  int choosed_col = 0;
  for (row = 1; row < rows; row++)
  {
    for (col = 0; col < cols; col++)
    {
      int parent = 0;
      double tmp_score;
      score = 0;

      for (int col0 = 0; col0 < cols; col0++)
      {

        double yaw_gap_raw = abs(M_yaw[row][col] - M_yaw[row - 1][col0]);
        double yaw_gap = (yaw_gap_raw > M_PI) ? (2 * M_PI - yaw_gap_raw) : yaw_gap_raw;
        // if (config.yaw_gap_max - yaw_gap < 0 && row > 1)
        //   continue;

        if ((M_vis_score[row][col]) == 0)
        {
          yaw_gap_raw = abs(M_yaw[row][col] - v_psi);
          yaw_gap = (yaw_gap_raw > M_PI) ? (2 * M_PI - yaw_gap_raw) : yaw_gap_raw;
        }
        // cout << "yaw gap:" << yaw_gap << " "<<yaw_gap_raw<< endl;

        tmp_score = M_vis_score[row][col] + M_score[row - 1][col0] + (config.yaw_gap_max - yaw_gap > 0) ? 0 : (config.yaw_gap_max - yaw_gap) * config.yaw_w;
        if (tmp_score > score)
        {
          parent = col0;
          score = tmp_score;
        }
        if (row == 1)
          break;
      }
      // cout << "parent:" <<parent <<" cols: "<<cols<<endl;
      M_parent[row][col] = parent;
      M_score[row][col] = score;

      if (row == rows - 1 && score > max_total_score)
      {
        max_total_score = score;
        choosed_col = col;
      }
    }
  }
  cout << "yaw plan finish--0" << endl;
  cout << "choosed_col:" <<choosed_col;
  for (row = rows - 1; row > 0; row--)
  {
    yaw_plan.push_back(M_yaw[row][choosed_col]);
    choosed_col = M_parent[row][choosed_col];
    // if (config.if_debug)
      std::cout << "yaw planned: " << yaw_plan.back() << " row:" << row << " choosed_col: " << choosed_col << " rows: " << rows << std::endl;
  }
  yaw_plan.push_back(M_yaw[row][choosed_col]);
  std::reverse(std::begin(yaw_plan), std::end(yaw_plan));
  cout << "yaw plan finish" << endl;
}

inline bool TrajectoryGenerator_fast::inFOV(Matrix<double, 3, 5> camera_vertex, Vector3d ct_center)
{
  return !((ct_center.head(3) - camera_vertex.col(0)).dot((camera_vertex.col(0) - camera_vertex.col(1)).cross(camera_vertex.col(0) - camera_vertex.col(2))) < 0 || (ct_center.head(3) - camera_vertex.col(0)).dot((camera_vertex.col(0) - camera_vertex.col(2)).cross(camera_vertex.col(0) - camera_vertex.col(3))) < 0 || (ct_center.head(3) - camera_vertex.col(0)).dot((camera_vertex.col(0) - camera_vertex.col(3)).cross(camera_vertex.col(0) - camera_vertex.col(4))) < 0 || (ct_center.head(3) - camera_vertex.col(0)).dot((camera_vertex.col(0) - camera_vertex.col(4)).cross(camera_vertex.col(0) - camera_vertex.col(1))) < 0 || (ct_center.head(3) - camera_vertex.col(1)).dot((camera_vertex.col(1) - camera_vertex.col(4)).cross(camera_vertex.col(1) - camera_vertex.col(2))) < 0);
}

inline double TrajectoryGenerator_fast::cal_yaw_with_2dVec(Vector2d p)
{
  Vector2d v1(1.0, 0.0);
  double psi = acos(v1.dot(p) / (p.norm()));
  if (p(1) < 0)
  {
    psi = -psi;
  }
  return psi;
}

inline bool TrajectoryGenerator_fast::HaveToCheckUnkownPoint(double traj_t_start, Vector2d &first_unkown_point_traj, mlmap::Ptr map)
{
  double speed = traj.getVel(traj_t_start).norm();
  Vector3d pos0 = traj.getPos(traj_t_start);
  for (double dt = 0.2; dt < speed / config.horizontalAccMax + 0.5; dt += 0.1)
  {
    Vector3d pos = traj.getPos(traj_t_start + dt);
    if (map->getOccupancy(pos) == mlmap::UNKNOWN)
    {
      first_unkown_point_traj = (pos - pos0).head(2);
      return true;
    }
  }
  return false;
}
Vector2d TrajectoryGenerator_fast::getYaw(double t, mlmap::Ptr map)
{
  uint indx = (t - yaw_plan_tm) / config.delta_t_yaw;
  Vector2d desire_yaw; // yaw,yaw_rate
  Vector2d first_unkown_point_traj;
  if (HaveToCheckUnkownPoint(t - plan_tm, first_unkown_point_traj, map))
  {
    // cout << "Turn camera to the unkown area1" << endl;
    desire_yaw(0) = cal_yaw_with_2dVec(first_unkown_point_traj);
    desire_yaw(1) = 0.0;
    return desire_yaw;
  }
  double yaw_gap;

  // std::cout<<"v2: "<<v2<<std::endl;
  double v_psi = cal_yaw_with_2dVec((traj.getPos(total_t) - traj.getPos(min((t - plan_tm), total_t - 0.2))).head(2));

  // std::cout<<"Set yaw begin "<< indx<<" "<<yaw_plan.size()<<"\n"
  // <<((yaw_plan.size()>0)?yaw_plan.back():0)<<" "<<dynobs_pointer->dyn_number<<std::endl;
  if ((indx + 2) > yaw_plan.size() || !config.yawplan || dynobs_pointer->dyn_number == 0)
  {
    // yaw_timeout = true;
    desire_yaw(0) = v_psi;
    desire_yaw(1) = 0;
  }
  else
  {
    try
    {

      yaw_gap = (yaw_plan[indx + 1] - yaw_plan[indx]);
      if (abs(yaw_gap) > M_PI)
      {
        yaw_gap = copysign((2 * M_PI - abs(yaw_gap)), drone_state.Euler(2));
      }
      desire_yaw(0) = yaw_plan[indx] + ((t - yaw_plan_tm) - indx * config.delta_t_yaw) / config.delta_t_yaw * yaw_gap;
      desire_yaw(1) = 0;
      if (abs(desire_yaw(0)) > M_PI)
        cout << "yaw out range! " << desire_yaw(0) << endl;
      // std::cout << "Set yaw-2:\n"
      // << desire_yaw(0) << " " << yaw_plan[indx + 1] << " " << yaw_plan[indx] << std::endl;
    }
    catch (...)
    { // exception should be caught by reference
      cout << "get yaw exception "
           << "\n";
      yaw_timeout = true;
      desire_yaw(0) = v_psi;
      desire_yaw(1) = 0;
    }
  }
  return desire_yaw;
}

Vector2d TrajectoryGenerator_fast::getYaw(double t)
{

  uint indx = (t - yaw_plan_tm) / config.delta_t_yaw;
  Vector2d desire_yaw; // yaw,yaw_rate
  double yaw_gap;
  Vector2d v1, v2;
  v1 << 1.0, 0.0;
  //  v2 = traj.getVel(min(t+1,total_t)).head(2);

  v2 = (traj.getPos(total_t) - traj.getPos(min((t - plan_tm), total_t - 0.2))).head(2);
  // std::cout<<"v2: "<<v2<<std::endl;
  double v_psi = acos(v1.dot(v2) / (v2.norm()));
  if (v2(1) < 0)
  {
    v_psi = -v_psi;
  }
  // std::cout<<"Set yaw begin "<< indx<<" "<<yaw_plan.size()<<"\n"<<((yaw_plan.size()>0)?yaw_plan.back():0)<<" "<<dynobs_pointer->dyn_number<< " "<<(yaw_plan.size() - 1)<<std::endl;
  if ((indx + 2) > yaw_plan.size() || !config.yawplan || dynobs_pointer->dyn_number == 0)
  {
    // yaw_timeout = true;
    desire_yaw(0) = v_psi;
    desire_yaw(1) = 0;
  }
  else
  {
    try
    {

      yaw_gap = (yaw_plan[indx + 1] - yaw_plan[indx]);
      if (abs(yaw_gap) > M_PI)
      {
        yaw_gap = copysign((2 * M_PI - abs(yaw_gap)), desire_yaw(0));
      }
      desire_yaw(0) = yaw_plan[indx] + ((t - yaw_plan_tm) - indx * config.delta_t_yaw) / config.delta_t_yaw * yaw_gap;
      desire_yaw(1) = 0;
      if (abs(desire_yaw(0)) > M_PI)
        cout << "yaw out range! " << desire_yaw(0) << endl;
      // std::cout << "Set yaw-2:\n"
      // << desire_yaw(0) << " " << yaw_plan[indx + 1] << " " << yaw_plan[indx] << std::endl;
    }
    catch (...)
    { // exception should be caught by reference
      cout << "get yaw exception "
           << "\n";
      yaw_timeout = true;
      desire_yaw(0) = v_psi;
      desire_yaw(1) = 0;
    }
  }
  return desire_yaw;
}


bool TrajectoryGenerator_fast::inGlobalBound(const Vector3d &pos)
{
  return (pos(0) > config.global_bound(0, 0) &&
          pos(0) < config.global_bound(0, 1) &&
          pos(1) > config.global_bound(1, 0) &&
          pos(1) < config.global_bound(1, 1) &&
          pos(2) > config.global_bound(2, 0) &&
          pos(2) < config.global_bound(2, 1));
}
bool TrajectoryGenerator_fast::check_traj_safe(const double plan_t, Vector3d &start, mlmap::Ptr map, dynobs_tmp *dynobs, double start_t, bool &if_time_out, bool pcl_update)
{

  dynobs_pointer = dynobs;
  // plan_tm = plan_t;
  double start_t1 = start_t - plan_t;
  //  check_sfc_ind = 0;
  Vector3d pos;
  double dt = 0.1;
  pt[0] = start(0);
  pt[1] = start(1);
  pt[2] = start(2);

  if (start_t1 > total_t * 0.5)
  {
    if_time_out = true;
    return false;
  }
  {
    if_time_out = false;
  }
  
  if (start_t - yaw_plan_tm > 0.2)
  yaw_timeout = true;

  if (config.yawplan && dynobs_pointer->dyn_number > 0 && yaw_timeout)
  {
    
    Yaw_plan(yaw_plan_tm);
    yaw_plan_tm = ros::Time::now().toSec();
    cout<<" Yaw plan time (ms): "<<(yaw_plan_tm - start_t) * 1000 <<endl;
    yaw_timeout = false;
  }
  if ((traj.getPos(total_t) - wPs.back()).norm() > config.horizon * 0.3)
  {
    cout << "Path replanned, so traj replan. Traj end and goal distance:" << (traj.getPos(total_t) - wPs.back()).norm() << endl;
    return false;
  }
  if (start_t1 + dt > total_t)
  {
    cout << "Safety check timestamp is close to traj end!" << endl;
    return true;
  }

  //  cout<< "mk2"<<endl;
  last_dyn_check_fail = false;
  for (double j = start_t1 + dt; j < total_t * 0.8; j += dt)
  {
    // cout<<"(in func) check traj safe:\n"<<j<<endl<<total_t<<endl;
    pos = traj.getPos(j);
    if (!dyn_safe_check(pos, j + start_t - start_t1))
    {
      last_dyn_check_fail = true;
      cout << "dyn check fail!  " << dynobs_pointer->dyn_number << endl; // dynobs_pointer->ballvel[0]<<endl
      return false;
    }

    if (map->getOccupancy(pos, config.safeMargin - 0.1) == mlmap::OCCUPIED || !inGlobalBound(pos))
    {
      cout << "map check to be collide or out side the global box!  " << pos.transpose() << endl; // dynobs_pointer->ballvel[0]<<endl
      return false;
    }
  }

  // cout << "old traj is safe for dyn obs, and all in new corridor! " << dynobs_pointer->dyn_number << " " << ros::Time::now().toSec() - dynobs_pointer->time_stamp << endl;
  return true;
}


inline bool TrajectoryGenerator_fast::dyn_safe_check(Vector3d pt, double check_t)
// { return true; // for rosbag tests!!!
{
  double t_base;
  Vector3d ct_center;
  Vector3d conv_vec;
  double obj_prop_conv;
  if (dynobs_pointer->dyn_number > 0)
  {
    t_base = check_t - dynobs_pointer->time_stamp;

    for (int j = 0; j < dynobs_pointer->dyn_number; j++)
    {
      ct_center = dynobs_pointer->centers[j] + t_base * dynobs_pointer->vels[j];
      obj_prop_conv = pow(dynobs_pointer->max_accs[j](1) + dynobs_pointer->max_accs[j](2) * t_base * t_base, 0.5);
      obj_prop_conv = obj_prop_conv > 0.15 ? 0.15 : obj_prop_conv;
      conv_vec = {obj_prop_conv, obj_prop_conv, 0.0};
      // cout << "ct_center:\n"<<ct_center <<endl<<"pt: \n"<<pt<<endl<<"obs size: \n"<<dynobs_pointer->obs_sizes[j]*0.5<<"\ntime gap: "<<t_base<<"\n pos gap:"<<(ct_center - pt).cwiseAbs()<<endl<<(((ct_center - pt).cwiseAbs() - dynobs_pointer->obs_sizes[j]*0.5).array()<0)<<endl;
      if ((((ct_center - pt).cwiseAbs() - dynobs_pointer->obs_sizes[j] * 0.5 - conv_vec).array() < config.safeMargin).all())
      {
        // cout<<"1111"<<endl;
        return false;
      }
    }
  }
  if (dynobs_pointer->ball_number > 0)
  {
    t_base = check_t - dynobs_pointer->ball_time_stamp;
    for (int j = 0; j < dynobs_pointer->ball_number; j++)
    {
      ct_center = dynobs_pointer->ballpos[j] + t_base * dynobs_pointer->ballvel[j] + 0.5 * t_base * t_base * dynobs_pointer->ballacc[j];
      if ((((ct_center - pt).cwiseAbs() - dynobs_pointer->ball_sizes[j] * 0.5).array() < config.safeMargin).all())
      {
        return false;
      }
    }
  }
  return true;
}

void TrajectoryGenerator_fast::Traj_opt_noPolyH(mlmap::Ptr map, const MatrixXd &iniState, const MatrixXd &finState, double plan_t)
{
  chrono::high_resolution_clock::time_point tic = chrono::high_resolution_clock::now();
  // Trajectory traj;
  //  cout << "received dynamic obs number:" << dynobs_pointer->dyn_number << endl << dynobs_pointer<<endl<<"dyn_timestamp:"<<dynobs_pointer->time_stamp<<endl<<plan_t<<endl;
  ROS_INFO("Begin to optimize the traj~");
  
  // TODO
  if (!MINCOOpt.setup(config.rho, config.totalT, iniState, finState, map, wPs, config.qdIntervals,
                      config.safeMargin, (config.velMax < 2.0 && last_dyn_check_fail) ? 2.0 : config.velMax, config.thrustAccMin, config.thrustAccMax,
                      config.bodyRateMax, config.gravAcc, config.penaltyPVTB, config.useC2Diffeo, plan_t, dynobs_pointer, config.global_bound))
  {
    ROS_INFO("gcopter initialize fail!");
    return;
  }
  double finalObj = MINCOOpt.optimize(traj, config.optRelTol);

  chrono::high_resolution_clock::time_point toc = chrono::high_resolution_clock::now();
  double compTime = chrono::duration_cast<chrono::microseconds>(toc - tic).count() * 1.0e-3;

  // printf("finished!!!\n");
  cout << "Optimization time usage: " << compTime << " ms" << endl;
  cout << "Final jerk cost: " << finalObj << endl;
  Max_traj_vel = traj.getMaxVelRate();
  cout << "Maximum Vel: " << Max_traj_vel << endl;
  cout << "Maximum Acc: " << traj.getMaxAccRate() << endl;
  cout << "Total traj Duration: " << traj.getTotalDuration() << endl;
  total_t = traj.getTotalDuration();
  plan_tm = ros::Time::now().toSec();
}

void TrajectoryGenerator_fast::Traj_opt(const MatrixXd &iniState, const MatrixXd &finState, double plan_t)

{
  chrono::high_resolution_clock::time_point tic = chrono::high_resolution_clock::now();
  // Trajectory traj;
  //  cout << "received dynamic obs number:" << dynobs_pointer->dyn_number << endl << dynobs_pointer<<endl<<"dyn_timestamp:"<<dynobs_pointer->time_stamp<<endl<<plan_t<<endl;
  ROS_INFO("Begin to optimize the traj~");

  if (!GCOpt.setup(config.rho, config.totalT, iniState, finState, hPolys, INFINITY,
                   config.qdIntervals, config.horizHalfLen, config.vertHalfLen,
                   config.safeMargin, (dynobs_pointer->dyn_number > 0 && config.velMax < 2.0) ? 2.0 : config.velMax, config.thrustAccMin, config.thrustAccMax,
                   config.bodyRateMax, config.gravAcc, config.penaltyPVTB, config.useC2Diffeo, plan_t, dynobs_pointer))
  {
    ROS_INFO("gcopter initialize fail!");
    return;
  }
  double finalObj = GCOpt.optimize(traj, config.optRelTol);

  chrono::high_resolution_clock::time_point toc = chrono::high_resolution_clock::now();
  double compTime = chrono::duration_cast<chrono::microseconds>(toc - tic).count() * 1.0e-3;

  printf("finished!!!\n");
  cout << "Optimization time usage: " << compTime << " ms" << endl;
  cout << "Final jerk cost: " << finalObj << endl;
  cout << "Maximum Vel: " << traj.getMaxVelRate() << endl;
  cout << "Maximum Acc: " << traj.getMaxAccRate() << endl;
  cout << "Total traj Duration: " << traj.getTotalDuration() << endl;
  total_t = traj.getTotalDuration();
  plan_tm = ros::Time::now().toSec();
}

void TrajectoryGenerator_fast::sort_vec(const VectorXd &vec, VectorXd &sorted_vec, VectorXi &ind)
{
  ind = VectorXi::LinSpaced(vec.size(), 0, vec.size() - 1); //[0 1 2 3 ... N-1]
  auto rule = [vec](int i, int j) -> bool
  {
    return vec(i) > vec(j);
  };
  sort(ind.data(), ind.data() + ind.size(), rule);

  sorted_vec.resize(vec.size());
  for (uint i = 0; i < vec.size(); i++)
  {
    sorted_vec(i) = vec(ind(i));
  }
}

// allocate time for waypoints
VectorXd TrajectoryGenerator_fast::allocateTime(const MatrixXd &wayPs,
                                                double vel,
                                                double acc)
{
  int N = (int)(wayPs.cols()) - 1;
  VectorXd durations(N);
  if (N > 0)
  {
    Eigen::Vector3d p0, p1;
    double dtxyz, D, acct, accd, dcct, dccd, t1, t2, t3;

    for (int k = 0; k < N; k++)
    {
      p0 = wayPs.col(k);
      p1 = wayPs.col(k + 1);
      D = (p1 - p0).norm(); // distance

      acct = vel / acc;               // accelerate time
      accd = (acc * acct * acct / 2); // accelerate distance
      dcct = vel / acc;               // de-accelerate time
      dccd = acc * dcct * dcct / 2;   // de-accelerate distance

      if (D < accd + dccd)
      {
        t1 = sqrt(acc * D) / acc;
        t2 = (acc * t1) / acc;
        dtxyz = t1 + t2;
      }
      else
      {
        t1 = acct;
        t2 = (D - accd - dccd) / vel;
        t3 = dcct;
        dtxyz = t1 + t2 + t3;
      }

      durations(k) = dtxyz;
    }
  }

  return durations;
}


void TrajectoryGenerator_fast::get_desire(double timee, Vector3d &p_d, Vector3d &v_d, Vector3d &a_d, Vector3d &p_d_yaw)
{
  p_d = traj.getPos(timee);
  v_d = traj.getVel(timee);
  a_d = traj.getAcc(timee);
  p_d_yaw = traj.getPos(clip(timee + 5.0, 0.0, total_t - 0.01));
  // p_d_yaw = traj.getPos(total_t-0.1);
  // cout<<"pd,vd,ad,p_d_yaw: \n"<<p_d<<"\n"<<v_d<<"\n"<<a_d<<"\n"<<p_d_yaw<<endl;
}

void TrajectoryGenerator_fast::get_traj_samples(MatrixXd &sp_pos, MatrixXd &sp_vel, MatrixXd &sp_acc, double start_t)
{
  // total_t = traj.getTotalDuration();
  // cout << "pub traj:" << total_t <<endl;
  start_t = min(start_t, total_t);
  double delta_t = 0.3;
  int num = (total_t - start_t) / delta_t;
  num = max(num, 1);
  sp_pos.resize(num + 1, 3);
  sp_vel.resize(num + 1, 3);
  sp_acc.resize(num + 1, 3);
  for (int i = 0; i < num; i++)
  {
    sp_pos.row(i) = traj.getPos(start_t + i * delta_t);
    sp_vel.row(i) = traj.getVel(min(start_t + i * delta_t + 0.1, total_t));
    sp_acc.row(i) = traj.getAcc(start_t + i * delta_t);
  }
  sp_pos.row(num) = traj.getPos(total_t);
  sp_vel.row(num) = traj.getVel(total_t);
  sp_acc.row(num) = traj.getAcc(total_t);
}
