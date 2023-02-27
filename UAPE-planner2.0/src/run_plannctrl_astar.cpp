#include <stdio.h>
#include <iostream>
#include <Eigen/Eigen>
#include <math.h>
#include <Tools.h>
#include <ros/ros.h>
#include <cstdlib>
#include <ros_missn.h>
#include <planner_fast.h>
#include <Tools.h>
#include <path_searching/include/path_searching/kinodynamic_astar.h>
#include <mlmap.h>
#include <time.h>

using namespace std;
using namespace Eigen;


double sign(double x)
{
  return (x > 0) - (x < 0);
}
unique_ptr<KinodynamicAstar> kino_path_finder_;
int main(int argc, char **argv)
{
  // Initialize ROS
  ros::init(argc, argv, "traj_planner");
  // Create a handle to this process node.
  ros::NodeHandle nh("~");

  // init the mapping kit
  mlmap::Ptr mapPtr;
  mapPtr.reset(new mlmap);
  mapPtr->init_map(nh);
  Vector3d ct_pos, ct_vel, ct_acc, home_pos;
  ct_pos.setZero();
  Vector3d p_d, v_d, a_d, p_d_yaw;
  a_d.setZero();
  MatrixXd sp_pos, sp_vel, sp_acc, waypoints_m;
  // double last_path_t = 0;

  Eigen::Vector3d end_state = Eigen::Vector3d::Zero(3);
  vector<double> goalp, gbbox_o, gbbox_l, timecosts;
  Matrix<double, 3, 5> camera_vertex, camera_vertex_b, camera_vertex_bv;
  vector<Eigen::Vector3d> waypoints, start_end_derivatives;
  int FE_status=0;
  double d2r = 3.14159265 / 180;
  double cam_depth = 10.0;
  double cam_depth_v = 2.0;
  double h_fov = 87; // in degree
  double v_fov = 58;
  double dis_goal, tem_dis_goal;
  // double sfck_t;
  bool ifMove, if_rand;
  int CtrlFreq;
  bool if_debug;
  bool path_replan;
  bool if_reach = false;
  bool last_if_reach = false;
  double gap;
  double singlestep_time;
  bool rc_goal = false;
  int return_home = 0;
  bool if_raw_pcl = false; // if directly use the raw point cloud from sensor
  bool if_depth_img = false;
  bool if_time_out = false;
  nh.getParam("goal", goalp);
  nh.getParam("search/horizon", dis_goal);
  nh.getParam("ifMove", ifMove);
  nh.getParam("cam_depth", cam_depth);
  nh.getParam("h_fov", h_fov);
  nh.getParam("v_fov", v_fov);
  nh.getParam("if_RandomGoal", if_rand);
  nh.getParam("GlobalBox_min", gbbox_o);
  nh.getParam("GlobalBox_size", gbbox_l);
  nh.getParam("CtrlFreq", CtrlFreq);
  nh.getParam("if_debug", if_debug);
  nh.getParam("UseRawPcl", if_raw_pcl);
  nh.getParam("UseRawDepth", if_depth_img);

  nh.getParam("ReturnHome", return_home);
  nh.getParam("UseRcGuide", rc_goal);

  ros::Rate loop_rate(CtrlFreq);
  camera_vertex_b.col(0) << 0, 0, 0;
  camera_vertex_b.col(1) << cam_depth, tan(h_fov / 2 * d2r) * cam_depth, tan(v_fov / 2 * d2r) * cam_depth;
  camera_vertex_b.col(2) << cam_depth, -tan(h_fov / 2 * d2r) * cam_depth, tan(v_fov / 2 * d2r) * cam_depth;
  camera_vertex_b.col(3) << cam_depth, -tan(h_fov / 2 * d2r) * cam_depth, -tan(v_fov / 2 * d2r) * cam_depth;
  camera_vertex_b.col(4) << cam_depth, tan(h_fov / 2 * d2r) * cam_depth, -tan(v_fov / 2 * d2r) * cam_depth;

  camera_vertex_bv.col(0) << 0, 0, 0;
  camera_vertex_bv.col(1) << cam_depth_v, tan(h_fov / 2 * d2r) * cam_depth_v, tan(v_fov / 2 * d2r) * cam_depth_v;
  camera_vertex_bv.col(2) << cam_depth_v, -tan(h_fov / 2 * d2r) * cam_depth_v, tan(v_fov / 2 * d2r) * cam_depth_v;
  camera_vertex_bv.col(3) << cam_depth_v, -tan(h_fov / 2 * d2r) * cam_depth_v, -tan(v_fov / 2 * d2r) * cam_depth_v;
  camera_vertex_bv.col(4) << cam_depth_v, tan(h_fov / 2 * d2r) * cam_depth_v, -tan(v_fov / 2 * d2r) * cam_depth_v;
  TrajectoryGenerator_fast planner(camera_vertex_b);
  // dis_goal = dis_goal_ini-0.5;
  Eigen::Vector3d g_goal = {goalp[0], goalp[1], goalp[2]};
  Eigen::Vector3d goal = g_goal;
  Eigen::Vector3d local_goal = {0, 0, 0};
  Eigen::Vector3d initial_goal = g_goal;
  bool if_initial = true;
  bool if_end = false;
  bool if_safe = true;
  bool ball_time_out = false;
  double ball_pass_time = 2.0, dyn_pass_time = 0.5;
  double min_dist2dynobs = 1e6;
  double tmp_dist, t_gap_ball;
  double timee = 0;
  ros::Time traj_last_t = ros::Time::now();
  chrono::high_resolution_clock::time_point last_traj_tic = chrono::high_resolution_clock::now();
  int rand_num = -1, rand_num_tmp;
  planner.read_param(&nh);
  States state;
  cout << "Traj node initialized!" << endl;
  RosClass RosMsg(&nh, CtrlFreq, if_raw_pcl, if_depth_img);
  do
  {
    state = RosMsg.get_state();
    ct_pos = state.P_E;
    ros::Duration(0.1).sleep();
  } while ((ct_pos.norm() < 1e-3) | isnan(ct_pos(0)) | isnan(ct_pos(1)) | isnan(ct_pos(2)));
  cout << "UAV message received!" << ct_pos << endl;
  planner.last_check_pos = ct_pos;
  home_pos = ct_pos;
  // RosMsg.dynobs_pointer->ball_number = 0;
  while (!mapPtr->has_data)
  {
    ros::spinOnce();
    ros::Duration(0.1).sleep();
  }
  cout << "Point cloud received!\n"
       << endl;

  // ------- start the safety check and re-planning loop--------------------------------------------------------------
  kino_path_finder_.reset(new KinodynamicAstar);
  kino_path_finder_->setParam(nh);
  kino_path_finder_->init();
  srand((unsigned)time(NULL));
  while (nh.ok()) // main loop
  {
    // dis_goal = dis_goal_ini-0.5;

    path_replan = false;
    if ((return_home > 0 && if_end) || (if_rand && (if_end || (rand_num < 0)))) // choose goal randomly at the global bounding box boundary
    {
      if (if_rand)
      {
        if (rand_num < 0)
          rand_num = rand() % 4;
        else
        {
          do
            rand_num_tmp = rand() % 4;
          while (rand_num_tmp == rand_num);
          rand_num = rand_num_tmp;
        }
        if (rand_num == 0)
        {
          g_goal(0) = gbbox_o[0] + gbbox_l[0] - 1.0;
          g_goal(1) = gbbox_o[1] + rand() % int(gbbox_l[1]) + 1.0;
        }
        else if (rand_num == 1)
        {
          g_goal(0) = gbbox_o[0] + rand() % int(gbbox_l[0]) + 1.0;
          g_goal(1) = gbbox_o[1] + gbbox_l[1] - 1.0;
        }
        else if (rand_num == 2)
        {
          g_goal(0) = gbbox_o[0] + 1.0;
          g_goal(1) = gbbox_o[1] + rand() % int(gbbox_l[1]) + 1.0;
        }
        else
        {
          g_goal(0) = gbbox_o[0] + rand() % int(gbbox_l[0]) + 1.0;
          g_goal(1) = gbbox_o[1] + 1.0;
        }
        g_goal(2) = 1.5;
        goal = g_goal;
      }
      else if (return_home > 1 && if_end)
      {
        if ((home_pos - ct_pos).norm() > 1)
          g_goal = home_pos;
        else
          g_goal = initial_goal;
        goal = g_goal;
        cout << "one goal: " << g_goal << endl;
        return_home--;
      }

      if_end = false;
      if_initial = true;
      if_safe = true;
      if_reach = false;
      FE_status = 0;
      waypoints.clear();
      RosMsg.dynobs_pointer->dyn_number = 0;
      Vector2d v2 = (g_goal - state.P_E).head(2);
      Vector2d v1;
      v1 << 1.0, 0.0;
      double desire_yaw = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
      if (v2(1) < 0)
      {
        desire_yaw = -desire_yaw;
      }
      ct_pos = state.P_E;
      double yaw_rate = abs(desire_yaw - state.Euler(2)) > 3.14159 ? -sign(desire_yaw - state.Euler(2)) * 0.6 : sign(desire_yaw - state.Euler(2)) * 0.6;
      ros::Duration(0.5).sleep();
      do
      {
        state = RosMsg.step(state.Euler(2) + yaw_rate * 0.5, yaw_rate, ct_pos, Vector3d::Zero(3), Vector3d::Zero(3), "pos_vel_acc_yaw_icuas");
        ros::Duration(0.05).sleep();
      } while (abs(state.Euler(2) - desire_yaw) > 0.3);
      ros::Duration(0.5).sleep();
    }
    else if (rc_goal)
    {
      // state = RosMsg.get_state();
      local_goal = {RosMsg.rc_data.ch[1], RosMsg.rc_data.ch[0], clip(RosMsg.rc_data.ch[3], -0.3, 0.3)};
      ct_pos = state.P_E;
      double yaw_fix = state.Euler(2);
      while ((local_goal * dis_goal).norm() < 0.5)
      {
        state = RosMsg.step(yaw_fix, 0, ct_pos, Vector3d::Zero(3), Vector3d::Zero(3), "pos_vel_acc_yaw_icuas");
        local_goal = {RosMsg.rc_data.ch[1], RosMsg.rc_data.ch[0], 0};
        ros::Duration(0.1).sleep();
        if_end = false;
        if_initial = true;
        if_reach = false;
        waypoints.clear();
        RosMsg.dynobs_pointer->dyn_number = 0;
        cout << "RC hover" << endl;
      }

      while ((local_goal * dis_goal).norm() < 2.0)
      {
        state = RosMsg.step(yaw_fix, 0, state.P_E + local_goal, Vector3d::Zero(3), Vector3d::Zero(3), "pos_vel_acc_yaw_icuas");
        local_goal = {RosMsg.rc_data.ch[1], RosMsg.rc_data.ch[0], clip(RosMsg.rc_data.ch[3], -0.3, 0.3)};
        ros::Duration(0.1).sleep();
        if_end = false;
        if_initial = true;
        if_reach = false;
        waypoints.clear();
        RosMsg.dynobs_pointer->dyn_number = 0;
        cout << "RC move" << endl;
      }

      g_goal = state.P_E + local_goal * dis_goal;
      goal = g_goal;
      Vector2d v2 = (g_goal - state.P_E).head(2);
      Vector2d v1;
      v1 << 1.0, 0.0;
      double desire_yaw = acos(v1.dot(v2) / (v1.norm() * v2.norm()));
      if (v2(1) < 0)
      {
        desire_yaw = -desire_yaw;
      }
      ct_pos = state.P_E;
      double yaw_rate = abs(desire_yaw - state.Euler(2)) > 3.14159 ? -sign(desire_yaw - state.Euler(2)) * 0.6 : sign(desire_yaw - state.Euler(2)) * 0.6;
      while (abs(state.Euler(2) - desire_yaw) > 0.3 && if_initial)

      {
        state = RosMsg.step(state.Euler(2) + yaw_rate * 0.5, yaw_rate, ct_pos, Vector3d::Zero(3), Vector3d::Zero(3), "pos_vel_acc_yaw_icuas");
        ros::Duration(0.05).sleep();
      }
      planner.config.velMax = clip(RosMsg.rc_data.ch[2] + 1, 0.3, 2.0);
    }

    // ros::Time t1 = ros::Time::now();
    chrono::high_resolution_clock::time_point t1 = chrono::high_resolution_clock::now();
    if (mapPtr->map_updated || RosMsg.dyn_update)
    {
      mapPtr->map_updated = false;
      if (RosMsg.dynobs_pointer->dyn_number > 0)
      {
        for (auto i = 0; i < RosMsg.dynobs_pointer->dyn_number; i++)
        {
          mapPtr->setFree_map_in_bound(RosMsg.dynobs_pointer->centers[i] - RosMsg.dynobs_pointer->obs_sizes[i]*0.6, 
          RosMsg.dynobs_pointer->centers[i] + RosMsg.dynobs_pointer->obs_sizes[i]*0.6);
        }
      }
      state = RosMsg.get_state();
      planner.drone_state = state;
      // RosMsg.set_cod_update(false);

      if (RosMsg.dynobs_pointer->ball_number > 0 && ros::Time::now().toSec() - RosMsg.dynobs_pointer->ball_time_stamp > ball_pass_time) // for bag sim// RosMsg.dynobs_pointer->ball_number>0 && (RosMsg.dynobs_pointer->ballvel[0](0) > -0.2)||
      {
        if (if_debug)
          cout << "dyn ball time out: " << RosMsg.dynobs_pointer->ballvel[0] << endl;
        RosMsg.dynobs_pointer->ball_number = 0;
        ball_time_out = true;
      }
      if (RosMsg.dynobs_pointer->dyn_number > 0 && ros::Time::now().toSec() - RosMsg.dynobs_pointer->time_stamp > dyn_pass_time) // for bag sim// RosMsg.dynobs_pointer->ball_number>0 && (RosMsg.dynobs_pointer->ballvel[0](0) > -0.2)||
      {
        RosMsg.dynobs_pointer->dyn_number = 0;
        cout << "dyn time out! " << endl;
      }
      min_dist2dynobs = 1e6;
      for (int bi = 0; bi < RosMsg.dynobs_pointer->dyn_number; bi++)

      {
        t_gap_ball = ros::Time::now().toSec() - RosMsg.dynobs_pointer->time_stamp;
        tmp_dist = (state.P_E - RosMsg.dynobs_pointer->centers[bi] - t_gap_ball * RosMsg.dynobs_pointer->vels[bi]).norm();
        if (tmp_dist < min_dist2dynobs)
        {
          min_dist2dynobs = tmp_dist;
          if (if_debug)
            cout << "min distance from objects to drone:" << min_dist2dynobs << endl
                 << state.P_E << endl
                 << RosMsg.dynobs_pointer->centers[bi] + t_gap_ball * RosMsg.dynobs_pointer->vels[bi];
        }
      }
      min_dist2dynobs = 1e6;
      for (int di = 0; di < RosMsg.dynobs_pointer->ball_number && RosMsg.dynobs_pointer->ballvel[0](0) < -0.4; di++)
      {
        t_gap_ball = ros::Time::now().toSec() - RosMsg.dynobs_pointer->ball_time_stamp;
        tmp_dist = (state.P_E - (RosMsg.dynobs_pointer->ballpos[di] + t_gap_ball * RosMsg.dynobs_pointer->ballvel[di] + 0.5 * t_gap_ball * t_gap_ball * RosMsg.dynobs_pointer->ballacc[di])).norm();
        if (tmp_dist < min_dist2dynobs)
        {
          min_dist2dynobs = tmp_dist;
          if (if_debug)
            cout << "min distance from ball to drone:" << min_dist2dynobs << endl;
        }
      }
      if (if_initial || !ifMove)
      {
        ct_pos = state.P_E;
        ct_vel.setZero();
        ct_acc.setZero();
      }
      else
      {
        ct_pos = p_d;
        ct_vel = v_d;
        ct_acc = a_d;
      }
      last_if_reach = if_reach;
      if (if_debug)
        cout << "track goal dist: " << (state.P_E - p_d).norm() << endl;
      // cout<<"if initial: "<<if_initial<<endl;
      camera_vertex = (state.Rota * camera_vertex_b).array().colwise() + state.P_E.array();
      double dis2goal = (g_goal - ct_pos).norm();

      if (!if_initial)
      {
        if (ifMove)
          if_safe = planner.check_traj_safe(traj_last_t.toSec(), ct_pos, mapPtr, RosMsg.dynobs_pointer, ros::Time::now().toSec(), if_time_out, RosMsg.pcl_update);
        else
          if_safe = planner.check_traj_safe(ros::Time::now().toSec(), ct_pos, mapPtr, RosMsg.dynobs_pointer, ros::Time::now().toSec(), if_time_out, RosMsg.pcl_update);
      }
      // cout << "if time out: " << if_time_out << endl;

      bool if_old_kinopath_collide = !kino_path_finder_->checkOldPath(waypoints, ct_pos, false);
      if (!if_time_out && !(!if_old_kinopath_collide == if_safe))
      {
        cout << "frond/back end safety check result is different: " << !if_old_kinopath_collide << "  " << if_safe << endl;
      }
      if ((if_time_out | (!if_safe | if_initial)) && FE_status != KinodynamicAstar::GOAL_CLOSE)
      {
        cout << "Ready to search" << endl;
        chrono::high_resolution_clock::time_point tic = chrono::high_resolution_clock::now();
        kino_path_finder_->setEnvironment(mapPtr, RosMsg.dynobs_pointer, camera_vertex, gbbox_o, gbbox_l);
        // kino_path_finder_->setEnvironment(RosMsg.obs_pointer, RosMsg.dynobs_pointer, camera_vertex, gbbox_o, gbbox_l);
        if (dis2goal > dis_goal)
        {
          goal = ct_pos + (g_goal - ct_pos) / dis2goal * dis_goal;
          if_reach = false;
          tem_dis_goal = dis_goal;
        }
        else
        {
          goal = g_goal;
          if_reach = true;
          tem_dis_goal = dis2goal;
        }
        kino_path_finder_->reset();
        FE_status = kino_path_finder_->search(ct_pos, ct_vel, ct_acc, goal, end_state, true);
        // last_path_t = ros::Time::now().toSec();
        while (FE_status == KinodynamicAstar::GOAL_OCC)
        {
          cout << "[kino replan]: Goal occluded:\n"
               << goal << " tem_dis_goal: " << tem_dis_goal << endl;
          tem_dis_goal -= 0.3;
          goal = ct_pos + (g_goal - ct_pos).normalized() * tem_dis_goal;
          if (if_reach)
            g_goal = goal;
          if ((goal - state.P_E).norm() < 0.5 || tem_dis_goal < 0.5)
          {
            g_goal = goal;
            if_reach = false;
            break;
          }
          kino_path_finder_->reset();
          FE_status = kino_path_finder_->search(ct_pos, ct_vel, ct_acc, goal, end_state, true);
          // if_reach = false;
        }
        if (FE_status == KinodynamicAstar::NO_PATH)
        {
          cout << "[kino replan]: kinodynamic search fail!" << endl;

          // retry searching with discontinuous initial state
          kino_path_finder_->reset();
          FE_status = kino_path_finder_->search(ct_pos, ct_vel, ct_acc, goal, end_state, false);
          // last_path_t = ros::Time::now().toSec();
          if (FE_status == KinodynamicAstar::NO_PATH)
          {
            if_reach = false;
            cout << "[kino replan]: Can't find path. Please restart" << endl;
            // return 0;
          }
          else
          {
            cout << "[kino replan]: retry search success." << endl;
          }
        }
        else if (FE_status == KinodynamicAstar::GOAL_OCC)
        {
          return 0;
        }

        double compTime = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - tic).count() * 1.0e-3;
        cout << "kino path planning finished! time cost (ms)： " << compTime << endl;
        if (FE_status != KinodynamicAstar::NO_PATH && FE_status != KinodynamicAstar::GOAL_OCC && FE_status != KinodynamicAstar::GOAL_CLOSE)
        {
          kino_path_finder_->getSamples(0.2, waypoints, start_end_derivatives);
        }
        else
        {
          waypoints.clear();
          waypoints.emplace_back(ct_pos);
          waypoints.emplace_back(goal);
        }

        path_replan = true;
      }
      else
      {
        // cout << "[kino replan]: Old kino path is safe" << endl;
      }

      waypoints_m = Map<MatrixXd>(waypoints[0].data(), 3, waypoints.size());

      if (if_debug)
        cout << "point cloud update, check safety result:" << if_safe << endl;
      if ((if_safe && !if_initial && !((!last_if_reach) && if_reach) && !ball_time_out) | 
           (FE_status == KinodynamicAstar::GOAL_CLOSE && !if_initial)) //&& (timee+sfck_t < planner.total_t || planner.total_t <sfck_t)
      {
        gap = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - last_traj_tic).count() * 1.0e-6;
        timee = gap; //+ 1 / CtrlFreq / 4;
      }
      else
      {
        cout << "Dynamic obj number: " << RosMsg.dynobs_pointer->dyn_number << endl;
        planner.replan_traj_noPolyH(mapPtr, ct_pos, ct_vel, ct_acc, waypoints_m,
                                      start_end_derivatives, RosMsg.dynobs_pointer,
                                      ros::Time::now().toSec(),
                                      true, if_reach);

        if (ball_time_out)
          ball_time_out = false;

        traj_last_t = ros::Time::now();
        gap = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - last_traj_tic).count() * 1.0e-6 - timee;
        last_traj_tic = chrono::high_resolution_clock::now();
        if (if_initial)
          timee = 0; // 1 / CtrlFreq / 4;
        else
          timee = gap;
        // cout<<"if reach: "<<last_if_reach<<"  "<<if_reach<<endl;
        if (if_initial)
          if_initial = false;
        if (if_debug)
          cout << "old traj is not safe, get new traj!" << endl;
        // if (planner.Max_traj_vel > max(2.0, planner.config.velMax) + 0.5)
        // {
        //   cout<<"traj maxvel is too large! not safe!"<<endl;
        //   return 0;
        // }

      }
      singlestep_time = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - t1).count() * 1.0e-3;
      // cout << "Total time cost (ms): " << singlestep_time << endl;
      timecosts.push_back(singlestep_time);
      if (timecosts.size() > 1000)
      {
        cout << "Average Total time cost (ms): " << std::accumulate(std::begin(timecosts), std::end(timecosts), 0.0) / timecosts.size() << endl;
        timecosts.clear();
      }
      //  <<ros::Time::now().toSec()<<endl<<t1<<endl<<traj_last_t<<endl<<RosMsg.dynobs_pointer->time_stamp<<endl;
      if (ifMove)
        planner.get_traj_samples(sp_pos, sp_vel, sp_acc, (ros::Time::now() - traj_last_t).toSec());
      else
        planner.get_traj_samples(sp_pos, sp_vel, sp_acc, 0.0);
      RosMsg.pub_traj(sp_pos, sp_vel, sp_acc, planner.fail_pt);
      // RosMsg.pub_fovlist (sp_pos, sp_vel, sp_acc,camera_vertex_bv, planner.yaw_plan);
      RosMsg.pub_path(waypoints);
    }

    else // if pcl has not been updated
    {

      gap = chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - last_traj_tic).count() * 1.0e-6;
      timee = gap; // + 1 / CtrlFreq / 4;
      // cout<<"control sample time gap (no pcl updated): "<<gap<<endl;
      state = RosMsg.get_state();
      planner.drone_state = state;

    }

    Vector2d desire_psi;
    if (timee < planner.total_t - 0.01)
    {
      planner.get_desire(timee, p_d, v_d, a_d, p_d_yaw);

      desire_psi = planner.getYaw(timee + traj_last_t.toSec(), mapPtr);
    }
    else
    {
      planner.get_desire(planner.total_t, p_d, v_d, a_d, p_d_yaw);
      cout << "goal reached: " << g_goal << endl;
      if_end = true;

    }

    if (ifMove)
    {
      state = RosMsg.step(desire_psi[0], desire_psi[1], p_d, v_d, a_d, "pos_vel_acc_yaw_icuas");

    }

    RosMsg.pub_fovshape(camera_vertex);
    if (RosMsg.dynobs_pointer->ball_number > 0)
    {
      RosMsg.pub_ballstates();
    }

    loop_rate.sleep();
    if (if_end && !if_rand && return_home <= 1 && !rc_goal)
    {
      break;
    }

  }

  return 0;
}
