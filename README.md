# uape_2.0_ws


## FIRST->Requirements: 
you will need to have the four packages in the workspace as shown in the figure below:

<img src="https://user-images.githubusercontent.com/59171742/221513025-44c8d0d6-4d1a-4ef2-9bd2-56c1b8d73cca.png"  width="500">

### Download other 3 packages,

#### 1.1 [MLMapping](https://github.com/chenhanpolyu/MLMapping) 

#### 1.2 [yolo_fdsst_piv_module](https://github.com/chenhanpolyu/AutoFly-demo/tree/master/src/yolo_fdsst_piv)

#### 1.3 [px4control](https://github.com/ZJU-FAST-Lab/Fast-Drone-250/tree/master/src/realflight_modules/px4ctrl)

and put them into the src/ folder in your workspace.

## SECOND->Installation

### step 1: install mlmapping dependencies:

$ cd uape_2.0_ws/src/MLMapping/3rdPartLib/

allen:~/uape_2.0_ws/src/MLMapping/3rdPartLib$ ./install3rdPartLib.sh

## step2: in the workspace direcotry:
~/uape_2.0_ws$ catkin_make

## THIRD->Launch
allen@allen:~/uape_2.0_ws$ source devel/setup.bash 

allen@allen:~/uape_2.0_ws$ roslaunch uape_planner_2 icuas_sim_traj_node.launch 

it runs with the docker-Gazebo simulation env provided by ICUAS-2023 competition.

Ros node is shown below:

<img src="https://user-images.githubusercontent.com/59171742/221510752-85ff7544-8dd2-4ff8-89fd-fc39fecf558c.png" width="800">

And the running scene is shown in below

![icuas_palnner](https://user-images.githubusercontent.com/59171742/221510875-22b40571-8f39-46e9-916a-b2bbe3230f32.png)
