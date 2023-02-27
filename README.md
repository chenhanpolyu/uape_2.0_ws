# uape_2.0_ws
you will need to have the four packages in the workspace as shown in the figure below:

![uape_2 0_ws](https://user-images.githubusercontent.com/59171742/221505002-1ec91827-7f4e-4658-a9af-7fadae330e67.png)

install steps:
## step 1: download other 3 packages
### 1.1 [MLMapping](https://github.com/chenhanpolyu/MLMapping) 
### 1.2 [yolo_module](https://github.com/chenhanpolyu/AutoFly-demo/tree/master/src/yolo_fdsst_piv)
### 1.3 [px4control](https://github.com/ZJU-FAST-Lab/Fast-Drone-250/tree/master/src/realflight_modules/px4ctrl)

## step 2: mlmapping dependencies:
$ cd uape_2.0_ws/src/MLMapping/3rdPartLib/
allen:~/uape_2.0_ws/src/MLMapping/3rdPartLib$ ./install3rdPartLib.sh

## step3: in the workspace direcotry:
~/uape_2.0_ws$ catkin_make
