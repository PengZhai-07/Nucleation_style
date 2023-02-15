#!/bin/bash
# from wozhi to current working director
scp -r pengzhai@wozhi.earth.lsa.umich.edu:/home/pengzhai/Spear-normal-stress/data/wholespace/$1 /Users/pengzhai/Desktop/Spear-normal-stress/data/wholespace/phase_diagram_L_b/ 


# from Geatlakes to current working director
scp -r pengzhai@greatlakes.arc-ts.umich.edu:/home/pengzhai/Spear-normal-stress/data/wholespace/phase_diagram_L_b/0_500_16_0.8_0.0_4_0.75_0.021429_0.0* /Users/pengzhai/Desktop/Spear-normal-stress/data/wholespace/phase_diagram_L_b/ 

# from turbo to current working director

scp -r pengzhai@greatlakes.arc-ts.umich.edu:/nfs/turbo/lsa-yiheh/yiheh-mistorage/pengz/data/wholespace/phase_diagram_L_b/0.75_32_100_0_0_1.0_0.0_4_0.6_0.04 /Users/pengzhai/Desktop/Spear-normal-stress/data/wholespace/phase_diagram_L_b/
