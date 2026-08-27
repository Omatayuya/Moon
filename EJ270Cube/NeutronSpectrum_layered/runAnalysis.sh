#!/bin/sh

#SHAPE="box_10x10x0.05m"
#ARRAY1=("H_3000ppm" "H_15000ppm" "H_30000ppm")
#ARRAY2=("depth_0.025m" "depth_0.075m" "depth_0.125m" "depth_0.225m" "depth_0.325m" "depth_0.525m" "depth_0.725m" "depth_1.025m" "depth_1.225m")

#SHAPE="box_10x10x0.1m"
#ARRAY1=("H_1500ppm" "H_7500ppm" "H_15000ppm")
#ARRAY2=("depth_0.05m" "depth_0.1m" "depth_0.15m" "depth_0.25m" "depth_0.35m" "depth_0.55m" "depth_0.75m" "depth_1.05m" "depth_1.25m")

#SHAPE="box_10x10x0.2m"
#ARRAY1=("H_750ppm" "H_3750ppm" "H_7500ppm")
#ARRAY2=("depth_0.1m" "depth_0.15m" "depth_0.2m" "depth_0.3m" "depth_0.4m" "depth_0.6m" "depth_0.8m" "depth_1.1m" "depth_1.3m")

#SHAPE="box_10x10x0.5m"
#ARRAY1=("H_300ppm" "H_1500ppm" "H_3000ppm")
#ARRAY2=("depth_0.25m" "depth_0.3m" "depth_0.35m" "depth_0.45m" "depth_0.55m" "depth_0.75m" "depth_0.95m" "depth_1.25m")

SHAPE="box_10x10x1m"
ARRAY1=("H_150ppm" "H_750ppm" "H_1500ppm")
ARRAY2=("depth_0.5m" "depth_0.55m" "depth_0.6m" "depth_0.7m" "depth_0.8m" "depth_1m")

#SHAPE="box_10x10x1.5m"
#ARRAY1=("H_100ppm" "H_500ppm" "H_1000ppm")
#ARRAY2=("depth_0.75m")


for i in "${!ARRAY1[@]}"
do
	for j in "${!ARRAY2[@]}"
	do
		root -l -q 'toRootNeutron_angle.C+("'${SHAPE}'","'${ARRAY1[$i]}'","'${ARRAY2[$j]}'")'
	done
done

