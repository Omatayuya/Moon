#!/bin/sh

ARRAY1=("H_0ppm" "H_10ppm" "H_20ppm" "H_50ppm" "H_100ppm" "H_200ppm" "H_500ppm" "H_1000ppm" "H_2000ppm" "H_5000ppm" "H_10000ppm")

for i in "${!ARRAY1[@]}"
do
	root -l -q 'toRootNeutron_angle.C+("'${ARRAY1[$i]}'")'
done

