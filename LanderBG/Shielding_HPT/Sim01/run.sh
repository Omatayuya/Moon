#!/bin/sh

DIR=$(cd $(dirname $0) && pwd)
PHI="W_0"
#PHI="W_118.5"

# ARRAY1=("01-H/$PHI" "02-He/$PHI")
# ARRAY1=("01-H/$PHI")
ARRAY1=("02-He/$PHI")
#ARRAY1=("Radioactive/K-40" "Radioactive/Th-232" "Radioactive/U-238")

echo "start" `date` >&2

for i in "${!ARRAY1[@]}"
do
	cd $DIR
	cd ./${ARRAY1[$i]}
	echo `pwd` `date` >&2
	../../../../bin/runMain
done

cd $DIR

echo "complete" `date` >&2

