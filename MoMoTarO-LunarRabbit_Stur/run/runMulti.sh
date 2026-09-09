#!/bin/bash

DIR=$(cd $(dirname $0) && pwd)

PPM_LIST=(0 10 20 50 100 200 500 1000 2000 5000 10000)

for ppm in "${PPM_LIST[@]}"; do
    # mkdir -p "$DIR/runs/${ppm}ppm"
    cd "$DIR/${ppm}ppm"

    echo `date` >&2
    echo `pwd` >&2

    bash ./run.sh
done

cd "$DIR"
echo `date` >&2
echo "completed" >&2
