#!/bin/bash

DIR=$(cd $(dirname $0) && pwd)

PPM_LIST=(0 10 20 50 100 200 500 1000 2000 5000 10000)

read -p "本当に各ディレクトリの results.root を削除しますか？ [y/N] " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "中止しました" >&2
    exit 0
fi

for ppm in "${PPM_LIST[@]}"; do
    target="$DIR/${ppm}ppm/results.root"

    if [ -f "$target" ]; then
        echo "削除: $target" >&2
        rm -f "$target"
    else
        echo "無し: $target" >&2
    fi
done

echo "completed" >&2
