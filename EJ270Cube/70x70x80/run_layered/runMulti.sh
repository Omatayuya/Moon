#!/bin/bash
# run/runMulti.sh を一般化したもの。
# ppm リストを直書きする代わりに、run.sh を持つサブディレクトリを自動で列挙して順に実行する
# (層構造データは 126 通りあり手書きリストは非現実的なため)。

DIR=$(cd $(dirname $0) && pwd)

for d in "$DIR"/*/; do
    name=$(basename "$d")
    [ -f "$d/run.sh" ] || continue

    cd "$d"

    echo `date` >&2
    echo `pwd` >&2

    bash ./run.sh
done

cd "$DIR"
echo `date` >&2
echo "completed" >&2
