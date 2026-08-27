#!/bin/bash
# run/delete.sh を一般化したもの。
# ppm リストを直書きする代わりに、results.root を持つサブディレクトリを自動で列挙して削除する。

DIR=$(cd $(dirname $0) && pwd)

read -p "本当に各ディレクトリの results.root を削除しますか？ [y/N] " confirm
if [[ "$confirm" != "y" && "$confirm" != "Y" ]]; then
    echo "中止しました" >&2
    exit 0
fi

for d in "$DIR"/*/; do
    target="$d/results.root"

    if [ -f "$target" ]; then
        echo "削除: $target" >&2
        rm -f "$target"
    else
        echo "無し: $target" >&2
    fi
done

echo "completed" >&2
