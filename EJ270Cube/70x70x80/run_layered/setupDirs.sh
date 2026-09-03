#!/bin/bash
# NeutronSpectrum_layered の root ファイル1本につき1サブディレクトリ (run.sh) を生成する。
# 既存の run/{ppm}ppm/run.sh を層構造データ向けに一般化したもの。

DIR=$(cd $(dirname $0) && pwd)
SRC_DIR="/home/yomata/work/Moon/EJ270Cube/70x70x80/NeutronSpectrum_layered"

mkdir -p "$DIR/analysis/macro" "$DIR/analysis/fig"
: > "$DIR/analysis/folders.list"

for f in "$SRC_DIR"/Neutron_*_angle.root; do
    fname=$(basename "$f")
    name=${fname#Neutron_}
    name=${name%_angle.root}
    name=${name/_380MV_/_}

    mkdir -p "$DIR/$name"

    cat > "$DIR/$name/run.sh" <<EOS
#!/bin/bash

export MOONNEUTRON_SPECTRUM_FILE="$f"

../../bin/exampleB2a \\
    ../../macro/beamOn.mac
EOS
    chmod +x "$DIR/$name/run.sh"

    echo "$name" >> "$DIR/analysis/folders.list"
done

sort -o "$DIR/analysis/folders.list" "$DIR/analysis/folders.list"

echo "completed: $(wc -l < "$DIR/analysis/folders.list") directories created" >&2
