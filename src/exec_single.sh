#!/bin/bash
# ------------------------------------------------------------
# Copyright (c) 2024, National Institute for Materials Science
#
# This software is released under the MIT License.
#
# Contributor:
#     Hiroshi Shinotsuka
# ------------------------------------------------------------

export LD_LIBRARY_PATH=/usr/local/lib:${LD_LIBRARY_PATH}

# if [ -n "${NSLOTS}" ]; then
#     # 1. NSLOTSがあれば最優先
#     export OMP_NUM_THREADS=${NSLOTS}
# elif [ -n "${OMP_NUM_THREADS}" ]; then
#     # 2. すでにOMP_NUM_THREADSが設定されていればそれを維持（何もしない）
#     :
# elif [ -n "$(nproc)" ]; then
#     # 3. なければ利用可能なCPU数を取得
#     export OMP_NUM_THREADS=$(nproc)
# else
#     # 4. それ以外（nprocが使えない環境など）はデフォルト4
#     export OMP_NUM_THREADS=4
# fi

# 強制的に1に設定
export OMP_NUM_THREADS=1

# SCRIPT_DIR="${PWD}"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "Start:" `TZ=JST-9 date "+%Y/%m/%d %H:%M:%S"`


if [ $# != 1 ]; then
    echo "Usage: $0 <csv_file>"
    exit 1
elif [ ! -s "$1" ]; then
    echo "Error. No data in $1"
    exit 1
fi

chmod u+x "${SCRIPT_DIR}/activeshirley_gcc.exe"

original_filename="$1"
original_dirname=`basename "${original_filename}" ".csv"`

temp_filename="_sample.csv"
cp "${original_filename}" "${temp_filename}"

work_dir=`basename "${temp_filename}" ".csv"`
rm -rf "${original_dirname}"
mkdir -p "${work_dir}"
cp "${temp_filename}" "${work_dir}"

# Use "recipe.txt" if it exists in $PWD; otherwise, not use it.
PWD_RECIPE="${SCRIPT_DIR}/recipe.txt"
if [ -f "${PWD_RECIPE}" ]; then
    cp "${PWD_RECIPE}" "${work_dir}"
fi

cd "${work_dir}"
    python3 ${SCRIPT_DIR}/automatic_xps_peak_separation_single.py -p "${temp_filename}" >stdout 2>stderr
    mv "${temp_filename}" "${original_filename}"
cd ..

mv "${work_dir}" "${original_dirname}"
rm -f "${temp_filename}"

# rm -rf __pycache__

echo "End:  " `TZ=JST-9 date "+%Y/%m/%d %H:%M:%S"`
