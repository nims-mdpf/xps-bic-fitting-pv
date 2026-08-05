# ------------------------------------------------------------
# Copyright (c) 2024, National Institute for Materials Science
#
# This software is released under the MIT License.
#
# Contributor:
#     Hiroshi Shinotsuka
# ------------------------------------------------------------
# # coding: utf-8

import os
import pathlib
import shutil
import subprocess
import sys
from glob import glob

import pandas as pd
from myfunc import get_value_from_settings

# Specify the directory which has "activeshirley_gcc.exe" and "plot_result.py"
software_dir = os.path.abspath(os.path.dirname(__file__))
activeshirley_exe_file = f"{software_dir}/activeshirley_gcc.exe"
plot_py_file = f"{software_dir}/plot_result_shirley_and_linear.py"
settingsfile = f"{software_dir}/settings.inp"
ppt_py_file = f"{software_dir}/make_result_figures_ppt2.py"
PYTHON = "python3"
toolname = "xps-ps-pv-20240510"

print(f"# toolname               = {toolname}")
print(f"# software_dir             = {software_dir}")
print(f"# activeshirley_exe_file   = {activeshirley_exe_file}")
print(f"# plot_py_file             = {plot_py_file}")
print(f"# ppt_py_file              = {ppt_py_file}")


def exec_command(cmd, targetdir="./", logfile=""):
    """
    Execute a command cmd assuming a Linux system.
    If specified, stdout and stderr will be output to logfile.
    Otherwise, the log is output to the standard console.

    >>> exec_command('echo hello')
    0
    """
    if logfile == "":
        p = subprocess.Popen(cmd, cwd=targetdir, shell=True)
    else:
        p = subprocess.Popen(
            cmd.split(), cwd=targetdir, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
        )
        with open(logfile, mode="w") as f:
            for line in p.stdout:
                line = line.decode("utf-8").rstrip()
                f.write(line + "\n")
    p.wait()
    return p.returncode


#
def get_num_model(kwd):
    # 最小のピーク半値半幅が0.1eVより大きいモデルの個数を返す
    # minFWHM_threshold = get_minFWHM_threshold(f"./{kwd}/settings.inp")
    minFWHM_threshold = get_value_from_settings(
        "minFWHM_threshold", 0.2, f"./{kwd}/settings.inp"
    )
    numModel = 0
    paramfileList = list(glob(f"./{kwd}/parameters_peak_search_mSG*.csv"))
    for paramfile in paramfileList:
        dfp = pd.read_csv(paramfile, usecols=["width"])
        min_HWHM = dfp.width[dfp.width.idxmin()]
        if min_HWHM > minFWHM_threshold / 2:
            numModel = numModel + 1
    return numModel


def main_core(fnm, kwd, devFlag=False):
    # pwd = os.getcwd()

    # kwdはshirleyまたはlinear
    # kwdディレクトリを作る
    targetdir = kwd
    os.makedirs(targetdir, exist_ok=True)

    # スペクトルデータファイルをコピーする
    shutil.copy2(fnm, targetdir)

    # デフォルトの settings.inp を targetdir にコピーする．
    # もしデフォルトの  settings.inp が存在しなければ， targetdir に空のファイルを作る．
    target_settings_file = f"{targetdir}/settings.inp"
    if os.path.exists(settingsfile):
        shutil.copy2(settingsfile, target_settings_file)
    else:
        pathlib.Path(target_settings_file).touch()

    # shirley の場合は settings.inp をそのまま使う．
    # linear の場合は settings.inp に backgroundType と allowedNegativeRatio の設定を追記する．
    if kwd == "shirley":
        os.path.exists(target_settings_file)
    elif kwd == "linear":
        with open(target_settings_file, "a") as f:
            print("backgroundType               Linear", file=f)
            print("allowedNegativeRatio         0.05", file=f)
    else:
        print("kwd must be either shirley or linear")
        sys.exit(1)

    # 自動ピーク分離を実行し，標準出力をlog_as_{kwd}.txtに書き出す
    cmd = f"{activeshirley_exe_file} {fnm}"
    print(cmd)
    exec_command(cmd, targetdir, f"log_as_{kwd}.txt")
    numModel = get_num_model(kwd)
    print(f"Number of valid models for {kwd} BG = {numModel}")
    return numModel


def main(fnm, bgType="auto", figType="-minfig", pptFlag=False, devFlag=False):
    if bgType == "auto":
        numModel = main_core(fnm, "shirley", devFlag)

        if numModel == 0:
            numModel = main_core(fnm, "linear", devFlag)

    elif bgType == "shirley":
        numModel = main_core(fnm, "shirley", devFlag)
    elif bgType == "linear":
        numModel = main_core(fnm, "linear", devFlag)

    # 事後処理を実行し，標準出力をlog2.txtに書き出す
    cmd = f"{PYTHON} {plot_py_file} {fnm} {figType}"
    print(cmd)
    exec_command(cmd, logfile="log2.txt")

    if numModel == 0:
        print("Warning. No result found.")
    elif pptFlag:
        # 結果をpptにまとめる
        cmd = f"{PYTHON} {ppt_py_file}"
        print(cmd)
        exec_command(cmd, logfile="log3.txt")

    # 無用の中間ファイルを削除
    if not devFlag:
        print(
            "remove linear/ shirley/ prefered_*.png py_*.png summary_evalFunc_minvalue.csv summary_evalFunc_vs_mSG.csv"
        )
        [shutil.rmtree(d) for d in ["linear", "shirley"] if os.path.exists(d)]
        [
            [os.remove(f) for f in glob(key)]
            for key in [
                "prefered_*.png",
                "py_*.png",
                "summary_evalFunc_minvalue.csv",
                "summary_evalFunc_vs_mSG.csv",
            ]
        ]

    return


def show_usage():
    print(
        f"Usage:\n        {sys.argv[0]} <csv file> [-d] [-s|-l] [-p] [-allfig|-minfig|-nofig]"
    )
    print(f"Example:\n        {sys.argv[0]} spectrumdata.csv")
    print("        <csv file> is a XPS data file in CSV-format.")
    print('        "-d" for developer mode')
    print('        "-s" for only Shirley background')
    print('        "-l" for only Linear background')
    print('        "-p" creats a summary ppt file')
    print('        "-allfig" will create all figures')
    print('        "-minfig" will create only minimum necessary figures (default)')
    print('        "-nofig" will not create any figures')
    sys.exit(1)


#
if __name__ == "__main__":
    args = sys.argv

    figType = "-minfig"

    # 引数に '-d' があれば開発者モード．全ての図を作る．中間ファイルを削除せずに残す．
    devFlag = "-d" in args
    if devFlag:
        figType = "-allfig"

    # 引数に '-p' があれば結果をまとめたpptファイルを作る．
    pptFlag = "-p" in args

    # デフォルトではまずShirley法を実施し，モデルが見つからない場合は直線法を実施する
    # 引数に '-l' があればバックグラウンドを直線法で近似する
    # 引数に '-s' があればバックグラウンドをShirley法で近似する
    # 引数に '-nofig' があれば図を作らない．時間短縮用
    # 引数に '-minfig' があれば必要最小限の図を作る．時間短縮用
    # 引数に '-allfig' があれば全ての図を作る．
    bgType = "auto"
    if "-s" in args:
        bgType = "shirley"
    elif "-l" in args:
        bgType = "linear"
    if "-nofig" in args:
        figType = "-nofig"
    elif "-minfig" in args:
        figType = "-minfig"
    elif "-allfig" in args:
        figType = "-allfig"

    args = [a for a in args if (not a.startswith("-"))]
    print(f"bgType = {bgType}")
    print(f"figType = {figType}")
    print(f"devFlag = {devFlag}")
    print(args)

    if len(args) <= 1:
        show_usage()

    main(args[1], bgType, figType, pptFlag, devFlag)
