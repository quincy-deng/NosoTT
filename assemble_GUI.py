from PyQt5.QtWidgets import QPushButton, QLabel, QVBoxLayout, QHBoxLayout, QWidget, QScrollArea, QFileDialog, \
    QTextEdit, QTableWidget, QTableWidgetItem
from PyQt5.QtCore import QThread, pyqtSignal

import pandas as pd
import glob
import subprocess
import os
import yaml


def replaceSuffix(filename):
    filename = os.path.basename(filename)
    suffixs = ["{}{}{}.{}.gz".format(i, j, n, k) for i in ["_", "."] for n in ["1", "2"] for j in ["R", ""] for k in
               ["fastq", "fq"]]
    for suffix in suffixs:
        filename = filename.replace(suffix, "")
    return filename


def produceFileConfig(files):
    filenames = [replaceSuffix(i) for i in files]
    file_dict = dict()
    for k, v in zip(filenames, files):
        if k not in file_dict:
            file_dict[k] = [v]
        else:
            file_dict[k].append(v)
    F_dict = {k: v for k, v in file_dict.items() if len(v) == 2}
    return F_dict


def checkInputDir(inputDir):
    types = ["[1-2].fq.gz", "[1-2].fastq.gz"]
    files = sum(
        [glob.glob("{input}/**/*{typen}".format(input=inputDir, typen=type1), recursive=True) for type1 in types],
        [])
    return produceFileConfig(files)


def commandExecute(command):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    result = process.communicate()
    process.wait()
    if process.returncode != 0:
        stderr = str(result[0], encoding='utf8')
        return stderr
    return None


def generate_config(samples):
    with open("config/spades.yaml", "w") as fl:
        yaml.dump(samples, fl)


class Assembly(QWidget):
    def __init__(self):
        super().__init__()
        self.inputBtn = QPushButton("选择原始数据文件夹")
        self.inputBtn.setFixedSize(225, 30)
        self.inputText = QLabel('未指定输入文件夹')
        self.outputBtn = QPushButton("选择输出文件夹")
        self.outputBtn.setFixedSize(225, 30)
        self.outputText = QLabel('output')
        self.executeBtn = QPushButton("执行")
        self.executeBtn.setFixedSize(225, 40)
        self.resTable = QTableWidget()
        self.logInfo = QTextEdit()
        self.set_layout()

        # 槽函数
        self.inputBtn.clicked.connect(lambda: self.file_dialog(self.inputText))
        self.outputBtn.clicked.connect(lambda: self.file_dialog(self.outputText))
        self.executeBtn.clicked.connect(self.execute)

    def set_layout(self):
        layout = QVBoxLayout()

        up_widget = QWidget()
        up_layout = QHBoxLayout(up_widget)
        up_left_widget = QWidget()
        up_left_layout = QVBoxLayout(up_left_widget)
        up_left_layout.addWidget(self.inputBtn)
        up_left_layout.addWidget(self.inputText)
        up_left_layout.addStretch(1)
        up_left_layout.addWidget(self.outputBtn)
        up_left_layout.addWidget(self.outputText)
        up_left_layout.addStretch(1)
        up_left_layout.addWidget(self.executeBtn)
        up_left_layout.addStretch(2)
        upRightWidget = QWidget()
        upRightlayout = QVBoxLayout(upRightWidget)
        upRightlayout.addWidget(self.logInfo)
        up_layout.addWidget(up_left_widget)
        up_layout.addWidget(upRightWidget)

        down_widet = QWidget()
        down_layout = QVBoxLayout(down_widet)
        down_layout.addWidget(QLabel("运行结果"))
        down_layout.addWidget(self.resTable)

        layout.addWidget(up_widget)
        layout.setStretchFactor(up_widget, 3)
        layout.addWidget(down_widet)
        layout.setStretchFactor(down_widet, 7)
        self.setLayout(layout)

    def file_dialog(self, line):
        selectDir = QFileDialog.getExistingDirectory(self, '选择文件夹', './')
        line.setText(selectDir)

    def check_fastq(self):
        return checkInputDir(self.inputText.text())

    def show_quast_result(self):
        df = pd.read_table("data/quast_results/report.tsv")
        row, col = df.shape
        print(df.shape)
        self.resTable.setRowCount(row)
        self.resTable.setColumnCount(col)
        self.resTable.setHorizontalHeaderLabels(df.columns.tolist())
        [self.resTable.setItem(i, j, QTableWidgetItem(str(df.iloc[i, j]))) for i in range(row) for j in range(col)]

    def show_log(self, msg):
        self.logInfo.append(msg)

    def execute(self):
        if not os.path.isdir(self.inputText.text()):
            self.logInfo.append("没有选择测序文件夹")
        else:
            samples = self.check_fastq()
            if samples:
                self.logInfo.append("找到如下样本:\n{}\n".format("\t".join(samples)))
                generate_config(samples)
                spades_command = ['conda', 'run', '-n', 'nosott', 'snakemake', '--core', '15', '-s', 'assembly.rules']
                self.logInfo.append("开始运行组装流程\n" + "-" * 50)
                self.assemble_subprocess = Execute(spades_command)
                self.assemble_subprocess.res_signal.connect(self.show_log)
                self.assemble_subprocess.show_signal.connect(self.show_quast_result)
                self.assemble_subprocess.start()
            else:
                self.logInfo.append("Error:文件夹下找不到测序文件")


class Execute(QThread):
    res_signal = pyqtSignal(str)
    show_signal = pyqtSignal(str)

    def __init__(self, command):
        super().__init__()
        self.command = command

    def run(self):
        res = commandExecute(self.command)
        if res:
            self.res_signal.emit("运行错误:\n" + res + "\n运行错误")
        else:
            self.res_signal.emit("运行成功\n" + "-" * 50)
            self.show_signal.emit("ok")
