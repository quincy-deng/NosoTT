import glob
import os
import datetime
import yaml
import subprocess

import numpy as np
import pandas as pd
from PyQt5.QtCore import QDate, QThread, pyqtSignal
from PyQt5.QtGui import QPixmap
from PyQt5.QtWidgets import (QPushButton, QLabel, QVBoxLayout, QHBoxLayout,
                             QWidget, QFileDialog, QTextEdit, QScrollArea,
                             QDialog, QTableWidget, QTableWidgetItem, QDateEdit, QMessageBox)


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


def checkFasta(select_items):
    items = {k: v for k, v in select_items.items() if
             os.path.isfile("data/fasta/{sample}.fasta".format(sample=k))}
    return items


def get_samples(names):
    fastas = names
    if os.path.isfile('data/mlst_database/STs.txt'):
        df = pd.read_csv('data/mlst_database/STs.txt', sep='\t',
                         names=['name', 'species', 'ST'] + ["gene" + str(i) for i in range(7)])
        fastas = df['name'][~df['name'].isin(names)].tolist()
    if len(fastas) > 0:
        os.system('conda run -n mlst mlst {} >> data/mlst_database/STs.txt'.format(' '.join(fastas)))
    df = pd.read_csv('data/mlst_database/STs.txt', sep='\t',
                     names=['name', 'species', 'ST'] + ["gene" + str(i) for i in range(7)])
    df = df[df['name'].isin(names)]
    df = df[df['ST'] != '-']
    df = df.drop_duplicates()
    items = {"{}_{}".format(sp, st): group['name'].tolist() for (sp, st), group in df.groupby(['species', 'ST']) if
             len(group) > 1}
    with open('config/nosott.yaml', 'w') as fl:
        yaml.dump(items, fl)
    return items


class NosoTT_UI(QWidget):
    def __init__(self):
        super().__init__()

        self.select_items = dict()

        # 控件
        self.btn_casedate = QPushButton('输入样本出入院/采样时间')
        self.text_output = QLabel('未选择文件夹,请指定')
        self.btn_output = QPushButton('选择文件夹')
        self.text_input = QLabel('未选择文件夹,请指定')
        self.btn_input = QPushButton('选择文件夹')
        self.btn_execute = QPushButton('执行')
        self.text_info = QTextEdit('')
        self.indirect_res = QScrollArea()
        self.indirect_widget = QWidget()
        self.indirect_layout = QVBoxLayout(self.indirect_widget)
        self.direct_res = QScrollArea()
        self.direct_widget = QWidget()
        self.direct_layout = QVBoxLayout(self.direct_widget)

        self.tab_nstUI()

        # 槽函数
        self.btn_input.clicked.connect(lambda: self.choose_dir(self.text_input))
        self.btn_output.clicked.connect(lambda: self.choose_dir(self.text_output))
        self.btn_casedate.clicked.connect(self.show_casedate)
        self.btn_execute.clicked.connect(self.execute)

    def tab_nstUI(self):
        layout = QHBoxLayout()

        left_widget = QWidget()
        left_panel = QVBoxLayout(left_widget)

        self.btn_input.setFixedSize(225, 30)
        self.btn_output.setFixedSize(225, 30)
        self.btn_casedate.setFixedSize(300, 30)

        left_panel.addWidget(QLabel('数据文件夹输入'))
        left_panel.addWidget(self.btn_input)
        left_panel.addWidget(self.text_input)
        left_panel.addStretch(1)
        left_panel.addWidget(self.btn_output)
        left_panel.addWidget(self.text_output)
        left_panel.addStretch(1)
        left_panel.addWidget(self.btn_casedate)
        left_panel.addStretch(2)
        left_panel.addWidget(self.btn_execute)
        left_panel.addStretch(1)
        left_panel.addWidget(QLabel("运行日志"))
        left_panel.addWidget(self.text_info)
        left_panel.addStretch(1)

        right_widget = QWidget()
        right_panel = QVBoxLayout(right_widget)
        right_panel.addWidget(QLabel('直接传播结果'))
        right_panel.addWidget(self.direct_res)
        right_panel.addWidget(QLabel('间接传播结果'))
        right_panel.addWidget(self.indirect_res)
        layout.addWidget(left_widget)
        layout.addWidget(right_widget)
        layout.setStretchFactor(left_widget, 3)
        layout.setStretchFactor(right_widget, 7)
        self.setLayout(layout)

    # NosoTT casedate窗口
    def show_casedate(self):
        if not os.path.isdir(self.text_input.text()):
            QMessageBox.information(self, "Error", "未指定输入文件夹")
        self.select_items = checkInputDir(self.text_input.text())
        if not self.select_items:
            QMessageBox.information(self, "Error", "请检查输入文件夹")
        casedate = CaseDate(self.select_items.keys())
        casedate.exec_()

    def execute(self):
        if not os.path.isdir(self.text_input.text()):
            QMessageBox.information(self, "Error", "未指定输入文件夹")
        elif not os.path.isdir(self.text_output.text()):
            QMessageBox.information(self, "Error", "未指定输出文件夹")
        elif not checkInputDir(self.text_input.text()):
            QMessageBox.information(self, "Error", "输入文件夹下没有测序文件")
        elif self.checkDateConfig():
            # 检查fasta文件是否存在,不存在组装,存在则跳过
            items = checkFasta(self.select_items)
            spades_command = []
            if not items:
                with open('config/spades.yaml', 'w') as fl:
                    yaml.dump(items, fl)
                spades_command = ['conda', 'run', '-n', 'nosott', 'snakemake', '-s', 'assembly.rules']
            nosott_command = ['conda', 'run', '-n', 'nosott', 'snakemake', '-s', 'NosoTT.rules']

            self.main_command = ExecuteCommand(spades_command, self.select_items, nosott_command)
            self.main_command.log_signal.connect(self.show_log)
            self.main_command.show_res.connect(self.show_results)
            self.main_command.start()

    def show_log(self, msg):
        self.text_info.append(msg)

    def show_results(self, msg):
        for i in msg:
            direct_res = QLabel()
            direct_res.setPixmap(QPixmap("{}/{}_direct_transmissions.jpg".format(self.text_output.text(),i)))
            print("{}/{}_direct_transmissions.jpg".format(self.text_output.text(),i))
            indirect_res = QLabel()
            indirect_res.setPixmap(QPixmap("{}/{}_indirect_transmissions.jpg".format(self.text_output.text(),i)))
            self.direct_layout.addWidget(QLabel(i))
            self.direct_layout.addWidget(direct_res)
            self.indirect_layout.addWidget(QLabel(i))
            self.indirect_layout.addWidget(indirect_res)
        self.direct_res.setWidget(self.direct_widget)
        self.indirect_res.setWidget(self.indirect_widget)

    def checkDateConfig(self):
        if not self.select_items:
            self.select_items = checkInputDir(self.text_input.text())
            if not self.select_items:
                QMessageBox.information(self, "Error", "输入文件夹下没有测序文件")
                return False
        db_df = pd.read_csv('data/casedate/casedatabase.csv')
        names = [i for i in self.select_items.keys() if i not in list(db_df['name'])]
        if names:
            QMessageBox.information(self, "Error", "缺少 {} 的出入院信息".format('\t'.join(names)))
            return False
        return True

    def choose_dir(self, line):
        file_text = QFileDialog.getExistingDirectory(self, '选择fastq/fasta数据文件夹', './')
        line.setText(file_text)


class ExecuteCommand(QThread):
    log_signal = pyqtSignal(str)
    show_res = pyqtSignal(list)

    def __init__(self, command, mlst_items, command2):
        super(ExecuteCommand, self).__init__()
        self.assemble_command = command
        self.items = mlst_items
        self.nosott_command = command2

    def run(self):
        if self.assemble_command:
            self.log_signal.emit("运行组装流程\n" + "-" * 50 + "\n")
            res1 = commandExecute(self.assemble_command)
            if res1:
                self.log_signal.emit("\n".join(["运行组装流程错误", "-" * 50, res1, "-" * 50]))
            self.log_signal.emit("-" * 50 + "\n" + "组装运行完毕\n" + "-" * 50 + "\n")
        else:
            self.log_signal.emit("样本已完成组装,无需重复运行\n" + "-" * 50 + "\n")
        items = self.items.keys()
        names = ["data/fasta/{}.fasta".format(i) for i in items]
        # 根据ST型结果,生成配置文件
        self.log_signal.emit("运行ST分型\n" + "-" * 50 + "\n")
        nosott_config = get_samples(names)
        if nosott_config:
            self.log_signal.emit("运行ST分型完毕\n" + "-" * 50 + "\n")
            self.log_signal.emit("运行传播路径分析流程\n" + "-" * 50 + "\n")
            res2 = commandExecute(self.nosott_command)
            if res2:
                self.log_signal.emit("\n".join(["运行传播分析错误", "-" * 50, res2, "-" * 50]))
            else:
                self.log_signal.emit("传播路径分析完成!\n" + "-" * 50 + "\n")
                self.show_res.emit(list(nosott_config.keys()))
        else:
            self.log_signal.emit("所选样本的ST型不一致,不存在传播关系")


class CaseDate(QDialog):
    def __init__(self, items):
        super().__init__()
        self.df = pd.read_csv('data/casedate/casedatabase.csv', parse_dates=['Admdate', 'Chgdate', 'Sampdate'])
        self.df = self.df[['name', 'item', 'Admdate', 'Sampdate', 'Chgdate']]
        dfn = pd.DataFrame({'name': list(items)})
        self.dfn = dfn.merge(self.df, how='left')
        self.tableWidge = QTableWidget()
        self.setFixedSize(720, 540)
        self.mainUI()
        self.show()

    def mainUI(self):
        layout = QHBoxLayout()
        scrollArea = QScrollArea()
        layout_sc = QHBoxLayout(scrollArea)
        row, col = self.dfn.shape
        self.tableWidge.setRowCount(row)
        self.tableWidge.setColumnCount(col)
        self.tableWidge.setHorizontalHeaderLabels(self.dfn.columns.tolist())
        for i in range(row):
            for j in range(2):
                text = ''
                if not pd.isna(self.dfn.iloc[i, j]):
                    text = self.dfn.iloc[i, j]
                self.tableWidge.setItem(i, j, QTableWidgetItem(text))
        for i in range(row):
            for j in range(2, 5):
                date = self.dfn.iloc[i, j]
                if pd.isna(date):
                    date = datetime.date.today()
                qdate = QDate(date.year, date.month, date.day)
                datedit = QDateEdit(qdate)
                datedit.setDisplayFormat('yyyy-MM-dd')
                datedit.setCalendarPopup(True)
                self.tableWidge.setCellWidget(i, j, datedit)
        layout_sc.addWidget(self.tableWidge)

        layout.addWidget(scrollArea)
        self.setLayout(layout)

    def checkCSV(self):
        rowLength = self.tableWidge.rowCount()
        names, items, Admdate, Sampdate, Chgdate = [], [], [], [], []
        caseDict = {'name': names, 'item': items, 'Admdate': Admdate, 'Sampdate': Sampdate, 'Chgdate': Chgdate}
        for i in range(rowLength):
            names.append(self.tableWidge.item(i, 0).text())
            items.append(self.tableWidge.item(i, 1).text())
            Admdate.append(self.tableWidge.cellWidget(i, 2).date().toString("yyyy-MM-dd"))
            Sampdate.append(self.tableWidge.cellWidget(i, 3).date().toString("yyyy-MM-dd"))
            Chgdate.append(self.tableWidge.cellWidget(i, 4).date().toString("yyyy-MM-dd"))
        caseDataframe = pd.DataFrame(caseDict)
        lack_name = caseDataframe[caseDataframe.isna().T.any()]['name'].to_list()
        if lack_name:
            QMessageBox.information(self, 'Error', '未指定名称{}'.format('\t'.join(lack_name)))
        time_error = caseDataframe[~((caseDataframe['Admdate'] < caseDataframe['Sampdate'])
                                     & (caseDataframe['Sampdate'] < caseDataframe['Chgdate']))]['name'].tolist()
        if time_error:
            QMessageBox.information(self, 'Error', '未指定名称{}\n'.format('\t'.join(time_error)))
        caseDateDF = pd.read_csv('data/casedate/casedatabase.csv')
        rm_samples = list(set(lack_name + time_error))
        caseDataframe = caseDataframe[~caseDataframe['name'].isin(rm_samples)]
        db = caseDataframe.merge(caseDateDF, how='outer')
        db = db.drop_duplicates()
        try:
            db.to_csv('data/casedate/casedatabase.csv', index=False)
        except IOError as e:
            QMessageBox(self, '警告', '无法保存流行病学信息文件{}'.format(e))

    def closeEvent(self, event):
        self.checkCSV()
        super(CaseDate, self).closeEvent(event)
