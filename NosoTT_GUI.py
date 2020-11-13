import glob
import os
import datetime

import numpy as np
import pandas as pd
from PyQt5.QtCore import QDate
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
        self.direct_res = QScrollArea()
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
        right_panel.addWidget(self.direct_res)
        right_panel.addWidget(self.indirect_res)
        # print('xx')
        right_panel.addWidget(QLabel('直接传播结果'))
        layout.addWidget(left_widget)
        right_panel.addWidget(QLabel('间接传播结果'))
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
        if not checkInputDir(self.text_input.text()):
            QMessageBox.information(self, "Error", "输入文件夹下没有测序文件")
        self.checkDateConfig()
        # todo:检查出入院时间表格是否完整,符合规范.
        # TODO:第一步调用组装流程
        # TODO: 第二步运行主程序

    def checkDateConfig(self):
        if not self.select_items:
            self.select_items = checkInputDir(self.text_input.text())
            if not self.select_items:
                QMessageBox.information(self, "Error", "输入文件夹下没有测序文件")
        db_df = pd.read_csv('data/casedate/casedatabase.csv')
        for i in self.select_items.keys():
            if i not in db_df['item']:
                QMessageBox.information(self, "Error", "缺少 {} 的出入院信息".format(i))

    def choose_dir(self, line):
        file_text = QFileDialog.getExistingDirectory(self, '选择fastq/fasta数据文件夹', './')
        line.setText(file_text)


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
        time_error = caseDataframe[~((caseDataframe['Admdate'] < caseDataframe['Sampdate'])
                                     & (caseDataframe['Sampdate'] < caseDataframe['Chgdate']))]['name'].tolist()
        QMessageBox.information(self,'Error','')
        caseDateDF = pd.read_csv('data/casedate/casedatabase.csv')
        db = pd.concat([caseDataframe, caseDateDF])
        db.to_csv('data/casedate/casedatabase.csv', index=False)
        return True

    def closeEvent(self, event):
        self.saveCSV()
        super(CaseDate, self).closeEvent(event)

    def saveCSV(self):
        if not self.checkCSV():
            QMessageBox.information(self, '警告', '无法保存流行病学信息')
