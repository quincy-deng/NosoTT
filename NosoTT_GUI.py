import pandas as pd
from PyQt5.QtCore import QDate
from PyQt5.QtWidgets import (QPushButton, QLabel, QVBoxLayout, QHBoxLayout,
                             QWidget, QFileDialog, QTextEdit, QScrollArea,
                             QDialog, QTableWidget, QTableWidgetItem, QDateEdit)


class NosoTT_UI(QWidget):
    def __init__(self):
        super().__init__()

        self.select_df = pd.read_csv('casedatabase.csv', parse_dates=['Admdate', 'Chgdate', 'Sampdate'])

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
        casedate = CaseDate(self.select_df)
        casedate.btn_savecsv.clicked.connect(self.saveCSV)
        casedate.exec_()

    # 子窗口槽函数
    def saveCSV(self):
        pass

    def checkCSV(self):
        pass

    def execute(self):
        pass

    def choose_dir(self, line):
        file_text = QFileDialog.getExistingDirectory(self, '选择fastq/fasta数据文件夹', './')
        line.setText(file_text)


class CaseDate(QDialog):
    def __init__(self, df):
        super().__init__()
        self.btn_savecsv = QPushButton('保存')
        self.tableWidge = QTableWidget()
        self.setFixedSize(720, 540)
        self.df = df[['name', 'item', 'Admdate', 'Sampdate', 'Chgdate']]
        self.mainUI()
        self.show()

    def mainUI(self):
        layout = QHBoxLayout()
        scrollArea = QScrollArea()
        layout_sc = QHBoxLayout(scrollArea)
        row, col = self.df.shape
        self.tableWidge.setRowCount(row)
        self.tableWidge.setColumnCount(col)
        self.tableWidge.setHorizontalHeaderLabels(self.df.columns.tolist())
        for i in range(row):
            for j in range(2):
                self.tableWidge.setItem(i, j, QTableWidgetItem(self.df.iloc[i, j]))
        for i in range(row):
            for j in range(2, 5):
                date = self.df.iloc[i, j]
                qdate = QDate(date.year, date.month, date.day)
                datedit = QDateEdit(qdate)
                datedit.setDisplayFormat('yyyy-MM-dd')
                datedit.setCalendarPopup(True)
                self.tableWidge.setCellWidget(i, j, datedit)
        layout_sc.addWidget(self.tableWidge)
        layout_sc.addWidget(self.btn_savecsv)

        layout.addWidget(scrollArea)
        self.setLayout(layout)