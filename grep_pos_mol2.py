import pandas as pd
import numpy as np
import glob,os,time
import xlrd
from PyQt5.QtWidgets import (QApplication,QFileDialog,QPushButton,QLineEdit,QTextEdit,QWidget,QHBoxLayout,
                             QVBoxLayout,QLabel,QDialog,QMessageBox,QGraphicsScene,QGraphicsView,QTableView,
                             QScrollArea,QGridLayout,QRadioButton,QComboBox)
from rdkit import Chem
from rdkit.Chem import Draw,AllChem
from rdkit import DataStructs
from PyQt5.QtPrintSupport import QPrintDialog,QPrinter
from PyQt5.QtGui import QIntValidator,QPixmap
from PyQt5.QtCore import QAbstractTableModel,Qt,QThread,pyqtSignal,QSize
import sys


class Demo(QWidget):
    def __init__(self):
        super().__init__()
        self.setFixedSize(800,600)
        self.loci = dict()
        # 阳性化合物库
        self.select_df = dict()
        self.db = pd.read_table('database/merge_db.tab')
        self.printer = QPrinter()
        self.mainUI()
        self.dir_input.clicked.connect(self.choose_dir)
        self.execute.clicked.connect(self.execut)
        self.btn3.clicked.connect(self.print_dialog)
        self.btn5.clicked.connect(self.show_mol)
        self.btn4.clicked.connect(self.set_output_filename)

    def mainUI(self):
        layout = QHBoxLayout()

        left_widget = QWidget()
        left_panel = QVBoxLayout(left_widget)
        self.dir_input = QPushButton("输入文件夹")
        self.dir_input.setFixedSize(100,30)
        self.qline = QLabel('J:/DQY')
        self.qline.setFixedSize(225,30)
        self.execute = QPushButton('执行')
        self.execute.setFixedSize(225,30)
        self.line2 = QLineEdit()
        self.line2.setValidator(QIntValidator(0,99)),self.line2.setText('50')
        self.btn3 = QPushButton("打印结果")
        self.btn3.setEnabled(False)
        self.btn3.setFixedSize(225,30)
        self.btn4 = QPushButton('导出阳性化合物库信息')
        self.btn4.setEnabled(False)
        self.btn4.setFixedSize(225,30)
        self.btn5 = QPushButton('显示阳性分子结构')
        self.btn5.setEnabled(False)
        self.btn5.setFixedSize(225,30)
        left_panel.addWidget(self.dir_input)
        left_panel.addWidget(self.qline)
        left_panel.addWidget(QLabel('设置抑菌比例(1-99),默认50%'))
        left_panel.addWidget(self.line2)
        left_panel.addStretch(1)
        left_panel.addWidget(self.execute)
        left_panel.addStretch(6)
        left_panel.addWidget(self.btn3)
        left_panel.addWidget(self.btn4)
        left_panel.addWidget(self.btn5)
        left_panel.addStretch(6)

        right_widget = QWidget()
        right_panel = QVBoxLayout(right_widget)
        res = QLabel('小于指定比例(如50%)酶活(DMSO平均值)的孔: ')
        self.res1 = QTextEdit()
        self.res1.setFontFamily('Consolas')
        self.res1.setFontPointSize(12)
        error = QLabel('xls文件有误: ')
        self.error1 = QTextEdit()
        self.error1.setFontFamily('Consolas')
        self.error1.setFontPointSize(12)
        right_panel.addWidget(res)
        right_panel.addWidget(self.res1)
        right_panel.addWidget(error)
        right_panel.addWidget(self.error1)

        layout.addWidget(left_widget)
        layout.setStretchFactor(left_widget,3)
        layout.addWidget(right_widget)
        layout.setStretchFactor(right_widget,7)
        self.setLayout(layout)

    def choose_dir(self):
        self.file_text =  QFileDialog.getExistingDirectory(self,'选择xls文件夹','./')
        self.qline.setText(self.file_text)

    def set_output_filename(self):
        self.directory = QFileDialog.getSaveFileName(self,"输出文件名","./","Excel Files  (*.xlsx)") 
        write = pd.ExcelWriter(self.directory[0])
        for sheet_name,df in self.select_df.items():
            df.to_excel(write,sheet_name=sheet_name,index=False)
        write.save()
        if os.path.isfile(self.directory[0]):
            QMessageBox.information(self,'','文件保存成功')

    def execut(self):
        self.btn3.setEnabled(False)
        self.btn4.setEnabled(False)
        self.btn5.setEnabled(False)
        self.error1.setText('')
        self.res1.setText('')
        path = self.qline.text()
        ratio = 1 - int(self.line2.text())/100
        if not os.path.exists(path):
            QMessageBox.information(self,'错误','文件夹未指定')
        xlss = [i for i in glob.glob('{}/**/*[0-9].xls'.format(path), recursive=True) if os.path.isfile(i)]
        if len(xlss) == 0:
            QMessageBox.information(self,'错误','所选文件夹下没有xls文件')
            
        # 子线程运行,如果放主线程容易卡死
        self.run_search = Search_Mol(xlss,ratio,self.db)
        self.run_search.error_signal.connect(self.show_error)
        self.run_search.res_singnal.connect(self.show_result)
        self.run_search.selected_df_signal.connect(self.get_select_df)
        self.run_search.start()

    def show_error(self,msg):
        self.error1.append(msg)

    def show_result(self,msg):
        self.res1.append(msg)

    def get_select_df(self,df):
        self.select_df = df
        self.btn3.setEnabled(True)
        self.btn4.setEnabled(True)
        self.btn5.setEnabled(True)

    def print_dialog(self):
        printdialog = QPrintDialog(self.printer,self)
        if QDialog.Accepted == printdialog.exec_():
            self.res1.print(self.printer)

    def show_mol(self):
        x = Show_mol_structure(self.select_df)
        x.show()
        x.exec_()


class Search_Mol(QThread):
    error_signal = pyqtSignal(str)
    res_singnal = pyqtSignal(str)
    selected_df_signal = pyqtSignal(pd.DataFrame)
    def __init__(self,xlss,ratio,db):
        super().__init__()
        self.xlss = sorted(xlss)
        self.ratio = ratio
        self.db = db

    def run(self):
        ind = [chr(i).upper() for i in range(97,123)]
        select_df = pd.DataFrame(columns=self.db.columns.tolist())
        n = 0
        for excel in self.xlss:
            filename = os.path.basename(excel).rstrip('.xls')
            try:
                wb = xlrd.open_workbook(excel,encoding_override="gb2312")
                df = pd.read_excel(wb,header=1)
                col_index = [str(i) for i in range(2,12)]
                dmso_mean = df['12'][:4].mean()*self.ratio
                dff = df[col_index]
                x = np.where(dff < dmso_mean)
                n += len(x[0])
                if len(x[0]) > 0:
                    for k,j in zip(x[0],x[1]):
                        db_1,p_num = filename.split('-')
                        p_num = int(p_num)
                        mol_name,temp_df = self.get_molname(db_1,p_num,ind[k],j+2)
                        self.res_singnal.emit('{} {}行,{:<2}列; 名称:{} 值:{:<5}.\n'.format(filename,k+1,j+2,mol_name,round(dff.iloc[k,j],3)))
                        select_df = select_df.append(temp_df,ignore_index=True)
                        # print(select_df)
            except:
                self.error_signal.emit("错误文件: "+filename)
        if select_df.shape[0] > 0:
            self.selected_df_signal.emit(select_df)
        self.res_singnal.emit('一共找到 {} 孔'.format(n))

    def get_molname(self,db,plate,row,col):
        x = self.db.loc[(self.db['db'] == db) & (self.db['Row'] == row )& (self.db['Col'] == col) & (self.db['plate'] == plate)]
        return [list(x['MOLENAME'])[0],x]

class Show_mol_structure(QDialog):
    def __init__(self,select_df):
        super().__init__()
        self.select_df = select_df
        self.setWindowTitle('分子结构图')
        self.col = 0
        self.row = 0
        self.max_col = 5

        self.mainUI()
        self.start_parse_mol()
        
    def mainUI(self):
        layout = QVBoxLayout()
        layout.addWidget(QLabel('根据SMIILES绘制的分子结构式如下:'))

        self.scrollarea = QScrollArea()
        self.scrollarea.setWidgetResizable(True)
        
        self.scroll_content = QWidget()
        self.gridlayout = QGridLayout(self.scroll_content)
        self.scrollarea.setWidget(self.scroll_content)
        desktop = QApplication.desktop()
        width,height = desktop.width(),desktop.height()
        self.scrollarea.setFixedSize(width-100,height-100)

        layout.addWidget(self.scrollarea)
        self.setLayout(layout)

    def start_parse_mol(self):
        for index,line in self.select_df[['db','plate','Row','Col','SMILES','MOLENAME']].iterrows():
            db,plate,row,col,smile,name = list(line)
            n = Chem.MolFromSmiles(smile)
            m = Draw.MolToQPixmap(n)
            pixmap = QPixmap(m)
            pixmap_id = '-'.join([db,str(plate),str(row),str(col),str(name)])
            item = MolImgItem(pixmap,pixmap_id,smile)
            self.addPixmap(item)
            QApplication.processEvents()

    def addPixmap(self,item):
        self.gridlayout.addWidget(item,self.row,self.col)
        if self.col < self.max_col-1:
            self.col += 1
        else:
            self.col = 0
            self.row += 1        

class Search_Similar_Mol(QThread):
    res_signal = pyqtSignal(pd.DataFrame)
    def __init__(self,db,select_mol):
        super().__init__()
        self.db = db
        self.select_mol = select_mol

    def get_fingerprint(self,smile):
        try:
            m = AllChem.MolFromSmiles(smile)
            mol_fps = AllChem.GetMorganFingerprint(m,2)
            return mol_fps
        except:
            return None
    
    def run(self):
        fp = self.get_fingerprint(self.select_mol)
        fps = [self.get_fingerprint(smile) if self.get_fingerprint(smile) is not None else 0 for smile in self.db['SMILES'] ]
        self.db['similar'] = [0 if i==0 else round(DataStructs.TanimotoSimilarity(fp, i),2) for i in fps]
        top20 = self.db.sort_values('similar',ascending=False).head(20)
        self.res_signal.emit(top20)

class MolImgItem(QWidget):
    def __init__(self,pixmap,imgID,smile,width=100,height=100):
        super().__init__()
        self.pixmap = pixmap
        self.imgID = imgID
        self.smile = smile
        self.width = width
        self.height = height
        self.mainUI()
        self.smilar.clicked.connect(self.show_similar_mol)
    
    def mainUI(self):
        self.resize(self.width,self.height)
        layout = QVBoxLayout()
        self.label1 = QLabel()
        self.label2 =QLabel()
        self.smilar = QPushButton('寻找类似物(前20)')
        self.smilar.setFixedSize(225,30)
        if self.pixmap:
            # pixmap = self.pixmap.scaled(QSize(self.width,self.height),Qt.KeepAspectRatio,Qt.SmoothTransformation)
            self.label1.setPixmap(self.pixmap)
            self.label1.setAlignment(Qt.AlignCenter)
            layout.addWidget(self.label1)
        if self.imgID:
            self.label2.setText(self.imgID)
            self.label2.setAlignment(Qt.AlignCenter)
            # self.label2.adjustSize()
            layout.addWidget(self.label2)
            layout.addWidget(self.smilar)
        self.setLayout(layout)
    
    def show_similar_mol(self):
        db = pd.read_table('database/merge_db.tab')
        self.run_search = Search_Similar_Mol(db,self.smile)
        self.run_search.res_signal.connect(self.show_result)
        self.run_search.start()

    def show_result(self,msg):
        x = Show_mol_structure(msg)
        x.show()
        x.exec_()


if __name__ == "__main__":
    app = QApplication(sys.argv)
    demo = Demo()
    demo.show()
    sys.exit(app.exec_())
