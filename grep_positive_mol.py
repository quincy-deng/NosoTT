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
# from rdkit.Chem.Draw import IPythonConsole
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
        self.db = self.get_database()
        self.printer = QPrinter()
        self.mainUI()
        self.dir_input.clicked.connect(self.choose_dir)
        self.execute.clicked.connect(self.execut)
        self.btn3.clicked.connect(self.print_dialog)
        self.btn5.clicked.connect(self.show_mol)
        self.btn4.clicked.connect(self.set_output_filename)
        self.btn6.clicked.connect(self.search_like_mol)

    def mainUI(self):
        layout = QHBoxLayout()

        left_widget = QWidget()
        left_panel = QVBoxLayout(left_widget)
        self.dir_input = QPushButton("输入文件夹")
        self.dir_input.setFixedSize(100,30)
        self.qline = QLabel('未选择文件夹,请指定')
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
        self.btn6 = QPushButton('寻找化合类似物')
        self.btn6.setEnabled(False)
        self.btn6.setFixedSize(225,30)
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
        left_panel.addWidget(self.btn6)
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

    def get_database(self):
        d7800 = r'./database/D7800-Bioactive_Compound_Library_Plus.xlsx'
        l4000 = r'./database/L4000-Bioactive_Compound_Library.xlsx'
        l2500 = r'./database/L2500-Human_Endogenous_Metabolite_Compound_Library.xlsx'

        df7800 = pd.read_excel(d7800,sheet_name='Compound Information')
        df4000 = pd.read_excel(l4000,sheet_name='Compound Information')
        df2500 = pd.read_excel(l2500,sheet_name='Compound Information')

        p1,p2,p3 = list(set(df2500['Plate'])),list(set(df4000['Plate'])),list(set(df7800['Plate']))
        x1 = {str(i+1):j for i,j in enumerate(sorted(p1))}
        x2 = {str(i+1):j for i,j in enumerate(sorted(p2))}
        x3 = {str(i+1):j for i,j in enumerate(sorted(p3))}
        x = {'L2500':x1,'L4000':x2,'D7800':x3}
        db = [x,{'L2500':df2500,'L4000':df4000,'D7800':df7800}]
        return db

    def execut(self):
        self.btn3.setEnabled(False)
        self.btn4.setEnabled(False)
        self.btn5.setEnabled(False)
        self.btn6.setEnabled(False)
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
        self.btn6.setEnabled(True)

    def print_dialog(self):
        printdialog = QPrintDialog(self.printer,self)
        if QDialog.Accepted == printdialog.exec_():
            self.res1.print(self.printer)

    def show_mol(self):
        x = Dialog(self.select_df)
        x.show()
        x.exec_()

    def search_like_mol(self):
        db = pd.read_table('database/merge_db.tab')
        x = Show_similar_mol(self.select_df,db)
        x.show()
        x.exec_()
                
class Search_Mol(QThread):
    error_signal = pyqtSignal(str)
    res_singnal = pyqtSignal(str)
    selected_df_signal = pyqtSignal(dict)
    def __init__(self,xlss,ratio,db):
        super().__init__()
        self.xlss = sorted(xlss)
        self.ratio = ratio
        self.db = db

    def run(self):
        ind = [chr(i).upper() for i in range(97,123)]
        res = {'L2500':[],'L4000':[],'D7800':[]}
        select_df = dict()
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
                        mol_name,line_df = self.get_molname(db_1,p_num,ind[k],j+2)
                        res[db_1].append(line_df)
                        self.res_singnal.emit('{} {}行,{:<2}列; 名称:{} 值:{:<5}.\n'.format(filename,k+1,j+2,mol_name,round(dff.iloc[k,j],3)))

            except:
                self.error_signal.emit("错误文件: "+filename)
        for  k,v in res.items():
            if len(v) > 0:
                df_temp = pd.concat(v)
                select_df[k] = df_temp
        if len(select_df) > 0:
            self.selected_df_signal.emit(select_df)
        self.res_singnal.emit('一共找到 {} 孔'.format(n))

    def get_molname(self,db_1,p_num,row,col):
        plate_dict,df_dict = self.db
        p_id = plate_dict[db_1][p_num]
        df = df_dict[db_1]
        x = df.loc[(df['Plate']==p_id) & (df['Row'] == row )& (df['Col'] == col)]
        return [list(x['MOLENAME'])[0],x]

class Dialog(QDialog):
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
        for database,df in self.select_df.items():
            for smile,col,row,name,plate in zip(df['SMILES'],df['Col'],df['Row'],df['MOLENAME'],df['Plate']):
                n = Chem.MolFromSmiles(smile)
                m = Draw.MolToQPixmap(n)
                pixmap = QPixmap(m)
                pixmap_id = '-'.join([database,plate,str(row),str(col),str(name)])
                item = MolImgItem(pixmap,pixmap_id)
                self.addPixmap(item)
                QApplication.processEvents()

    def addPixmap(self,item):
        self.gridlayout.addWidget(item,self.row,self.col)
        if self.col < self.max_col-1:
            self.col += 1
        else:
            self.col = 0
            self.row += 1

class Show_similar_mol(QDialog):
    def __init__(self,select_df,db):
        super().__init__()
        self.select_df = select_df
        self.db = db
        self.dff = {database+'\t'+molname:smile for database,df in self.select_df.items() for molname,smile in zip(df['MOLENAME'],df['SMILES'])}
        self.setWindowTitle('查找类似物')
        desktop = QApplication.desktop()
        width,height = desktop.width(),desktop.height()
        self.setFixedSize(width-200,height-200)
        self.mainUI()
        self.btn.clicked.connect(self.show_res)
    
    def mainUI(self):
        layout = QHBoxLayout()
        left_panel = QWidget()
        left_layout = QVBoxLayout(left_panel)
        left_layout.addWidget(QLabel('选择化合物'))
        self.combox = QComboBox()
        self.combox.addItems(self.dff.keys())
        self.combox.setCurrentIndex(1)
        left_layout.addWidget(self.combox)
        self.btn = QPushButton('执行(默认前20)')
        self.btn.setFixedSize(225,30)
        left_layout.addWidget(self.btn)

        right_panel =QWidget()
        right_layout = QHBoxLayout(right_panel)
        self.res = QTextEdit()
        right_layout.addWidget(self.res)

        layout.addWidget(left_panel)
        layout.setStretchFactor(left_panel,2)
        layout.addWidget(right_panel)
        layout.setStretchFactor(right_panel,8)
        self.setLayout(layout)
    
    def show_res(self):
        self.select_mol = self.combox.currentText()
        # 子线程运行,如果放主线程容易卡死
        self.run_search = Search_Similar_Mol(self.select_df,self.db,self.dff,self.select_mol)
        self.run_search.res_signal.connect(self.show_result)
        self.run_search.start()
    
    def show_result(self,msg):
        self.res.setText(msg)

class Search_Similar_Mol(QThread):
    res_signal = pyqtSignal(str)
    def __init__(self,select_df,db,dff,select_mol):
        super().__init__()
        self.dff = dff
        self.select_mol = ''
        self.select_df = select_df
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
        # if self.combox.currentText() == self.select_mol:
        #     QMessageBox.information(self,'错误','运行已完成')
        
        smile = self.dff[self.select_mol]
        fp = self.get_fingerprint(smile)
        fps = [self.get_fingerprint(smile) if self.get_fingerprint(smile) is not None else 0 for smile in self.db['SMILES'] ]
        self.db['similar'] = [0 if i==0 else round(DataStructs.TanimotoSimilarity(fp, i),2) for i in fps]
        top20 = self.db[['db','plate','Col','Row','similar','MOLENAME']].sort_values('similar',ascending=False).head(20)
        # top20 = sorted(x,key = lambda i: i[1], reverse = True)[:20]
        # df = self.db[['db','plate','Col','Row','MOLENAME']].iloc[[i for i,j in top20],]
        res = '\n'.join(['\t'.join([str(i) for i in line]) for index,line in top20.iterrows()])
        self.res_signal.emit(res)


class MolImgItem(QWidget):
    def __init__(self,pixmap,imgID,width=100,height=100):
        super().__init__()
        self.pixmap = pixmap
        self.imgID = imgID
        self.width = width
        self.height = height
        self.mainUI()
    
    def mainUI(self):
        self.resize(self.width,self.height)
        layout = QVBoxLayout()
        self.label1 = QLabel()
        self.label2 =QLabel()
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
        self.setLayout(layout)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    demo = Demo()
    demo.show()
    sys.exit(app.exec_())
