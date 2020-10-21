from PyQt5.QtWidgets import (QApplication, QWidget, QTabWidget, QSplashScreen, )
from PyQt5.QtCore import Qt
from PyQt5.QtGui import QPixmap
import time
from NosoTT_GUI import NosoTT_UI
from MutiDetect_GUI import MultiDetect
import sys


class Main(QTabWidget):
    def __init__(self):
        super().__init__()
        self.setFixedSize(720, 540)

        self.multi_tab = QWidget()
        self.assemble_tab = MultiDetect()
        self.nst_tab = NosoTT_UI()

        self.addTab(self.nst_tab, 'NosoTT')
        self.addTab(self.assemble_tab, 'Assemble RawData')
        self.addTab(self.multi_tab, 'Multi Type Detect')


class SplashScreen(QSplashScreen):
    def __init__(self):
        super(SplashScreen, self).__init__(QPixmap("qidong.jpg"), flags=Qt.WindowStaysOnTopHint)

    def effect(self):
        self.setWindowOpacity(0)
        t = 0
        while t <= 50:
            newOpacity = self.windowOpacity() + 0.02
            if newOpacity > 1:
                break
            self.setWindowOpacity(newOpacity)
            self.show()
            t -= 1
            time.sleep(0.05)

        time.sleep(1)
        t = 0
        while t <= 50:
            newOpacity = self.windowOpacity() - 0.02
            if newOpacity < 0:
                break
            self.setWindowOpacity(newOpacity)
            self.show()
            t += 1
            time.sleep(0.05)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    splash = SplashScreen()
    splash.show()
    splash.showMessage('启动中...', Qt.AlignBottom | Qt.AlignCenter, Qt.black)
    time.sleep(1.5)
    demo = Main()
    splash.finish(demo)
    demo.show()
    sys.exit(app.exec_())
