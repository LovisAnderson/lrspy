from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit, QFileDialog,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget, )
from PyQt5.QtCore import QRect,Qt, QVariant,QModelIndex
from crisscross import CrissCross
from reader import reader

import pandas as pd


class WidgetGallery(QDialog):
    def __init__(self, parent=None):
        super(WidgetGallery, self).__init__(parent)


        self.canvasWidth = 300
        self.canvasHeight = 300
        self.createCanvas(width=self.canvasWidth, height=self.canvasHeight)
        self.createControls()
        self.createHyperplanes()
        self.createMatrix()

        mainLayout = QGridLayout()
        mainLayout.addLayout(self.controls, 0, 0, 1, 2)
        mainLayout.addLayout(self.canvas, 1, 0, 2, 2)
        mainLayout.addWidget(self.matrixDisplay, 1, 1)
        mainLayout.addWidget(self.hyperplaneDisplay, 2, 1)
        self.setLayout(mainLayout)
        self.setWindowTitle("Lrs")


    def createCanvas(self, width=300, height=300):
        layout = QVBoxLayout()
        layout.setGeometry(QRect(300, 300, width, height))
        self.canvas = layout

    def createHyperplanes(self):
        self.hyperplaneDisplay = QGroupBox("Hyperplanes")
        layout = QVBoxLayout()
        self.hyperplaneDisplay.setLayout(layout)

    def createControls(self):
        self.controls = QGroupBox("Controls")

        pivotButton = QPushButton("Pivot")
        pivotButton.setDefault(True)

        searchStepButton = QPushButton("Search Step")
        self.button.clicked.connect(self.searchButton)
        searchStepButton.setDefault(True)

        fileButton = QPushButton("Open File")
        fileButton.setDefault(True)
        fileButton.clicked.connect(self.openFile)

        layout = QVBoxLayout()
        layout.addWidget(pivotButton)
        layout.addWidget(searchStepButton)
        layout.addWidget(fileButton)
        self.controls = layout

    def searchButton(self):
        self.lrs.search()

    def createMatrix(self):
        self.matrixDisplay = QGroupBox("Matrix")
        layout = QVBoxLayout()
        self.matrixDisplay.setLayout(layout)

    def openFile(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Ine Files (*.ine)",
                                                  options=options)
        if fileName:
            self.lrs = CrissCross(*reader(fileName))

    def drawHyperplanes(self):
        pass


class PandasModel(QtCore.QAbstractTableModel):
    def __init__(self, df = pd.DataFrame(), parent=None):
        QtCore.QAbstractTableModel.__init__(self, parent=parent)
        self._df = df

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant()

        if orientation == QtCore.Qt.Horizontal:
            try:
                return self._df.columns.tolist()[section]
            except (IndexError, ):
                return QVariant()
        elif orientation == Qt.Vertical:
            try:
                # return self.df.index.tolist()
                return self._df.index.tolist()[section]
            except (IndexError, ):
                return QVariant()

    def data(self, index, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant()

        if not index.isValid():
            return QVariant()

        return QVariant(str(self._df.ix[index.row(), index.column()]))

    def setData(self, index, value, role):
        row = self._df.index[index.row()]
        col = self._df.columns[index.column()]
        if hasattr(value, 'toPyObject'):
            # PyQt4 gets a QVariant
            value = value.toPyObject()
        else:
            # PySide gets an unicode
            dtype = self._df[col].dtype
            if dtype != object:
                value = None if value == '' else dtype.type(value)
        self._df.set_value(row, col, value)
        return True

    def rowCount(self, parent=QModelIndex()):
        return len(self._df.index)

    def columnCount(self, parent=QModelIndex()):
        return len(self._df.columns)

    def sort(self, column, order):
        colname = self._df.columns.tolist()[column]
        self.layoutAboutToBeChanged.emit()
        self._df.sort_values(colname, ascending= order == Qt.AscendingOrder, inplace=True)
        self._df.reset_index(inplace=True, drop=True)
        self.layoutChanged.emit()

import sys
app = QApplication(sys.argv)
gallery = WidgetGallery()
gallery.show()
sys.exit(app.exec_())