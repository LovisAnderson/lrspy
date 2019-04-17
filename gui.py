from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit, QFileDialog,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget, )
from PyQt5.QtCore import QRect, Qt, QVariant,QModelIndex, QAbstractTableModel

import pandas as pd
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from plot import plot_arrangement
from crisscross import CrissCross
from reader import reader
from lrs import SearchStatus


class WidgetGallery(QDialog):
    def __init__(self, parent=None):
        super(WidgetGallery, self).__init__(parent)

        self.search_status = SearchStatus.NONE
        self.canvasWidth = 300
        self.canvasHeight = 300
        self.createCanvas(width=self.canvasWidth, height=self.canvasHeight)
        self.createControls()
        self.createHyperplanes()
        self.createMatrix()
        self.first_basis_found = False
        mainLayout = QGridLayout()
        mainLayout.addLayout(self.controls, 0, 0, 1, 2)
        mainLayout.addWidget(self.canvas, 1, 0)
        mainLayout.addWidget(self.matrixDisplay, 1, 1)
        mainLayout.addWidget(self.hyperplaneDisplay, 2, 1)
        self.setLayout(mainLayout)
        self.setWindowTitle("Lrs")


    def createCanvas(self, width=300, height=300):
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

    def createHyperplanes(self):

        self.hyperplaneDisplay = QGroupBox("Hyperplanes")
        layout = QVBoxLayout()
        self.hyperplaneDisplay.setLayout(layout)

    def createControls(self):
        self.controls = QGroupBox("Controls")

        self.pivotButton = QPushButton("First Basis")
        self.pivotButton.clicked.connect(self.pivot)
        self.pivotButton.setDefault(True)

        self.searchStepButton = QPushButton("First Basis")
        self.searchStepButton.clicked.connect(self.searchButton)
        self.searchStepButton.setDefault(True)

        fileButton = QPushButton("Open File")
        fileButton.setDefault(True)
        fileButton.clicked.connect(self.openFile)

        layout = QVBoxLayout()
        layout.addWidget(self.pivotButton)
        layout.addWidget(self.searchStepButton)
        layout.addWidget(fileButton)
        self.controls = layout

    def plot(self):
        self.figure.clear()
        # create an axis
        ax = self.figure.add_subplot(111)

        if self.first_basis_found:
            point = self.lrs.get_vertex()
        else:
            point = None
        # plot data
        plot_arrangement(self.lrs.hyperplanes, ax=ax, point=point)
        # refresh canvas
        self.canvas.draw()

    def pivot(self):
        if not self.first_basis_found:
            self.first_basis_found = True
            self.lrs.first_basis()
            self.plot()
            self.start_search()

        elif self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            while self.search_status not in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED, SearchStatus.DONE]:
                self.search_status = self.search.__next__()
            self.plot()

    def start_search(self):
        self.lrs.set_objective()
        self.searchStepButton.setText('Search Step')
        self.pivotButton.setText('Pivot')
        self.plot()
        self.search = self.lrs.search()

    def searchButton(self):
        if not self.first_basis_found:
            self.lrs.first_basis()
            self.first_basis_found = True
            self.plot()
            self.start_search()
        elif self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            if self.search_status in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED]:
                self.plot()

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
            self.lrs.augment_matrix_with_objective()
            self.lrs.init_dicts()
            self.plot()


class PandasModel(QAbstractTableModel):
    def __init__(self, df = pd.DataFrame(), parent=None):
        QAbstractTableModel.__init__(self, parent=parent)
        self._df = df

    def headerData(self, section, orientation, role=Qt.DisplayRole):
        if role != Qt.DisplayRole:
            return QVariant()

        if orientation == Qt.Horizontal:
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