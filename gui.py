from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit, QFileDialog,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget, )
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from PyQt5.QtCore import QRect
from crisscross import CrissCross
from reader import reader
from plot import plot_arrangement


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

        pivotButton = QPushButton("Pivot")
        pivotButton.setDefault(True)

        searchStepButton = QPushButton("Search Step")
        searchStepButton.setDefault(True)

        fileButton = QPushButton("Open File")
        fileButton.setDefault(True)
        fileButton.clicked.connect(self.openFile)

        layout = QVBoxLayout()
        layout.addWidget(pivotButton)
        layout.addWidget(searchStepButton)
        layout.addWidget(fileButton)
        self.controls = layout

    def plot(self):
        self.figure.clear()
        # create an axis
        ax = self.figure.add_subplot(111)
        # plot data
        plot_arrangement(self.lrs.hyperplanes, ax=ax)
        # refresh canvas
        self.canvas.draw()


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
            self.plot()

    def drawHyperplanes(self):
        pass


def computeCoordinate(hyperplane_array, x):
    # hyperplane array [b, a0, a1] corresponds to b + ao*x +a1*y = 0
    # Therefore y = (-b -a0*x) / a1
    if hyperplane_array[2] == 0:
        return ValueError
    else:
        return (-hyperplane_array[0] - hyperplane_array[1] * x) / hyperplane_array[2]


import sys
app = QApplication(sys.argv)
gallery = WidgetGallery()
gallery.show()
sys.exit(app.exec_())