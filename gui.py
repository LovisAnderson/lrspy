from PyQt5.QtWidgets import (QApplication, QCheckBox, QComboBox, QDateTimeEdit, QFileDialog,
        QDial, QDialog, QGridLayout, QGroupBox, QHBoxLayout, QLabel, QLineEdit,
        QProgressBar, QPushButton, QRadioButton, QScrollBar, QSizePolicy,
        QSlider, QSpinBox, QStyleFactory, QTableWidget, QTabWidget, QTextEdit,
        QVBoxLayout, QWidget, )
from PyQt5.QtCore import QRect, Qt, QVariant,QModelIndex, QAbstractTableModel
from PyQt5 import QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from plot import plot_arrangement, well_distinguishable_colors
from crisscross import CrissCross
from reader import reader
from lrs import SearchStatus, hyperplane_string

LabelFont = QtGui.QFont('SansSerif', 12)


class WidgetGallery(QDialog):
    def __init__(self, parent=None):
        super(WidgetGallery, self).__init__(parent)

        self.set_background_color()
        self.search_status = SearchStatus.NONE
        self.create_canvas()
        self.create_controls()
        self.create_hyperplane_display()
        self.create_matrix_dispay()
        self.create_status_display()
        self.first_basis_found = False
        mainLayout = QGridLayout()
        mainLayout.addLayout(self.controls, 0, 0, 1, 2)
        self.display_layout = QVBoxLayout()
        self.display_layout.addWidget(self.canvas)
        mainLayout.addLayout(self.display_layout, 1, 0, 3, 1)
        mainLayout.addWidget(self.matrixDisplay, 1, 1)
        mainLayout.addLayout(self.hyperplaneDisplay, 2, 1, 2, 1)
        mainLayout.addWidget(self.status_display, 3, 1)
        self.setLayout(mainLayout)
        self.setWindowTitle("Lrs")

    def set_background_color(self, color=Qt.white):

        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), color)
        self.setPalette(p)

    def create_canvas(self):
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

    def create_controls(self):

        self.pivot_button = QPushButton("First Basis")
        self.pivot_button.clicked.connect(self.pivot_step)
        self.pivot_button.setDefault(True)

        self.search_step_button = QPushButton("First Basis")
        self.search_step_button.clicked.connect(self.search_step)
        self.search_step_button.setDefault(True)

        fileButton = QPushButton("Open File")
        fileButton.setDefault(True)
        fileButton.clicked.connect(self.open_file)

        layout = QVBoxLayout()
        layout.addWidget(self.pivot_button)
        layout.addWidget(self.search_step_button)
        layout.addWidget(fileButton)
        self.controls = layout

    def create_hyperplane_display(self):
        self.hyperplaneDisplay = QVBoxLayout()
        self.hyperplaneDisplay.setContentsMargins(0, 0, 0, 0)

    def create_matrix_dispay(self):
        self.matrixDisplay = QLabel()

        self.matrixDisplay.setFont(LabelFont)

    def create_status_display(self):
        self.status_display = QLabel()
        self.matrixDisplay.setFont(LabelFont)

    def clear_layout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def write_hyperplanes(self):
        self.clear_layout(self.hyperplaneDisplay)
        colors = well_distinguishable_colors(len(self.lrs.hyperplanes) + 1)

        def get_label_stylesheet(i):
            return "QLabel {{ color : rgb({}, {}, {}); }}".format(*[int(c*255) for c in colors[i]])
        vars = self.lrs.hyperplane_variables()

        for i, hyperplane in enumerate(self.lrs.hyperplanes):
            hyperplane_label = QLabel()
            hyperplane_label.setAlignment(Qt.AlignCenter)
            var_text = 'Var: {}; Hyperplane: '.format(vars[i])
            hyperplane_label.setText(var_text + hyperplane_string(hyperplane))
            hyperplane_label.setStyleSheet(
                get_label_stylesheet(i)
            )
            hyperplane_label.setFont(LabelFont)
            self.hyperplaneDisplay.addWidget(hyperplane_label)
        self.hyperplaneDisplay.addStretch(1)

    def update(self, update_hyperplanes=False):
        self.plot()
        if update_hyperplanes:
            self.write_hyperplanes()
        self.matrixDisplay.setText(self.lrs.info_string())

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

    def pivot_step(self):
        if not self.first_basis_found:
            self.first_basis_found = True
            self.lrs.first_basis()
            self.start_search()

        elif self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            while self.search_status not in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED, SearchStatus.DONE]:
                self.search_status = self.search.__next__()
            self.update()

    def start_search(self):
        self.lrs.set_objective()
        self.search_step_button.setText('Search Step')
        self.pivot_button.setText('Pivot')
        self.update(update_hyperplanes=True)
        self.search = self.lrs.search()

    def search_step(self):
        if not self.first_basis_found:
            self.lrs.first_basis()
            self.first_basis_found = True
            self.plot()
            self.start_search()
        elif self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            if self.search_status in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED]:
                self.update()

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Ine Files (*.ine)",
                                                  options=options)
        if fileName:
            self.lrs = CrissCross(*reader(fileName))
            self.lrs.augment_matrix_with_objective()
            self.lrs.init_dicts()
            self.update(update_hyperplanes=True)


import sys
app = QApplication(sys.argv)
gallery = WidgetGallery()
gallery.show()
sys.exit(app.exec_())