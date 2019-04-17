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
        self.create_coordinate_controls()
        self.create_hyperplane_display()
        self.create_status_display()
        self.first_basis_found = False
        mainLayout = QGridLayout()
        mainLayout.addLayout(self.controls, 0, 0, 1, 3)
        self.display_layout = QVBoxLayout()
        self.display_layout.addWidget(self.canvas)
        status_layout = QGridLayout()
        status_layout.addWidget(self.matrixDisplay, 0, 0)
        status_layout.addWidget(self.status_display, 1, 0)
        status_layout.addLayout(self.hyperplaneDisplay, 2, 0, 1, 1)
        mainLayout.addLayout(self.coordinate_controls, 4, 0, 1, 3)
        mainLayout.addLayout(self.display_layout, 1, 0, 3, 1)
        mainLayout.addLayout(status_layout, 1, 1, 3, 1)
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

    def create_coordinate_controls(self):
        self.coordinate_controls = QGridLayout()

        # min x text field
        self.min_x_label = QLabel()
        self.min_x_label.setFont(LabelFont)
        self.min_x_label.setText('x min:')

        self.min_x_box = QLineEdit(self)
        self.min_x_box.setText('0')

        # max x text field
        self.max_x_label = QLabel()
        self.max_x_label.setFont(LabelFont)
        self.max_x_label.setText('x max:')

        self.max_x_box = QLineEdit(self)
        self.max_x_box.setText('10')

        # min y text field
        self.min_y_label = QLabel()
        self.min_y_label.setFont(LabelFont)
        self.min_y_label.setText('y min:')

        self.min_y_box = QLineEdit(self)
        self.min_y_box.setText('0')

        # max y text field
        self.max_y_label = QLabel()
        self.max_y_label.setFont(LabelFont)
        self.max_y_label.setText('y max:')

        self.max_y_box = QLineEdit(self)
        self.max_y_box.setText('10')

        # button for updating bounds
        self.set_coordinates = QPushButton("Set Coordinates")
        self.set_coordinates.clicked.connect(self.plot)
        self.set_coordinates.setDefault(True)

        # adding everything

        self.coordinate_controls.addWidget(self.min_x_label, 0, 0)
        self.coordinate_controls.addWidget(self.min_x_box, 0, 1)

        self.coordinate_controls.addWidget(self.max_x_label, 0, 2)
        self.coordinate_controls.addWidget(self.max_x_box, 0, 3)

        self.coordinate_controls.addWidget(self.min_y_label, 1, 0)
        self.coordinate_controls.addWidget(self.min_y_box, 1, 1)

        self.coordinate_controls.addWidget(self.max_y_label, 1, 2)
        self.coordinate_controls.addWidget(self.max_y_box, 1, 3)

        button_help_layout = QVBoxLayout()
        button_help_layout.addWidget(self.set_coordinates)
        self.coordinate_controls.addLayout(button_help_layout, 2, 0, 2, 4)

    def update_display_bounds(self):

        self.plot(min_x=x_min, max_x=x_max, min_y=y_min, max_y=y_max)

    def create_hyperplane_display(self):
        self.hyperplaneDisplay = QVBoxLayout()
        self.hyperplaneDisplay.setContentsMargins(0, 0, 0, 0)

    def create_status_display(self):
        self.status_display = QLabel()
        self.status_display.setFont(LabelFont)
        self.matrixDisplay = QLabel()
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
            color = 'rgb({}, {}, {})'.format(
                *[int(c*255) for c in colors[i]]
            )
            text_align = 'left'
            stylesheet = "QLabel {{ color : {}; text-align: {};}}".format(
                color,
                text_align
            )
            return stylesheet

        vars = self.lrs.hyperplane_variables()

        for i, hyperplane in enumerate(self.lrs.hyperplanes):
            hyperplane_label = QLabel()
            hyperplane_label.setAlignment(Qt.AlignLeft)
            var_text = 'Var: {}; Hyperplane: '.format(vars[i])
            hyperplane_label.setText(var_text + hyperplane_string(hyperplane))
            hyperplane_label.setStyleSheet(
                get_label_stylesheet(i)
            )
            hyperplane_label.setFont(LabelFont)
            self.hyperplaneDisplay.addWidget(hyperplane_label)
        self.hyperplaneDisplay.addStretch(1)

    def write_status(self):
        status_string = self.search_status.value + '\n'
        status_string += 'i: {} \n'.format(self.lrs.i)
        status_string += 'j: {} \n'.format(self.lrs.j)
        self.status_display.setText(status_string)

    def update(self, update_hyperplanes=False):
        self.plot()
        if update_hyperplanes:
            self.write_hyperplanes()
        self.matrixDisplay.setText(self.lrs.info_string())
        self.write_status()

    def plot(self):
        x_min = int(self.min_x_box.text())
        x_max = int(self.max_x_box.text())

        y_min = int(self.min_y_box.text())
        y_max = int(self.max_y_box.text())

        self.figure.clear()

        # create an axis
        ax = self.figure.add_subplot(111)

        if self.first_basis_found:
            point = self.lrs.get_vertex()
        else:
            point = None
        # plot data
        plot_arrangement(
            self.lrs.hyperplanes,
            ax=ax, point=point, x_limits=(x_min, x_max),  y_limits=(y_min, y_max)
        )
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
            self.write_status()

    def start_search(self):
        self.lrs.set_objective()
        self.search_step_button.setText('Search Step')
        self.pivot_button.setText('Pivot')
        self.update(update_hyperplanes=True)
        self.search = self.lrs.search()

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Ine Files (*.ine)",
                                                  options=options)
        if fileName:
            self.search_status = SearchStatus.NONE
            self.lrs = CrissCross(*reader(fileName))
            self.lrs.augment_matrix_with_objective()
            self.lrs.init_dicts()
            self.update(update_hyperplanes=True)


import sys
app = QApplication(sys.argv)
gallery = WidgetGallery()
gallery.show()
sys.exit(app.exec_())