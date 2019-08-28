from PyQt5.QtWidgets import (QApplication, QFileDialog, QDialog, QGridLayout, QLabel, QLineEdit,
        QPushButton, QVBoxLayout, QSpacerItem, QSizePolicy)
from PyQt5.QtCore import Qt
from PyQt5 import QtGui
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from plot import plot_arrangement, well_distinguishable_colors
from bland import Bland
from crisscross import CrissCross
from reader import reader
from lrs import SearchStatus, hyperplane_string, PrettyInfo

LabelFont = QtGui.QFont('SansSerif', 12)


class WidgetGallery(QDialog):
    def __init__(self, parent=None):
        super(WidgetGallery, self).__init__(parent)

        self.set_background_color()

        self.first_basis_found = False
        self.lrs = None
        self.pivot_rule = 'CrissCross'
        self.search_status = SearchStatus.NONE

        self.create_layout()

        self.setLayout(self.mainLayout)

        self.setWindowTitle("Lrs")

    def set_background_color(self, color=Qt.white):

        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), color)
        self.setPalette(p)

    def create_layout(self):
        self.create_canvas_layout()
        self.create_controls()
        self.create_hyperplane_display()
        self.create_status_display()
        self.create_status_layout()

        self.create_main_layout()

    def create_canvas_layout(self):
        # a figure instance to plot on
        self.figure = plt.figure()

        # this is the Canvas Widget that displays the `figure`
        # it takes the `figure` instance as a parameter to __init__
        self.canvas = FigureCanvas(self.figure)

        self.canvas_layout = QVBoxLayout()
        self.canvas_layout.addWidget(self.canvas)

    def create_controls(self):

        self.next_pivot_button = QPushButton("Next Pivot")
        self.next_pivot_button.clicked.connect(self.next_pivot)
        self.next_pivot_button.setDefault(True)

        self.search_step_button = QPushButton("Search Step")
        self.search_step_button.clicked.connect(self.search_step)
        self.search_step_button.setDefault(True)

        self.first_basis_button = QPushButton("First Basis")
        self.first_basis_button.setDisabled(True)
        self.first_basis_button.clicked.connect(self.first_basis)
        self.first_basis_button.setDefault(True)

        self.create_pivot_rule_button()

        fileButton = QPushButton("Open File")
        fileButton.setDefault(True)
        fileButton.clicked.connect(self.open_file)

        layout = QGridLayout()
        layout.addWidget(self.first_basis_button, 1, 0)
        layout.addWidget(fileButton, 2, 0)
        layout.addWidget(self.pivot_rule_button, 3, 0)

        self.controls = layout

    def create_pivot_rule_button(self):
        self.pivot_rule_button = QPushButton("Change Pivot Rule to Bland")
        self.pivot_rule_button.clicked.connect(self.change_pivot_rule)
        self.pivot_rule_button.setDefault(True)

    def create_coordinate_controls(self):
        coordinate_controls = QGridLayout()

        # min x text field
        min_x_label = QLabel()
        min_x_label.setFont(LabelFont)
        min_x_label.setText('x min:')

        self.min_x_box = QLineEdit(self)
        self.min_x_box.setText('-1')

        # max x text field
        max_x_label = QLabel()
        max_x_label.setFont(LabelFont)
        max_x_label.setText('x max:')

        self.max_x_box = QLineEdit(self)
        self.max_x_box.setText('3')

        # min y text field
        min_y_label = QLabel()
        min_y_label.setFont(LabelFont)
        min_y_label.setText('y min:')

        self.min_y_box = QLineEdit(self)
        self.min_y_box.setText('-1')

        # max y text field
        max_y_label = QLabel()
        max_y_label.setFont(LabelFont)
        max_y_label.setText('y max:')

        self.max_y_box = QLineEdit(self)
        self.max_y_box.setText('3')

        # button for updating bounds
        self.set_coordinates = QPushButton("Set Coordinates")
        self.set_coordinates.clicked.connect(self.plot)
        self.set_coordinates.setDefault(True)

        # adding everything

        coordinate_controls.addWidget(min_x_label, 0, 0)
        coordinate_controls.addWidget(self.min_x_box, 0, 1)

        coordinate_controls.addWidget(max_x_label, 0, 2)
        coordinate_controls.addWidget(self.max_x_box, 0, 3)

        coordinate_controls.addWidget(min_y_label, 1, 0)
        coordinate_controls.addWidget(self.min_y_box, 1, 1)

        coordinate_controls.addWidget(max_y_label, 1, 2)
        coordinate_controls.addWidget(self.max_y_box, 1, 3)

        button_help_layout = QVBoxLayout()
        button_help_layout.addWidget(self.set_coordinates)
        coordinate_controls.addLayout(button_help_layout, 2, 0, 2, 4)
        self.mainLayout.addLayout(coordinate_controls, 0, 1, 1, 1)

    def create_main_layout(self):
        self.mainLayout = QGridLayout()
        self.mainLayout.addLayout(self.controls, 0, 0, 1, 1)
        self.mainLayout.addLayout(self.canvas_layout, 1, 0, 3, 1)
        self.mainLayout.addLayout(self.status_layout, 1, 1, 3, 1)

    def create_status_layout(self):
        self.status_layout = QGridLayout()
        self.status_layout.addItem(QSpacerItem(0, 50, QSizePolicy.Minimum, QSizePolicy.Minimum), 0, 0)
        self.status_layout.addWidget(self.matrixDisplay, 1, 0)

        self.status_layout.addWidget(self.status_display, 2, 0)

        def pivotLayout():
            pivotHelperLayout = QGridLayout()
            pivotHelperLayout.addWidget(self.pivot_label_i, 0, 0)
            pivotHelperLayout.addWidget(self.pivot_box_i, 0, 1)

            pivotHelperLayout.addWidget(self.pivot_label_j, 1, 0)
            pivotHelperLayout.addWidget(self.pivot_box_j, 1, 1)

            pivotHelperLayout.addWidget(self.pivot_button, 1, 2)
            pivotHelperLayout.addWidget(self.select_pivot_button, 0, 2)
            return pivotHelperLayout

        self.status_layout.addLayout(pivotLayout(), 3, 0, 1, 1)
        self.status_layout.addLayout(self.hyperplaneDisplay, 4, 0, 1, 1)

    def create_hyperplane_display(self):
        self.hyperplaneDisplay = QVBoxLayout()
        self.hyperplaneDisplay.setContentsMargins(0, 0, 0, 0)

    def create_status_display(self):
        self.status_display = QLabel()
        self.status_display.setFont(LabelFont)

        self.pivot_label_i = QLabel()
        self.pivot_label_i.setText('i:')
        self.pivot_label_i.setFont(LabelFont)
        self.pivot_box_i = QLineEdit(self)

        self.pivot_label_j = QLabel()
        self.pivot_label_j.setText('j:')
        self.pivot_label_j.setFont(LabelFont)
        self.pivot_box_j = QLineEdit(self)

        self.pivot_button = QPushButton("Pivot")
        self.pivot_button.clicked.connect(self.pivot)
        self.pivot_button.setDefault(True)

        self.select_pivot_button = QPushButton("Select Pivot")
        self.select_pivot_button.clicked.connect(self.select_pivot)
        self.select_pivot_button.setDefault(True)

        self.matrixDisplay = QLabel()
        self.matrixDisplay.setFont(LabelFont)
        self.matrixDisplay.setAlignment(Qt.AlignVCenter)

    def change_pivot_rule(self):
        if self.pivot_rule == 'CrissCross':
            self.pivot_rule_button.setText('Change pivot rule to CrissCross')
            self.pivot_rule = 'Bland'
        elif self.pivot_rule == 'Bland':
            self.pivot_rule_button.setText('Change pivot rule to CrissCross')
            self.pivot_rule = 'CrissCross'

    def pivot(self):
        self.lrs.i = int(self.pivot_box_i.text())
        self.lrs.j = int(self.pivot_box_j.text())
        self.lrs.pivot()
        self.update()

    def select_pivot(self):
        i, j = self.lrs.select_pivot()
        self.pivot_box_i.setText(str(i))
        self.pivot_box_j.setText(str(j))

    def clear_layout(self, layout):
        while layout.count():
            child = layout.takeAt(0)
            if child.widget():
                child.widget().deleteLater()

    def write_hyperplanes(self):
        self.clear_layout(self.hyperplaneDisplay)
        colors = well_distinguishable_colors(
            len(self.lrs.hyperplanes) + len(self.lrs.bounding_box) + 1
        )

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

        for i, hyperplane in enumerate(self.lrs.hyperplanes + self.lrs.bounding_box):
            hyperplane_label = QLabel()
            hyperplane_label.setAlignment(Qt.AlignLeft)
            if vars[i].box_variable:
                var_text = 'Var: {}; Box: '.format(vars[i])
            else:
                var_text = 'Var: {}; Hyperplane: '.format(vars[i])
            hyperplane_label.setText(var_text + hyperplane_string(hyperplane))
            hyperplane_label.setStyleSheet(
                get_label_stylesheet(i)
            )
            hyperplane_label.setFont(LabelFont)
            self.hyperplaneDisplay.addWidget(hyperplane_label)
        self.hyperplaneDisplay.addStretch(1)

    def write_status(self):
        status_string = 'Search Status: {} \n'.format(self.search_status.value)
        self.status_display.setText(status_string)
        self.pivot_box_i.setText(str(self.lrs.i))
        self.pivot_box_j.setText(str(self.lrs.j))

    def update(self, update_hyperplanes=False):
        self.plot()
        if update_hyperplanes:
            self.write_hyperplanes()
        self.matrixDisplay.setText(PrettyInfo(self.lrs))
        self.write_status()

    def plot(self):
        x_min = float(self.min_x_box.text())
        x_max = float(self.max_x_box.text())

        y_min = float(self.min_y_box.text())
        y_max = float(self.max_y_box.text())

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
            self.lrs.bounding_box,
            ax=ax, point=point, x_limits=(x_min, x_max),  y_limits=(y_min, y_max)
        )
        # refresh canvas
        self.canvas.draw()

    def next_pivot(self):

        if self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            while self.search_status not in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED, SearchStatus.DONE]:
                self.search_status = self.search.__next__()
            self.update()

    def reset_controls(self):
        self.controls.addWidget(self.first_basis_button, 1, 0)
        self.create_pivot_rule_button()
        self.controls.addWidget(self.pivot_rule_button, 3, 0)
        self.next_pivot_button.setParent(None)
        self.search_step_button.setParent(None)

    def search_step(self):
        if self.search_status != SearchStatus.DONE:
            self.search_status = self.search.__next__()
            if self.search_status in [SearchStatus.NEWBASIS, SearchStatus.BACKTRACKED]:
                self.update()
            self.write_status()

    def first_basis(self):
        self.lrs.first_basis()
        self.first_basis_found = True
        self.update(update_hyperplanes=True)
        self.first_basis_button.setParent(None)
        self.controls.addWidget(self.next_pivot_button, 0, 0)
        self.controls.addWidget(self.search_step_button, 1, 0)
        self.start_search()

    def start_search(self):
        self.lrs.set_objective()
        self.update(update_hyperplanes=True)
        self.search = self.lrs.search()

    def open_file(self):
        options = QFileDialog.Options()
        options |= QFileDialog.DontUseNativeDialog
        fileName, _ = QFileDialog.getOpenFileName(self, "QFileDialog.getOpenFileName()", "",
                                                  "All Files (*);;Ine Files (*.ine)",
                                                  options=options)
        if fileName:
            if self.lrs is not None:
                self.reset_controls()
                self.first_basis_found = False
            else:
                self.create_coordinate_controls()
            self.pivot_rule_button.deleteLater()
            self.search_status = SearchStatus.NONE
            if self.pivot_rule == 'CrissCross':
                self.lrs = CrissCross(*reader(fileName))
            elif self.pivot_rule == 'Bland':
                self.lrs = Bland(*reader(fileName))
            self.lrs.augment_matrix_with_objective()
            self.lrs.init_dicts()
            if len(self.lrs.bounding_box) > 0:
                self.lrs.add_box_constraints(self.lrs.bounding_box)
            self.first_basis_button.setDisabled(False)
            self.update(update_hyperplanes=True)

import sys
app = QApplication(sys.argv)
gallery = WidgetGallery()
gallery.show()
sys.exit(app.exec_())