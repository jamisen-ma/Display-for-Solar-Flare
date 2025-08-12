from __future__ import annotations

from typing import Tuple

import matplotlib.ticker as mticker
from PyQt5 import QtWidgets
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width: float = 2, height: float = 4, dpi: int = 100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        super().__init__(self.fig)


class MathTextSciFormatter(mticker.Formatter):
    def __init__(self, fmt: str = "%1.2e"):
        self.fmt = fmt

    def __call__(self, x, pos=None):
        s = self.fmt % x
        decimal_point = '.'
        positive_sign = '+'
        significand, exp = s.split('e')
        significand = significand.rstrip(decimal_point)
        sign = exp[0].replace(positive_sign, '')
        exponent = exp[1:].lstrip('0')
        if exponent:
            exponent = f"10^{{{sign}{exponent}}}"
        s = rf"{significand}{{\times}}{exponent}" if (significand and exponent) else rf"{significand}{exponent}"
        return f"${s}$"


class MatplotlibWidget(QtWidgets.QWidget):
    def __init__(self, size: Tuple[float, float] = (2, 2), dpi: int = 100):
        super().__init__()
        self.fig = Figure(size, dpi=dpi)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self)
        # toolbar can be added here if desired

        self.vbox = QtWidgets.QVBoxLayout()
        self.vbox.addWidget(self.canvas)
        self.setLayout(self.vbox)

    def getFigure(self) -> Figure:
        return self.fig

    def draw(self) -> None:
        self.canvas.draw()