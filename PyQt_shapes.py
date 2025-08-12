from __future__ import annotations

from PyQt5 import QtCore, QtWidgets, QtGui


class RectItem(QtWidgets.QGraphicsRectItem):
    def paint(self, painter, option, widget=None):
        super().paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        painter.restore()


class EclipseItem(QtWidgets.QGraphicsEllipseItem):
    def paint(self, painter, option, widget=None):
        super().paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        painter.restore()