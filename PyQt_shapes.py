from PyQt5 import QtCore, QtWidgets, QtGui

class RectItem(QtWidgets.QGraphicsRectItem):
    def paint(self, painter, option, widget=None):
        super(RectItem, self).paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        # painter.drawEllipse(option.rect)
        painter.restore()


class EclipseItem(QtWidgets.QGraphicsEllipseItem):
    def paint(self, painter, option, widget=None):
        super(EclipseItem, self).paint(painter, option, widget)
        painter.save()
        painter.setRenderHint(QtGui.QPainter.Antialiasing)
        painter.setBrush(QtCore.Qt.red)
        painter.setPen(QtCore.Qt.green)
        # painter.drawEllipse(option.rect)
        painter.restore()