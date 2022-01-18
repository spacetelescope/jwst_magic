"""
This module is for any extra widgets/objects needed for the QT windows to keep mainGUI.py from
getting too messy

"""
import sys

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QComboBox, QListView
from PyQt5.QtGui import QStandardItemModel


class CheckableComboBox(QComboBox):
    def __init__(self, *args, **kwargs):
        super(CheckableComboBox, self).__init__(*args, **kwargs)
        self.setView(QListView(self))
        self.view().pressed.connect(self.handleItemPressed)
        self.setModel(QStandardItemModel(self))

    def handleItemPressed(self, index):
        item = self.model().itemFromIndex(index)
        if item.checkState() == Qt.Checked:
            item.setCheckState(Qt.Unchecked)
        else:
            item.setCheckState(Qt.Checked)

    def checkedItems(self):
        checkedItems = []
        for index in range(self.count()):
            item = self.model().item(index)
            if item.checkState() == Qt.Checked:
                checkedItems.append(item)
        return checkedItems
