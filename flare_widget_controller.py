from __future__ import annotations

import csv
import copy
import os
import subprocess
import sys
import threading
from typing import List

import pandas as pd
import pyqtgraph as pg
import pyqtgraph.ptime as ptime
from PyQt5 import QtCore, QtWidgets, QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *  # noqa: F401,F403 - used by Qt model setup below

from flare_display import Ui_MainWindow
from logging_setup import get_logger

logger = get_logger(__name__)


class MainWindow:

    def __init__(self):

        """MainWindow constructor.
        This widget will be our main window.
        We'll define all the UI components in here.
        """
        file_name = 'flare_list.csv'
        self.main_window = QMainWindow()
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self.main_window)
        self.model = self.ui.model


        self.csv_data = pd.read_csv(file_name)
        self.loadCsv(file_name)
        self.flare_info = self.ui.flare_info
        self.run_widget_button = self.ui.run_widget_button
        self.run_widget_button.clicked.connect(self.run_widget)
        self.flare_table = self.ui.flare_table
        #self.flare_table.cellClicked.connect(self.cell_was_clicked)
        self.flare_table.setSortingEnabled(True)
        self.flare_table.clicked.connect(self.cell_was_clicked)
        # Optional buttons may not exist on this UI; guard connections.
        if hasattr(self.ui, 'update_data_button'):
            self.ui.update_data_button.clicked.connect(self.update_data)
        if hasattr(self.ui, 'run_auto_data_collection_button'):
            self.ui.run_auto_data_collection_button.clicked.connect(self.run_auto_data_collection)

    def run_auto_data_collection(self):
        processThread = threading.Thread(target=self.thread_auto_data_collection)
        processThread.start()

    def are_checkboxes_checked(self, list1):
        list2 = []
        for checkbox in list1:
            if checkbox.checkState() == QtCore.Qt.Checked:
                list2.append(QtCore.Qt.Checked)
            else:
                list2.append(False)
        return list2

    def update_download(self):
        data_index_list = list(pd.read_csv('flare_pickle_links.csv', index_col='index').index)
        logger.info("Updating download state for %d entries", len(data_index_list))
        for num in data_index_list:
            self.checkbox_list[num].setCheckState(QtCore.Qt.Checked)
        list2 = range(0, 15171)
        list3 = set(list2) - set(data_index_list)
        for num in list3:
            self.checkbox_list[num].setCheckState(False)


    def update_data(self):
        #self.update_download()
        validation_list = self.are_checkboxes_checked(self.validation_checkbox_list)
        flagged_list = self.are_checkboxes_checked(self.flagged_checkbox_list)
        pd.to_pickle([validation_list, flagged_list],'validation_flagged_list.pickle')

    def cell_was_clicked(self):
        index = (self.flare_table.selectionModel().currentIndex().row())+1
        logger.debug("Row clicked: %s", index)
        self.index = index
        self.get_metadata()



    def get_metadata(self):
        self.flare_metadata = self.csv_data.iloc[int(self.index)-1]
        self.flare_metadata = self.flare_metadata[1:]
        flare_metadata = copy.deepcopy(self.flare_metadata)
        flare_metadata.drop('rhessi_image_link', inplace=True)
        self.flare_info.setText(str(flare_metadata))



    def show(self):
        self.main_window.show()
        logger.debug("Main window shown")

    def thread_auto_data_collection(self):
        subprocess.call([sys.executable, "automatic_data_collection.py"])  # nosec B603

    def thread_second(self):
        subprocess.call([sys.executable, "flare_display.py"])  # nosec B603

    def run_widget(self):
        event_peak_time1 = self.flare_metadata['GOES Peak Time']
        flare_class = self.flare_metadata['Flare Class']
        file_name = "Flare_Pickle_Files/"+ flare_class+" " + event_peak_time1[:16]+'.pickle'
        flare_links = pd.read_csv('flare_pickle_links.csv')['filename'].to_list()
        if not file_name in flare_links:
            pd.to_pickle(dict(self.flare_metadata), 'current_flare_info.pickle')
        logger.info("Running data collection for selected flare")
        subprocess.call([sys.executable, 'data_collection.py'])  # nosec B603
            flare_pickle_links_dataframe = pd.read_csv('flare_pickle_links.csv', index_col='index')
            flare_pickle_links_dataframe.loc[self.index] = file_name
            flare_pickle_links_dataframe.to_csv('flare_pickle_links.csv')

        pd.to_pickle(file_name,'pickle_file_name.pickle')
        processThread = threading.Thread(target=self.thread_second)  # <- note extra ','
        processThread.start()



    def loadCsv(self, fileName):
        with open(fileName, "r") as fileInput:
            for row in csv.reader(fileInput):
                items = [
                    QtGui.QStandardItem(field)
                    for field in row
                ]
                self.model.appendRow(items)
        self.model.removeRow(0)
        self.model.insertColumn(0)
        self.model.insertColumn(1)

        self.checkbox_list = []
        self.validation_checkbox_list = []
        self.flagged_checkbox_list = []
        validation_list,flagged_list = pd.read_pickle('validation_flagged_list.pickle')
        for column in range(15171):
            item_checked = QStandardItem()
            item_checked.setCheckable(True)
            item_checked.setCheckState(False)
            self.checkbox_list.append(item_checked)
            self.model.setItem(column, 0, item_checked)

        for column in range(15171):
            item_checked = QStandardItem()
            item_checked.setCheckable(True)
            item_checked.setCheckState(validation_list[column])
            self.validation_checkbox_list.append(item_checked)
            self.model.setItem(column, 1, item_checked)

        for column in range(15171):
            item_checked = QStandardItem()
            item_checked.setCheckable(True)
            item_checked.setCheckState(flagged_list[column])
            self.flagged_checkbox_list.append(item_checked)
            self.model.setItem(column, 2, item_checked)


        #self.update_download()


        self.model.setHorizontalHeaderLabels(['Data?','Validated?','Flagged','RHESSI Start Time','RHESSI Peak Time','RHESSI End Time','RHESSI List X center','RHESSI List Y center','Active Region','Flare Class','GOES Start Time','GOES Peak Time','GOES End Time','rhessi_header_x','rhessi_header_y','rhessi_conversion','rhessi_image_link','CME','Speed','Angle'])

if __name__ == '__main__':
    app = QApplication(sys.argv)
    # it's required to save a reference to MainWindow.
    # if it goes out of scope, it will be destroyed.
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec())
