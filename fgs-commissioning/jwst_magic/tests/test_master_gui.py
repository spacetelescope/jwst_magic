"""Collection of unit tests to verify the correct function of the
master GUI module.

Authors
-------
    - Lauren Chambers

Use
---
    ::
        pytest test_master_gui.py
"""
import sys

from PyQt5 import QtCore
from PyQt5.QtWidgets import QApplication, QDialogButtonBox
# import qtbot

from jwst_magic import masterGUI

def test_no_input_image(qtbot):
    # Run GUI
    app  = QApplication(sys.argv)
    mg = masterGUI.MasterGui(app=app)
    mg.show()

    # # Schedule press of "Run" button
    # QtCore.QTimer.singleShot(0, mg.pushButton_run.clicked)
    qtbot.mouseClick(mg.pushButton_run, QtCore.Qt.LeftButton)
    assert hasattr(mg, 'no_inputImage_dialog_box')

    # Schedule press of "Ok" button on no_input_image_dialog
    print()
    ok_button = mg.no_input_im_dialog_box.buttonBox.button(QDialogButtonBox.Ok)
    # QtCore.QTimer.singleShot(0, .clicked)
    qtbot.mouseClick(ok_button, QtCore.Qt.LeftButton)

    # app.exec()

    # Make sure the dialog box was spawned


