#!/usr/bin/env python3

"""GUIllelic: a graphical front-end interface for pyllelic."""

import PySimpleGUI as sg
import methyl_level
import os
import threading
import time

# Gui Interface:
sg.theme("Dark Blue 13")

layout = [
    [sg.Text("Pyllelic GUI:")],
    [
        sg.Text("Search for file: ", size=(14, 2)),
        sg.InputText(key="_FILES_"),
        sg.FileBrowse(),
    ],
    [
        sg.Button("Add"),
    ],
    [sg.Button("Submit")],
    [sg.Text("Your selected .BAM files: "), sg.Text(size=(40, 2), key="-OUTPUT-")],
]

window = sg.Window("Pyllelic 1.00", layout)

file_set = []
mutant_set = []

# Will be True until closing windows.
while True:
    event, values = window.read()

    if event == sg.WIN_CLOSED or event == "Cancel":
        window.close()
        break

    # Can add multiple cell lines at once.
    # But needs work before running multiples without errors.

    if event == "Add":

        # Divides each of the chosen file names and obtains only the mutant cell names.

        file_set.append(values["_FILES_"].split("/")[-1])
        mut = values["_FILES_"].split("/")[-1].split("_")[1]
        mutant_set.append(mut)
        window["-OUTPUT-"].update(mutant_set)

    if event == "Submit":

        # File browser to choose base location of both tables
        # (not ready yet) and figures.

        file_save = sg.popup_get_folder(
            "Please select location for your raw data to be saved: "
        )

        from pathlib import Path  # filesystem pathways module

        # can change to different folder
        py_path = Path.cwd() / "pyllelic"

        with open(f"{py_path}/config.py", "a") as configure:

            # Adds a new entry to the configure file--added code to pyllelic
            # code to communicate with guillelic
            configure.write("\n" + f"file_save: str ='{file_save}'")

        print(f"Your file is saved in: {file_save}")

        # Takes the file location defined above and appends mutant cell type.
        # Takes the file location defined above and appends mutant cell type.

        image_name = f"{file_save}/{mut}.png"

        """To avoid overwriting multiple files saved with the same cell type name,
           code autogenerates new file name for each run."""

        if not os.path.isfile(image_name):
            fig_save = f"{file_save}/{mut}.png"
            excel_file = f"{file_save}/{mut}.xlsx"

            from pathlib import Path  # filesystem pathways module

            py_path = Path.cwd() / "pyllelic"

            with open(f"{py_path}/test.py", "a") as configure:
                configure.write("\n" + f"fig_save: str ='{fig_save}'")

                # This is where raw data can be found--new file name for each run.
                print(f"Image is saved under: {fig_save}")

        # If file name already exists, below codes for the same mutant type,
        # just one sequential number above the last file.
        else:
            n = 1
            while os.path.isfile(image_name):
                n = n + 1
                image_name = f"{file_save}/{mut}_{n}.png"
                excel_file = f"{file_save}/{mut}_{n}.xlsx"

                if not os.path.isfile(image_name):
                    fig_save = image_name

                    from pathlib import Path  # filesystem pathways module

                    py_path = Path.cwd() / "pyllelic"

                    with open(f"{py_path}/test.py", "a") as configure:
                        configure.write("\n" + f"fig_save: str ='{fig_save}'")
                        print(f"Image is saved under: {fig_save}")

        # Window that contains output data:

        layout = [
            [sg.Output(size=(100, 20), key="-OUTPUT-")],
            [sg.Button("Run Analysis")],
            [sg.Button("Exit")],
        ]

        window = sg.Window("Pyllelic", layout)

        def methylation():
            print("Analysis Running... This can take a few minutes. \n \n ")

            time.sleep(2)
            methyl_level.main(file_set, file_save, fig_save, excel_file, mut)
            time.sleep(2)
            window.refresh()

            print(f"Figure saved in the following path: {fig_save}")

            return ()

        time.sleep(2)

        # Uploads graph to desktop and opens it in a new tab.
        def pic():
            from PIL import Image

            im = Image.open(f"{fig_save}")

            return im.show()

        # The following is for the output window:
        while True:
            # window read outputs data to user-interface
            event, values = window.read()

            # This is the actual methyl analysis
            if event == "Run Analysis":
                x = threading.Thread(target=methylation())

                x.start()
                # Tells CPU to give 2 second rest so GUI can continue to run.
                time.sleep(2)

                y = threading.Thread(target=pic)
                y.start()

                time.sleep(2)

                break

    if event == "Exit" or sg.WIN_CLOSED:
        window.close()


#    methyl_level.config_reset()

""" This needs some work. It was working fine until recently it started deleting
    the entire .config file for no apparent reason program still runs without it...
    but prefer to clean it up.
    This will delete the configure entries so each time it runs will start from
    the original .config file."""


#             if event == 'Plot':

#                 window2 = sg.Window('Pylleic 1.00--Mean vs. Mode Figure:', layout2)
#                 event, values = window2.read()
#                 sg.popup(pyllelic.histogram(individual_data, 'SW1710', '1295089'))
#                 if event in (sg.WIN_CLOSED, 'Cancel'):
#                     break

#             if event == 'Data_Table':
#                 pyllelic.histogram(individual_data, 'SW1710', '1295089')
#                 window2 = sg.Window('Pylleic 1.00--Mean vs. Mode Figure:', layout2)
#                 event, values = window2.read()
#                 sg.popup(print(individual_data))
#                 if event in (sg.WIN_CLOSED, 'Cancel'):
#                     break

window.close()

#           elif event == 'Popup':
#             sg.popup('Yes, your application is still running')while True:


# change the "output" element to be the value of "input" element
#         window['-OUTPUT-'].update(values['-IN-'])

# set up your disk location:
# base_path should be the directory we'll do our work in
# make a sub-directory under base_path with a folder named "test"
# and put the .bam and .bai files in "test"
