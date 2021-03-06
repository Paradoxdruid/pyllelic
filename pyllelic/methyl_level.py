# !/usr/bin/env python3
import os

def methyl_levels(files,file_saver,fig_saver,x_file,mutation):
    import pyllelic
    import PySimpleGUI as sg
    import pandas as pd 
    import numpy as np
    import pysam
    import os #allows interaction with operating system
    import errno
    from skbio import DNA
    from skbio.alignment import StripedSmithWaterman
    import plotly.express as px
    import subprocess #https://docs.python.org/2/library/subprocess.html
    from subprocess import check_output
    from pathlib import Path #filesystem pathways module
    from io import StringIO
    import sys
    import statistics as stat
    import cv2
    
    
    pyllelic.set_up_env_variables(
        base_path= '/home/dylan/research/methyl/tert',
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
        location=file_saver, #Saves the file to a predetermined location. Sequentially different name as to not overwrite data.
        image = fig_saver #Saves the image---automated
        
        
    )
    
    print("Loading. . .")
    
    positions = pyllelic.index_and_fetch(files)
    
#   pyllelic.genome_parsing()
    cell_types = pyllelic.extract_cell_types(files)
    
    df_list = pyllelic.run_quma_and_compile_list_of_df(cell_types, "tester5.xlsx") # to skip quma: , run_quma=False)
    df_list.keys()
    
    means = pyllelic.process_means(df_list, positions, files)
    means
    modes = pyllelic.process_modes(df_list, positions, files)
    modes
    
    diff = pyllelic.find_diffs(means, modes)
    
    #Converts mean/mode difference into an excel file and saves to unique filename
    def excel_file(x_file):
#         data_path = os.getcwd() + f'/pyllelic_raw_data {x_file}'
        diff.to_excel(x_file, index = False)
        
        pd.set_option('display.max_columns', None)
        excel_sheet = pd.read_excel(x_file)
        excel_sheet
        
    excel_file(x_file)

    pyllelic.write_means_modes_diffs(means, modes, diff, "Test5")
    
    #creates the folder that the excel sheets will be stored in.
    
    #OPTIONAL:
    
#     if not os.path.exists("pyllelic_raw_data"):
#         os.mkdir("pyllelic_raw_data")
    
    
    final_data = pyllelic.pd.read_excel(pyllelic.config.base_directory.joinpath("Test5_diff.xlsx"),
        dtype=str,
        index_col=0,
    )
    
    individual_data = pyllelic.return_individual_data(df_list, positions, files)
    
    print('Raw, Individual Data: ' + '\n' + '\n')
    print(individual_data)
    print('\n')
    print('Difference between means vs modes: ' + '\n' + '\n')
    print(diff)
    
    import matplotlib
    matplotlib.use('tkAgg')
    import matplotlib.pyplot as plt

    hist = pyllelic.histogram(individual_data, mutation, '1295089')
    

    #CAUTION!!!!!
#Work in progress--was working fine, now deletes .config file...
#Currently commented out in code.

def config_reset(): 
    py_path = './pyllelic'
    with open(f"{py_path}/config.py", "r") as config:
        lines = f"{py_path}/config.py"
        with open(lines,'r') as f:
            ff = f.read()
            h = "\n".join(ff.split("\n")[:-2])

            with open(lines, 'w') as trunc:
                trunc.write(h) 
                
