# !/usr/bin/env python3

"""Wrapper module for running pyllelic inside gui."""

import pyllelic
from typing import List
import matplotlib

matplotlib.use("tkAgg")  # for pillow image show


def main(
    files: List[str], file_saver: str, fig_saver: str, x_file: str, cell_line: str
) -> None:
    """Wrapper to run pyllelic for the gui frontend.

    Args:
        files (List[str]): list of bam files to process
        file_saver (str): name of file to save
        fig_saver (str): name of figure to save
        x_file (str): name of excel file to save
        cell_line (str): name of cell_line being analyzed
    """
    # Configure Pyllelic
    pyllelic.set_up_env_variables(
        base_path="/home/dylan/research/methyl/tert",
        prom_file="TERT-promoter-genomic-sequence.txt",
        prom_start="1293000",
        prom_end="1296000",
        chrom="5",
        location=file_saver,
        image_location=fig_saver,
        offset=1298163,
    )

    print("Loading. . .")

    # Run core pyllelic processing steps
    positions = pyllelic.index_and_fetch(files)
    #   pyllelic.genome_parsing()  # skip for time during testing
    cell_types = pyllelic.extract_cell_types(files)
    df_list = pyllelic.run_quma_and_compile_list_of_df(
        cell_types, "tester5.xlsx"  # modify appropriately?
    )  # to skip quma: , run_quma=False)
    df_list.keys()
    means = pyllelic.process_means(df_list, positions, files)
    modes = pyllelic.process_modes(df_list, positions, files)
    diff = pyllelic.find_diffs(means, modes)

    diff.to_exel(x_file, index=False)  # Unneeded with command below

    pyllelic.write_means_modes_diffs(means, modes, diff, "Test5")

    # creates the folder that the excel sheets will be stored in.
    # OPTIONAL:
    #     if not os.path.exists("pyllelic_raw_data"):
    #         os.mkdir("pyllelic_raw_data")

    # final_data = pyllelic.pd.read_excel(
    #     pyllelic.config.base_directory.joinpath("Test5_diff.xlsx"),
    #     dtype=str,
    #     index_col=0,
    # )
    individual_data = pyllelic.return_individual_data(df_list, positions, files)

    print("Raw, Individual Data: " + "\n" + "\n")
    print(individual_data)
    print("\n")
    print("Difference between means vs modes: " + "\n" + "\n")
    print(diff)

    pyllelic.histogram(individual_data, cell_line, "1295089")

    # CAUTION!!!!!


# Work in progress--was working fine, now deletes .config file...
# Currently commented out in code.


# def config_reset():
#     py_path = "./pyllelic"
#     with open(f"{py_path}/config.py", "r"):
#         lines = f"{py_path}/config.py"
#         with open(lines, "r") as f:
#             ff = f.read()
#             h = "\n".join(ff.split("\n")[:-2])

#             with open(lines, "w") as trunc:
#                 trunc.write(h)
