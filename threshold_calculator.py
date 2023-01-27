"""
Threshold Calculator
====================

This scripts takes one or several fasta **SILVA alignment** files and
 calculates the percentage of positions matched out of all non-gap
 positions with valid base value.

This script would be given names of fasta files as arguments and then
 it reads the files create list of uniq sequence headers and makes a
 combination out of them.

The list of combination would be...

"""
import csv
import logging
import argparse
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from collections import defaultdict as ddict
from concurrent.futures import ProcessPoolExecutor, as_completed

# New imports
from sys import version_info
if version_info[0] < 3:
    from pathlib2 import Path  # pip2 install pathlib2
else:
    from pathlib import Path


"""
Logger
======
"""
db_feeder_dir = Path(__file__).parent
logger = logging.getLogger("Threshold Calculator")
# set the logging level
logger.setLevel(logging.DEBUG)
# create a console and file handler
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
fh = logging.FileHandler(db_feeder_dir.joinpath("log_thershold_calculator.txt"))
fh.setLevel(logging.WARNING)
# create a formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
# set the formatter to the console handler
ch.setFormatter(formatter)
fh.setFormatter(formatter)
# add the console handler to the logger
logger.addHandler(ch)
logger.addHandler(fh)


# Setting global paths
global graph_path, excel_path, seq_tax_map_path
graph_path = Path(__file__).parent.joinpath("similarity_heatmap.png")
excel_path = Path(__file__).parent.joinpath("similarity_matrix.xlsx")
seq_tax_map_path = Path(__file__).parent.joinpath("seqid_taxonomy_map.tsv")


def defval_str():
    return ""


def sort_tax_dict(tax_dict: tuple, sort_level: int):
    """
    Tries to get the sort_level item of second element of
     tax_dict items which is a list of taxonomy levels in
     order to sort it.
    """
    for item_tuple in tax_dict.items():
        target_tax = ""
        try:
            if sort_level not in list(range(6)):
                raise ValueError("sort_level must be between 0 to 5")
            sta_ind = sort_level - 1
            end_ind = sort_level
            tax_list = item_tuple[1]
            if isinstance(tax_list, list):
                target_slice = tax_list[sta_ind:end_ind]
                target_tax = target_slice[0] if target_slice else ""
        except Exception as exc:
            logger.warning(f"{exc}. Invalid taxonomy level for {str(item_tuple)}")

    return target_tax


def parse_fasta_headers(file_list):
    headers = set()
    for file in file_list:
        with open(file, 'r') as f:
            for line in f:
                if line.startswith(">"):
                    headers.add(line[1:].strip())
    return headers


def create_pairwise_combinations(headers):
    return list(itertools.combinations(headers, 2))


def get_header_seq_dict(file_path,
                        head_sep=";",
                        ind_from=0,
                        ind_to=1,
                        tax_sep=";",
                        ):
    """
    Parses fasta file to dictionary.

    **Arguments:**
    =========

    - file_path: the path to fasta file
    - head_sep: the seperator to cut the header for shorter header creationg
    - ind_from: the index to slice from the cut header
    - ind_to: the index to slice to the cut header
    - tax_sep: seperator of taxonomy levels
    """
    fasta_dict = ddict(defval_str)
    tax_dict = ddict(defval_str)
    file_path = Path(file_path).resolve()
    logger.info(f"Reading and parsing file: {file_path.name}")
    # should_sort = sort_level in list(range(6))
    with open(file_path, "r") as infile:
        line = infile.readline().strip()
        curr_header = ""
        while line:
            if line.startswith(">"):
                line_list = line[1:].split(head_sep)  # [1:] excluding '>'
                header = line_list[ind_from:ind_to][0]
                curr_header = header
                # Creating tax_dict
                tax_list = str(str(line.split("tax=")[1]).split(";;")[0]).split(tax_sep)
                tax_dict[header] = tax_list

                line = infile.readline().strip()
                continue
            if curr_header:
                fasta_dict[curr_header] += line
            line = infile.readline().strip()

    return fasta_dict, tax_dict


def parse_fastas_to_dict(file_list, sort_level=-1):
    """
    Takes a list of files and transform the fasta to dictionary
    **Arguments:**
    ==============
    - file_list: a list of fasta files
    - sort_level: sort key for headers. Any value rathar than 0-5
         means no sorting
    """
    sort_level_dict = {
        0: "kingdom",
        1: "phylum",
        2: "class",
        3: "order",
        4: "family",
        5: "genus",
        6: "species"
    }
    fastas_dict = {}
    taxs_dict = {}
    readers_tasks = {}
    should_sort = sort_level in list(range(7))
    sort_level_name = sort_level_dict.get(sort_level, "")
    with ProcessPoolExecutor() as executer:
        for filo in tqdm(file_list,
                         total=int(len(file_list)),
                         desc="Fasta Files Parsing",
                         unit="Files"):
            readers_tasks.update({
                executer.submit(get_header_seq_dict, filo, ";", 0, 1, ";"): filo
            })
    completed_tasks = as_completed(readers_tasks)
    for read_task in completed_tasks:
        try:
            fasta_dict, tax_dict = read_task.result()
            fastas_dict.update(fasta_dict)
            taxs_dict.update(tax_dict)
        except Exception as exc:
            logger.warning(f"Could not parse file {str(exc)}")
    # sorting tax_dict
    if should_sort:  # It means sorting
        taxs_dict = dict(
            sorted(
                taxs_dict.items(),
                key=lambda x: x[1][sort_level:sort_level + 1]
            )
        )
        # Updating the seq_tax_map_path path to include sortlevel
        global seq_tax_map_path
        new_seqtaxmap = seq_tax_map_path.parent.joinpath(
            f"{seq_tax_map_path.stem}_{sort_level_name}.tsv")
        # Write the seqid_tax map and order fastas_dict
        with open(new_seqtaxmap, "w") as fileo:
            csvwriter = csv.writer(fileo, delimiter="\t")
            # csvheader
            csvheader = [
                "SeqID",
                "Kingdom",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus",
                "Species"
            ]
            csvwriter.writerow(csvheader)
            for header, tax_list in taxs_dict.items():
                # writing to taxs_dict file
                write_list = [""] * 7  # 7 levels of taxonomy
                for i, tax in enumerate(tax_list):
                    write_list[i] = tax
                write_list = [header] + write_list
                csvwriter.writerow(write_list)
                taxs_dict[header] = fastas_dict[header]
        fastas_dict = taxs_dict
    # return
    return fastas_dict, sort_level_name


def calculat_similarity(seq_one_t: tuple, seq_two_t: tuple):
    """
    Calculates the percentage of positions matched out of all non-gap
        positions with valid base value.
    """
    ind_f = 0
    ind_t = 1
    score = 0
    sim_score = 0.0
    total_match = 0
    h1, seq_one = seq_one_t
    h2, seq_two = seq_two_t
    min_len = min((len(seq_one), len(seq_two)))
    while ind_f < min_len:
        seq_o_pos = seq_one[ind_f:ind_t]
        seq_t_pos = seq_two[ind_f:ind_t]
        nucs = ["a", "c", "t", "g", "n", "u"]
        one_is_valid_nuc = str(seq_o_pos[0]).lower() in nucs if len(seq_o_pos) == 1 else False
        two_is_valid_nuc = str(seq_t_pos[0]).lower() in nucs if len(seq_t_pos) == 1 else False
        if one_is_valid_nuc and two_is_valid_nuc:
            total_match += 1
            if str(seq_o_pos[0]).lower() == str(seq_t_pos[0]).lower():
                score += 1
        ind_f += 1
        ind_t += 1
    try:
        sim_score = float(score / total_match)
    except Exception as exc:
        logger.warning(f"{exc}")

    return h1, h2, sim_score


def create_sim_matrix(files_list, sort_level=-1):
    """
    **Arguments:**
    ==============
    - file_list: a list of fasta files
    - sort_level: taxonomy level with which the headers are sorted
    """
    fastas_dict, sort_level_name = parse_fastas_to_dict(files_list, sort_level=sort_level)
    headers_list = list(fastas_dict.keys())
    init_zero_matrix = {el: [0.0] * len(headers_list) for el in headers_list}
    matrix = pd.DataFrame(init_zero_matrix, index=headers_list, columns=headers_list)
    headers_comb = create_pairwise_combinations(headers_list)
    with ProcessPoolExecutor() as executer:
        sim_tasks = {}
        for h1, h2 in tqdm(headers_comb,
                           total=int(len(headers_comb)),
                           desc="Sequence Similarity",
                           unit="Sequence Pair"):
            seq_1_t = (h1, fastas_dict[h1])
            seq_2_t = (h2, fastas_dict[h2])
            sim_tasks.update({
                executer.submit(calculat_similarity, seq_1_t, seq_2_t): f"{h1}_{h2}"
            })
    completed_tasks = as_completed(sim_tasks)
    for sim_task in completed_tasks:
        try:
            h1, h2, sim = sim_task.result()
            sim = round(sim, 4)
            matrix.at[h1, h2] = sim
            matrix.at[h2, h1] = sim
        except Exception:
            logger.warning(f"Could not calculate similarity {str((h1, h2))}")
    # TODO for lengthy sequence headers the plot should be adjusted visually
    sns_plot = sns.heatmap(
        matrix,
        vmin=0,
        vmax=1,
        cmap="YlGnBu",
        yticklabels=False,
        xticklabels=False,
    )
    fig = sns_plot.get_figure()
    # Updating the output path according to sort level
    global excel_path, graph_path
    new_excel_path = excel_path
    new_graph_path = graph_path
    if sort_level_name:
        new_excel_path = excel_path.parent.joinpath(
            f"{excel_path.stem}_{sort_level_name}.xlsx")
        new_graph_path = graph_path.parent.joinpath(
            f"{graph_path.stem}_{sort_level_name}.png")
    # Writing to the graph and matrix to files
    matrix.to_excel(new_excel_path)
    fig.savefig(new_graph_path)
    # Clears the Axes to avoid multiple colorbar
    plt.clf()


def main(files_list, sort_level_range: list = range(4, 6)):
    """
    **Arguments:**
    - files_list:
    - sort_level_range: list of sorting keys according to sort_level_dict.
        matrices, mapping files and heatmaps gets produced for each level
    """
    for sort_l in sort_level_range:
        create_sim_matrix(files_list, sort_level=sort_l)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs='+', help="list of fasta files")
    args = parser.parse_args()
    main(args.files, range(4, 7))
