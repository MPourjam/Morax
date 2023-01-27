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


# Setting global paths
global SEQ_TAX_MAP_PATH, EXCEL_PATH, GRAPH_PATH, OUTPUT_DIR
OUTPUT_DIR = Path(__file__).parent.joinpath("threshold_calculator_output")
GRAPH_PATH = OUTPUT_DIR.joinpath("similarity_heatmap.png")
EXCEL_PATH = OUTPUT_DIR.joinpath("similarity_matrix.xlsx")
SEQ_TAX_MAP_PATH = OUTPUT_DIR.joinpath("seqid_taxonomy_map.tsv")


def gimmelogger():
    """
    Logger
    ======
    """
    logger = logging.getLogger("Threshold Calculator")
    # set the logging level
    logger.setLevel(logging.DEBUG)
    # create a console and file handler
    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)
    fh = logging.FileHandler(OUTPUT_DIR.joinpath("log_thershold_calculator.txt"))
    fh.setLevel(logging.WARNING)
    # create a formatter
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    # set the formatter to the console handler
    ch.setFormatter(formatter)
    fh.setFormatter(formatter)
    # add the console handler to the logger
    logger.addHandler(ch)
    logger.addHandler(fh)
    return logger


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
    fasta_dict = ddict(str)
    tax_dict = ddict(str)
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
        sort_factor_d = ddict(list)
        for ke, vl in taxs_dict.items():
            tax_level_factor = vl[sort_level:sort_level + 1][0] if vl[sort_level:sort_level + 1] else ""
            sort_factor_d[tax_level_factor] += [ke]
        # Updating the SEQ_TAX_MAP_PATH path to include sortlevel
        global SEQ_TAX_MAP_PATH
        new_seqtaxmap = OUTPUT_DIR.joinpath(
            f"{SEQ_TAX_MAP_PATH.stem}_{sort_level_name}.tsv")
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
    return fastas_dict, sort_level_name, sort_factor_d


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


def calculate_group_sim_avg(matrix, seq_headers: list, tax_level):
    """
    It gets sub-matrix of of matrix by rows and columns
    then calculates average similarity score and returns it.
    **Arguments:**
    ==============
    - matrix: a dataframe containing the rows and columns as sequnce headers
    - seq_headers: a list of sequence headers we want to calculate their combination averge
        similairy
    - tax_level: the group member of sorted level tax set
    """
    cum_sim = 0.0
    tot_count = 0
    if len(seq_headers) > 1:
        seq_h_comb = itertools.combinations(seq_headers, 2)
    elif len(seq_headers) == 1:
        seq_h_comb = [(seq_headers[0], seq_headers[0])]
    else:
        raise ValueError("seq_heades should be a list containg sequnce headers. It's empty!")

    for comb in seq_h_comb:
        cum_sim += matrix.loc[comb[0], comb[1]]
        tot_count += 1
    avg_sim = round(cum_sim / tot_count, 4)
    return seq_headers, tax_level, avg_sim


def write_avg_sim_tsv(taxlevel_seqheader: dict, matrix, sort_level_name=""):
    """
    Takes a taxlevel_seqheader dictionary and a similarity matrix and writes a tsv file
    containing the sequence headers as first column, tax_level value to which sort has
     happend as the second and the average similarity of the level.
     **Arguments:**
     - taxlevel_seqheader: a dictionary like {"uniq_val_sorted_tax_level": [seqheader1, seqheader2, ...]}
     - matrix: a panadas dataframe containing the pairwise similarities of sequences
    """
    global SEQ_TAX_MAP_PATH
    avg_sim_dir = OUTPUT_DIR
    # Calculating averages
    group_sim_tasks = {}
    with ProcessPoolExecutor() as executer:
        for tax_level, seqlist in tqdm(taxlevel_seqheader.items(),
                                       total=int(len(taxlevel_seqheader)),
                                       desc="Calulating Group Similarities",
                                       unit="Groups"):
            group_sim_tasks.update({
                executer.submit(calculate_group_sim_avg, matrix, seqlist, tax_level): tax_level
            })
    completed_tasks = as_completed(group_sim_tasks)
    rows_to_write = []
    for ct in completed_tasks:
        try:
            seq_headers, tax_level, avg_sim = ct.result()
        except Exception as exc:
            logger.warning(f"Group similarity calculation failed for {tax_level}. {str(exc)}")
        else:
            for seq_h in seq_headers:
                rows_to_write.append([seq_h, tax_level, avg_sim])
    # Writing all to a file
    sort_part = f"_{sort_level_name}" if sort_level_name else ""
    avg_sim_file = avg_sim_dir.joinpath(f"avg_sims{sort_part}.tsv")
    # avg_sim_dir.mkdir(parents=True, exist_ok=True)
    # logger.info(f"Created directory: {avg_sim_dir}")
    if rows_to_write:
        with open(avg_sim_file, "w") as filew:
            csvwriter = csv.writer(filew, delimiter="\t")
            csv_header = ["SeqID", "Group", "Average_Similarity"]
            csvwriter.writerow(csv_header)
            for map_row in rows_to_write:
                csvwriter.writerow(map_row)
    else:
        avg_sim_file = ""
        try:
            # if directory is empty it removes it in case of failure
            avg_sim_dir.rmdir()
        except Exception:
            pass
        logger.warning(f"No values to write to {avg_sim_file}!")

    return avg_sim_file


def create_sim_matrix(files_list, sort_level=-1):
    """
    **Arguments:**
    ==============
    - file_list: a list of fasta files
    - sort_level: taxonomy level with which the headers are sorted
    """
    parse_fastas_result = parse_fastas_to_dict(files_list, sort_level=sort_level)
    fastas_dict, sort_level_name, sortedtaxlevel_seqh_dict = parse_fastas_result
    # Preparing headers
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
    global EXCEL_PATH, GRAPH_PATH
    new_excel_path = EXCEL_PATH
    new_graph_path = GRAPH_PATH
    if sort_level_name:
        new_excel_path = OUTPUT_DIR.joinpath(
            f"{new_excel_path.stem}_{sort_level_name}.xlsx")
        new_graph_path = OUTPUT_DIR.joinpath(
            f"{new_graph_path.stem}_{sort_level_name}.png")
    # Writing to the graph and matrix to files
    matrix.to_excel(new_excel_path)
    fig.savefig(new_graph_path)
    # Clears the Axes to avoid multiple colorbar
    plt.clf()
    # TODO Calculate similarity averages using sortedtaxlevel_seqh_dict
    avg_sim_file_path = write_avg_sim_tsv(
        sortedtaxlevel_seqh_dict,
        matrix,
        sort_level_name
    )
    logger.info("Averge similarities written to {avg_sim_file_path}")


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
    # global OUTPUT_DIR
    parser = argparse.ArgumentParser()
    parser.add_argument('-o', '--output-dir',
                        type=str,
                        default=OUTPUT_DIR,
                        required=False,
                        help="The dirctory to put the result")
    parser.add_argument("files", nargs='+', help="list of fasta files")
    args = parser.parse_args()
    OUTPUT_DIR = Path(args.output_dir)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    # setting logger
    global logger
    logger = gimmelogger()
    main(args.files, range(4, 6))
