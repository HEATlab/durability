#!/usr/bin/env python

import argparse
from dataparser import DataParser

HEADER_DICT = {"agt": "agent_count",
               "tp": "timepoint_count",
               "ia": "interagent_constraint_density",
               "ctg": "contingent_constraint_density",
               "run": "run_index",
               "try": "attempts_to_gen_valid_stn",
               "early": "early_first_robustness",
               "stat": "srea_robustness",
               "n_stat": "new_SREA_robustness",
               "dyn": "DREA_robustness",
               "dc": "dynamic_controllability_robustness",
               "ct#": "contingent_edge_count",
               "ia#": "interagent_edge_count",
               "samp": "samples",
               "ps": "threads",
               "gnt": "gnt",
               "n_gnt": "n_gnt",
               "early_t": "early_first_computation_time_sec",
               "n_stat_t": "new_SREA_computation_time_sec",
               "stat_t": "SREA_computation_time_sec",
               "dyn_t": "DREA_computation_time_sec",
               "dc_t": "dynamic_controllability_computation_time_sec",
               "std": "stddev_of_contingent_edges_sec",
               "tght": "interagent_tightness_millisec"}


##
# \fn main
# \brief
def main():
    args = parse_args()
    filepaths = args.f
    for fp in filepaths:
        convert_file(fp)


##
# \fn convert_file
# \brief Converts a multiComparsion output file to a CSV file.
#
# @param fp Filepath of the multiComparison output file.
def convert_file(fp):
    dp = DataParser()
    dp.load_matrix(fp)
    dp = dp.rename_headers(HEADER_DICT)
    new_name = fp
    if fp[-4:] == ".dat":
        new_name = fp[:-4] + ".csv"
    else:
        new_name += "_converted"
    dp.to_csv(outpath=new_name)


##
# \fn parse_args
# \brief Parse command line arguments.
def parse_args():
    # Handle argument parsing
    parser = argparse.ArgumentParser(prog='multiconverter',
                                     description='Converts a multiComparison'
                                     + 'data output to a CSV.')
    parser.add_argument('f',
                        metavar='DATAFILE',
                        type=str,
                        nargs='+',
                        help='File from multiComparison.py to convert.')
    return parser.parse_args()

if __name__ == "__main__":
    main()
