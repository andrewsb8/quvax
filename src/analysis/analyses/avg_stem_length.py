#Script to take in a single .ct file and gather total stem length and number of stems  

import numpy as np
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Gets Average Length of Stems')
parser.add_argument('--path', type=str, required=True, help='Path to the .ct file')
args = parser.parse_args()

def main():

        connect_table_df = pd.read_csv(args.path, sep=r'\s+', skiprows=1, header=None)
        connect_table_df.columns =["Index","Nucleotide","Previous","Next","Paired With","Counter"]

        sequence_len = len(connect_table_df["Index"])
        stems = []
        i=0

        while i < sequence_len:
            # move on from unpaired bases and don't double count base pairs
                if (
                        connect_table_df["Paired With"].iloc[i]
                        < connect_table_df["Index"].iloc[i]
                        or connect_table_df["Paired With"].iloc[i] == 0
                ):
                        i += 1
                # if base pair is found, need to define the length of stem
                elif connect_table_df["Paired With"].iloc[i] != 0:
                # loop through connect table until nonsequential base pair is
                # found and then append to stems
                        for j in range(i + 1, sequence_len - 1):
                    # if nonsequential base pair ordering or a unpaired base is found, use information to generate stem
                                if (
                                        connect_table_df["Paired With"].iloc[j]
                                        != connect_table_df["Paired With"].iloc[j - 1] - 1
                                        or connect_table_df["Paired With"].iloc[j] == 0
                                ):
                                        stems.append(
                                                (
                                                        connect_table_df["Index"].iloc[i],
                                                        connect_table_df["Paired With"].iloc[i],
                                                        j - i,
                                                )
                                        )
                                        i += j-i
                                        break
        num_stems = len(stems)
        tot_length = sum(inner_list[2] for inner_list in stems)
        #want to scale the total length with sequence length
        scaled_avg_length = tot_length/(num_stems)
        print(scaled_avg_length)

if __name__ == "__main__":
        main()

