import numpy as np
import pandas as pd
import os
import gzip
import matplotlib.pyplot as plt
import seaborn as sns
# from snakemake.script import snakemake

def get_entropy(x):
    # remove all barcode fractions which are zero
    x = x[x>0]
    # calculate the shannon entropy of the barcode classifications
    return abs(-np.sum(x * np.log2(x)))

def perform_test(df, ref_entropy, iterations):
    # get entropy of classifications
    entropy = ref_entropy
    # convert dataframe to array
    my_arr = df.values
    # iterate over array, shuffling the classifications for each iterion and calculating an
    # entropy for the shuffled classifications
    for i in range(iterations):
        for col in range(my_arr.shape[1]):
            np.random.shuffle(my_arr[:, col])
        pdist = [np.unique(row, return_counts=True)[1] / len(row) for row in my_arr]
        shuffled_entropy = [get_entropy(arr) for arr in pdist]
        entropy = np.column_stack((entropy, shuffled_entropy))
    # get the fraction of shuffled entropies that are higher than the entropy calculated
    # from the classifications. This is our p-value
    p_vals = np.apply_along_axis(lambda row: sum(row[1:] <= row[0]) / iterations, axis=1, arr=entropy)

    return p_vals

def perform_test_remove_na(df, ref_entropy, iterations):
    # get entropy of classifications
    entropy = ref_entropy
    # convert dataframe to array
    my_arr = df.values
    # iterate over array, shuffling the classifications for each iterion and calculating an
    # entropy for the shuffled classifications
    for i in range(iterations):
        for col in range(my_arr.shape[1]):
            np.random.shuffle(my_arr[:, col])
        # additonnally apply filter to remove NA values from probability distribution calculations
        # this is neccessary since we removed "unclassified" classifications.    
        pdist = [np.unique(row[~pd.isna(row)], return_counts=True)[1] / len(row[~pd.isna(row)]) for row in my_arr]
        shuffled_entropy = [get_entropy(arr) for arr in pdist]
        entropy = np.column_stack((entropy, shuffled_entropy))
    # get the fraction of shuffled entropies that are higher than the entropy calculated
    # from the classifications. This is our p-value
    p_vals = np.apply_along_axis(lambda row: sum(row[1:] <= row[0]) / iterations, axis=1, arr=entropy)

    return p_vals

def split_fastq_by_classification(fastq_file, classification_df, output_dir = "classify_reads_output", compressed = False):
    """
    Splits a FASTQ file into separate files based on classifications provided in a TSV file.

    Args:
        fastq_file (str): Path to the input FASTQ file.
        tsv_file (str): Path to the TSV file with read IDs and classifications.
        output_dir (str): Directory to save the classified FASTQ files.
        compressed (bool): Whether the FASTQ file is gzipped.
    """
    # initialize a dictionary of reads and their barcode classifications
    classifications = {}
    for index, value in zip(classification_df["read_id"], classification_df["barcode"]):
        read_id, classification = index, value 
        classifications[read_id] = classification

    # create dicitonary for file handles for each classification
    barcode_files = {}

    # open FASTQ file
    open_func = gzip.open if compressed else open
    with open_func(fastq_file, 'rt') as fastq:
        while True:
            # read four lines, equivilant to a single read, at a time
            header = fastq.readline()
            if not header:
                break  # end of file
            sequence = fastq.readline()
            plus = fastq.readline()
            quality = fastq.readline()

            # extract read-id (part of header before first whitespace)
            read_id = header.split()[0][1:] # not including @ at start of the line

            # create output folder
            os.makedirs(output_dir, exist_ok=True)

            # determine classification
            classification = classifications.get(read_id)
            if classification not in barcode_files:
                # create a new file for the classification
                output_file =  "/".join([output_dir, f"{classification}.fastq"])
                barcode_files[classification] = open(output_file, 'w')
            # write the read to the appropriate file
            barcode_files[classification].write(f"{header}{sequence}{plus}{quality}")

    # close all files
    for handle in barcode_files.values():
        handle.close()


classifications = pd.read_csv(snakemake.input[1], sep="\t")
barcode_assignment = pd.read_csv(snakemake.input[2], sep="\t")

# remove read_id column and column for dorado_rear tool with very poor performance
no_read_id = classifications.drop(columns=["read_id", "dorado_rear"])

# get p-values
p_vals = perform_test_remove_na(no_read_id, barcode_assignment["entropy"].to_numpy(), 100)

barcode_assignment["p_val"] = p_vals

# split reads into files equal to their barcode assignment, if p-val from test
# is less than 0.05
split_fastq_by_classification(snakemake.input[0], 
                              barcode_assignment[barcode_assignment["p_val"] <= 0.05], 
                              os.path.join(snakemake.params.output_dir, "classified_reads"),
                              snakemake.params.compressed)

# create abrplot fo read assignments
order = ["unclassified", "barcode02", "barcode03", "barcode04", 
         "barcode05", "barcode07", "barcode08", "barcode10", "barcode11"]
barcode_frac = sns.countplot(data=barcode_assignment[barcode_assignment["p_val"] <= 0.05], 
                             x="barcode", 
                             edgecolor='black', 
                             linewidth=1,
                             order=order)
barcode_frac.set_title("Number of reads assigned to each barcode")
plt.xticks(rotation=45, ha="right")
plt.legend([f"Total Count: {len(barcode_assignment)}"], loc="upper right")
plt.tight_layout()
plt.savefig(f"{snakemake.params[0]}/assignments_plot.png")
plt.close()

barcode_assignment.to_csv(f"{snakemake.params[0]}/assignments_table.txt", sep="\t", index=False)
