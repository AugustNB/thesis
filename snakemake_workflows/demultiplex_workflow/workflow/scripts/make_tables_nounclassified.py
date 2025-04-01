import numpy as np
import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
# from snakemake.script import snakemake

def extract_classification(file_path):
    # extract information from file path
    split_path = file_path.split(os.path.normcase("/"))
    tool = "_".join(split_path[-3:-1])
    classification = split_path[-1]
    classification = classification.split("_")[-1]
    classification = classification.split(".")[0]
    read_id_list = []
    # parse over fastq file
    with open(file_path, 'r') as file:
        while True:
            # read the four lines of each FASTQ entry
            # extract first part of read-ID, which is constant between
            # dorado and cutadapt classified files
            split_id = re.split(r'[ \t]+', file.readline())
            modified_id = split_id[0][1:]
            read_id = modified_id.strip()
            if not read_id:  # if no line, file is finished
                break
            # iterate over last three lines. sequence and quality are not extarcted
            sequence = file.readline().strip() #
            _ = file.readline().strip()  # skip line with "+" sign
            quality = file.readline().strip()

            # extract read-id
            read_id_list.append(read_id)
    # create dataframe with read-id, classification and tool for each classification made
    df = pd.DataFrame({"classification": classification, 
                      "read_id": read_id_list, 
                      "tool": tool})

    return df

def get_entropy(x):
    # remove all barcode fractions which are zero
    x = x[x>0]
    # calculate the shannon entropy of the barcode classifications
    return abs(-np.sum(x * np.log2(x)))

# remove unclassified suffix. The unclassified file is set as target, because it will always
# be created. But we are only interested in the file path lead to the "unclassified.fastq" file
dorado_base_dirs = [path.removesuffix("unclassified.fastq") for path in snakemake.input["dorado_input"]]

# extract file paths to classifications made by dorado
dorado_path_full = [
    path
    for base_dir in dorado_base_dirs
    for path in glob.glob(os.path.join(base_dir, "*.fastq"))
]

# grab cutadapt file paths from snakemake input, while concatenating the just extracted
# dorado paths
list_of_paths = snakemake.input["cutadapt_input"] + dorado_path_full

# now we have file paths for all our classification files and we can extract the
# actual classifications into a list of dataframes
classifications = []
for path in list_of_paths:
    classifications.append(extract_classification(path))

# concatenate dataframes along rows
df = pd.concat(classifications, axis = 0)

# merge dataframe on read-id
reshaped_df = df.pivot(index="read_id", columns="tool", values="classification")
# replace unknown with unclassified to have coherence
reshaped_df.replace("unknown", "unclassified", inplace=True)
# calculate probability distribution for each read
pdist = reshaped_df.apply(lambda x: x.value_counts() / len(x), axis=1)
# get entropy for each read
entropy = np.array(pdist.apply(get_entropy, axis = 1))
# assign the barcode with the highest probability
barcode_assignment = pdist.idxmax(axis=1)
# convert to dataframe
barcode_assignment = barcode_assignment.to_frame("barcode")
# same entropy
barcode_assignment["entropy"] = entropy

# remove unclassified barcodes
df_nounknown = df[df["classification"]!="unknown"]
df_nounclassified = df_nounknown[df_nounknown["classification"]!="unclassified"]

# perform same steps as above for the dataframe without "unclassified" classificatons
reshaped_df_nounclassified = df_nounclassified.pivot(index="read_id", columns="tool", values="classification")
pdist_nounclassified = reshaped_df_nounclassified.apply(lambda x: x.value_counts() / len(x[x.notna()]), axis=1)
entropy_nounclassified = np.array(pdist_nounclassified.apply(get_entropy, axis = 1))
barcode_assignment_nounclassified = pdist_nounclassified.idxmax(axis=1)
barcode_assignment_nounclassified = barcode_assignment_nounclassified.to_frame("barcode")
barcode_assignment_nounclassified["entropy"] = entropy_nounclassified

# set order of plot labels
order = ["unclassified", "barcode02", "barcode03", "barcode04", 
         "barcode05", "barcode07", "barcode08", "barcode10", "barcode11"]

# save assignments to file
barcode_assignment.to_csv(snakemake.output[0], sep='\t', index=True)

# save file with all classifications (not just assigned)
reshaped_df.to_csv(snakemake.output[3], sep='\t', index=True)

# create different plots to investigate summary statistics of the classifications

# create histogram of read entropies
sns.set_theme()
entropy_hist = sns.histplot(data=barcode_assignment, x="entropy", edgecolor='black', 
                            linewidth=1, bins=12)
entropy_hist.set_title("Entropy distribution")
plt.savefig(snakemake.output[1])
plt.close()

# create barplot of barcode classifications
barcode_frac = sns.countplot(data=barcode_assignment, x="barcode", edgecolor='black', 
                             linewidth=1, order=order)
barcode_frac.set_title("Tootal number of classifications for each barcode")
plt.xticks(rotation=45, ha="right")
plt.legend([f"Total Count: {len(barcode_assignment)}"], loc="upper right")
plt.tight_layout()
plt.savefig(snakemake.output[2])
plt.close()

dataframes = []
for barcode in pd.unique(reshaped_df.values.ravel("K")):
    df = reshaped_df.apply(lambda col: (col.str.contains(barcode)).mean())
    df.name = barcode
    dataframes.append(df)

tool_stats = pd.concat(dataframes, axis=1)

# save statistics about the classifications made
tool_stats.to_csv(snakemake.output[4], sep='\t', index=True)

# create heatmap og classifications by each tool
plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(tool_stats, annot=True, fmt=".2f", cmap="crest", cbar=True, 
            linewidths=0.5, linecolor="black", edgecolor="black")

cbar = heatmap.collections[0].colorbar
cbar.outline.set_visible(True)
cbar.outline.set_edgecolor('black')
cbar.outline.set_linewidth(1)

plt.title("Classification rates for different tools", fontsize=16)
plt.xlabel("Barcode", fontsize=15)
plt.ylabel("Tool", fontsize=15)

plt.xticks(fontsize=13, rotation=45, ha='right')
plt.yticks(fontsize=13, rotation=0)

plt.tight_layout()
plt.savefig(snakemake.output[5])
plt.close()

# create plots for the classifications made without "unclassified" reads
barcode_assignment_nounclassified.to_csv(snakemake.output[6], sep='\t', index=True)

sns.set_theme()
entropy_hist_nounclassified = sns.histplot(data=barcode_assignment_nounclassified, x="entropy", edgecolor='black', 
                            linewidth=1, bins=12)
entropy_hist_nounclassified.set_title("Entropy distribution")
plt.savefig(snakemake.output[7])
plt.close()

barcode_frac_nounclassified = sns.countplot(data=barcode_assignment_nounclassified, x="barcode", edgecolor='black', 
                             linewidth=1)
barcode_frac_nounclassified.set_title("Barcode fractions")
plt.xticks(rotation=45, ha="right")
plt.legend([f"Total Count: {len(barcode_assignment_nounclassified)}"], loc="upper right")
plt.tight_layout()
plt.savefig(snakemake.output[8])
plt.close()

reshaped_df_nounclassified.to_csv(snakemake.output[9], sep='\t', index=True)