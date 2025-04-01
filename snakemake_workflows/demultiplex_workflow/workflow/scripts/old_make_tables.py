import numpy as np
import pandas as pd
import glob
import os
import re
import matplotlib.pyplot as plt
import seaborn as sns
from snakemake.script import snakemake

def extract_classification(file_path):
    split_path = file_path.split(os.path.normcase("/"))
    tool = "_".join(split_path[-3:-1])
    classification = split_path[-1]
    classification = classification.split("_")[-1]
    classification = classification.split(".")[0]
    read_id_list = []
    # parse over fastq file
    with open(file_path, 'r') as file:
        while True:
            # Read the four lines of each FASTQ entry
            split_id = re.split(r'[ \t]+', file.readline())
            modified_id = split_id[0][1:]
            read_id = modified_id.strip()
            if not read_id:  # If no line, we're done
                break
            sequence = file.readline().strip()
            _ = file.readline().strip()  # Skip the plus line
            quality = file.readline().strip()

            # Extract relevant information
            read_id_list.append(read_id)

    df = pd.DataFrame({"classification": classification, 
                      "read_id": read_id_list, 
                      "tool": tool})

    return df

def get_entropy(x):
    x = x[x>0]
    return abs(-np.sum(x * np.log2(x)))


dorado_base_dirs = [path.removesuffix("unclassified.fastq") for path in snakemake.input["dorado_input"]]

dorado_path_full = [
    path
    for base_dir in dorado_base_dirs
    for path in glob.glob(os.path.join(base_dir, "*.fastq"))
]

list_of_paths = snakemake.input["cutadapt_input"] + dorado_path_full

# list_of_paths = glob.glob(os.path.normcase(
#     "/data/August/speciale_kode/demultiplex_workflow/results/*/*/*"))

classifications = []
for path in list_of_paths:
    classifications.append(extract_classification(path))

# concatenate dataframes along rows
df = pd.concat(classifications, axis = 0)
# reshape df
reshaped_df = df.pivot(index="read_id", columns="tool", values="classification")
reshaped_df.replace("unknown", "unclassified", inplace=True)
pdist = reshaped_df.apply(lambda x: x.value_counts() / len(x), axis=1)
pdist.apply(get_entropy, axis = 1)
entropy = np.array(pdist.apply(get_entropy, axis = 1))
barcode_assignment = pdist.idxmax(axis=1)
barcode_assignment = barcode_assignment.to_frame("barcode")
barcode_assignment["entropy"] = entropy
# p_vals, adj_p_vals = perform_test(reshaped_df, entropy, 1000)
# classification["p-value"] = p_vals
# classification["bonf. adj. p-value"] = adj_p_vals
# classification.to_csv("table.txt", sep='\t', index=False)

barcode_assignment.to_csv(snakemake.output[0], sep='\t', index=True)

sns.set_theme()
entropy_hist = sns.histplot(data=barcode_assignment, x="entropy", edgecolor='black', 
                            linewidth=1, bins=12)
entropy_hist.set_title("Entropy distribution")
plt.savefig(snakemake.output[1])
plt.close()

barcode_frac = sns.countplot(data=barcode_assignment, x="barcode", edgecolor='black', 
                             linewidth=1)
barcode_frac.set_title("Barcode fractions")
plt.xticks(rotation=45, ha="right")
plt.legend([f"Total Count: {len(barcode_assignment)}"], loc="upper right")
plt.tight_layout()
plt.savefig(snakemake.output[2])
plt.close()

reshaped_df.to_csv(snakemake.output[3], sep='\t', index=True)

dataframes = []
for barcode in pd.unique(reshaped_df.values.ravel("K")):
    df = reshaped_df.apply(lambda col: (col.str.contains(barcode)).mean())
    df.name = barcode
    dataframes.append(df)

tool_stats = pd.concat(dataframes, axis=1)

tool_stats.to_csv(snakemake.output[4], sep='\t', index=True)

plt.figure(figsize=(10, 8))
heatmap = sns.heatmap(tool_stats, annot=True, fmt=".2f", cmap="crest", cbar=True, 
            linewidths=0.5, linecolor="black", edgecolor="black")

# Add a border to the colorbar
cbar = heatmap.collections[0].colorbar
cbar.outline.set_visible(True)  # Make the border of the colorbar visible
cbar.outline.set_edgecolor('black')  # Set the color of the colorbar border
cbar.outline.set_linewidth(1)  # Set the width of the colorbar border

# Add labels and title
plt.title("Classification rates for different tools", fontsize=16)
plt.xlabel("Barcode", fontsize=15)
plt.ylabel("Tool", fontsize=15)

# Ensure row and column labels are visible
plt.xticks(fontsize=13, rotation=45, ha='right')
plt.yticks(fontsize=13, rotation=0)

# Show the plot
plt.tight_layout()
plt.savefig(snakemake.output[5])
plt.close()

# split_fastq_by_classification(fastq_file=),
#                                   classification_df=classification,
#                                   output_dir=),
#                                   compressed=)


# pdist = df.apply(lambda x: x.value_counts() / len(x), axis=1)
# barcode_assignment = pdist.idxmax(axis=1)
# entropy = np.array(pdist.apply(get_entropy, axis = 1))
# barcode_assignment = barcode_assignment.to_frame("barcode")
# barcode_assignment["entropy"] = entropy

# p_vals, adj_p_vals = perform_test(df, 1000)

# classification["p-value"] = p_vals
# classification["bonf. adj. p-value"] = adj_p_vals

# classification.sort_values(by="p-value", ascending=True)