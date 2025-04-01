# Thesis
Repository for August Bodilsen / kfl916 KU thesis

Code for plotting and performing statistical tests on derived length and methylation data is available in the scripts folder.

Pipeline applied to each sequencing run/sample is available in the bash_code.md file

Snakemake workflows are available in the snakemake_workflows folders.

demultiplex_workflow:
  - Basefolder for developed demultiplex pipeline integrated as a snakemake workflow.

alignment_workflow:
  - Snakemake workflow for performing filtering, alignment and other actions applied to all samples.

methylation_workflow:
  - Small seperate snakemake workflow for perfoming methylation analysis of each sample.
