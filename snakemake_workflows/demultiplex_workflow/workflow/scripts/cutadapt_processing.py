from snakemake.shell import shell
# script which performs cutadapt classifications for each configuration with snakemake
name = "{{name}}"

if "front" in snakemake.wildcards.location:
    shell_command = f"""
    cutadapt -g file:{snakemake.input.barcodes} -O {snakemake.params.O} -e {snakemake.params.E} \
    -j 0 --quiet \
    --action none \
    --buffer-size=16000000 \
    -o results/dorado{snakemake.wildcards.score}_O{snakemake.wildcards.O}_E{snakemake.wildcards.E}/cutadapt/{snakemake.wildcards.location}/{name} \
    {snakemake.input.fastq}
    """
else:
    shell_command = f"""
    cutadapt -a file:{snakemake.input.barcodes} -O {snakemake.params.O} -e {snakemake.params.E} \
    -j 0 --quiet \
    --action none \
    --buffer-size=16000000 \
    -o results/dorado{snakemake.wildcards.score}_O{snakemake.wildcards.O}_E{snakemake.wildcards.E}/cutadapt/{snakemake.wildcards.location}/{name} \
    {snakemake.input.fastq}
    """

shell(shell_command)
