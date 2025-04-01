from snakemake.shell import shell
# script which performs dorado classifications for each configuration with snakemake
if ("front" in snakemake.wildcards.location or "double" in snakemake.wildcards.location):
            shell_command = f'''
            dorado demux \
            --barcode-arrangement resources/custom_arrangements/toml_files/{snakemake.params.score}_{snakemake.wildcards.location}.toml \
            --barcode-sequences {snakemake.input.barcodes} \
            --emit-fastq \
            --output-dir results/dorado{snakemake.wildcards.score}_O{snakemake.wildcards.O}_E{snakemake.wildcards.E}/dorado/{snakemake.wildcards.location} \
            {snakemake.input.fastq}
            '''
else:
    shell_command = f'''
            dorado demux \
            --barcode-arrangement resources/custom_arrangements/toml_files/{snakemake.params.score}_{snakemake.wildcards.location}.toml \
            --barcode-sequences {snakemake.input.barcodes_reversed} \
            --emit-fastq \
            --output-dir results/dorado{snakemake.wildcards.score}_O{snakemake.wildcards.O}_E{snakemake.wildcards.E}/dorado/{snakemake.wildcards.location} \
            {snakemake.input.fastq}
            '''

shell(shell_command)