#convert fast5 files from sequencing to pod5 files
```bash
pod5 convert fast5 file.fast5 --output OUTDIR
```

#basecall pod5 files
```bash
dorado basecaller --modified-bases 5mCG_5hmCG --no-trim \
hac@v4.1.0 output.pod5 > basecalled.bam
```

#convert unaligned bams to fastq and filter out reads shorter than 80 bases
```bash
samtools fastq -T'*' basecalled.bam > basecalled.fastq
```

#get qc statistics of non-filtered basecalled reads
```r
nanoQC basecalled.fastq -o OUTDIR
```
```r
NanoPlot --fastq basecalled.fastq -o OUTDIR \
--format png --plots dot --loglength
```

#convert unaligned bams to fastq and filter out reads shorter than 80 bases
```bash
seqtk seq -L 80 basecalled.fastq > basecalled_minlength80.fastq
```

#perform demultiplexing with either developed workflow or dorado
```bash
/workflow_demultiplex

#OR

dorado demux \
--barcode-arrangement /data/August/speciale_kode/snakemake_workflows/demultiplex_workflow/resources/custom_arrangements/toml_files/2_double.toml \
--barcode-sequences /data/August/speciale_kode/snakemake_workflows/demultiplex_workflow/resources/barcodes.fasta \
--emit-fastq \
--output-dir classified_reads \
--no-trim \
basecalled_minlength80.fastq.gz
```

#trimming from inner flanks and barcode with cutadapt
```bash
cutadapt -a ACAGACGACTACAAACGGAATCGACAGCACCT...TTAACCTTAGCAATTCGATTCCGTTTGTAGTCGTCTGT -O 10 -e 0.3 -o trimmed_reads/barcode02.fastq -j 0 --buffer-size=16000000 classified_reads/barcode02.fastq && \
cutadapt -a CCTGGTAACTGGGACACAAGACTCCAGCACCT...TTAACCTTAGCAATGAGTCTTGTGTCCCAGTTACCAGG -O 10 -e 0.3 -o trimmed_reads/barcode03.fastq -j 0 --buffer-size=16000000 classified_reads/barcode03.fastq && \
cutadapt -a TAGGGAAACACGATAGAATCCGAACAGCACCT...TTAACCTTAGCAATTTCGGATTCTATCGTGTTTCCCTA -O 10 -e 0.3 -o trimmed_reads/barcode04.fastq -j 0 --buffer-size=16000000 classified_reads/barcode04.fastq && \
cutadapt -a AAGGTTACACAAACCCTGGACAAGCAGCACCT...TTAACCTTAGCAATCTTGTCCAGGGTTTGTGTAACCTT -O 10 -e 0.3 -o trimmed_reads/barcode05.fastq -j 0 --buffer-size=16000000 classified_reads/barcode05.fastq && \
cutadapt -a AAGGATTCATTCCCACGGTAACACCAGCACCT...TTAACCTTAGCAATGTGTTACCGTGGGAATGAATCCTT -O 10 -e 0.3 -o trimmed_reads/barcode07.fastq -j 0 --buffer-size=16000000 classified_reads/barcode07.fastq && \
cutadapt -a ACGTAACTTGGTTTGTTCCCTGAACAGCACCT...TTAACCTTAGCAATTTCAGGGAACAAACCAAGTTACGT -O 10 -e 0.3 -o trimmed_reads/barcode08.fastq -j 0 --buffer-size=16000000 classified_reads/barcode08.fastq && \
cutadapt -a GAGAGGACAAAGGTTTCAACGCTTCAGCACCT...TTAACCTTAGCAATAAGCGTTGAAACCTTTGTCCTCTC -O 10 -e 0.3 -o trimmed_reads/barcode10.fastq -j 0 --buffer-size=16000000 classified_reads/barcode10.fastq && \
cutadapt -a TCCATTCCCTCCGATAGATGAAACCAGCACCT...TTAACCTTAGCAATGTTTCATCTATCGGAGGGAATGGA -O 10 -e 0.3 -o trimmed_reads/barcode11.fastq -j 0 --buffer-size=16000000 classified_reads/barcode11.fastq

```

#run alignment workflow
```bash
/workflow_alignment
```

#run methylaton workflow
```bash
/workflow_methylation
# calculate methylation statistics per sample
```
