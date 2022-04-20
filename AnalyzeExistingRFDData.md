

```bash
SRAIDS=( $(cut -d ',' -f1 sra_explorer_metadata.csv ) )
unset SRAIDS[0]
```

```bash
ml ncbi_sra

mkdir fastqs
cd fastqs
for t in ${SRAIDS[@]}; do
  sbatch --cpus-per-task 12 --wrap="fasterq-dump --threads 12 $t"
  done
```

```bash
OKSeqHeLaRep1=( $( grep 'HeLa EdU Rep1' sra_explorer_metadata.csv | cut -d ',' -f12 ) )
OKSeqHeLaRep2=( $( grep 'HeLa EdU Rep2' sra_explorer_metadata.csv | cut -d ',' -f12 ) )
OKSeqGM06990Rep1=( $( grep 'GM06990 EdU Rep1' sra_explorer_metadata.csv | cut -d ',' -f12 ) )
OKSeqGM06990Rep2=( $( grep 'GM06990 EdU Rep2' sra_explorer_metadata.csv | cut -d ',' -f12 ) )

sbatch --cpus-per-task 2 --wrap="cat SRR2913037.fastq SRR2913042.fastq SRR2913039.fastq SRR2913041.fastq SRR2913038.fastq SRR2913040.fastq SRR2913036.fastq > OKSeqHeLaRep1.fastq"
sbatch --cpus-per-task 2 --wrap="cat SRR2913058.fastq SRR2913055.fastq SRR2913054.fastq SRR2913057.fastq SRR2913056.fastq SRR2913051.fastq SRR2913053.fastq > OKSeqHeLaRep2.fastq"
sbatch --cpus-per-task 2 --wrap="cat SRR2913061.fastq SRR2913063.fastq SRR2913062.fastq SRR2913060.fastq > OKSeqGM06990Rep1.fastq"
sbatch --cpus-per-task 2 --wrap="cat SRR2913064.fastq SRR2913066.fastq SRR2913069.fastq SRR2913065.fastq SRR2913067.fastq SRR2913068.fastq SRR2913070.fastq > OKSeqGM06990Rep2.fastq"
```