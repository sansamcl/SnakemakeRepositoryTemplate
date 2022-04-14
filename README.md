# SnakemakeRepositoryTemplate

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm
ml miniconda
```
#### 3A.  FIRST TIME ONLY:  Setup conda environment with snakemake
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/SnakemakeEnv.yml -p /s/sansam-lab/SnakemakeEnv 
```

#### 3B.  Activate conda environment with snakemake
```bash
conda activate /s/sansam-lab/SnakemakeEnv
```


```bash
git clone https://github.com/sansamcl/SnakemakeRepositoryTemplate.git
```

#### 7A. Use conda environments
If conda is to be used for rule-specific environments, you may find it useful to create the environments first. Running 'snakemake' with the '--conda-create-envs-only' option will create the environments without running the pipeline. The '--conda-prefix' option is used to set a directory in which the ‘conda’ and ‘conda-archive’ directories are created. This directory may be changed to a stable or shared location.
```bash
sbatch --mem 32G \
--wrap="\
snakemake \
--cores all \
--use-conda \
--conda-prefix ../condEnvs/ \
--conda-create-envs-only \
--conda-frontend conda"
```

Once the environments are setup, you may execute pipeline with conda environments using the following command:
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-conda \
--conda-prefix ../condEnvs/ \
--conda-frontend conda \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```

#### 7B. Use environment modules.
Rather than using conda environments, you may prefer to use modules installed on your computing cluster. These modules are defined for each rule in 'workflow/Snakefile'. This must be customized for your environment, and you must modify the Snakefile yourself.

To execute the pipeline with environment modules, enter the following:
```bash
sbatch --constraint=westmere \
--wrap="\
snakemake \
-R \
-j 999 \
--use-envmodules \
--latency-wait 100 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output}'"
```


# HiC-SnakeMake
<img src="https://github.com/SansamLab/HiC-SnakeMake/blob/main/HiCSnakemake2.png?raw=true" alt="" width="300"/>

## Project Description:

## Table of contents:
* [Description of individual steps in pipeline](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#description-of-individual-steps-in-pipeline)
  * [1.  run_bwa_mem](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#1--run_bwa_mem)
  * [2.  make_pairs_with_pairtools_parse](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#2--make_pairs_with_pairtools_parse)
  * [3.  sort_pairs_with_pairtools_sort](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--sort_pairs_with_pairtools_sort)
  * [4.  mark_duplicates_with_pairtools_dedup](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--mark_duplicates_with_pairtools_dedup)
  * [5.  filter_pairs](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--filter_pairs)
  * [6.  add_frag2Pairs](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--add_frag2pairs)
  * [7.  run_cooler](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#7--run_cooler)
* [Step-by-step instructions on running Snakemake pipeline:](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#step-by-step-instructions-on-running-snakemake-pipeline)
  * [1.  Load slurm and miniconda](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#1--load-slurm-and-miniconda)
  * [2.  Clone repository](https://github.com/SansamLab/Process_HiC_SnakeMake#2--clone-repository)
  * [3.  Start the conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3--start-the-conda-environment)
    * [3A.  FIRST TIME ONLY:  Setup conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3a--first-time-only--setup-conda-environment)
    * [3B.  Activate conda environment](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#3b--activate-conda-environment)
  * [4.  Modify the job-specific configuration files.](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--modify-the-job-specific-coniguration-files)
    * [4A.  Modify the config/config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4a--modify-the-configconfigyml-file)
    * [4B.  Modify the config/samples.csv file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4b--modify-the-configsamplescsv-file)
    * [4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4c--if-slurm-resource-changes-are-needed-modify-the-configcluster_configyml-file)
  * [5.  Do a dry run](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#4--do-a-dry-run)
  * [6.  Make a DAG diagram](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#5--make-a-dag-diagram)
  * [7.  Run on cluster with slurm](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/README.md#6--run-on-cluster-with-slurm)

## Description of individual steps in pipeline:
![DAG of Pipeline](dag.svg)

### 1.  run_bwa_mem
```bash
# align pair of .fastq files to the genome and convert the .sam output to .bam
bwa mem \
 -t {params.bwaThreads} \`            # number of threads to use for alignment`
 -SP5M \`                             # -5 for split alignment, take the alignment with the smallest coordinate as primary
                                      # -S skip mate rescue
                                      # -P skip pairing
                                      # -M mark shorter split hits as secondary`
 {params.bwaIndex} \`                 # path to bwa indexed genome`
 {input.fq1} \`                       # path to fastq 1`
 {input.fq2} | \`                     # path to fastq 2`
 samtools view -Shb - > {output.bam}` # convert bwa mem output to .bam file with header`
```
### 2.  make_pairs_with_pairtools_parse
```bash
pairtools parse \
  -c {params.chrom_sizes} \
  --drop-sam \
  --add-columns mapq \
  --output {output.pairs} \
  {input.bam}
```

### 3.  sort_pairs_with_pairtools_sort
```bash
# if temporary directory is not available, make it
[ -d {params.tempdir} ] || mkdir {params.tempdir}
# sort .pairs file created in step 2
pairtools sort \
  --nproc {params.nproc} \
  --memory {params.memory} \
  --tmpdir {params.tempdir} \
  --output {output.sorted_pairs} \
  {output.pairs}
# remove temporary directory after sorting
rm -rf {params.tempdir}
```

### 4.  mark_duplicates_with_pairtools_dedup
```bash
# Find and mark PCR/optical duplicates in .pairs file from step 3.
pairtools dedup \
 --mark-dups \`                    # duplicate pairs are marked as DD in pair_type and as a duplicate in the sam entries.`
 --output-dups - \`                # output duplicates together with deduped pairs`
 --output-unmapped - \`            # output unmapped pairs together with deduped pairs`
 --output {output.marked_pairs} \` # name of output .pairs file with duplicates marked`
 {input.sorted_pairs}`             # .pairsam file input from step 3`
 
# index the .pairsam
pairix {output.marked_pairs}
```
### 5.  filter_pairs
```bash
## Generate lossless bam from the pairsam file
pairtools split \
 --output-sam {output.lossless_bam} \`# name of .bam file produced`
 {input.marked_pairs}`                # .pairsam file input from step 3`
 
# Select UU, UR, RU reads
 ## UU = unique-unique, both alignments are unique
 ## UR or RU = unique-rescued, one alignment unique, the other rescued
pairtools select \
 '(pair_type == "UU") or (pair_type == "UR") or (pair_type == "RU")' \
 --output-rest {output.unmapped_sam} \ `# name of file with the remainder of read pairs`
 --output {params.temp_file} \`         # temporary output file with the selected read pairs`
 {input.marked_pairs}`                  # .pairs file input from step 4`
 
# Generate .pairs file from the UU, UR, and RU reads selected above
pairtools split \
 --output-pairs {params.temp_file1} \`  # temporary .pairs output file`
 {params.temp_file}`                    # input .pairsam file generated with pairtools select above`
 
# Make a .pairs file with only pairs in chromosomes of interest
pairtools select 'True' \
 --chrom-subset {params.chrom_sizes} \` # path to a chromosomes file containing a chromosome subset of interest.`
 -o {output.dedup_pairs} \`             # ouput path and filename for .pairs file in chromosomes of interest`
 {params.temp_file1}`                   # input .pairsam file generated with pairtools split above`

# index the .pairs file
pairix {output.dedup_pairs}            
```
### 6.  add_frag2Pairs
The pearl script used here (fragment_4dnpairs.pl) is from this [GitHub](https://github.com/aidenlab/juicer.git) (Release 1.6).
Github - aidenlab/juicer: A one-click system for analyzing loop-resolution hi-c experiments. (n.d.). GitHub. Retrieved March 10, 2022, from https://github.com/aidenlab/juicer. This requres a restriction site file. See [these instructions](https://github.com/SansamLab/Process_HiC_SnakeMake/blob/main/RE_File_Instructions.md) for generating this file outside of the snakemake pipeline.

#### 6A.  FIRST TIME ONLY. Generate restriction site file.

```bash
python generate_site_positions.py 'HindIII' 'ecoli' '/net/qlotsam.lan.omrf.org/qlotsam/data/sansam/hpc-nobackup/scripts/CutAndRun_2020Aug1/ecoliBowtie2/GCF_000005845.2_ASM584v2_genomic.fna'"
```


```bash
# convert to fragment map
gunzip -ck {input.dedup_pairs} | \`   # unzip dsthe .pairs file generated in step 4`
 workflow/scripts/fragment_4dnpairs.pl \
  -a - `                              # allows replacing existing frag1/frag2 columns`\
  {params.frag2_pairs_basename} \ `   # filename for .pairs file generated`
  {params.restriction_file}`          # a restriction site file, which lists on each line, the sorted locations of the enzyme restriction sites.`

# Compress the .pairs file generated
bgzip -f {params.frag2_pairs_basename}

# Index the .pairs file generated
pairix -f {output.frag2_pairs}
```

### 7.  run_cooler

```bash
# use for all chromosomes and contigs
cp {params.chrom_sizes} {params.cooler_tempchrsize}
            
# the cload command requires the chrom size file to exist besides the chrom size bin file.
cooler cload pairix \
 -p {params.cooler_n_cores} \
 -s {params.cooler_max_split} \
 {params.cooler_tempchrsize}:{params.cooler_bin_size} \
 {input.frag2_pairs} \
 {output.cooler}
```

## Step-by-step instructions on running Snakemake pipeline:

### 1.  Load slurm and miniconda
Note. The commands to do this will be different on your machine. These commands are specific to an HPC using slurm with these modules installed.

```bash
ml slurm
ml miniconda
```
### 2.  Clone repository
```bash
git clone https://github.com/SansamLab/Process_HiC_SnakeMake.git
# rename folder with project name
mv Process_HiC_SnakeMake/ My_HiC_Project_Folder/
# change directory into root of your project folder
cd My_HiC_Project_Folder
```
### 3.  Start the conda environment
### 3A.  FIRST TIME ONLY:  Setup conda environment
```bash
# -f is the location of the environment .yml file. 
## The relative path assumes that you are in the root directory of this repository.
# -p is the path where you want to install this environment
conda env create -f workflow/envs/HiCSnakemake.yml -p /s/sansam-lab/HiC_Conda_Environment 
```

### 3B.  Activate conda environment
```bash
conda activate /s/sansam-lab/HiC_Conda_Environment
```

### 4.  Modify the job-specific configuration files.
#### 4A.  Modify the config/config.yml file

You must enter paths to the following:
* bwa_genome:
  * location of bwa indexed genome for the alignment
* chrom_sizes
  * chromosome sizes file
* juicer_RE_file
  * restriction enzyme file generated with juicer

#### 4B.  Modify the config/samples.csv file

The samples.csv file in the config folder has paths to the test fastq files. You must replace those paths with those for your own fastq files. The first column of each row is the sample name. This name will be used for all output files.

#### 4C.  IF SLURM RESOURCE CHANGES ARE NEEDED. Modify the config/cluster_config.yml file

CPU and memory requests for each rule in the pipeline are detailed in this file. If you are using SLURM, you may need to alter this file to fit your needs/system.

### 5.  Do a dry run.
A dry run produces a text output showing exactly what commands will be executed. Look this over carefully before submitting the full job. It is normal to see warnings about changes made to the code, input, and params.
```bash
snakemake -npr
```

### 6.  Make a DAG diagram.
```bash
snakemake --dag | dot -Tpdf > dag.pdf
```

### 7.  Run on cluster with slurm.
This snakemake pipeline could be executed without slurm, but if an hpc with slurm is used, the following will start the pipeline with the parameters defined in the config/cluster_config.yml file.
```bash
sbatch --wrap="\
snakemake \
-R \
-j 999 \
--cluster-config config/cluster_config.yml \
--cluster '\
sbatch \
-A {cluster.account} \
-p {cluster.partition} \
--cpus-per-task {cluster.cpus-per-task}  \
--mem {cluster.mem} \
--output {cluster.output} \
--constraint=westmere'"
```

### 8.  Check results, and when finished, exit environment.
The results will be saved to the "results" folder. Look over log files generated in either the logs/ or logs/snakelogs folders (depending on whether slurm was used).
```bash
conda deactivate
```

## References:

Goloborodko, A., Nezar Abdennur, Venev, S., Hbbrandao, & Gfudenberg. (2019). mirnylab/pairtools v0.3.0 (v0.3.0) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.2649383

Abdennur, N., Goloborodko, A., Imakaev, M., Kerpedjiev, P., Fudenberg, G., Oullette, S., Lee, S., Strobelt, H., Gehlenborg, N., & Mirny, L. (2021). open2c/cooler: v0.8.11 (v0.8.11) [Computer software]. Zenodo. https://doi.org/10.5281/ZENODO.4655850

Lee, S., Bakker, C. R., Vitzthum, C., Alver, B. H., & Park, P. J. (2022). Pairs and Pairix: a file format and a tool for efficient storage and retrieval for Hi-C read pairs. In C. Alkan (Ed.), Bioinformatics (Vol. 38, Issue 6, pp. 1729–1731). Oxford University Press (OUP). https://doi.org/10.1093/bioinformatics/btab870

Durand, N. C., Shamim, M. S., Machol, I., Rao, S. S. P., Huntley, M. H., Lander, E. S., & Aiden, E. L. (2016). Juicer provides a one-click system for analyzing loop-resolution hi-c experiments. Cell Systems, 3(1), 95–98. https://doi.org/10.1016/j.cels.2016.07.002

Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., Whitwham, A., Keane, T., McCarthy, S. A., Davies, R. M., & Li, H. (2021). Twelve years of samtools and bcftools. GigaScience, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008

Li, H. (2013). Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. ArXiv:1303.3997 [q-Bio]. http://arxiv.org/abs/1303.3997

Köster, J., & Rahmann, S. (2012). Snakemake—A scalable bioinformatics workflow engine. Bioinformatics (Oxford, England), 28(19), 2520–2522. https://doi.org/10.1093/bioinformatics/bts480

Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc. Retrieved from https://docs.anaconda.com/

