### Mapping to a scaffolded genome assembly

We additionally examined gene expression using an existing genomic resource. [Rodriguez et al.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6719-5) scaffolded the Rogers et al. genome assembly. We used the predicted transcriptome from their study in our analyses.

```
#!/bin/bash
#SBATCH --job-name="kallisto_quant"
#SBATCH --output=kallisto_log.txt
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=append
#SBATCH --partition=macmanes,shared
# echo commands to stdout
set -x

module purge

conda activate base

# global variables
DIR=$(pwd)
TRANSCRIPTOME=$"Rodriguez_genomic_resources/SuperDuper.fasta"

### Kallisto quantification ###
# list samples
samples=$(ls raw_data/*R1_001.fastq.gz | sed "s/1_001.fastq.gz//g" | grep -v Undetermined)

# index transcriptome
kallisto index -i pumilioRodriguez.idx $TRANSCRIPTOME

# make directories
mkdir kallisto_quants_Rodriguez
cd kallisto_quants_Rodriguez
mkdir raw_data/
for i in $samples; do mkdir $i; done
cd ..

# perform pseudo-quantification
for sample in $samples
do
kallisto quant -i pumilioRodriguez.idx -t 24 -o kallisto_quants_Rodriguez/$sample -b 100 \
${sample}1_001.fastq.gz ${sample}2_001.fastq.gz
done
```
