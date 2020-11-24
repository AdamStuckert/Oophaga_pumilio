### Mapping to a scaffolded genome assembly

We additionally examined gene expression using an existing genomic resource. [Rodriguez et al.](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-020-6719-5) scaffolded the [Rogers et al. genome assembly](https://academic.oup.com/mbe/article/35/12/2913/5106668). We used the predicted transcriptome from their study in our analyses.

```
#!/bin/bash
#SBATCH --job-name=salmonRod
#SBATCH --output=SalmonRod.txt
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=append
#SBATCH --partition=macmanes,shared
# echo commands to stdout
set -x

# activate OPR environment
#source activate orp_v2

# global variables
DIR=$(pwd)
TRANSCRIPTOME=$"Rodriguez_genomic_resources/SuperDuper.fasta"
assembly=$(basename $TRANSCRIPTOME)


## salmon index
salmon index -t $TRANSCRIPTOME -i $assembly.idx


### Salmon quantification ###
# list samples
samples=$(ls raw_data/*R1_001.fastq.gz | sed "s/1_001.fastq.gz//g" | grep -v Undetermined)

# do quantification
for sample in $samples
do
echo quantifying $sample
salmon quant -i $assembly.idx  \
        --no-version-check --seqBias --gcBias -l a \
         -1 ${sample}1_001.fastq.gz \
         -2 ${sample}2_001.fastq.gz \
         -p 24 -o quants/${sample}


done

# what is the mapping rate?
for sample in $samples
do
rate=$(grep "Mapping rate" ${DIR}/quants/${sample}/logs/salmon_quant.log | cut -d "=" -f 2 | sed "s/\s//g")
printf '%s \t %s \n' "$sample" "$rate" >> ${DIR}/quants/Rodriguezmappingrates.txt
done

```
