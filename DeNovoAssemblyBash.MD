## This details producing our de novo assembly

First we randomly chose individuals from each treatment to use to build an assembly. We then used the Oyster River Protocol to assemble a transcriptome.

```#!/bin/bash
#SBATCH -N 1
#SBATCH --job-name="pumilio_assembly"
#SBATCH --output=pumilioassemblylog.txt
#SBATCH --open-mode=append
# echo commands to stdout
set -x


DIR=$(pwd)

# concatenate randomly chosen individuals, forward reads
cd raw_data/
zcat GS26_2_D_S478_L007_R1_001.fastq.gz GS26_2_V_S479_L007_R1_001.fastq.gz GS37_2_D_S440_L007_R1_001.fastq.gz GS37_2_V_S456_L007_R1_001.fastq.gz GS42_1_D_S430_L007_R1_001.fastq.gz GS42_1_V_S427_L007_R1_001.fastq.gz GS46_3_D_S449_L007_R1_001.fastq.gz GS46_3_V_S435_L007_R1_001.fastq.gz AD_4_D2_S444_L007_R1_001.fastq.gz AD_4_V2_S434_L007_R1_001.fastq.gz > pumilio4assembly.R1.fastq

# concatenate randomly chosen individuals, forward reads
zcat GS26_2_D_S478_L007_R2_001.fastq.gz GS26_2_V_S479_L007_R2_001.fastq.gz GS37_2_D_S440_L007_R2_001.fastq.gz GS37_2_V_S456_L007_R2_001.fastq.gz GS42_1_D_S430_L007_R2_001.fastq.gz GS42_1_V_S427_L007_R2_001.fastq.gz GS46_3_D_S449_L007_R2_001.fastq.gz GS46_3_V_S435_L007_R2_001.fastq.gz AD_4_D2_S444_L007_R2_001.fastq.gz AD_4_V2_S434_L007_R2_001.fastq.gz > pumilio4assembly.R2.fastq

# move these reads for the assembly
cd $DIR
cp $DIR/raw_data/pumilio4assembly.R1.fastq $DIR/
cp $DIR/raw_data/pumilio4assembly.R2.fastq $DIR/

# transcriptome assembly
source activate orp_v2

/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/oyster.mk main \
MEM=120 \
CPU=24 \
READ1=pumilio4assembly.R1.fastq \
READ2=pumilio4assembly.R2.fastq \
TPM_FILT=1 \
RUNOUT=O.pumilio
```

After this, we annotated the assembly using diamond and the uniprot database.

```
#!/bin/bash
#SBATCH --job-name=annotate
#SBATCH --output=annotate.log
#SBATCH --cpus-per-task=24
#SBATCH --partition=macmanes,shared

source activate orp_v2

diamond blastx -p 24 -e 1e-8 --top 0.1 -q O.pumilio.ORP.fasta -d /mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/software/diamond/swissprot -o pumilio.all.annotation.hits.txt

sort pumilio.all.annotation.hits.txt -k 1,1 -k11,11g | sort -u -k 1,1 --merge > pumilio.annotation.txt
sed "s/|/\t/g" pumilio.annotation.txt > tmp
mv tmp pumilio.annotation.txt
```

We then ran BUSCO and quantified reads for each sample relative to our assembly using salmon.

```
#!/bin/bash
#SBATCH --job-name="busco_quant"
#SBATCH --output=busco_quant_log.txt
#SBATCH --cpus-per-task=24
#SBATCH --open-mode=append
#SBATCH --partition=macmanes,shared
# echo commands to stdout
set -x

# activate OPR environment
source activate orp_v2

# global variables
DIR=$(pwd)
TRANSCRIPTOME=$"O.pumilio.ORP.fasta"
BUSCO=$(which run_BUSCO.py)
TETRAPODA="/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/busco_dbs/tetrapoda_odb9/"


# move transcriptome into the correct location
cp ${DIR}/assemblies/${TRANSCRIPTOME} $DIR

### BUSCO on tetrapoda database ###
python ${BUSCO} -i ${DIR}/${TRANSCRIPTOME} -o pumilio_tetrapoda -l $TETRAPODA -m transcriptome -c 24


### Salmon quantification ###
# list samples
samples=$(ls raw_data/*R1_001.fastq.gz | sed "s/1_001.fastq.gz//g" | grep -v Undetermined)

# do quantification
for sample in $samples
do
/mnt/lustre/macmaneslab/ams1236/Oyster_River_Protocol/quant.mk all \
MEM=500 CPU=24 \
SAMPLE=$sample \
SUFFIX=_001.fastq.gz \
TRANSCRIPTOME=$TRANSCRIPTOME
done

# what is the mapping rate?
for sample in $samples
do
rate=$(grep "Mapping rate" ${DIR}/quants/${sample}/logs/salmon_quant.log | cut -d "=" -f 2 | sed "s/\s//g")
printf '%s \t %s \n' "$sample" "$rate" >> ${DIR}/quants/mappingrates.txt
done

```

