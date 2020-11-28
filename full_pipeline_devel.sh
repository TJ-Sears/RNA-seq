#!/bin/bash
#
#SBATCH --job-name=RNAseq
#
#SBATCH --time=04-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=8G
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tsears
#SBATCH --partition=sfgf

module load java
module load biology samtools
module load pigz
module load biology fastx_toolkit
module load biology salmon
module load R
module load python/3.6
ml pigz

export fastq_dir="/oak/stanford/groups/emoding/sequencing/fastq/hiseq5"
export working_dir="/oak/stanford/groups/emoding/analysis/rnaseq/demux_outs/test5/"
export scripts_dir="/oak/stanford/groups/emoding/scripts/rnaseq_scripts"
export index="/oak/stanford/groups/emoding/analysis/rnaseq/index/SalmonIndex/gencodeIndex27/"
export R1="$fastq_dir/hiseq5_R1.fastq.gz"
export R2="$fastq_dir/hiseq5_R2.fastq.gz"
export txt_dir="/oak/stanford/groups/emoding/analysis/rnaseq/txtfiles"

### Demux CAPPseq outputs
cd $working_dir

#zcat $R1 > ${R1%.*}
#zcat $R2 > ${R2%.*} 

demultiplex demux -m 2 $txt_dir/hiseq5_2.txt ${R1%.*} ${R2%.*}

### Trim 3 bases

for FILE in *R2*; do fastx_trimmer -Q33 -f 4 -i $(echo $FILE) -o $(echo $FILE | sed 's/.fastq/.trimmed.fastq/g'); done

### Run Salmon Quant
# for file in * ; do
#     mv -v "$file" "${file#*_}"
# done

# for fn in *trimmed.fq;
# do
# base=`basename ${fn}`
# samp=$(echo "$base" | cut -f 2 -d '_')
# echo "Processing sample ${samp}"
# salmon quant -i $index -l A \
#          -1 ${samp}.fastq \
#          -2 ${samp}.trimmed.fastq \
#          -p 8 -o ../../quants/${samp}_quant
# done

# #Run tximport
# cd ../../quants/
# Rscript $scripts_dir/tximport-HiSeq1342.R
