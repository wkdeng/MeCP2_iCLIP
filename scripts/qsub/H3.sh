mkdir -p /home/dengw1/workspace/MeCP2_iCLIP2/cutadapt
cutadapt -a AGATCGGAAGAGCGGTTCAG --quality-cutoff 6 -m 15 -o /home/dengw1/workspace/MeCP2_iCLIP2/cutadapt/H3.fq /home/dengw1/workspace/MeCP2_iCLIP2/raw_fq/H3.fq
fastqc -o /home/dengw1/workspace/MeCP2_iCLIP2/cutadapt/ /home/dengw1/workspace/MeCP2_iCLIP2/cutadapt/H3.fq
mkdir -p /home/dengw1/workspace/MeCP2_iCLIP2/star/H3/
STAR --genomeDir /home/dengw1/workspace/genome/star/hg19 --readFilesIn /home/dengw1/workspace/MeCP2_iCLIP2/cutadapt/H3.fq  --outSAMtype BAM Unsorted --outFileNamePrefix /home/dengw1/workspace/MeCP2_iCLIP2/star/H3/ --outFilterMultimapNmax 100 --runThreadN 20 --alignEndsProtrude 15 ConcordantPair --twopassMode Basic --limitOutSJcollapsed 2000000 --outStd Log > /home/dengw1/workspace/MeCP2_iCLIP2/star/H3/log.txt 2>&1
bedtools intersect -f 0.90 -abam /home/dengw1/workspace/MeCP2_iCLIP2/star/H3/Aligned.out.bam -b /home/dengw1/workspace/genome/raw/hg19/rRNA_tRNA.bed -v > /home/dengw1/workspace/MeCP2_iCLIP2/star/H3/Aligned.out.mask_rRNA.bam
