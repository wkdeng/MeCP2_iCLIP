## Parameter setting, please edit it to your own path
script_folder=~/snakerun/MeCP2_iCLIP/scripts
data_folder=~/workspace/MeCP2_iCLIP
genome_folder=~/workspace/genome
## Demultiplex
python $script_folder/demultiplex_read1.py $data_folder/raw_fq/Undetermined_S0_L003_R1_001.fastq $data_folder/raw_fq

## Cut adapter, QC, genome mapping, rRNA/tRNA masking, deduplication, 
mkdir -p $data_folder/cutadapt

for group in M1 M2 M3 H1 H2 H3
do
genome=hg19
if [$group=M*]
then
    genome=mm10
fi
echo $group
# Cut adapter
mkdir -p $data_folder/star/$group/
cutadapt -a AGATCGGAAGAGCGGTTCAG --quality-cutoff 6 -m 15 -o $data_folder/cutadapt/$group.fq $data_folder/raw_fq/$group.fq
# QC
fastqc -o $data_folder/MeCP2_iCLIP/cutadapt/ $data_folder/cutadapt/$group.fq
# STAR mapping
STAR --genomeDir $genome_folder/star/$genome \
--readFilesIn $data_folder/cutadapt/$group.fq  --outSAMtype BAM Unsorted \
--outFileNamePrefix $data_folder/star/$group/ \
--outFilterMultimapNmax 100 \
--runThreadN 20 \
--alignEndsProtrude 15 ConcordantPair \
--twopassMode Basic --limitOutSJcollapsed 2000000 \
--outStd Log > $data_folder/star/$group/log.txt 2>&1

# Mapping statistic
python $data_folder/mapping_stat.py $data_folder/star/$group/Aligned.out.bam $data_folder/cutadapt/${group}_fastqc.zip > $data_folder/star/$group/mapping_stat.txt

# rRNA/tRNA masking
bedtools intersect -f 0.90 -abam $data_folder/star/$group/Aligned.out.bam -b $genome_folder/raw/$genome/rRNA_tRNA.bed -v > $data_folder/star/$group/Aligned.out.mask_rRNA.bam

## PCR duplication removal
python $script_folder/collapse_dup.py -b $data_folder/star/$group/Aligned.out.mask_rRNA.bam -o $data_folder/star/$group/Aligned.out.mask_rRNA.collapsed_dup.bam -m $data_folder/star/$group/dup_removal.txt

# CLAM preparation 
mkdir -p $data_folder/clam/$group
CLAM preprocessor -i $data_folder/star/$group/Aligned.out.mask_rRNA.collapsed_dup.bam -o $data_folder/clam/$group --read-tagger-method start --lib-type sense

    # To save some time, binw of 100 and 150 can be removed as 200bp peak width gave us better results
    for binw in '100' '150' '200'
    do
    mkdir -p $data_folder/clam/${group}_${binw}
    mkdir -p $data_folder/clam/peaks_${binw}/$group
    mkdir -p $data_folder/homer/$binw/$group
    mkdir -p $data_folder/bigwig/${binw}

    # CLAM realign, test bin width of 100, 150, 200 
    CLAM realigner -i $data_folder/clam/$group -o $data_folder/clam/${group}_${binw} --winsize $binw --max-tags -1 --read-tagger-method start --lib-type sense

    python $script_folder/make_bw.py --ub $data_folder/clam/${group}/unique.sorted.bam --rb \
    $data_folder/clam/${group}_${binw}/realigned.sorted.bam -g $genome_folder/raw/hg19/hg19.gtf \
    -u -t 20 -o $data_folder/bigwig/${binw} -s ${group} -c ${genome_folder}/raw/${genome}/chr_size.txt

    # Peakcalling, using permutation instead of negative binomial test for a better result without control
    CLAM permutation_callpeak -i $data_folder/clam/${group}/unique.sorted.bam $data_folder/clam/${group}_${binw}/realigned.sorted.bam \
    -o $data_folder/clam/peaks_${binw}/${group} -p 30 \
    --gtf $genome_folder/raw/${genome}/${genome}.gtf --qval-cutoff 1 --merge-size ${binw}

    # Filter peaks by qval<=0.005
    mv $data_folder/clam/peaks_${binw}/${group}/narrow_peak.permutation.bed $data_folder/clam/peaks_${binw}/${group}/narrow_peak.permutation.bed.all
    awk '{split($4,a,":");if(a[3]<=0.005){print}}' $data_folder/clam/peaks_${binw}/${group}/narrow_peak.permutation.bed.all > $data_folder/clam/peaks_${binw}/${group}/narrow_peak.permutation.bed

    # Search for motif
    findMotifsGenome.pl $data_folder/clam/peaks_${binw}/${group}/narrow_peak.permutation.bed ${genome} $data_folder/homer/${binw}/${group} \
    -rna -len 5,6,7 \
    -p 20 -size 100 -S 10
    done
done
