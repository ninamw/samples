#!/bin/bash 
#QC
fastqc *fastq.gz 
mkdir fastqc_results
mv *fastqc.* fastqc_results/ 
#Align
for fq in *_PE1*rep*.fastq.gz
do
        name=$(echo $fq | awk -F"/" '{print $NF}' | awk -F"_PE" '{print $1}')
        suffix=$(echo $fq | awk -F"/" '{print $NF}' | awk -F"_rep" '{print $2}')
        end=$(echo $fq | awk -F"_rep" '{print $NF}' | awk -F".fastq" '{print $1}')
        name1=$(echo ${name}_PE1_rep${suffix})
        name2=$(echo ${name}_PE2_rep${suffix})
        #echo $name1
        #echo $name2
	#echo $name
        bowtie2 -p 7 -x /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38 -1 $name1 -2 $name2 | samtools view -bS - -@ 1 | samtools sort -m 2G - -@ 7 -o ${name}_rep${end}.sorted.bam
	samtools rmdup -S ${name}_rep${end}.sorted.bam ${name}_rep${end}.sorted.rmdup.bam
done
#seqOutBias for in silico trimming and to get SNR bigWig and bed
cp /media/ninadmw/sgate4/06052017-T1DGC_ATACSeq-40229196/hg38_38.4.2.2.tbl .
cp /media/ninadmw/sgate4/06052017-T1DGC_ATACSeq-40229196/hg38.tal_38.gtTxt.gz .
for i in *rmdup*bam
do
name=$(echo $i | awk -F"." '{print $1}')
echo $name
echo $i
echo $name.SNR20_130.bigWig
echo $name.SNR20_130.bed
seqOutBias /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.fa $i --no-scale --bw=$name.SNR20_130.bigWig --shift-counts --bed=$name.SNR20_130.bed --pdist=20:130 --only-paired --read-size=38
#sort -k1,1 -k2,2n -o $i $i
echo $name.SNR20_500.bigWig
echo $name.SNR20_500.bed
seqOutBias /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.fa $i --no-scale --bw=$name.SNR20_500.bigWig --shift-counts --bed=$name.SNR20_500.bed --pdist=20:500 --only-paired --read-size=38
done
rm -rf hg38_38.4.2.2.tbl
rm -rf hg38.tal_38.gtTxt.gz

for i in *.bed
do
sort -k1,1 -k2,2n -o $i $i
done

#smooth bigWig and bedGhraphs

for i in *.rmdup.bam
do
name=$(echo $i | awk -F".sorted" '{print $1}')
echo $name
echo $i
awk '{$2-=20}1''{$3+=20}1' $name.SNR20_130.bed| awk '{ print $1"\t"$2"\t"$3"\t"$5 }'> $name.SNR130.bed.smooth
bedClip $name.SNR130.bed.smooth /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes $name.SNR130.bed.new
reads=$(samtools view -c -F 4 $i)
norm=$(echo 10000000/$reads | bc -l | xargs printf "%.*f\n" 3)
echo $reads
echo $norm
genomeCoverageBed -bg -scale $norm -trackline -trackopts name=$name.130.w41.bedGraph -i $name.SNR130.bed.new -g /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes > $name.130.w41.bedGraph
rm -rf $name.SNR130.bed.smooth
rm -rf $name.SNR130.bed.new
gzip $name.130.w41.bedGraph
wigToBigWig $name.130.w41.bedGraph.gz /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes $name.130.w41.bigWig 
awk '{$2-=20}1''{$3+=20}1' $name.SNR20_500.bed| awk '{ print $1"\t"$2"\t"$3"\t"$5 }'> $name.SNR500.bed.smooth
bedClip $name.SNR500.bed.smooth /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes $name.SNR500.bed.new
genomeCoverageBed -bg -scale $norm -trackline -trackopts name=$name.500.w41.bedGraph -i $name.SNR500.bed.new -g /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes > $name.500.w41.bedGraph
rm -rf $name.SNR500.bed.new
rm -rf $name.SNR500.bed.smooth
gzip $name.500.w41.bedGraph
wigToBigWig $name.500.w41.bedGraph.gz /media/ninadmw/OS_Install/ATAC/hg38_indices/hg38.chrom.sizes $name.500.w41.bigWig 
done
dir=$(pwd)
#peak calling with hotspot
foldername=$(date +%Y%m%d)
#foldername=$(date +%Y%m%d%H%M%S)
mkdir -p /home/ninadmw/Documents/"hotspot_$foldername"
cp *.rmdup.bam /home/ninadmw/Documents/"hotspot_$foldername"/
cd /home/ninadmw/Documents/"hotspot_$foldername"
cp /media/ninadmw/OS_Install/ATAC/hotspot_scripts/runhotspot_all_sh.sh .
cp /media/ninadmw/OS_Install/ATAC/hotspot_scripts/token_shell.sh .
chmod +x runhotspot_all_sh.sh 
chmod +x token_shell.sh
./runhotspot_all_sh.sh
chmod +x *runhotspot
./token_shell.sh
for file in *.runhotspot
do
        ./$file
done
rm -rf *tok *txt *.sh *bam *runhotspot
mv ../hs_out .
wc -l *final/*fdr0.01.hot.bed > hotspot.log
head hs_out/*spot.out >> hotspot.log
dir2=$(pwd)
cd ..
mv echo "$dir2" echo "$dir" 2>/dev/null

cd $dir

echo 'name totalreads mappingrate duprate mitorate filteredreads frip'> project.stats
for i in *sorted.bam
do 
  name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".sorted" '{print $1}')
  samtools index $name.sorted.bam
  samtools index $name.sorted.rmdup.bam
  totalreads=$(samtools view -c $name.sorted.bam)  # reads off the sequencer
  mappedtotal=$(samtools view -c -F 4 $name.sorted.bam) # total mapped reads, F (not matching the flag), 4 (read unmapped (0x4)
  mappingrate=$(echo "$totalreads" "$mappedtotal" | awk '{print 100-(($1-$2) * 100/ $1)}'|xargs printf "%.*f\n" 2)
  nonmitoreads=$(samtools view -c -F 4 $name.sorted.bam chrM)
  mitorate=$(echo "$mappedtotal" "$nonmitoreads" | awk '{print 100-(($1-$2) * 100/ $1)}'|xargs printf "%.*f\n" 2)
  mappedtotal2=$(samtools view -c -F 0x904 $name.sorted.bam) # this is to simulate samtools rmdup log. flag '0x904' is for read unmapped (0x4), not primary alignment (0x100) and supplementary alignment (0x800)
  mappeddeduped=$(samtools view -F 0x904 -c $name.sorted.rmdup.bam)
  duprate=$(echo "$mappedtotal2" "$mappeddeduped" | awk '{print ($1-$2) * 100/ $1}'|xargs printf "%.*f\n" 2)
  filteredreads=$(samtools view -c -F 4 $name.sorted.rmdup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY)
  RiP=$(samtools view -F 4 -cL CD19_hotSpot_liftedhg19tohg38.bed  $name.sorted.rmdup.bam)
  totaldeduped=$(samtools view -c $name.sorted.rmdup.bam)
  FRiP=$(echo "$RiP" "$totaldeduped" | awk '{OFS="\t";} {print ($1/$2) * 100}'|xargs printf "%.*f\n" 2)  
echo $name $totalreads $mappingrate $duprate $mitorate $filteredreads $FRiP >> project.stats
done

Rscript plotting.R
mkdir LibraryStats
mv *pdf LibraryStats/ 

for i in *sorted.bam
do 
  name=$(echo $i | awk -F"/" '{print $NF}' | awk -F".sorted" '{print $1}')
samtools view $name.sorted.rmdup.bam | awk '$9>0' | cut -f 9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \t]*//' > Insert_count.txt
  Rscript insertsize.R
  mv insertsize.pdf $name.insertsize.pdf
  rm -rf Insert_count.txt insertsize.pdf
done
mkdir LibraryStats/Insertsize
mv *insertsize.pdf LibraryStats/Insertsize/ 

#Align

#samtools view -F 4 -cL CD19_hotSpot_liftedhg19tohg38.bed  nmw_test_rep1.sorted.rmdup.bam

#samtools view -c -F 4 nmw_test_rep1.sorted.rmdup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY

#samtools view -c nmw_test_rep1.sorted.rmdup.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
