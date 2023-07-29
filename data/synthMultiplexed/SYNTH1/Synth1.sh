

targetReadsPerCell=15000
TSVbarcodesPath="/Sample5/0B_SynthCreate/1_PureLinesAlign/"
sing='singularity exec  singularity.sif'
BamPath="/Sample5/0B_SynthCreate/1_PureLinesAlign/"
tag="SYNTH5"
seed=0
samples=(Sample1 Sample2 Sample3 Sample4 Sample5)



zcat ${TSVbarcodesPath}/Sample_S202*/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | sort | uniq -c | awk '$1 > 1' | awk '{print $2}' > blacklist.txt

for i in "${samples[@]}"
do
zcat ${TSVbarcodesPath}/${i}/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | grep -v -f blacklist.txt | sort -R | head -n 2000 >  $i.goodBarcodes.tsv
done
rm blacklist.txt


for i in "${samples[@]}"
do
$sing \
samtools view -@10 -m3G ${BamPath}/${i}/outs/possorted_genome_bam.bam | fgrep -o -w -f $i.goodBarcodes.tsv | sort | uniq -c  > ${i}.readCount.txt
done


for i in "${samples[@]}"
do
awk '{print $1}' ${i}.readCount.txt | sort -n | awk -f ~/scripts/functions/median.awk |  awk '{printf "%.0f\n", $1}' > ${i}.median
done





for i in "${samples[@]}"
do
  median=`cat ${i}.median`
  if (( ${median} > ${targetReadsPerCell} )); then
    Ratio=`bc -l <<< "((${targetReadsPerCell} / ${median}))"`
    Subsampling=`echo $Ratio | awk -v S = $seed '{print "S"$1}' | awk '{printf "%0.3f\n", $1}'`
    ($sing samtools view -m 4 -@ 15 -H ${BamPath}/${i}/outs/possorted_genome_bam.bam ; \
    $sing samtools view -m 4 -@ 15 ${BamPath}/${i}/outs/possorted_genome_bam.bam | \
    fgrep -f ${i}.goodBarcodes.tsv ) | $sing \
    samtools view -m 4 -@ 15 -h -s ${Subsampling} | $sing samtools view -Sb > ${i}.SS.bam
  else
    ($sing samtools view -m 4 -@ 15 -H ${BamPath}/${i}/outs/possorted_genome_bam.bam ; \
    $sing samtools view -m 4 -@ 15 ${BamPath}/${i}/outs/possorted_genome_bam.bam | \
    fgrep -f ${i}.goodBarcodes.tsv ) | $sing samtools view -Sb > ${i}.SS.bam
  fi
done && rm *readCount.txt && rm *.median



#Create tsv with barcode-id info for demultiplexing
for id in "${samples[@]}"
do
cat $id.goodBarcodes.tsv | awk -v OFS="\t" -v experiment="$id" '{print $1, experiment}' > $id.goodBarcodes.ID.tsv
done

cat *.goodBarcodes.ID.tsv > Barcodes-IDs.tsv && rm *.goodBarcodes.ID.tsv && rm *goodBarcodes.tsv


array=( "${samples[@]/#/-I=}" )
array2=( "${array[@]/%/.SS.bam}" )
function join_by { local IFS="$1"; shift; echo "$*"; }
SSbamList=`join_by , "${array2[@]}" | sed 's/,/ /g'`



$sing gatk MergeSamFiles \
$SSbamList \
-O=SilicoMultiplexed.bam


$sing gatk AddOrReplaceReadGroups \
-I SilicoMultiplexed.bam  \
-O $tag.bam  \
-SO coordinate \
--RGID $tag \
--RGLB $tag \
--RGPL ILLUMINA \
--RGPU machine \
--RGSM $tag


$sing samtools index $tag.bam


$sing samtools view -@10 -m3G SilicoMultiplexed.bam | grep -oP '(?<=CB:Z:)\w+' | sort | uniq -c | awk '{print $1}' | sort | awk -f ~/scripts/functions/median.awk |  awk '{printf "%.0f\n", $1}' > medianReads.SilicoMultiplexed.txt
