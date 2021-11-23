#java -Xmx40G -Djava.io.tmpdir=./tmp -jar /Software/Packages/picard.jar MarkDuplicates \
# REMOVE_DUPLICATES=true \
# I=/nfs01data1/Groups/Liuningning/Base_marrow_project/Whole_genome_sequencing/SC_WGS_DATA/1_Month_WT_Terc/H101SC20041250_20200520/RawData/Cleandata/processing/MDA-Terc-1.bam \
# O=MDA-Terc-1.dedup.bam \
# METRICS_FILE=MDA-Terc-1.dedup.metrics.txt \
# PROGRAM_RECORD_ID= MarkDuplicates PROGRAM_GROUP_VERSION=null \
# PROGRAM_GROUP_NAME=MarkDuplicates

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-Terc-1.dedup.bam

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-Terc-3.dedup.bam

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-Terc-4.dedup.bam

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-WT-2.dedup.bam

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-WT-4.dedup.bam

java -jar /Software/Packages/picard.jar BuildBamIndex \
 I=MDA-WT-5.dedup.bam
