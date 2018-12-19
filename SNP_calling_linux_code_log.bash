#Linux and R code log for Schnepp et al paper

#using the hg38 human build from the NCBI

##Alignment to hg38 genome build

java -jar $PICARD_JARS/picard.jar FastqToSam \
 FASTQ=combined.fastq.gz O=unmapped.bam \
 SAMPLE_NAME=pseudo-bulk QUALITY_FORMAT=Standard \
 READ_GROUP_NAME=pseudo-bulk LIBRARY_NAME=Illumina_Nextra_HiSeq PLATFORM=Illumina PLATFORM_UNIT=unit1
java -jar $PICARD_JARS/picard.jar MarkIlluminaAdapters INPUT=unmapped.bam OUTPUT=markadapters.bam M=markadapters_metrics.txt ADAPTERS=FLUIDIGM
java -jar $PICARD_JARS/picard.jar SamToFastq INPUT=markadapters.bam FASTQ=foruse.fastq CLIPPING_ATTRIBUTE=XT CLIPPING_ACTION=2 NON_PF=true
bwa mem -t 8 ../BWAIndex/genome.fa foruse.fastq > aligned_reads.sam
java -jar $PICARD_JARS/picard.jar SortSam INPUT=aligned_reads.sam OUTPUT=align_sorted.bam SORT_ORDER=queryname
java -jar $PICARD_JARS/picard.jar SortSam INPUT=unmapped.bam OUTPUT=unmapped_sorted.bam SORT_ORDER=queryname
java -jar $PICARD_JARS/picard.jar MergeBamAlignment UNMAPPED_BAM=unmapped_sorted.bam ALIGNED_BAM=align_sorted.bam OUTPUT=aligned_reads.merged.bam REFERENCE_SEQUENCE=../../BWAIndex/genome.fa PAIRED_RUN=true CREATE_INDEX=true ADD_MATE_CIGAR=true MAX_INSERTIONS_OR_DELETIONS=-1 PRIMARY_ALIGNMENT_STRATEGY=MostDistant ATTRIBUTES_TO_RETAIN=XS
java -jar $PICARD_JARS/picard.jar SortSam INPUT=aligned_reads.merged.bam OUTPUT=aligned_sorted.bam SORT_ORDER=coordinate
java -jar $PICARD_JARS/picard.jar MarkDuplicates INPUT=aligned_sorted.bam OUTPUT=dedup_reads.bam METRICS_FILE=metrics.txt 

java -jar $PICARD_JARS/picard.jar ValidateSamFile I=dedup_reads.bam MODE=SUMMARY
java -jar $PICARD_JARS/picard.jar BuildBamIndex INPUT=dedup_reads.bam
java -jar $PICARD_JARS/picard.jar SortVcf I=mgp.v5.merged.snps_all.dbSNP142.sorted.vcf O=mgp.v5.merged.snps_all.dbSNP142.updated.vcf SEQUENCE_DICTIONARY=genome.dict
java -jar $PICARD_JARS/picard.jar UpdateVcfSequenceDictionary INPUT=mgp.v5.merged.snps_all.dbSNP142.vcf OUTPUT=mgp.v5.merged.snps_all.dbSNP142.sorted.vcf SEQUENCE_DICTIONARY=genome.dict
java -jar $PICARD_JARS/picard.jar ReorderSam I=dedup_reads.bam O=dedup_reads.sorted.bam R=../genome.fa CREATE_INDEX=TRUE

#GATK specific code

java -jar $GATK_JARS/GenomeAnalysisTK.jar -T BaseRecalibrator -R ../../BWAIndex/genome.fa -I dedup_reads.bam -knownSites mouse.dbsnp.vcf -o recal_data.table
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T PrintReads -R ../../BWAIndex/genome.fa -I dedup_reads.bam -BQSR recal_data.table -o recal_reads.bam
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../../BWAIndex/genome.fa -I recal_reads.bam --genotyping_mode DISCOVERY -stand_emit_conf 10 -stand_call_conf 30 -o raw_variants.vcf
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T GenotypeGVCFs -R genome.fa -V bulk/raw_variants.g.vcf -V ERR1211176/raw_variants.g.vcf -V ERR1211177/raw_variants.g.vcf -V ERR1211178/raw_variants.g.vcf -V ERR1211179/raw_variants.g.vcf -V ERR1211180/raw_variants.g.vcf -V ERR1211181/raw_variants.g.vcf -V ERR1211182/raw_variants.g.vcf -V ERR1211183/raw_variants.g.vcf -V ERR1211184/raw_variants.g.vcf -V ERR1211185/raw_variants.g.vcf -V ERR1211186/raw_variants.g.vcf -V ERR1211187/raw_variants.g.vcf -V ERR1211188/raw_variants.g.vcf -V ERR1211189/raw_variants.g.vcf -V ERR1211190/raw_variants.g.vcf -V ERR1211191/raw_variants.g.vcf -V ERR1211192/raw_variants.g.vcf -V ERR1211193/raw_variants.g.vcf -V ERR1211194/raw_variants.g.vcf -V ERR1211195/raw_variants.g.vcf -V ERR1211196/raw_variants.g.vcf -V ERR1211197/raw_variants.g.vcf -V ERR1211198/raw_variants.g.vcf -V ERR1211199/raw_variants.g.vcf -V ERR1211200/raw_variants.g.vcf -V ERR1211201/raw_variants.g.vcf -V ERR1211202/raw_variants.g.vcf -V ERR1211203/raw_variants.g.vcf -V ERR1211204/raw_variants.g.vcf -V ERR1211205/raw_variants.g.vcf -V ERR1211206/raw_variants.g.vcf -V ERR1211207/raw_variants.g.vcf -V ERR1211208/raw_variants.g.vcf -V ERR1211209/raw_variants.g.vcf -V ERR1211210/raw_variants.g.vcf -V ERR1211211/raw_variants.g.vcf -V ERR1211212/raw_variants.g.vcf -V ERR1211213/raw_variants.g.vcf -V ERR1211214/raw_variants.g.vcf -V ERR1211215/raw_variants.g.vcf -V ERR1211216/raw_variants.g.vcf -V ERR1211217/raw_variants.g.vcf -V ERR1211218/raw_variants.g.vcf -V ERR1211219/raw_variants.g.vcf -V ERR1211220/raw_variants.g.vcf -V ERR1211221/raw_variants.g.vcf -V ERR1211222/raw_variants.g.vcf -V ERR1211223/raw_variants.g.vcf -V ERR1211224/raw_variants.g.vcf -V ERR1211225/raw_variants.g.vcf -V ERR1211226/raw_variants.g.vcf -V ERR1211227/raw_variants.g.vcf -V ERR1211228/raw_variants.g.vcf -V ERR1211229/raw_variants.g.vcf -V ERR1211230/raw_variants.g.vcf -V ERR1211231/raw_variants.g.vcf -V ERR1211232/raw_variants.g.vcf -V ERR1211233/raw_variants.g.vcf -V ERR1211234/raw_variants.g.vcf -V ERR1211235/raw_variants.g.vcf -V ERR1211236/raw_variants.g.vcf -V ERR1211237/raw_variants.g.vcf -V ERR1211238/raw_variants.g.vcf -V ERR1211239/raw_variants.g.vcf -V ERR1211240/raw_variants.g.vcf -o combined_variants.g.vcf
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T ValidateVariants -R ../BWAIndex/genome.fa -V combined_variants.g.vcf --validateGVCF
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T VariantRecalibrator -R ../../BWAIndex/genome.fa -input bulk_one_sample.g.vcf -resource:hapmap,known=false,training=true,truth=true,prior=15.0 ../../hapmap_3.3.hg38.vcf -resource:1000G,known=false,training=true,truth=false,prior=10.0 ../../V1000G_phase1.snps.high_confidence.hg38.vcf -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ../../dbsnp_146.hg38.vcf -an DP -an QD -an FS -an SOR -an MQ -mode SNP -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -rscriptFile recalibrate_SNP_plots.READ_GROUP_NAME
java -jar $GATK_JARS/GenomeAnalysisTK.jar -T ApplyRecalibration -R ../../BWAIndex/genome.fa -input bulk_one_sample.g.vcf -mode SNP --ts_filter_level 99.0 -recalFile recalibrate_SNP.recal -tranchesFile recalibrate_SNP.tranches -o recalibrated_snps_raw_indels.vcf

#Monovar specific code

nano filenames.txt

#start of filenames file
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211176/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211177/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211178/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211179/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211180/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211181/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211182/dedup_reads.bam
/nfs/turbo/etkeller-fluxod/pschnepp/SNP_calling/ERR1211183/dedup_reads.bam
#end of filenames file

java -jar $PICARD_JARS/SortSam.jar INPUT=star_Aligned.out.bam OUTPUT=reduced_aligned_sorted.bam SORT_ORDER=coordinate
samtools mpileup -BQ0 -d10000 -q 40 -f /home/pschnepp/turbo/Libraries/Homo_sapiens/UCSC/hg38/Sequence/StarIndex/genome.fa -b filenames.txt | python /home/pschnepp/Programs/monovar/src/monovar.py -f /home/pschnepp/turbo/Libraries/Homo_sapiens/UCSC/hg38/Sequence/StarIndex/genome.fa -p 0.002 -a 0.2 -t 0.05 -m 2 -b filenames.txt -o NA19098-r3_bulk_star_monovar.vcf


##Annovar to identify the genomic locations and overlapping SNPs

awk 'length($4) == 1' recalibrated_snps_raw_indels.vcf | awk 'length($5) == 1' > inter.vcf
head -n 20 recalibrated_snps_raw_indels.vcf | cat - inter.vcf > sample_clean.vcf

convert2annovar.pl -format vcf sample_clean.vcf -outfile sample_clean.avinput
annotate_variation.pl -dbtype refGene -buildver hg38 sample_clean.avinput ../../../../humandb/

