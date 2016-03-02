# Genome-SNP
Discovery of SNPs from genome
Tools: bwa-mem, Picard(V1.141), GATK(V3.5), bedtools(V2.17), samtools(V1.2)

**A. Mapping with BWA**

<br>1. Index the reference genome using bwa-index
       bwa index papaya.fasta

<br>2. Align the paired reads to reference genome using bwa mem
       bwa mem -M -t 8 /reference/japonica/reference.fa /reads/filename_1.fq.gz /reads/filename_2.fq.gz 
       > out.sam

<br>3. Sort SAM file and output as BAM file
       java -jar picard.jar SortSam INPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sam 
      OUTPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.bam 
      VALIDATION_STRINGENCY=LENIENT 
      CREATE_INDEX=TRUE
      SO=coordinate 

<br>4. Alignment statistics with samtools flagstat and picard
       samtools flagstat bwa_aln.sorted.bam > aln.stats
       java -Xmx10000m -jar picard-tools-1.58/BamIndexStats.jar I=bwa_aln.sorted.bam > <sorted.stats>

<br>5. FixMate Information
       java -jar picard.jar FixMateInformation 
       INPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.bam 
       OUTPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.bam 
       VALIDATION_STRINGENCY=LENIENT 
       SO=coordinate 
       CREATE_INDEX=TRUE


<br>6. Mark duplicate reads
       java -Xmx8g -jar picard.jar MarkDuplicates 
       INPUT=/output/bwa_aln.sorted.fxmt.bam 
       OUTPUT=/output/bwa_aln.sorted.fxmt.mdup.bam 
       METRICS_FILE=/output/filename.metrics 
       VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=TRUE


<br>7. Add or replace read groups
       java -jar picard.jar AddOrReplaceReadGroups 
       INPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.mdup.bam 
       OUTPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.mdup.arg.bam 
       RGID=MD2 
       PL=illumina 
       SM=MD2sample 
       RGLB=MD2project 
       VALIDATION_STRINGENCY=LENIENT 
       SO=coordinate 
       CREATE_INDEX=TRUE 
       RGPU=none
       ##notes RGID = ReadGroup ID - the Unique ID of the Sample isolate
               RGLB = ReadGroup Library - the name of library (eg. 200bp, RNA)
               RGPL = ReadGroup Platform (e.g. illumina, solid)
               RGPU = ReadGroup Platform Unit (eg. run barcode)
               RGSM = ReadGroup Sample - the Sample name
       
<br>9. Sort Bam file
       java -jar picard.jar BuildBamIndex INPUT=/home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.mdup.addrg.bam
       
<br>10.Create index and dictionary for reference genome
       java -Xmx8g -jar picard.jar CreateSequenceDictionary 
       REFERENCE=/reference/japonica/reference.fa 
       OUTPUT=/reference/reference.dict

<br>11.Realign Target
      java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator 
      -I /home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.mdup.addrg.bam  
      -R /home/cmdv/adra/ananas/SNP/Acomosus_v3.fa 
      -o /home/cmdv/adra/ananas/SNP/bwa.intervals
      -fixMisencodedQuals 
      -nt 16

<br>12.Indel Realignment
       java -jar GenomeAnalysisTK.jar -T IndelRealigner 
       -R /home/cmdv/adra/ananas/SNP/Acomosus_v3.fa 
       -I /home/cmdv/adra/ananas/SNP/bwa_aln.sorted.fxmt.mdup.arg.bam 
       -targetIntervals /home/cmdv/adra/ananas/SNP/bwa.intervals 
       -o bwa_realigned.bam

<br>13.SNP calling
       java -jar GenomeAnalysisTK.jar
       -T HaplotypeCaller -R /home/cmdv/adra/ananas/SNP/Acomosus_v3.fa 
       -I bwa_realigned.bam
       -o gatk.raw.vcf
       -nct 32
       --genotyping_mode DISCOVERY
       -stand_call_conf 30
       -stand_emit_conf 10
       
<br>14.Select Variant
       java -Xmx1g -jar GenomeAnalysisTK.jar
       -T SelectVariants
       -R /home/cmdv/adra/ananas/SNP/Acomosus_v3.fa 
       -V gatk.raw.vcf
       -selectType SNP
       -o snp_gatk_raw.vcf
       
       java -Xmx1g -jar GenomeAnalysisTK.jar
       -T SelectVariants
       -R /home/cmdv/adra/ananas/SNP/Acomosus_v3.fa 
       -V indel_gatk_raw.vcf 
       -selectType INDEL
       -o gatk.snp.raw.vcf


<br>15.Variant Filtration
       java -jar ../../../bioinfo_sware/GATK/GenomeAnalysisTK.jar 
       -T VariantFiltration -R Acomosus_v3.fa 
       -V snp_gatk_raw.vcf -o snp_gatk_hardfilter.vcf
       --filterExpression "QD < 2.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
       --filterName "ananasSNP"

       java -jar ../../../bioinfo_sware/GATK/GenomeAnalysisTK.jar 
       -T VariantFiltration -R Acomosus_v3.fa 
       -V indel_gatk_raw.vcf 
       --filterExpression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" --filterName "indel_gatk" 
       -o indel_gatk_hardfilter.vcf


<br>18. SNP annotation
        java -jar snpEff.jar ann -v Acomosus_v3 ./data/Acomosus_v3/snp_gatk_hardfilter.vcf 
        -s annot_summary.html >  ./data/Acomosus_v3/annot.vcf

