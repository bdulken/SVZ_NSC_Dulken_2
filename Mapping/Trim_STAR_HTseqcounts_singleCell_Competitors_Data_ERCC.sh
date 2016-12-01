#!/bin/bash

#This program takes as input a folder containing unzipped fastq files and performs
#trimming (trim_galore), alignment (STAR), and post processing of the aligned reads
#(Convert to bam and sort). The resulting sorted bam file can then be used with cufflinks
#to determine differential expression.

#Input management
	if [[ "$#" -lt 1 ]]
	then
		echo "$(basename $0) [Dir] "  1>&2
		echo "   [Dir]: directory containing files" 1>&2
		exit 1
	fi
	
Dir=$(echo $1 | sed 's:/$::g')

PICARD_PATH="/Volumes/guacamole/Software/ATAC_Analysis/picard-tools-1.74"
ANNOTATE_PATH="/Volumes/guacamole/Software/Annotations/mouseRefFlat.txt"	

# Makes output directories to write the trimmed fastq files and star aligned bam files
#into.

	[[ ! -d "${Dir}/STAR" ]] && mkdir "${Dir}/STAR"
	[[ ! -d "${Dir}/TRIMMED" ]] && mkdir "${Dir}/TRIMMED"
	
	source ~/.bash_profile

#1. run quality trimming
#Adjust trim_galore for read length
	filePath="${Dir}"
	outPath="${Dir}/TRIMMED"
	for f in $(find "$filePath" -name '*_1.fastq')
	do
		base2=$(basename "${f}" | sed 's/1\.fastq/2\.fastq/g');
		f2="${Dir}/${base2}"
		
		trim_galore -q 15 --phred33 --paired --length 45 -a CTGTCTCTTATACACATCT -a2 CTGTCTCTTATACACATCT --stringency 3 $f $f2
		mv *_val_*.fq $outPath
		mv *txt $outPath
		gzip $f
		gzip $f2
	done
	echo "Finished quality trimming\n"


#2. run STAR mapping into BAM format
	filePath="${outPath}"
	for f in $(find "$filePath" -name '*_val_1.fq')
	do
		file2=$(basename "${f}" | sed 's/1_val_1\.fq/2_val_2\.fq/g');
		f2="${filePath}/${file2}"
	
		folderName=$(basename "${f}" | sed 's/1_val_1\.fq//g');
		oPath="${Dir}/STAR"
		oFname="${oPath}/${folderName}"
		
		mkdir "$oFname"
		cd "$oFname"
		/Volumes/guacamole/Software/STAR_2.3.0/STAR --genomeDir /Volumes/guacamole/Software/STAR_2.3.0/Genome_annotated_ERCC --readFilesIn $f $f2 --outSAMstrandField intronMotif --runThreadN 4	
			
		oFname2="${oFname}/Aligned.out.sam"
		fileName=$(basename "${f}" | sed 's/_val_1\.fq/_STAR\.sam/g');
		oFname3="${oPath}/${fileName}"
	
		ln -s $oFname2 $oFname3
		rm -f $f
		rm -f $f2
	done
	echo "Finished STAR Mapping\n"

#3. convert to BAM
	echo "Now converting to BAM" 
	
	mkdir "${Dir}/STAR/BAM"
	mkdir "${Dir}/STAR/Sorted"
	
	cd "${Dir}/STAR/BAM"
	newDir="${Dir}/STAR"
	for f in $(find "${newDir}" -name "*.sam" -maxdepth 1)
	do
		fileName=$(basename "${f}" | sed 's/1_STAR\.sam/_STAR\.bam/g')
		samtools view -b -S $f > $fileName
	done
	
	echo "Done converting to bam"
	
	for f in $(find "${newDir}" -name "*.sam")
	do
		rm -f $f
	done
	
	echo "Sam files deleted"

#4. Sort Bam file
	filePath2="${newDir}/BAM"
	cd "${Dir}/STAR/Sorted"
	
	for f in $(find "$filePath2" -name "*.bam")
	do
		fileName=$(basename "${f}" | sed 's/\.bam/_sorted/g')
		samtools sort $f $fileName
		rm -f $f
	done
	
	echo "Done sorting bam"

#5. Count reads
ANNOTATIONS='/Volumes/guacamole/Software/Annotations/Katje_mm9_annotations_ERCC92.gtf'
BAM_DIR="${Dir}/STAR/Sorted"

[[ ! -d "${BAM_DIR}/HTseqcounts" ]] && mkdir "${BAM_DIR}/HTseqcounts"

for RAW_BAM_FILE in $(find "$BAM_DIR" -name '*.bam')
do
	OFPREFIX=$(basename "${RAW_BAM_FILE}" | sed 's/\.bam//g')  
	COUNT_BAM="${BAM_DIR}/HTseqcounts/${OFPREFIX}.counts.txt" 

    python -m HTSeq.scripts.count -f bam -r pos -s no ${RAW_BAM_FILE} ${ANNOTATIONS} > ${COUNT_BAM}
    gzip ${RAW_BAM_FILE}
done
	
echo "Done determining counts"

[[ ! -d "${BAM_DIR}/Duplication" ]] && mkdir "${BAM_DIR}/Duplication"
[[ ! -d "${BAM_DIR}/InsertSize" ]] && mkdir "${BAM_DIR}/InsertSize"
[[ ! -d "${BAM_DIR}/RNAseqMetrics" ]] && mkdir "${BAM_DIR}/RNAseqMetrics"

for RAW_BAM_FILE in $(find "$BAM_DIR" -name '*.bam')
do
	OFPREFIX=$(basename "${RAW_BAM_FILE}" | sed 's/\.bam//g')    
	
	FILT_BAM_PREFIX="${OFPREFIX}.filt.srt" # <name>.filt.srt
    TMP_FILT_BAM_FILE="${BAM_DIR}/Duplication/${FILT_BAM_PREFIX}.dupmark.bam" # <name>.filt.srt.dupmark.bam
    DUP_FILE_QC="${BAM_DIR}/Duplication/${FILT_BAM_PREFIX}.dup.qc" # QC file <name>.filt.srt.sup.qc

    java -Xmx2g -jar "$PICARD_PATH/MarkDuplicates.jar" INPUT=${RAW_BAM_FILE} OUTPUT=${TMP_FILT_BAM_FILE} METRICS_FILE=${DUP_FILE_QC} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false
	rm $TMP_FILT_BAM_FILE
	
    INSERT_TXT_FILE="${BAM_DIR}/InsertSize/${FILT_BAM_PREFIX}.insertSizes.txt"
	INSERT_HISTO_FILE="${BAM_DIR}/InsertSize/${FILT_BAM_PREFIX}.insertSizes.pdf"
	
	java -Xmx2g -jar "$PICARD_PATH/CollectInsertSizeMetrics.jar" METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=${INSERT_TXT_FILE} HISTOGRAM_FILE=${INSERT_HISTO_FILE} INPUT=${RAW_BAM_FILE}

	RNAMETRICS_TXT_FILE="${BAM_DIR}/RNAseqMetrics/${FILT_BAM_PREFIX}.RNAseqMetrics.txt"
	RNA_PLOT_FILE="${BAM_DIR}/RNAseqMetrics/${FILT_BAM_PREFIX}.RNAseqMETRICS_plot.pdf"

	java -Xmx2g -jar "$PICARD_PATH/CollectRnaSeqMetrics.jar" STRAND_SPECIFICITY=NONE METRIC_ACCUMULATION_LEVEL=ALL_READS REF_FLAT=${ANNOTATE_PATH} OUTPUT=${RNAMETRICS_TXT_FILE} CHART_OUTPUT=${RNA_PLOT_FILE} INPUT=${RAW_BAM_FILE}
done