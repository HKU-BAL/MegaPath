#!/usr/bin/env bash
set -e
set -m
set -o pipefail

SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})

# configurations
MEGAHIT=megahit
SAMTOOLS=samtools
BEDTOOLS=bedtools

BBMAP_DIR=${SCRIPT_PATH}/bbmap
BBDUK=${BBMAP_DIR}/bbduk2.sh

DB=${SCRIPT_PATH}/db
SCRIPT_DIR=${SCRIPT_PATH}/scripts

export PATH=${SCRIPT_PATH}:${PATH}

# default parameters
PREFIX="megapath-amplicon"
THREADS=45
BBDUK_MAQ=10 #default 10
BBDUK_TRIMQ=10 #default
MIN_LEN=150
ENTORPY=0.75
AS_OVER_LEN_RATIO_THRES=1
MAPQ_TO_TARGET_THRES=10
AS_TO_TARGET_THRES=150

while getopts "p:1:2:t:d:" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		1) READ1=${OPTARG};;
		2) READ2=${OPTARG};;
		t) THREADS=${OPTARG};; 
		d) DB=${OPTARG};;
		*) exit 1;;
	esac
done


if [ -z "${READ1}" ] || [ -z "${READ2}" ]; then
   echo "Usage: $0 -1 <read1.fq> -2 <read2.fq> [options]"
   echo "    -p  output prefix [megapath-amplicon]"
   echo "    -t  number of threads [45]"
   echo "    -L  max read length [250]"
   echo "    -d  database directory [${SCRIPT_PATH}/db]"
   exit 1
fi

TARGET_IDX=$DB/amplicon/Mycobacterium_tuberculosis_H37Rv_genome_v3.fasta
HUMAN_IDX=$DB/amplicon/GCF_000001405.39_GRCh38.p13_genomic.fna.gz
ORAL_MICROBIOME_IDX=$DB/amplicon/ORAL_MICROBIOME_genomic.no_Myc.fna
REFSEQ_IDX=$DB/amplicon/refseq.fna.no_blastid85_tb.nomyc.onlytb.gz
AMP_REGION_BED=$DB/amplicon/MTBC_H37RV_19_target_regions_V2.std.bed
TARGET_SEQ_ID=`head -1 $AMP_REGION_BED|cut -f1`

run_bwa_mem(){
if [ ! -f $1.bwt ];then
    bwa index $1
fi
BWA_MEM_OUT=`basename $2`
bwa mem $1 $2 $3 -t $THREADS $4| samtools view -bh | samtools sort > $BWA_MEM_OUT.bam
samtools index $BWA_MEM_OUT.bam
}

run_minimap2(){
MINIMAP2_PREFIX=`basename $2`
minimap2 -ax map-ont $1 $2 | samtools sort -o $MINIMAP2_PREFIX.bam 
samtools index $MINIMAP2_PREFIX.bam
}

AS_over_len_ratio_filter(){
run_bwa_mem $1 $2 $3 
FILTER_PREFIX=`basename ${2}`
python $SCRIPT_DIR/amplicon/filter_bam.py --AS_and_AS-LEN_ratio_filter $FILTER_PREFIX.bam --AS_LEN_ratio $AS_OVER_LEN_RATIO_THRES --AS_threshold 0
python $SCRIPT_DIR/amplicon/process_fastx_by_readid.py --filter ${FILTER_PREFIX}.bam.AS*_AS-LEN_ratio*.bam.ASandAS-LEN_ratio_filter.readid $2 $3 --suffix $4
}


# 0. preprocessing
if [ -e ${PREFIX}.prep.done ]; then
	echo "Skipping quality control";
else
	echo "[TIMESTAMP] $(date) Running BBDuk to preprocess..."
	STARTTIME=$(date +%s)
	BBDUK_THREADS=$((${THREADS} / 2))
	if [ ${BBDUK_THREADS} -lt 1 ]; then
		BBDUK_THREADS=1
	fi

	${BBDUK} kmask=N qtrim=rl maq=${BBDUK_MAQ} trimq=${BBDUK_TRIMQ} threads=${BBDUK_THREADS} minlength=${MIN_LEN} in=${READ1} in2=${READ2} out=stdout.fq ref=${BBMAP_DIR}/resources/adapters.fa hdist=1  | ${BBDUK} entropy=${ENTORPY} in=stdin.fq outm=${PREFIX}.low_compl.fq.gz threads=${BBDUK_THREADS} out=${PREFIX}.bbduk_1.fq.gz out2=${PREFIX}.bbduk_2.fq.gz overwrite=t

	echo "[TIMESTAMP] $(date) Running BBDuk to preprocess... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] BBDuk took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.prep.done
fi

# 1. assembly filtering
if [ -e ${PREFIX}.assem.done ]; then
    echo "Skipping assembly filter";
else
    echo "[TIMESTAMP] $(date) Running assembly filter..."
    STARTTIME=$(date +%s)
    ASSEM_INPUT_1=${PREFIX}.bbduk_1.fq.gz
    ASSEM_INPUT_2=${PREFIX}.bbduk_2.fq.gz
	OUT_DIR=assem_$ASSEM_INPUT_1
    mkdir -p $OUT_DIR
    cd $OUT_DIR

    run_bwa_mem $TARGET_IDX ../$ASSEM_INPUT_1 ../$ASSEM_INPUT_2 -a

    FULL_BAM=$ASSEM_INPUT_1.bam
    export FULL_BAM TARGET_SEQ_ID TARGET_IDX THREADS SCRIPT_DIR
	export -f run_bwa_mem run_minimap2

    cut -f2,3 $AMP_REGION_BED |parallel --colsep '\t' 'samtools view -h $FULL_BAM $TARGET_SEQ_ID:{1}-{2} | samtools sort -n |samtools bam2fq - -1 $FULL_BAM.{1}-{2}.1.fq -2 $FULL_BAM.{1}-{2}.2.fq -0 /dev/null -s /dev/null -n'
    parallel -n2 "megahit -t 4 --presets meta-sensitive -1 {1} -2 {2} -o {1/}.megahit; \
        cd {1/}.megahit/; \
        run_minimap2 $TARGET_IDX final.contigs.fa; \
        run_bwa_mem final.contigs.fa ../{1} ../{2} -a" ::: *.fq || true
    cut -f2,3 $AMP_REGION_BED |parallel --colsep '\t' "cd $ASSEM_INPUT_1.bam.{1}-{2}.1.fq.megahit;\
        python $SCRIPT_DIR/amplicon/filter_contigs.py final.contigs.fa.bam {1} {2} ../*{1}*.fq"
    sort */assembly_filter.retain_readid |uniq> assembly_filter.retain_readid.sortuniq
    parallel 'seqtk subseq {} assembly_filter.retain_readid.sortuniq > {}.assembly_filtered' ::: ../$ASSEM_INPUT_1 ../$ASSEM_INPUT_2
    cd ..

    echo "[TIMESTAMP] $(date) Running assembly filter... Done"
    ENDTIME=$(date +%s)
    echo "[TIMER] Assembly filter took $(($ENDTIME - $STARTTIME)) sec."

    touch ${PREFIX}.assem.done
fi


# 2. human and decoy filtering
if [ -e ${PREFIX}.human_decoy.done ]; then
	echo "Skipping human and decoy filtering";
else

	STARTTIME=$(date +%s)
	HUMAN_DECOY_INPUT_1=${PREFIX}.bbduk_1.fq.gz.assembly_filtered
	HUMAN_DECOY_INPUT_2=${PREFIX}.bbduk_2.fq.gz.assembly_filtered
	echo "[TIMESTAMP] $(date) Mapping reads to human and decoy..."
	AS_over_len_ratio_filter $HUMAN_IDX $HUMAN_DECOY_INPUT_1 $HUMAN_DECOY_INPUT_2 human
	AS_over_len_ratio_filter $ORAL_MICROBIOME_IDX $HUMAN_DECOY_INPUT_1.humanfiltered $HUMAN_DECOY_INPUT_2.humanfiltered decoy
	echo "[TIMESTAMP] $(date) Mapping reads to human and decoy... Done"


	ENDTIME=$(date +%s)
	echo "[TIMER] Human and decoy Filtering took $(($ENDTIME - $STARTTIME)) sec."


	touch ${PREFIX}.human_decoy.done
fi

# 3. taxon filtering
if [ -e ${PREFIX}.taxon.done ]; then
	echo "Skipping taxon filtering";
else
	STARTTIME=$(date +%s)
	TAXON_INPUT_1=${PREFIX}.bbduk_1.fq.gz.assembly_filtered.humanfiltered.decoyfiltered	
	TAXON_INPUT_2=${PREFIX}.bbduk_2.fq.gz.assembly_filtered.humanfiltered.decoyfiltered	

	echo "[TIMESTAMP] $(date) Mapping reads to refseq DB..."
	run_bwa_mem $REFSEQ_IDX $TAXON_INPUT_1 $TAXON_INPUT_2 -a
	python $SCRIPT_DIR/amplicon/get_highestAS_read_match_target.py $TAXON_INPUT_1.bam
	python $SCRIPT_DIR/amplicon/process_fastx_by_readid.py --retain $TAXON_INPUT_1.bam.highestAS_read_match_target.list $TAXON_INPUT_1 $TAXON_INPUT_2 --suffix taxon

	echo "[TIMESTAMP] $(date) Mapping reads to refseq DB... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Refseq alignment took $(($ENDTIME - $STARTTIME)) sec."
	touch ${PREFIX}.taxon.done
fi

# 4. alignment filtering
if [ -e ${PREFIX}.alm_filter.done ]; then
	echo "Skipping alignment filtering";
else


	STARTTIME=$(date +%s)
	ALM_FILTER_INPUT_1=${PREFIX}.bbduk_1.fq.gz.assembly_filtered.humanfiltered.decoyfiltered.taxonfiltered
	ALM_FILTER_INPUT_2=${PREFIX}.bbduk_2.fq.gz.assembly_filtered.humanfiltered.decoyfiltered.taxonfiltered
	echo "[TIMESTAMP] $(date) Alignment filtering..."
	run_bwa_mem $TARGET_IDX $ALM_FILTER_INPUT_1 $ALM_FILTER_INPUT_2
	python $SCRIPT_DIR/amplicon/filter_bam.py --AS_and_MAPQ_filter $ALM_FILTER_INPUT_1.bam --AS_threshold $AS_TO_TARGET_THRES --MAPQ_threshold $MAPQ_TO_TARGET_THRES
	echo "[TIMESTAMP] $(date) Alignment filtering... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Alignment filtering took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.alm_filter.done
fi


# 5 Variant calling with GATK
if [ -e ${PREFIX}.gatk.done ]; then
	echo "Skipping variant calling with GATK";
else
	STARTTIME=$(date +%s)
	GATK_INPUT=${PREFIX}.bbduk_1.fq.gz.assembly_filtered.humanfiltered.decoyfiltered.taxonfiltered.bam.alignmentfiltered.bam

	echo "[TIMESTAMP] $(date) Variant calling with GATK..."
    gatk AddOrReplaceReadGroups \
    I=$GATK_INPUT \
    O=$GATK_INPUT.gatk_rg.tmp.bam \
    RGLB=lib_dummy \
    RGPL=illumina \
    RGPU=unit_dummy \
    RGSM=$TARGET_SEQ_ID
	samtools index $GATK_INPUT.gatk_rg.tmp.bam
	
	if [ ! -f ${TARGET_IDX%.*}.dict ];then
    gatk CreateSequenceDictionary -R $TARGET_IDX
	fi

	gatk HaplotypeCaller \
			-R $TARGET_IDX \
			-I $GATK_INPUT.gatk_rg.tmp.bam \
			-O $GATK_INPUT.vcf \
			--read-filter PrimaryLineReadFilter \
			--max-reads-per-alignment-start 0 
	rm -f $GATK_INPUT.gatk_rg.tmp.bam $GATK_INPUT.gatk_rg.tmp.bam.bai


	echo "[TIMESTAMP] $(date) Variant calling with GATK... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Variant calling with GATK took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.gatk.done
fi

# 6 Realignment
if [ -e ${PREFIX}.realignment.done ]; then
	echo "Skipping realignment";
else
	STARTTIME=$(date +%s)
	REALIGNMENT_INPUT_BAM=${PREFIX}.bbduk_1.fq.gz.assembly_filtered.humanfiltered.decoyfiltered.taxonfiltered.bam.alignmentfiltered.bam
	REALIGNMENT_INPUT_VCF=$REALIGNMENT_INPUT_BAM.vcf

	echo "[TIMESTAMP] $(date) Realignment..."
    bash $SCRIPT_DIR/realignment/realignment.sh \
        -b $REALIGNMENT_INPUT_BAM \
        -r $TARGET_IDX \
        -v $REALIGNMENT_INPUT_VCF \
        -c $TARGET_SEQ_ID \
        -p "illumina" \
        -t $THREADS \
        -o `readlink -f ${REALIGNMENT_INPUT_VCF}`_realigned

	echo "[TIMESTAMP] $(date) Realignment... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Realignment took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.realignment.done
fi
