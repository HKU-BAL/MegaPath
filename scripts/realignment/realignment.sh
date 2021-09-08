#!/usr/bin/env bash
while getopts b:v:r:c:p:o:t: option
do
    case "${option}"
        in
        b) BAM_PATH=("`readlink -f ${OPTARG}`");;
        v) VCF_PATH=("`readlink -f ${OPTARG}`");;
        r) REF_PATH=`readlink -f ${OPTARG}`;;
        c) CONTIGS=${OPTARG};;
        p) PLATFORM=${OPTARG};;
		o) OUT_FOLDER=${OPTARG};;
        t) THREADS=${OPTARG};;
    esac
done

SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})
DATE_TIME=`date "+%Y%m%d_%H%M%S"`
ILLUMINA_REALIGN_BAM_FOLDER=${OUT_FOLDER}/illumina_realignment
mkdir -p ${OUT_FOLDER}
mkdir -p ${OUT_FOLDER}/tmp
mkdir -p ${ILLUMINA_REALIGN_BAM_FOLDER}/
cd ${OUT_FOLDER}

if [ "${PLATFORM}" = "illumina" ]
then
    echo "[INFO] Realign Illumina BAM"
    time parallel --joblog ./illumina_reads_realignment.log -j${THREADS} \
    "pypy ${SCRIPT_PATH}/realign_illumina_reads.py \
        --bam_fn {1} \
        --ref_fn ${REF_PATH} \
        --read_fn ${ILLUMINA_REALIGN_BAM_FOLDER}/{1/} \
        --samtools samtools \
        --ctgName ${CONTIGS}" ::: ${BAM_PATH[@]}

    ls ${ILLUMINA_REALIGN_BAM_FOLDER}/*.bam | parallel -j20 samtools index {}
fi

for i in ${!VCF_PATH[@]}
do
BAM_FILE=${BAM_PATH[i]}
gzip -fdc ${VCF_PATH[i]} | grep -v '#' | cut -f2 > ${OUT_FOLDER}/tmp/${i}_all_pos
readarray ALL_POS < ${OUT_FOLDER}/tmp/${i}_all_pos
echo "[INFO] Extract Candidates"
time parallel --joblog  ./${i}_create_tensor_pileup.log -j${THREADS} \
"pypy ${SCRIPT_PATH}/local_realignment.py \
--bam_fn ${BAM_FILE} \
--ref_fn ${REF_PATH} \
--pos {1} \
--platform ${PLATFORM} \
--samtools `which samtools` \
--ctgName ${CONTIGS} \
--ilmn_realign_folder ${ILLUMINA_REALIGN_BAM_FOLDER} \
--output_fn ${OUT_FOLDER}/tmp/${i}_output_{#}" ::: ${ALL_POS[@]} |& tee  ./${i}_CTP.log

cat ${OUT_FOLDER}/tmp/${i}_output_* > ${OUT_FOLDER}/tmp/${i}_all_alt_info

pypy ${SCRIPT_PATH}/extract_vcf_position.py \
--vcf_fn ${VCF_PATH[i]} \
--alt_fn ${OUT_FOLDER}/tmp/${i}_all_alt_info \
--output_fn ${OUT_FOLDER}/${i}_realign.vcf |& tee  ./${i}_EVP.log
done
