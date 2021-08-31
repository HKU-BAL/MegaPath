#!/usr/bin/env bash

set -e
set -m
set -o pipefail

SCRIPT=$(readlink -f $0)
SCRIPT_PATH=$(dirname ${SCRIPT})

# configurations
MEGAHIT=megahit
ACD=${SCRIPT_PATH}/ac-diamond-0.1-beta-linux64/ac-diamond
SOAP4=${SCRIPT_PATH}/soap4/soap4
SOAP4_BUILD=${SCRIPT_PATH}/soap4/2bwt-lib/2bwt-builder
HG_INI=${SCRIPT_PATH}/soap4/soap4.ini
NT_INI=${SCRIPT_PATH}/soap4/soap4-nt2.ini
HC_INI=${SCRIPT_PATH}/soap4/soap4-nt2.ini
RIBO_INI=${SCRIPT_PATH}/soap4/soap4.ini

ACD_TMP=/tmp/
SAMTOOLS=samtools
BEDTOOLS=bedtools
TAX_LOOKUP=${SCRIPT_PATH}/cc/taxLookupAcc
DEINTERLEAVE=${SCRIPT_PATH}/cc/deinterleave
FASTQ2LSAM=${SCRIPT_PATH}/cc/fastq2lsam
GEN_COUNT_TB=${SCRIPT_PATH}/cc/genKrakenReport
M8_TO_LSAM=${SCRIPT_PATH}/m8_to_lsam.pl
SAM2CFQ=${SCRIPT_PATH}/cc/sam2cfq
reassign=${SCRIPT_PATH}/cc/reassign
FILTER_CROSS_FAMILY=${SCRIPT_PATH}/cc/filterCrossFamilyReads
EXTRACT_FROM_LSAM=${SCRIPT_PATH}/extractFromLSAM.pl
FQ2FA=${SCRIPT_PATH}/fastq2fasta.pl
CLEANUP=${SCRIPT_PATH}/cc/cleanup
GENOME_COV_FILTER=${SCRIPT_PATH}/cc/genomeCovFilter
LSAM_FILTER=${SCRIPT_PATH}/cc/lsamReadFilter

BBMAP_DIR=${SCRIPT_PATH}/bbmap
BBNORM=${BBMAP_DIR}/bbnorm.sh
BBDUK=${BBMAP_DIR}/bbduk2.sh

DB=${SCRIPT_PATH}/db

export PATH=${SCRIPT_PATH}:${PATH}

# default parameters
PREFIX="megapath"
THREADS=24
NT_CUT_OFF=40
READ_LEN=150
MIN_LEN=50
ENTORPY=0.75
SPIKE_STDEV=60
SPIKE_OVERLAP=0.5

while getopts "p:1:2:t:c:s:o:L:d:M:SHA" option; do
	case "${option}" in
		p) PREFIX=${OPTARG};;
		1) READ1=${OPTARG};;
		2) READ2=${OPTARG};;
		t) THREADS=${OPTARG};; # not working with 4 & BBduk
		c) NT_CUT_OFF=${OPTARG};;
        s) SPIKE_STDEV=${OPTARG};;
        o) SPIKE_OVERLAP=${OPTARG};;
		L) READ_LEN=${OPTARG};;
		d) DB=${OPTARG};;
		S) NO_RIBO=true;;
		H) NO_HG=true;;
		A) DO_ASM=true;;
		*) exit 1;;
	esac
done

if [ ${READ_LEN} -le 50 ]; then
	MIN_LEN=30
	ENTORPY=0
fi

let "HG_CUT_OFF=READ_LEN*3/5";

if [ ${READ_LEN} -le 120 ]; then
	READ_LEN=121 # to activate SOAP4 long read mode
fi

if [ -z "${READ1}" ] || [ -z "${READ2}" ]; then
   echo "Usage: $0 -1 <read1.fq> -2 <read2.fq> [options]"
   echo "    -p  output prefix [megapath]"
   echo "    -t  number of threads [24]"
   echo "    -c  NT alignment score cutoff [40]"
   echo "    -s  SPIKE filter number of stdev [60]"
   echo "    -o  SPIKE overlap [0.5]"
   echo "    -L  max read length [150]"
   echo "    -d  database directory [${SCRIPT_PATH}/db]"
   echo "    -S  skip ribosome filtering"
   echo "    -H  skip human filtering"
   echo "    -A  Perform assembly & protein alignment"
   exit 1
fi

SOAP_HG_IDX=${DB}/soap4-hg/human.maskViral.fna.index
SOAP_NT_IDX=${DB}/refseq/refseq_plus_complete_genome_fungi_and_protozoa_and_fungiUnite_masked
CLEANUP_IDX=${DB}/refseq/hc.ref.index
SOAP_16S_IDX=${DB}/ribo/SILVA.fasta.index
ACD_NR=${DB}/ac-diamond-nr-2/nr
NR_FAI=${DB}/ac-diamond-nr-2/nr.fai
ACC2TID=${DB}/tax/abhv_complete_genome_fungi_and_protozoa_and_fungiUnite.accession2taxid


# 0. preprocessing
if [ -e ${PREFIX}.prep.done ]; then
	echo "Skipping preprocessing";
else
	echo "[TIMESTAMP] $(date) Running BBDuk to preprocess..."
	STARTTIME=$(date +%s)
	BBDUK_THREADS=$((${THREADS} / 2))
	if [ ${BBDUK_THREADS} -lt 1 ]; then
		BBDUK_THREADS=1
	fi

	${BBDUK} kmask=N qtrim=rl trimq=10 threads=${BBDUK_THREADS} minlength=${MIN_LEN} in=${READ1} in2=${READ2} out=stdout.fq ref=${BBMAP_DIR}/resources/adapters.fa hdist=1  | ${BBDUK} entropy=${ENTORPY} in=stdin.fq outm=${PREFIX}.low_compl.fq.gz threads=${BBDUK_THREADS} out=${PREFIX}.bbduk_1.fq.gz out2=${PREFIX}.bbduk_2.fq.gz overwrite=t

	echo "[TIMESTAMP] $(date) Running BBDuk to preprocess... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] BBDuk took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.prep.done
fi

# 1. host filtering
if [ -e ${PREFIX}.host.done ]; then
	echo "Skipping host filtering";
else
	if [ -z ${NO_HG} ]; then

		STARTTIME=$(date +%s)
		echo "[TIMESTAMP] $(date) Mapping reads to host..."
		/usr/bin/time -v  ${SOAP4} pair ${SOAP_HG_IDX} ${PREFIX}.bbduk_1.fq.gz ${PREFIX}.bbduk_2.fq.gz -o ${PREFIX}.dummy -C ${HG_INI} -L ${READ_LEN} -T ${THREADS} -u 750 -F -nc | tee hostfilterresult | ${FASTQ2LSAM} 1 | ${EXTRACT_FROM_LSAM} -t ${HG_CUT_OFF} - | ${DEINTERLEAVE} ${PREFIX}.non_hg
		echo "[TIMESTAMP] $(date) Mapping reads to host... Done"


		ENDTIME=$(date +%s)
		echo "[TIMER] Human Genome Filtering took $(($ENDTIME - $STARTTIME)) sec."

		if [ ! -s ${PREFIX}.non_hg.pe_1.fq ]; then
			echo "No reads remained after host filtering" >&2;
			exit 1;
		fi
	else
		ln -s ${PREFIX}.bbduk_1.fq.gz ${PREFIX}.non_hg.pe_1.fq # WARNING: mind the gz
		ln -s ${PREFIX}.bbduk_2.fq.gz ${PREFIX}.non_hg.pe_2.fq
	fi

	touch ${PREFIX}.host.done
fi

# 1.5 map to ribo
if [ -e ${PREFIX}.ribosome.done ]; then
	echo "Skipping Ribo cleaning";
else
	if [ -z ${NO_RIBO} ]; then
		STARTTIME=$(date +%s)
		echo "[TIMESTAMP] $(date) Mapping reads to 16S DB..."
		/usr/bin/time -v ${SOAP4} pair ${SOAP_16S_IDX} ${PREFIX}.non_hg.pe_1.fq ${PREFIX}.non_hg.pe_2.fq -o ${PREFIX}.dummy -C ${RIBO_INI} -L ${READ_LEN} -T ${THREADS} -u 750 -F -P -nc -top 100 | ${FASTQ2LSAM} | ${EXTRACT_FROM_LSAM} -t 0.95 - | ${DEINTERLEAVE} ${PREFIX}.non_ribo
		echo "[TIMESTAMP] $(date) Mapping reads to 16S DB... Done"

		ENDTIME=$(date +%s)
		echo "[TIMER] Ribosome Filtering took $(($ENDTIME - $STARTTIME)) sec."
	fi
	touch ${PREFIX}.ribosome.done
fi

# 2. map to NT
if [ -e ${PREFIX}.nt.done ]; then
	echo "Skipping Nucleotide mapping";
else
	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Mapping reads to Nucleotide DB..."
	i=0
	rm -f ${PREFIX}.nt.tmp_out.pe_1.fq ${PREFIX}.nt.tmp_out.pe_2.fq

	if [ ${NO_RIBO} ]; then
		NT_INPUT_1=${PREFIX}.non_hg.pe_1.fq
		NT_INPUT_2=${PREFIX}.non_hg.pe_2.fq
	else
		NT_INPUT_1=${PREFIX}.non_ribo.pe_1.fq
		NT_INPUT_2=${PREFIX}.non_ribo.pe_2.fq
	fi

	ln -s $(readlink -f ${NT_INPUT_1}) ${PREFIX}.nt.tmp_out.pe_1.fq
	ln -s $(readlink -f ${NT_INPUT_2}) ${PREFIX}.nt.tmp_out.pe_2.fq

	while [ -e ${SOAP_NT_IDX}.$i.index.pac ]; do

		let "j=i+1"

		rm -f ${PREFIX}.nt.tmp_in_1.fq ${PREFIX}.nt.tmp_in_2.fq
		mv -f ${PREFIX}.nt.tmp_out.pe_1.fq ${PREFIX}.nt.tmp_in_1.fq
		mv -f ${PREFIX}.nt.tmp_out.pe_2.fq ${PREFIX}.nt.tmp_in_2.fq

		OPT="-b -o ${PREFIX}.nt.bam -L ${READ_LEN} -T ${THREADS} -u 750 -F -C ${NT_INI} -top 95"

		if [ $i -eq 0 ]; then
			OPT=${OPT}" -nc"
		fi

		if [ -e ${SOAP_NT_IDX}.$j ]; then
			/usr/bin/time -v ${SOAP4} pair ${SOAP_NT_IDX}.${i}.index ${PREFIX}.nt.tmp_in_1.fq ${PREFIX}.nt.tmp_in_2.fq ${OPT} | ${DEINTERLEAVE} ${PREFIX}.nt.tmp_out
		else
			/usr/bin/time -v ${SOAP4} pair ${SOAP_NT_IDX}.${i}.index ${PREFIX}.nt.tmp_in_1.fq ${PREFIX}.nt.tmp_in_2.fq ${OPT} | ${FASTQ2LSAM} 1 | ${TAX_LOOKUP} ${ACC2TID} ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp - | gzip -1 > ${PREFIX}.nt.lsam.id.gz
			#merge nt.bam
			#the following before fi (end of if statement) are spike-filter
			/usr/bin/time -v ${SAMTOOLS} merge ${PREFIX}.nt.unsorted.bam ${PREFIX}.nt.bam.*[0-9]
			/usr/bin/time -v ${SAMTOOLS} sort ${PREFIX}.nt.unsorted.bam ${PREFIX}.nt
			rm -f ${PREFIX}.nt.unsorted.bam
			/usr/bin/time -v ${SAMTOOLS} view -H ${PREFIX}.nt.bam | awk -F'\t|:| ' '{if ($1=="@SQ") print $3"\t"$5 }' > ${PREFIX}.genome
			/usr/bin/time -v ${BEDTOOLS} bamtobed -i ${PREFIX}.nt.bam > ${PREFIX}.nt.bed
			/usr/bin/time -v ${BEDTOOLS} genomecov -bga -ibam ${PREFIX}.nt.bam -g ${PREFIX}.genome > ${PREFIX}.genomecov
			/usr/bin/time -v ${GENOME_COV_FILTER} ${PREFIX}.genome ${PREFIX}.genomecov ${SPIKE_STDEV} > ${PREFIX}.genomecovFilter
			/usr/bin/time -v ${BEDTOOLS} annotate -i ${PREFIX}.nt.bed -files ${PREFIX}.genomecovFilter | awk -F'\t' -v threshold="${SPIKE_OVERLAP}" '{if ($7 >= threshold) print $4}' > ${PREFIX}.readFilter
			/usr/bin/time -v mv ${PREFIX}.nt.lsam.id.gz ${PREFIX}.nt.lsam.id.original.gz
			/usr/bin/time -v gunzip -fk ${PREFIX}.nt.lsam.id.original.gz
			/usr/bin/time -v ${LSAM_FILTER} ${PREFIX}.readFilter ${PREFIX}.nt.lsam.id.original > ${PREFIX}.nt.lsam.id
			/usr/bin/time -v gzip -1fk ${PREFIX}.nt.lsam.id
			
		fi

		i=$j
	done

	echo "[TIMESTAMP] $(date) Mapping reads to NT... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Nucleotide alignment took $(($ENDTIME - $STARTTIME)) sec."
	touch ${PREFIX}.nt.done
fi

# 3. Generate count-tables
if [ -e ${PREFIX}.count.done ]; then
	echo "Skipping NT counting";
else

	echo "[TIMESTAMP] $(date) Reassign..."
	# ${CLEANUP} ${PREFIX}.nt.lsam.id.gz | gzip -1 > ${PREFIX}.nt.cleanup.lsam.id.gz
	STARTTIME=$(date +%s)
    if [ ! -e ${PREFIX}.nt.lsam.id.gz ]; then
        echo "${PREFIX}.nt.lsam.id.gz not found!"
        exit 1
    fi
	${reassign} -p ${THREADS} -t ${NT_CUT_OFF} ${PREFIX}.nt.lsam.id.gz > ${PREFIX}.nt.ra.lsam.id
	echo "[TIMESTAMP] $(date) Reassign... Done"
	ENDTIME=$(date +%s)
	echo "[TIMER] Reassignment took $(($ENDTIME - $STARTTIME)) sec."


	STARTTIME=$(date +%s)
	echo "[TIMESTAMP] $(date) Generating count table..."
	${GEN_COUNT_TB} ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp ${PREFIX}.nt.lsam.id.gz ${NT_CUT_OFF} > ${PREFIX}.nt.report &
	${GEN_COUNT_TB} ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp ${PREFIX}.nt.ra.lsam.id ${NT_CUT_OFF} > ${PREFIX}.nt.ra.report &
	wait
	# ${GEN_COUNT_TB} ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp ${PREFIX}.nt.cleanup.lsam.id.gz ${NT_CUT_OFF} > ${PREFIX}.nt.cleanup.report
	echo "[TIMESTAMP] $(date) Generating count table... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Generate count table took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.count.done
fi

# 4. Assembly of viral & unmapped reads
if [ -e ${PREFIX}.asm.done ] || [ -z ${DO_ASM} ]; then
	echo "Skipping Megahit assembly";
else
	STARTTIME=$(date +%s)

	echo "[TIMESTAMP] $(date) Extracting viral & ummapped reads..."
	${EXTRACT_FROM_LSAM} -t ${NT_CUT_OFF} -v ${PREFIX}.nt.lsam.id.gz | gzip -1 > ${PREFIX}.nt.viral.and.unmapped.fq.gz
	cat ${PREFIX}.low_compl.fq.gz >> ${PREFIX}.nt.viral.and.unmapped.fq.gz
	echo "[TIMESTAMP] $(date) Extracting viral & ummapped reads... Done"

	echo "[TIMESTAMP] $(date) BBNorm";
	${BBNORM} -Xmx20G interleaved=true in=${PREFIX}.nt.viral.and.unmapped.fq.gz out=${PREFIX}.nt.viral.and.unmapped.bbnorm.fq.gz target=70 mindepth=1 threads=${THREADS} overwrite=t
	echo "[TIMESTAMP] $(date) BBNorm Done";

	echo "[TIMESTAMP] $(date) Assembly..."
	${MEGAHIT} -t ${THREADS} --presets meta-sensitive --12 ${PREFIX}.nt.viral.and.unmapped.bbnorm.fq.gz -o ${PREFIX}.asm
	echo "[TIMESTAMP] $(date) Assembly... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Assembly overall took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.asm.done
fi

# 4.1 map to contigs to NR
if [ -e ${PREFIX}.remap.done ] || [ -z ${DO_ASM} ]; then
	echo "Skipping remap";
else
	STARTTIME=$(date +%s)

	echo "[TIMESTAMP] $(date) Mapping to contigs..."
	if [ ! -s ${PREFIX}.asm/final.contigs.fa ]; then
		echo ">1" > ${PREFIX}.asm/final.contigs.fa
		echo "NNNNN" >> ${PREFIX}.asm/final.contigs.fa
	fi
	${SOAP4_BUILD} ${PREFIX}.asm/final.contigs.fa
	${EXTRACT_FROM_LSAM} -t ${NT_CUT_OFF} -i ${PREFIX}.nt.lsam.id.gz | ${DEINTERLEAVE} ${PREFIX}.nt.unmapped
	${SOAP4} pair ${PREFIX}.asm/final.contigs.fa.index ${PREFIX}.nt.unmapped.pe_1.fq ${PREFIX}.nt.unmapped.pe_2.fq -o ${PREFIX}.dummy -C ${NT_INI} -L ${READ_LEN} -T ${THREADS} -u 750 -F | ${FASTQ2LSAM} 1 | gzip -1 > ${PREFIX}.r2c.lsam.gz
	rm -f ${PREFIX}.dummy.*
	${EXTRACT_FROM_LSAM} -t ${NT_CUT_OFF} -s -g ${PREFIX}.r2c.lsam.gz | ${FQ2FA} > ${PREFIX}.contig.unmap.fa
	sed 's/^>/>contig_/' ${PREFIX}.asm/final.contigs.fa >> ${PREFIX}.contig.unmap.fa
	echo "[TIMESTAMP] $(date) Mapping to contigs... Done"

	echo "[TIMESTAMP] $(date) AC-DIAMOND to NR DB..."
	${ACD} blastx -p ${THREADS} -q ${PREFIX}.contig.unmap.fa -d ${ACD_NR} -a ${PREFIX}.nr -t ${ACD_TMP} --log
	${ACD} view -a ${PREFIX}.nr.daa -o ${PREFIX}.nr.m8 -t ${ACD_TMP}
	${ACD} view -a ${PREFIX}.nr.daa -f sam -t ${ACD_TMP} | tee ${PREFIX}.nr.sam | grep -v "^@" | cut -f3 | sort | uniq > aligned
	grep -w -F -f aligned ${NR_FAI} > aligned_nr.fai
	${SAMTOOLS} view -bS ${PREFIX}.nr.sam -t aligned_nr.fai -o ${PREFIX}.nr.bam
	echo "[TIMESTAMP] $(date) AC-DIAMOND to NR DB... Done"

	echo "[TIMESTAMP] $(date) Taxa lookup..."
	${M8_TO_LSAM} ${PREFIX}.nr.m8 | gzip -1 > ${PREFIX}.nr.lsam.gz
	${TAX_LOOKUP} ${DB}/tax/prot.accession2taxid.gz ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp ${PREFIX}.nr.lsam.gz | gzip -1 > ${PREFIX}.nr.lsam.id.gz
	${SCRIPT_PATH}/r2c_to_r2g.pl ${PREFIX}.r2c.lsam.gz ${PREFIX}.nr.lsam.id.gz | gzip -1 > ${PREFIX}.nt.unmap.r2g.lsam.id.gz
	zcat ${PREFIX}.nr.lsam.id.gz ${PREFIX}.nt.unmap.r2g.lsam.id.gz | grep -v '^contig_' | ${GEN_COUNT_TB} ${DB}/tax/nodes.dmp ${DB}/tax/names.dmp - ${NT_CUT_OFF} > ${PREFIX}.nr.report
	echo "[TIMESTAMP] $(date) Taxa lookup... Done"

	ENDTIME=$(date +%s)
	echo "[TIMER] Protein alignment overall took $(($ENDTIME - $STARTTIME)) sec."

	touch ${PREFIX}.remap.done
fi
