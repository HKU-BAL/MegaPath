import pysam
import argparse
import statistics
import sys
import os
from lib.get_list_from_file import get_list_from_file
def get_fastq_name(fastq_files,file_out):
    read_retain_set=set()
    for fastq_fn in fastq_files:
        with pysam.FastxFile(fastq_fn) as fin:
            for entry in fin:
                read_retain_set.add(entry.name)
    file_out.write('\n'.join(read_retain_set))

def filter_contigs(tb_bam_fn,start,end,fastq_files,out_prefix='filter_contigs',mean_mapq_thres=10):
    SAMTOOLS_EXCL_FLAGS=1796
    #extract start and end pos from bed
    with pysam.AlignmentFile(tb_bam_fn, "rb" ) as tb_bam_in, open(f'{out_prefix}.contig_filter.log','w') as log, open(f'{out_prefix}.contig_filter','w') as out, open('assembly_filter.retain_readid','w') as fo_retain:
        passed_contig=set()
        #filter contig unaligned to tb or aligned to non-target regions
        contigs_aligned_tb=set()
        for entry in tb_bam_in.fetch(start=start,end=end):
            if entry.flag&4==0:
                contigs_aligned_tb.add(entry.qname)
        #  empty
        if not contigs_aligned_tb:
            get_fastq_name(fastq_files,fo_retain)
            return

        #filter by meanmaq
        with pysam.AlignmentFile(f'{out_prefix}.bam', "rb" ) as bam_in:
            for ref in contigs_aligned_tb:
                mapq_per_read=[]
                for entry in bam_in.fetch(ref):
                    if entry.flag&SAMTOOLS_EXCL_FLAGS==0:
                        mapq_per_read.append(entry.mapping_quality)
                if mapq_per_read:
                    mean_mapq=statistics.mean(mapq_per_read)
                else:
                    mean_mapq=0
                log_text='\t'.join([ref,str(round(mean_mapq,1))])
                log.write(f'{log_text}\n')
                
                if  mean_mapq>=mean_mapq_thres:
                    passed_contig.add(ref)
        out_text='\n'.join(passed_contig)
        out.write(f'{out_text}\n')

def get_readid_aligned_seqid(in_bam,contig_name_fn):
    contig_name=get_list_from_file(contig_name_fn)
    if not contig_name:
        return
    read_retain_set=set()
    with pysam.AlignmentFile(in_bam) as fin, open('assembly_filter.retain_readid','w') as fo_retain:
        for entry in fin:
            if entry.reference_name in contig_name or entry.flag & 4 == 4:
                read_retain_set.add(entry.qname)
        fo_retain.write('\n'.join(read_retain_set))


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('tb_bam_fn', metavar='tb_bam_fn', type=str,)
    parser.add_argument('start', metavar='start', type=int,)
    parser.add_argument('end', metavar='end', type=int,)
    parser.add_argument('fastq_files', metavar='fastq_files', type=str, nargs='+',)
    parser.add_argument('--mean_mapq_thres', default=10, help='The threshold of mean mapping quality')
    parser.add_argument('--out_prefix', default='', help='Output prefix')
    args = parser.parse_args()
    if args.out_prefix=='':
        args.out_prefix=os.path.basename(args.fastq_files[0])
    filter_contigs(tb_bam_fn=args.tb_bam_fn,
            start=args.start,
            end=args.end,
            fastq_files=args.fastq_files,
            mean_mapq_thres=args.mean_mapq_thres,
            out_prefix=args.out_prefix)
    get_readid_aligned_seqid(in_bam=f'{args.out_prefix}.bam',contig_name_fn=f'{args.out_prefix}.contig_filter')
