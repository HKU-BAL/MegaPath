import pysam
import argparse
import sys
import os
def filter_bam(bam,AS_threshold=0,MAPQ_threshold=0,AS_LEN_ratio=0,if_AS_and_MAPQ_filter=True):
    suffix=''
    if if_AS_and_MAPQ_filter:
        suffix=f'alignmentfiltered'
    else:
        suffix=f'AS{AS_threshold}_AS-LEN_ratio{AS_LEN_ratio}'
    bam_path=f'{os.path.basename(bam)}.{suffix}.bam'
    with pysam.AlignmentFile(bam, "rb") as samfile, pysam.AlignmentFile(bam_path,'wb',template=samfile) as fo:
        for read in samfile.fetch():
            if if_AS_and_MAPQ_filter:
                if read.get_tag('AS')>=AS_threshold and read.mapping_quality>=MAPQ_threshold:
                    fo.write(read)
            else:
                #  AS_and_AS-LEN_ratio_filter
                with open(f'{bam_path}.ASandAS-LEN_ratio_filter.readid','w') as out_list:
                    readid_set=set()
                    for read in samfile.fetch():
                        if read.get_forward_sequence() == None:
                            continue
                        if read.get_tag('AS')>=AS_threshold and read.get_tag('AS')/len(read.get_forward_sequence())>=AS_LEN_ratio:
                            fo.write(read)
                            readid_set.add(read.qname)
                    out_list.write('\n'.join(readid_set))
    os.system(f'samtools index {bam_path}')

if __name__=="__main__":
    parser = argparse.ArgumentParser()
    group_action = parser.add_mutually_exclusive_group(required=True)
    group_action.add_argument('--AS_and_AS-LEN_ratio_filter', dest='if_AS_and_MAPQ_filter', action='store_false')
    group_action.add_argument('--AS_and_MAPQ_filter', dest='if_AS_and_MAPQ_filter', action='store_true')
    parser.add_argument('bam', metavar='bam_fn', type=str)
    parser.add_argument('--AS_threshold', type=int,required=True)
    parser.add_argument('--MAPQ_threshold', type=int,required='--AS_and_MAPQ_filter' in sys.argv)
    parser.add_argument('--AS_LEN_ratio', type=float,required='--AS_and_AS-LEN_ratio_filter' in sys.argv)
    args = parser.parse_args()
    filter_bam(bam=args.bam,AS_threshold=args.AS_threshold,MAPQ_threshold=args.MAPQ_threshold,AS_LEN_ratio=args.AS_LEN_ratio,if_AS_and_MAPQ_filter=args.if_AS_and_MAPQ_filter)
