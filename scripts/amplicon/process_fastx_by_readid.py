import pysam
import argparse
import os
import sys
from lib.get_list_from_file import get_list_from_file
#import gzip
def process_fastx_by_readid(readid_list,if_filter,read_fn,prefix='',suffix='',readid=''):
    if readid=='readid':
        prefix=f'readid_{prefix}'

    for filename in read_fn:
        out_filename=f'{filename}.{suffix}filtered'
        with pysam.FastxFile(filename) as fin, open(f'{prefix}{os.path.basename(out_filename)}', mode='w') as fout:
            for entry in fin:
                if if_filter:
                    if entry.name in readid_list: 
                        continue
                    else:
                        fout.write(''.join([str(entry),'\n']))
                else:

                    if entry.name in readid_list:
                        if readid=='readid':
                            fout.write(''.join([str(entry.name),'\n']))
                        else:
                            fout.write(''.join([str(entry),'\n']))
                    else:
                        continue
        if readid=='readid':
            break

if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('readid_fn', metavar='readid_fn', type=str, 
                        )
    parser.add_argument('read_fn', metavar='read_fn', type=str, nargs='+',
                        )
    group_action = parser.add_mutually_exclusive_group(required=True)
    group_action.add_argument('--retain', dest='if_filter', action='store_false')
    group_action.add_argument('--filter', dest='if_filter', action='store_true')
    parser.add_argument('--prefix', dest='prefix', action='store_const',
                        const='prefix', default='',
                        )
    parser.add_argument('--suffix', dest='suffix',
                         default='',
                        )
    parser.add_argument('--readid', dest='readid',
                         default='',
                        )

    args = parser.parse_args()
    readid_list=get_list_from_file(args.readid_fn)
    process_fastx_by_readid(readid_list=readid_list,if_filter=args.if_filter,read_fn=args.read_fn,prefix=args.prefix,suffix=args.suffix,readid=args.readid)
