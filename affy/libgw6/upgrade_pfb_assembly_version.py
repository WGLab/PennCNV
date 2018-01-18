#!/usr/bin/env python

import os
import sys

tab  = '\t'
endl = '\n'
arg = sys.argv[1:]

usage = 'python ' + __file__ + ' ' + '<in.pfb> '
argc  = 1 

def main():
    if len(arg) < argc:
        print usage
        sys.exit()

    in_pfb_file = arg.pop(0)

    in_pfb_bed_file = in_pfb_file + '.bed' 
    in_pfb_fp = open(in_pfb_file, 'r')
    in_pfb_bed_fp = open(in_pfb_bed_file, 'w') 
    while 1:
        line = in_pfb_fp.readline()
        if not line: break
        if 'Name' in line: continue

        line = line.strip().split(tab)

    



if __name__ == '__main__':
    main()
