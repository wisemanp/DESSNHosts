# -*- coding: utf-8 -*-

import os
import argparse

from utils.cigale_tools import cigale_config, run_cigale

def parser():
    parser=argparse.ArgumentParser()
    parser.add_argument('-n','--fname',help='Name of Cigale input data file',default = None)
    return parser.parse_args()

def main():
    args = parser()

    cigale_config(args.fname)
    run_cigale(args.fname.split(',')[0])
    print ('Done!')

if __name__=="__main__":
    main()
