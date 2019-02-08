# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import subprocess
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import logging
import os
import argparse
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u
from astroquery.irsa_dust import IrsaDust
from shutil import copyfile

from utils.cigale_tools import prep_cigale_data, cigale_config, run_cigale

def parser():
    parser=argparse.ArgumentParser()
    parser.add_argument('-l','--namelist',help='List of SN names as a .txt or .csv file',default = None)
    parser.add_argument('-d','--dered',help='Deredden the input magnitudes?',action='store_true')
    parser.add_argument('-ci','--docig',help='Do Cigale run?',action='store_true')
    parser.add_argument('-ml','--ml',help ='Do ML-based cut?',action='store_true')
    return parser.parse_args()

def main():
    args = parser()
    prep_cigale_data(sn_name_fn = args.namelist,dered=args.dered,ml=args.ml)
    if args.docig:
        dered_suffix=''
        if args.dered:
            dered_suffix='_dered'
        cigale_config('%s.dat'%os.path.split(args.namelist)[-1].split('.')[0]+dered_suffix)
        run_cigale(os.path.split(args.namelist)[-1].split('.')[0])
    print ('Done!')

if __name__=="__main__":
    main()
