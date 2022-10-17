#!/bin/env python

import sys

#sys.path.append("..")
import BE_endo_Pack.main.File_input_main  as PIP_F
import argparse


if __name__=='__main__':
    ap=argparse.ArgumentParser()
    #ap.add_argument("-M","--ModelType",help="Select model type (Sequence/Combination/Expression/RNAPoll/CTCF/H3K4me3/Dnase/ATAC/Methylation/H3K4me1/H3K27ac/H3K36me3)",type=str,required=True)
    ap.add_argument("-T","--Editor",help="BASE editor(ABE/CBE)",type=str,required=True)
    ap.add_argument("-W","--Wordir",help="The directory of BE_endo_smart",type=str,required=True)
    ap.add_argument("-B","--SgRNABED",help="Output sgRNA BED file path",type=str,required=True)
    ap.add_argument("-O","--Outfile",help="Outfile path to store effiency and propotion result",type=str,required=True)
    ap.add_argument("Input_file",help="40bp(10bp upstream + 20bp sgRNA+3bp PAM + 7bp downstream)(Tab)chr:start-end",type=str)
    args=ap.parse_args()
    PIP_F.pip_file(args.Input_file,args.Editor,args.SgRNABED,args.Wordir,args.Outfile)
    
    
    

    
