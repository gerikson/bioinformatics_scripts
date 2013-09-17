#!/bin/bash

cd $WORKINGDIR
export CNVHOME=/gpfs/home/gerikson/CNV_pipeline/

python $CNVHOME/CNV_parser.py $WORKINGDIR $INTERNALERROR $INPUT_FILE 
