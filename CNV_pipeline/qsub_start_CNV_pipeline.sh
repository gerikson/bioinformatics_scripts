#!/bin/bash

#setenv CNVHOME
export CNVHOME=/gpfs/home/gerikson/CNV_pipeline/
source /gpfs/home/gerikson/CNV_pipeline/cnv_environ.sh
#
#   $1 == working directory we are processing
#   $2 == internal error file name
#   $3 == input variants
#

#WORKINGDIR=$1
#INTERNALERROR=$2
#INPUT_FILE=$3

#echo $INPUT_FILE

qsub -v WORKINGDIR=$1,INTERNALERROR=$2,INPUT_FILE=$3 -q workq -M gerikson@scripps.edu -l mem=8G -l cput=36:00:00 -l walltime=72:00:00 $CNVHOME/start_cnv.sh 

 
#qsub -q workq -M gerikson@scripps.edu -l mem=8G -l cput=36:00:00 -l walltime=72:00:00 $CNVHOME/start_pipeline.qsub $WORKINGDIR $INTERNALERROR $INPUT_FILE


