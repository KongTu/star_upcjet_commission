#!/bin/bash

input=$1
dir=`pwd`
filelist=$dir/$input
scriptdir=`pwd`/scripts
mkdir -p $scriptdir
outdir=/star/u/starkong/trigger-testing/jobtest/batch_3
echo $filelist
mkdir -p $outdir
mkdir -p $outdir/log
mkdir -p $outdir/err
echo $outdir
star-submit-template -template Mudst.xml -entities filelist=$filelist,outdir=$outdir,scriptdir=$scriptdir
