#!/bin/env bash

Input=$1
Wdir=$2

bedtools intersect -loj -wa -wb -a ${Input} -b ${Wdir}/BE_endo_smart/model_data/HEK293_bed/endo_factors.bed >${Wdir}/BE_endo_smart/main/Intersection_temp/result.bed


if [ $? -eq 0 ];then
echo "bedtools intersecting Success!"

else
echo "bedtools interscting Failure!" >&2

fi

 
