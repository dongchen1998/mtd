#!/bin/bash
# 根据Mt的GSE号，批量下载GSE与相关GSM的实验信息原始数据
# 一般下的并不全，但是绝大部分都能下载到

Mt='/home/dongjc/mtd/Obtaining_experimental_information/Mt'
Mt_file='/home/dongjc/mtd/Obtaining_experimental_information/Mt/File'
Mt_GSE_soft='/home/dongjc/mtd/Obtaining_experimental_information/Mt/GSE_soft'
Mt_GSM_soft='/home/dongjc/mtd/Obtaining_experimental_information/Mt/GSM_soft'

mkdir $Mt
mkdir $Mt_file
mkdir $Mt_GSE_soft
mkdir $Mt_GSM_soft

cat Mt_gse.txt | while read gse; do geofetch -i $gse --processed -n $Mt_file --just-metadata ; done

cd $Mt_file

mv GSE*_GSE.soft $Mt/GSE_soft
mv GSE*_GSM.soft $Mt/GSM_soft