#!/bin/bash

lImgName="92659 95852 94119 94130 95068 95290 95404 95548 95745 96151 102509 102295 101504 101499 101268 96275 96395 96442 96748 96952 97053 97266 97498 97532 101098 101111"
cd /home/genouest/irset/archives/projects/ChemPSy/processed_data/
mkdir ../exclude_processed_data/
for elmt in $lImgName
	do
		filename="qc_image_"$elmt".CEL.gz.png"
		for i in $(find . -name "$filename")
			do
				IFS=’/’ read -ra NAMES <<< "$i"                 
				name=${NAMES[1]}
				mv $name"/" ../exclude_processed_data/
			done
	done