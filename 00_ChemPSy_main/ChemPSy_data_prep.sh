#!/bin/bash

#################################
#     Source .ini file          #
#################################

echo "##############################       ChemPSy       ##############################"


echo "--1-- Checking config file"
source /home/genouest/irset/tdarde/projects/ChemPSy/20171204/scripts/ChemPsy_processing/ChemPSy.ini
echo "Reading config...." >&2
echo "Reading scriptPath: $scriptPath" >&2
echo "Reading config: $Rscript " >&2
echo "Reading dataPath: $dataPath" >&2
echo "Reading config: $processedPath " >&2


#################################
#     Source env                #
#################################
#/local/env/envR-3.1.0.sh


#################################
#           Usage               #
#################################

function usage {
echo -e "\nUsage : $(basename $0) takes 1 :\n-p (--path) : Input path\n-t (--tissue)"
}


###############################################
#           Step 1 - process data             #
###############################################

function step_1 {
	echo "--2-- STEP_1 process_data"
	for tissue in $tissues
		do
			echo $tissue
			outputT=$processedPath$tissue"/"
			mkdir -p $outputT
			
			for gse in $gsePath
				do
					path=$dataPath$tissue"/"$gse"/"
					if [ -d $path ]
						then
							output=$processedPath$tissue"/"$gse"/Experimental_conditions/"
							mkdir -p $output
							scriptA=$scriptPath'Rlauncher.sh'
							$scriptA -p $path -t $tissue -o $output -c $cdfpath
					fi
			done
	done
	while [ $(qstat | grep "ChemPSy_" | wc -l) -ne 0 ]
	do
		echo "Running --2-- STEP_1 process_data"
		sleep 7
	done
	echo "--2-- STEP_1 process_data finish"
}

######################################################
#           Step 2 - listing conditions              #
######################################################

function step_2 {
	echo "--3-- STEP_1.5 Conditions listing"
	for tissue in $tissues
		do
			echo $tissue
			for gse in $gsePath
				do
					path=$processedPath$tissue"/"$gse"/"				
					if [ -d $path ]
						then
							echo "--3-- STEP_1.5 Conditions WORKING"
							scriptB=$scriptPath'listDir.sh' 
							$scriptB -p $path -t $tissue 
					fi
			done
	done
	echo "--3-- STEP_1.5 Conditions listing finish"
}

##############################################
#           Step 3 - Merge data              #
##############################################
function step_3 {
	echo "--4-- STEP_2 Merge data"
	for tissue in $tissues
		do
			echo $tissue
			path=$processedPath$tissue"/"$tissue".txt"	
			if [ -f $path ]
				then
					outdir=$processedPath$tissue"/2_merge_data/"
					mkdir -p $outdir
					logout=$logpath"STEP_2_Merge_data.out"
					scriptWrap=$scriptPath'wrapperR.sh'
					scriptC=$Rscript'2_merge_data.R'
					qsub -N "ChemPSy_STEP_2" -j y -o $logout $scriptWrap $scriptC $path $outdir
			fi
	done
	while [ $(qstat | grep "ChemPSy_" | wc -l) -ne 0 ]
	do
		echo "Running --4-- STEP_2 Merge data"
		sleep 7
	done
	echo "--4-- STEP_2 Merge data finish"
}

###############################################
#           Step 4 - Reduce data              #
###############################################
function step_4 {
	echo "--5-- STEP_3 Reduce Data"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/2_merge_data/data.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/3_reduce_data/"
					mkdir -p $outdir
					logout=$logpath"STEP_3_reduce_data.out"
					scriptWrap=$scriptPath'wrapperR.sh'
					scriptD=$Rscript'3_reduce_data.R'
					qsub -N "ChemPSy_STEP_3" -j y -o $logout $scriptWrap $scriptD $Rdata $outdir
			fi
	done
	while [ $(qstat | grep "ChemPSy_" | wc -l) -ne 0 ]
	do
		echo "Running --5-- STEP_3 Reduce Data"
		sleep 7
	done
	echo "--5-- STEP_3 Reduce Data finish"
}


################################################
#           Step 5 - PCA analysis              #
################################################
function step_5 {
	echo "--6-- STEP_4 PCA Analysis"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/3_reduce_data/data.reduced.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/4_pca_analysis/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath"STEP_4_pca_analysis.out"
					script=$Rscript'4_pca_analysis.R'
					jobName="ChemPSy_"$tissue"_STEP_4"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
	while [ $(qstat | grep "ChemPSy_" | wc -l) -ne 0 ]
	do
		echo "Running --6-- STEP_4 PCA Analysis"
		sleep 7
	done
	echo "--6-- STEP_4 PCA Analysis finish"
}

#######################################################
#           Step 6.1 - Mclust Classification          #
#######################################################
function step_6 {
	echo "--7.1-- STEP_5.1 Mclust Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/mclust/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath$tissue"STEP_5_classification_mclust.out"
					script=$Rscript'5_classification_mclust.R'
					jobName="ChemPSy_"$tissue"_STEP_5.1_mclust"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}

#######################################################
#           Step 7.1 - RecMclust Classification       #
#######################################################
function step_7 {
	echo "--8.1-- STEP_6.1 Mclust Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/recMclust/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath$tissue"STEP_5_classification_recmclust.out"
					script=$Rscript'5_classification_recmclust.R'
					jobName="ChemPSy_"$tissue"_STEP_6.1_recmclust"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}

#######################################################
#           Step 8.1 - dynamictreecut Classification  #
#######################################################
function step_8 {
	echo "--9.1-- STEP_7.1 cutreeDynamic Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/dynamicTreeCut/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					script=$Rscript'5_classification_cutreeDynamic.R'
					logout=$logpath$tissue"STEP_5_classification_cutreeDynamic.out"
					jobName="ChemPSy_"$tissue"_STEP_7.1cutreeDynamic"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}
##########################################################
#           Step 9.1 - recdynamictreecut Classification  #
##########################################################
function step_9 {
	echo "--10.1-- STEP_8.1 reccutreeDynamic Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/recDynamicTreeCut/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath$tissue"STEP_5_classification_reccutreeDynamic.out"
					script=$Rscript'5_classification_reccutreeDynamic.R'
					jobName="ChemPSy_"$tissue"_STEP_8.1reccutreeDynamic"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}
##########################################################
#           Step 10 - hopach Classification             #
##########################################################
function step_10 {
	echo "--11-- STEP_9.1 hopach Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/hopach/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath$tissue"STEP_5_classification_hopach.out"
					script=$Rscript'5_classification_hopach.R'
					jobName="ChemPSy_"$tissue"_STEP_9.1_hopach"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}
##########################################################
#           Step 11 - rechopach Classification          #
##########################################################
function step_11 {
	echo "--11-- STEP_10.1 rechopach Classification"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/recHopach/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath$tissue"_STEP_5_classification_rechopach.out"
					script=$Rscript'5_classification_rechopach.R'
					jobName="ChemPSy_"$tissue"_STEP_10.1_rechopach"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $Rdata $outdir
			fi
	done
}
##########################################################
#           Step 11.5 - Images                           #
##########################################################
function step_115 {

	while [ $(qstat | grep "ChemPSy_" | wc -l) -ne 0 ]
	do
		sleep 10
	done
	echo "--11.5-- STEP_10.5 Image"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					outdir=$processedPath$tissue"/5_classification/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR3_1.sh'
					logout=$logpath"STEP_5.5_"$tissue"_Heatmap.out"
					jobNameImg="ChemPSy_"$tissue"IMG_STEP_10.5_HeatMap"
					scriptimg=$Rscript'5_figures.R'
					qsub -N "$jobNameImg" -j y -o $logout $scriptWrap $scriptimg $outdir $Rdata $outdir
			fi
	done
}

##########################################################
#           Step 12 - Signature                          #
##########################################################
function step_12 {
	echo "--12-- STEP_11 signature"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					indir=$processedPath$tissue"/5_classification/"
					outdir=$processedPath$tissue"/6_signature/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath"STEP_6_signature.out"
					script=$Rscript'6_signatures.R'
					jobName="CheSig_"$tissue"_STEP_11_signature"
					qsub -N "$jobName" -j y -o $logout $scriptWrap $script $indir $Rdata $outdir

			fi
	done
}

##########################################################
#           Step 13 - Enrichissement                     #
##########################################################
function step_13 {
	while [ $(qstat | grep "CheSig_" | wc -l) -ne 0 ]
	do
		sleep 10
	done
	echo "--13-- STEP_12 Enrichisement"
	for tissue in $tissues
		do
			echo $tissue
			Rdata=$processedPath$tissue"/4_pca_analysis/PCA.RData"
			if [ -f $Rdata ]
				then
					echo $Rdata
					indir=$processedPath$tissue"/5_classification/"
					sindir=$processedPath$tissue"/6_signature/"
					CTD="/home/genouest/irset/archives/projects/ChemPSy/CTD_data/CTD_diseases_data.RData"
					outdir=$processedPath$tissue"/7_enrichissement/"
					mkdir -p $outdir
					scriptWrap=$scriptPath'wrapperR.sh'
					logout=$logpath"STEP_7_enrichissement.out"
					script=$Rscript'7_association_enrichment.R'
					jobName="ChemPSy_"$tissue"_STEP_12_enrichissement"
					qsub -N "$jobName" -q "sihp.q" -j y -o $logout $scriptWrap $script $indir $sindir $CTD $Rdata $outdir

			fi
	done
}



###################################
#           Main function         #
###################################

#step_1
step_2
step_3
#step_4
step_5
#step_6
#step_7
#step_8
#step_9
#step_10
#step_11
#step_115
#step_12
#step_13
exit

