#!/bin/bash

source /home/genouest/irset/tdarde/projects/ChemPSy/20171204/scripts/ChemPsy_processing/ChemPSy.ini

scriptWrap=$scriptPath'wrapperR.sh'
Rscript='/groups/irset/archives/projects/ChemPSy/scripts/1_process_data_gcrma.R'
function usage {
echo -e "\nUsage : $(basename $0) takes 1 :\n-p (--path) : Input path\n-t (--tissue): Tissue\n-o (--output): Output\n-c (--cdf)CDF file path"
}

if [ "$1" == "-h" ] || [ "$1" == "--help" ]
        then
                usage
                exit
fi


while (( $# > 1 ))
	do key="$1"
	case $key in
		-p|--path)
		path="$2"
		shift
		;;
		-t|--tissue)
		tissue="$2"
		shift
		;;
		-o|--output)
		outDir="$2"
		shift
		;;
		-c|--cdf)
		cdfpath="$2"
		shift
		;;
		*)
		echo "unknown option"
		usage
		exit
		;;
	esac
	shift
done

if [ -z $path ]
	then
		echo "Path not indicated"
		usage
		exit #exit the script if the genome of the input file was not indicated
fi
	

if [ ! -d $path ] #check if the path exists
	then
		echo "the path does not exists"
		exit
fi

if [[ "${path: -1}" != "/" ]] #add a slash in the end of the path if it doesn't exists yet (the lack of slash causes errors when setting the path later on)
        then
                path=$path"/"
fi

echo "RLauncher script"
Cellpath=$path"Individual_experiments"
Expath=$path"Experimental_conditions/"
#cdfpath="/home/genouest/irset/archives/projects/ChemPSy/scripts/CDF/Rat2302_Rn_ENTREZG.cdf.RData"
#outDir="/home/genouest/irset/archives/projects/ChemPSy/processed_data/Rattus_norvegicus/"$tissue"/"

z=1 #index of the current qsub job
index=1

for i in $(find $Expath -type d)
	do
		
		IFS=’/’ read -ra NAMES <<< "$i" 
		LENGTH=${#NAMES[@]} # Get the length.                                          
		LAST_POSITION=$((LENGTH - 1)) # Subtract 1 from the length.                   
		name=${NAMES[${LAST_POSITION}]}
		if [[ "$name" != "Experimental_conditions" ]]
			then
				Treatfile=$Expath$name"/treatment.info"
				finalOutDir=$outDir$name"/"
				let "z = $z + 1"
				if (( $z % 25 == 0))
					then
						let "index = $z - 1"
				fi
				nomSubmit="ChemPSy_"$name$z
				logout=$logpath"STEP_1_processed_data.out"
				mkdir $finalOutDir
				if (( $index == 1 ))
					then
						qsub -N "$nomSubmit" -o $logout -q "sihp.q" /home/genouest/irset/tdarde/projects/ChemPSy/20171204/scripts/ChemPsy_processing/scripts/wrapperRChempsy.sh $Treatfile $Cellpath $cdfpath $finalOutDir
					else
						lastJob="ChemPSy_"$name$index
						qsub -N "$nomSubmit" -j y -o $logout -hold_jid $lastJob -q "sihp.q" /home/genouest/irset/tdarde/projects/ChemPSy/20171204/scripts/ChemPsy_processing/scripts/wrapperRChempsy.sh $Treatfile $Cellpath $cdfpath $finalOutDir
				fi
		fi
	done

sleep 1

echo "All jobs submitted"
exit 0
