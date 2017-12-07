##!/bin/bash

##!/bin/bash

function usage {
echo -e "\nUsage : $(basename $0) takes 1 :\n-p (--path) : Input path"
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


cd $path"Experimental_conditions/"
lCellFileRemove="003017906019 003017906007 003017906016 003017906018 003017906006 003017906011 003017906008 003017906003 003017906001 003017922030_2 003017906013 003017905028 003017923005 003017906021 003017923001 003017906014 003017906020 003017906005 003017906015 003017906004 003017911029 003017906002 003017923004 003017906009 003017905029 003017906010 003017905030 003017906012 003017024023 003017048027 003017048028 003017048029 003017048030 003017049001 003017049002 003017049003 003017104017 003017108014 003017145016 003017165001 003017196002 003017242025 003017281012 003017281014 003017281029 003017391013 003017392005 003017419004 003017484016 003017500005 003017510020 003017513015 003017541008 003017609008 003017662027 003017687017 003017687029 003017746019 003017746025"
for elmt in $lCellFileRemove
	do
		elmtName=$elmt".CEL.gz"
		for i in $(grep -lR $elmtName *)
			do
				echo $i
				IFS=’/’ read -ra NAMES <<< "$i"                                                            
				name=${NAMES[0]}
				final=$name"/treatment.info"
				echo "name :"$final
				echo "element :"$elmtName
				sed -i.BAK "/$elmtName/d" $i
		done		
done
