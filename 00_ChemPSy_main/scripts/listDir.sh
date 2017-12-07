 #!/bin/bash
function usage {
echo -e "\nUsage : $(basename $0) takes 2 :\n-p (--path) : Input path\n-t (--tissue)"
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

ExPath=$path"Experimental_conditions/"
echo $ExPath
cd $ExPath
#remove directories with no filtration.txt files
for i in $(find $Expath -type d)
do
	IFS=’/’ read -ra NAMES <<< "$i" 
	LENGTH=${#NAMES[@]} # Get the length.                                          
	LAST_POSITION=$((LENGTH - 1)) # Subtract 1 from the length.                   
	name=${NAMES[${LAST_POSITION}]}
	filtrationFile=$Expath$name"/filtration.txt"
	if [[ "$filtrationFile" != "./filtration.txt" ]]
		then
			if [ -f $filtrationFile ]
				then
					continue
			else
				echo $filtrationFile" not exists"
				echo "Remove directory : "$path"Experimental_conditions/"$name"/"
				dirPath=$path"Experimental_conditions/"$name"/"
				echo $dirPath
				rm -rf $dirPath
			fi
	fi
done

ls $ExPath > $path"directories.txt"
sed -i 's/$/\/filtration.txt/' $path"directories.txt"
sed -i 's#^#'"${ExPath}"'#' $path"directories.txt"
cat $path"directories.txt" >> $path"../${tissue}.txt"
cat $path"directories.txt" >> $path"../../ALL/ALL.txt"
