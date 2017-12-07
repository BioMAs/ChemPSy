#!/bin/bash

pythonScript="/home/genouest/irset/tdarde/Projets/ChemPSy/script/Open-TG-GAte_process.py"
pathGate="/home/genouest/irset/tdarde/Projets/ChemPSy/TGGATE/Rat/in_vivo/Kidney/Repeat/ /home/genouest/irset/tdarde/Projets/ChemPSy/TGGATE/Rat/in_vivo/Kidney/Single/ /home/genouest/irset/tdarde/Projets/ChemPSy/TGGATE/Rat/in_vivo/Liver/Repeat/ /home/genouest/irset/tdarde/Projets/ChemPSy/TGGATE/Rat/in_vivo/Liver/Single/ /home/genouest/irset/tdarde/Projets/ChemPSy/TGGATE/Rat/in_vitro/"

for pathZ in $pathGate:
do
    echo $pathZ
    for i in $(find $pathZ -type d)
    do
        echo $i
        attributFile=$i"/Attribut.tsv"
        outdir=$pathZ$i
        echo $attributFile
        echo $outdir
    done
done