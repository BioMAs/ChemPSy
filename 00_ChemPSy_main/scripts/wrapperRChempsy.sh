#!/bin/bash

. /local/env/envR-3.1.0.sh 

Rscript /groups/irset/archives/projects/ChemPSy/scripts/1_process_data_gcrma.R $1 $2 $3 $4
