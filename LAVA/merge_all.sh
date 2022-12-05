#!/bin/bash
#
FILES=($(ls *.log))		  # assuming the output format is [phenotype]_rg.log
N=$(echo ${#FILES[@]})  		# and that all combinations of penotypes have been analysed
for I in ${FILES[@]}; do
        PHEN=$(echo $I | sed 's/.log//')

        # subset log files to relevant output
        tail -n$(($N+4)) $I | head -$((N+1)) > $PHEN.rg 	# (adapt as necessary)

        # add to single data set
        if [[ $I == ${FILES[0]} ]]; then
        	cat $PHEN.rg > all.rg		# only including the header for the first phenotypes
        else
        	cat $PHEN.rg | sed '1d' >> all.rg
        fi
done
