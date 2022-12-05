#!/bin/bash
#
# This script will reformat bed files for LAVA.

# Add rownumbers

for i in *.bed; do
	
	awk 'NR>1 {print $1, $2, $3}' $i > tmp && mv tmp ${i} 
	nl -n ln $i > tmp && mv tmp ${i}
	awk '{gsub("chr", "", $2); print}' $i > tmp && mv tmp ${i}
	echo -e "LOC CHR START STOP" | cat - $i > tmp && mv tmp ${i};
	
done
