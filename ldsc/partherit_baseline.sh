#!/bin/sh

ftrait=$1
fannot=$2
foutput=$3

trait=$(basename "$ftrait")
annot=$(basename "$fannot")

echo "*** Beginning LDSC partitioned heritability..."
echo "Trait: $trait"
echo "Annot: $annot"
echo "Output: $foutput"

module load python/2.7.15

python /home/gokala/programs/ldsc/ldsc.py  \
--h2 ${ftrait} \
--out ${foutput} \
--frqfile-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_Phase3_frq/1000G.EUR.QC. \
--overlap-annot \
--ref-ld-chr /data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/new_annotations/${fannot}/${annot}.,/data/clusterfs/lag/users/gokala/beat-dyslexiaevol/resources/baselineLD/baselineLD.  \
--w-ld-chr /data/workspaces/lag/shared_spaces/Resource_DB/LDscores/Phase3/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
--print-coefficients

echo "Finished!"
