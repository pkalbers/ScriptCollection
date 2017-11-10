#!/bin/bash

PREFIX="history" # history file prefix

cat > error_profile_list.txt <<EOT
~/Research/GenotypeErrorProfiles/1000G/_stat.prop_matrix.1000G.platinum.typed_TRUE.RData
~/Research/GenotypeErrorProfiles/Affymetrix_Axiom/_stat.prop_matrix.affymetrix_axiom.platinum.typed_TRUE.RData
~/Research/GenotypeErrorProfiles/Illumina_HiSeq4000_TruSeq/_stat.prop_matrix.illumina_hiseq4000_truseq.platinum.typed_FALSE.RData
~/Research/GenotypeErrorProfiles/Illumina_HiSeq2000/_stat.prop_matrix.illumina_hiseq2000.platinum.typed_FALSE.RData
~/Research/GenotypeErrorProfiles/10XGenomics/_stat.prop_matrix.10xgenomics.platinum.typed_FALSE.RData
EOT


while read EFILE; do

  ENAME=`echo $EFILE | sed 's|.*matrix\.\(.*\)\.RData|\1|'`
  ENAME="generror_${ENAME%%.*}"

  echo ""
  echo "Error file: ${EFILE}"
  echo "Error name: ${ENAME}"
  echo ""

  Rscript error.R ${PREFIX} ${EFILE} ${ENAME}

  echo ""

done < error_profile_list.txt





