conda activate quibl_env

# GENUS LEVEL
cat ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/*.treefile > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk

grep -l "ERR2040368_nuclear" ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/*.treefile | xargs cat > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees1.nwk
mv ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees1.nwk ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk
pxrmt -t ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk -n ERR3487350_nuclear > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees1.nwk
mv ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees1.nwk ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk

grep -oE "[a-zA-Z0-9_.-]+" ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk | grep -vE "^[0-9.-]+$" | sort | uniq -c | awk '$1 > 180 {print $2}' > ./QuIBL/core_taxa.txt
grep -oE "[a-zA-Z0-9_.-]+" ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk | grep -vE "^[0-9.-]+$" | sort | uniq -c | awk '$1 <= 180 {print $2}' > ./QuIBL/remove_taxa.txt

cat ./QuIBL/remove_taxa.txt

pxrmt -t ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees.nwk -n AXPJ_nuclear,D06_nuclear,D09_nuclear,D14_nuclear,D08_nuclear,D05_nuclear > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190.nwk

cat ./QuIBL/remove_sameSPP.txt

pxrmt -t ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190.nwk -n D10_nuclear,D11_nuclear,D12_nuclear > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190_prunedSPP.nwk

TOTAL_TAXA=$(grep -oE "[a-zA-Z0-9_.-]+" ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190_prunedSPP.nwk | grep -vE "^[0-9.-]+$" | sort -u | wc -l | xargs)

echo $TOTAL_TAXA

TARGET=$TOTAL_TAXA; awk -F, -v t="$TARGET" 'NF == t' ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190_prunedSPP.nwk > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190_prunedSPP_completeCases.nwk

nano QuIBL/myparameters.txt

mkdir ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/qibl_trees
mv ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/allbestrees_pruned190_prunedSPP_completeCases.nwk ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned190_prunedSPP_completeCases.nwk

python QuIBL/QuIBL.py ./QuIBL/myparameters.txt


# LINUM BIENNE NODE
## W/ ALL
grep -l "perenne" ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases.nwk | xargs cat > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode.nwk
pxrmt -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode.nwk -n perenne,D07,C96,D03,D02,D13 -c > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode1.nwk
pxrr -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode1.nwk -g perenne > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode1R.nwk

## W/ NARBONENSE, BIENNE, HOLOGYNUM
pxrmt -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode.nwk -n perenne,D07,D02,D13 -c > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode2.nwk
pxrr -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode2.nwk -g perenne > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode2R.nwk

## W/ DECUMBENS+GRANDIFLORUM, BIENNE, HOLOGYNUM
pxrmt -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode.nwk -n perenne,C96,D03,D02,D13 -c > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode3.nwk
pxrr -t ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode3.nwk -g perenne > ./data/Angiosperm353/iqtrees/gtr_g/Astral/gene_bestrees_phyxed/qibl_trees/allbestrees_pruned210_prunedSPP_completeCases_lbNode3R.nwk


nano QuIBL/myparameters.txt #save changes to QuIBL/myparameters_lbNode1.txt
nano QuIBL/myparameters.txt #save changes to QuIBL/myparameters_lbNode2.txt
nano QuIBL/myparameters.txt #save changes to QuIBL/myparameters_lbNode3.txt

python QuIBL/QuIBL.py ./QuIBL/myparameters_lbNode1.txt
python QuIBL/QuIBL.py ./QuIBL/myparameters_lbNode2.txt
python QuIBL/QuIBL.py ./QuIBL/myparameters_lbNode3.txt