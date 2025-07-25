# ANG353 LINUM PHYLOGENY BASED on SUPERMATRIX or SINGLE GENE TREES & ASTRAL

##NB1: all trees are made unrooted and then rooted with phyx using both phyllanthus and euphorbia as roots
###NB2: because the space character in 'Shared drives' creates problems in iqtree (does not find files) I run trees in alignment subfolder, then I move trees to output folder, so I circumnavigate problem with space in path
####it would be easier to change folder name and remove space, but this cannot be done with gDrive default folders



## PREPARE OUTPUT DIRECTORIES
#mkdir ./data/Angiosperm353/iqtrees
#mkdir ./data/Angiosperm353/iqtrees/gtr_g #in the end I did not use this model
#mkdir ./data/Angiosperm353/iqtrees/gtr_i_g
#mkdir ./data/Angiosperm353/iqtrees/gtr_i_g/Astral
#mkdir ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts



## SUPERMATRIX TREE
### MAKE TREE
iqtree -nt 5 -s ./data/Angiosperm353/aln/supermx353.mft.tall.nexus -m GTR+I+G -bb 1000 -wbtl -pre ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353 

### ROOT IT
pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.contree -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.phyxed.contree



## GENE TREES & ASTRAL
### MAKE GENE TREES
for aln in `ls ./data/Angiosperm353/aln/*.tall`

do

gene_name=$(basename -s .mft.tall $aln)
echo ${gene_name}

iqtree -nt 5 -s $aln -m GTR+I+G -bb 1000 -wbtl -pre "./data/Angiosperm353/iqtrees/gtr_i_g/"${gene_name}

done


### MAKE CONSENSUS TREE WITH ASTRAL & ROOT IT
cat ./data/Angiosperm353/iqtrees/gtr_i_g/*.contree >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre

head -n1 ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre

java -jar ./Astral/astral.5.7.8.jar -i ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.tre 2> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.log #/home/bland/Programs/Astral/astral.5.7.8.jar
pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.tre -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.phyxed.tre


### ROOT SINGLE GENE TREES 
#### Overall, I use euphorbia to root trees. If euphorbia does not appear, I use phyllanthus to root trees. If none appear, I discard tree (save gene tee name to file).
####NB: There is only one gene without any outgroups, I re-do astral tree without that gene later on


for f in `ls ./data/Angiosperm353/iqtrees/gtr_i_g/*.contree`

do

tree_name=$(basename -s .contree $f)
output_name=${tree_name}_phyxed.contree

if grep -q 'euphorbia' $f; then 
	outgroup_name="euphorbia"
	echo ${tree_name} 
	echo ${outgroup_name}
	pxrr -t $f -g ${outgroup_name} -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/gene_contrees_phyxed/${output_name}

else

if grep -q 'phyllanthus' $f && ! grep -q 'euphorbia' $f; then 
	outgroup_name="phyllanthus"
	echo ${tree_name} 
	echo ${outgroup_name}
	pxrr -t $f -g ${outgroup_name} -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/gene_contrees_phyxed/${output_name}

else

if ! grep -q 'phyllanthus' $f && ! grep -q 'euphorbia' $f; then 
	outgroup_name="NO OUTGROUP found for this gene"
	echo ${tree_name} 
	echo ${outgroup_name}
	echo ${tree_name} >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/genes2remove.txt

fi
fi
fi

done


### EXTRACT DATA & CREATE PLOT TO SHOW AGREEMENT/DISAGREEMENT BETWEEN GENE TREES
#### Re-run Astral tree without genes without outgroups (take genes to remove from ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/genes2remove.txt)
ls ./data/Angiosperm353/iqtrees/gtr_i_g/*.contree | cat $(grep -Ev "6792") >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre

head -n1 ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre

java -jar ./Astral/astral.5.7.8.jar -i ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.tre 2> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.log

pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.tre -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre

#### Extract data on gene node support
java -jar ./phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -d ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/gene_contrees_phyxed/ -m ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre -a 1 -v -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts #/home/bland/Programs/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar

#### Run piecharts tool (if I turn off screen it produces plot, although in current folder, not in output folder?) - https://github.com/mossmatters/MJPythonNotebooks/blob/master/phypartspiecharts.py

conda activate ete3
QT_QPA_PLATFORM=offscreen python3 phypartspiecharts.py ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts/phyparts 253
conda deactivate
