# ANG353 LINUM PHYLOGENY BASED on SUPERMATRIX or SINGLE GENE TREES & ASTRAL
##NB1: all ML trees are made unrooted and then rooted with phyx using both phyllanthus and euphorbia as roots

## PREPARE OUTPUT DIRECTORIES
mkdir -p ./data/Angiosperm353/iqtrees
mkdir -p ./data/Angiosperm353/iqtrees/nuclear_gtr_g 
mkdir -p ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral
mkdir -p ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/phyparts
#mkdir -p ./data/Angiosperm353/iqtrees/gtr_i_g #in the end I did not use this model
#mkdir -p ./data/Angiosperm353/iqtrees/gtr_i_g/Astral
#mkdir -p ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts



# GTR+G model ###########################################################################################################################
conda activate phylo

## SUPERMATRIX TREE
### MAKE TREE
iqtree -nt 5 -s ./data/Angiosperm353/aln_nuclear/supermx353.mft.tall.nexus -m GTR+G -bb 1000 -wbtl -pre ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/supermx353 

### ROOT IT
pxrr -t ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/supermx353.treefile -g ERR2040368_nuclear -o ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/supermx353.phyxed.treefile


## GENE TREES & ASTRAL
### MAKE GENE TREES
for aln in `ls ./data/Angiosperm353/aln_nuclear/*.tall`

do

gene_name=$(basename -s .mft.tall $aln)
echo ${gene_name}

iqtree -nt 5 -s $aln -m GTR+G -bb 1000 -wbtl -pre "./data/Angiosperm353/iqtrees/nuclear_gtr_g/"${gene_name}

done


### MAKE CONSENSUS TREE WITH ASTRAL & ROOT IT
cat ./data/Angiosperm353/iqtrees/nuclear_gtr_g/*.treefile >> ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral_input.tre

head -n1 ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral_input.tre

java -jar ./Astral/astral.5.7.8.jar -i ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral_input.tre -o ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.tre 2> ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.log #/home/bland/Programs/Astral/astral.5.7.8.jar
pxrr -t ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.tre -g ERR2040368_nuclear -o ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.phyxed.tre


### ROOT SINGLE GENE TREES 
#### Overall, I use euphorbia to root trees. If euphorbia does not appear, I use phyllanthus to root trees. If none appear, I discard tree (save gene tree name to file).
####NB: There is only one gene without any outgroups, I re-do astral tree without that gene later on


for f in ./data/Angiosperm353/iqtrees/nuclear_gtr_g/*.treefile; do
    # Skip the loop if no files are found
    [ -e "$f" ] || continue

    tree_name=$(basename -s .treefile "$f")
    output_name="${tree_name}_phyxed.treefile"
    out_dir="./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed"

    if grep -q 'ERR2040368_nuclear' "$f"; then
        outgroup_name="ERR2040368_nuclear"
        echo "${tree_name}: ${outgroup_name}"
        pxrr -t "$f" -g "${outgroup_name}" -o "${out_dir}/${output_name}"

    elif grep -q 'ERR3487350_nuclear' "$f"; then
        outgroup_name="ERR3487350_nuclear"
        echo "${tree_name}: ${outgroup_name}"
        pxrr -t "$f" -g "${outgroup_name}" -o "${out_dir}/${output_name}"

    else
        echo "${tree_name}: NO OUTGROUP found for this gene"
        echo "${tree_name}" >> ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/genes2remove.txt
    fi
done


### EXTRACT DATA & CREATE PLOT TO SHOW AGREEMENT/DISAGREEMENT BETWEEN GENE TREES
#### Extract data on gene node support
java -jar ./phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -d ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/gene_bestrees_phyxed/ -m ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.phyxed.tre -a 1 -v -o ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/phyparts/phyparts #/home/bland/Programs/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar

conda deactivate

#### Run piecharts tool (if I turn off screen it produces plot, although in current folder, not in output folder?) - https://github.com/mossmatters/MJPythonNotebooks/blob/master/phypartspiecharts.py

conda activate ete3
QT_QPA_PLATFORM=offscreen python3 phypartspiecharts.py ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/astral353.phyxed.tre ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/phyparts/phyparts 208
conda deactivate



# GTR+I+G model ###########################################################################################################################
#conda activate phylo

## SUPERMATRIX TREE
### MAKE TREE
#iqtree -nt 5 -s ./data/Angiosperm353/aln/supermx353.mft.tall.nexus -m GTR+I+G -bb 1000 -wbtl -pre ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353 

### ROOT IT
#pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.treefile -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/supermx353.phyxed.treefile



## GENE TREES & ASTRAL
### MAKE GENE TREES
#for aln in `ls ./data/Angiosperm353/aln/*.tall`

#do

#gene_name=$(basename -s .mft.tall $aln)
#echo ${gene_name}

#iqtree -nt 5 -s $aln -m GTR+I+G -bb 1000 -wbtl -pre "./data/Angiosperm353/iqtrees/gtr_i_g/"${gene_name}

#done


### MAKE CONSENSUS TREE WITH ASTRAL & ROOT IT
#cat ./data/Angiosperm353/iqtrees/gtr_i_g/*.treefile >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre

#head -n1 ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre

#java -jar ./Astral/astral.5.7.8.jar -i ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input.tre -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.tre 2> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.log #/home/bland/Programs/Astral/astral.5.7.8.jar
#pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.tre -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353.phyxed.tre


### ROOT SINGLE GENE TREES 
#### Overall, I use euphorbia to root trees. If euphorbia does not appear, I use phyllanthus to root trees. If none appear, I discard tree (save gene tee name to file).
####NB: There is only one gene without any outgroups, I re-do astral tree without that gene later on


#for f in ./data/Angiosperm353/iqtrees/gtr_i_g/*.treefile; do
    # Skip the loop if no files are found
    #[ -e "$f" ] || continue

    #tree_name=$(basename -s .treefile "$f")
    #output_name="${tree_name}_phyxed.treefile"
    #out_dir="./data/Angiosperm353/iqtrees/gtr_i_g/Astral/gene_bestrees_phyxed"

    #if grep -q 'euphorbia' "$f"; then
        #outgroup_name="euphorbia"
        #echo "${tree_name}: ${outgroup_name}"
        #pxrr -t "$f" -g "${outgroup_name}" -o "${out_dir}/${output_name}"

    #elif grep -q 'phyllanthus' "$f"; then
        #outgroup_name="phyllanthus"
        #echo "${tree_name}: ${outgroup_name}"
        #pxrr -t "$f" -g "${outgroup_name}" -o "${out_dir}/${output_name}"

    #else
        #echo "${tree_name}: NO OUTGROUP found for this gene"
        #echo "${tree_name}" >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/genes2remove.txt
    #fi
#done


### EXTRACT DATA & CREATE PLOT TO SHOW AGREEMENT/DISAGREEMENT BETWEEN GENE TREES
#### Re-run Astral tree without genes without outgroups (take genes to remove from ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/genes2remove.txt)
#ls ./data/Angiosperm353/iqtrees/gtr_i_g/*.treefile | cat $(grep -Ev "6792") >> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre

#head -n1 ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre

#java -jar ./Astral/astral.5.7.8.jar -i ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral_input1.tre -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.tre 2> ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.log

#pxrr -t ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.tre -g euphorbia -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre

#### Extract data on gene node support
#java -jar ./phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar -d ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/gene_bestrees_phyxed/ -m ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre -a 1 -v -o ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts #/home/bland/Programs/phyparts/target/phyparts-0.0.1-SNAPSHOT-jar-with-dependencies.jar
#conda deactivate

#### Run piecharts tool (if I turn off screen it produces plot, although in current folder, not in output folder?) - https://github.com/mossmatters/MJPythonNotebooks/blob/master/phypartspiecharts.py
#conda activate ete3
#QT_QPA_PLATFORM=offscreen python3 phypartspiecharts.py ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/astral353_1.phyxed.tre ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/phyparts/phyparts 253
#conda deactivate
