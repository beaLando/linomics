# treePL

## NB1: convergence is not reached when calibrating tree with linum fossil + either (or both) of the older calibration points (linoids + phyllantoids and crown node - euphorbioids)
## convergence is reached with linum calibration alone or the other two calibration points together
## NB2: configuration files should be in data folder, but I added them to treePL scripts folder for upload to github

#mkdir ./output/Angiosperm353/treePL/

#GTR+I+G
## PRIMING
### SUPERMX ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN.log


## DATING
### SUPERMX ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN_primed.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN_primed.log




## PRIMING
### ASTRAL ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN.log


## DATING
### ASTRAL ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN_primed.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN_primed.log





#GTR+G
## PRIMING
### SUPERMX ROOTED
treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s0_linumCN.log #/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL

treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s100_linumCN.log


## DATING
### SUPERMX ROOTED
treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s0_linumCN_primed.log

treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/supermx/config_supermx_s100_linumCN_primed.log




## PRIMING
### ASTRAL ROOTED
treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s0_linumCN.log

treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s100_linumCN.log


## DATING
### ASTRAL ROOTED
treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s0_linumCN_primed.log

treepl ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/nuclear_gtr_g/Astral/config_astral_s100_linumCN_primed.log





