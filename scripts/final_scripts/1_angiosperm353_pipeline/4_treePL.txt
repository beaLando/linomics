# treePL

## convergence is not reached when calibrating tree with linum fossil + either (or both) of the older calibration points (linoids + phyllantoids and crown node - euphorbioids)
## convergence is reached with linum calibration alone or the other two calibration points together
## NB: configuration files should be in data folder, but I added them to treePL scripts folder for upload to github

#mkdir ./output/Angiosperm353/treePL/


#PRIMING
##SUPERMX ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN.log


#DATING
##SUPERMX ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s0_linumCN_primed.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/supermx/config_supermx_s100_linumCN_primed.log




#PRIMING
##ASTRAL ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN.log


#DATING
##ASTRAL ROOTED
/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s0_linumCN_primed.log

/home/linuxbrew/.linuxbrew/Cellar/treepl/2022.04.06/bin/treePL ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN_primed.txt > ./data/Angiosperm353/iqtrees/gtr_i_g/Astral/config_astral_s100_linumCN_primed.log





