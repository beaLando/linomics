# LINUM BIENNE GENOMICS

# OVERVIEW & DATA AVAILABILITY

This is a repository to conduct phylogenomic and phylogeographic analyses in Linum and Linum bienne based on nuclear and plastome data in relation to the following publication: XXXX.
Note that another repository (XXXX) deals with demographic analyses based on nuclear data for the same publication.

Most genomic data have been produced specifically for this study and have been uploaded under project XXXXX of the European Nucleotide Archive (ENA).
A subset of Linum species were available online and data can be found at XXXX.
Only a subset of genomic data for 7 L. bienne accessions was previously uploaded to NCBI (project XXXX).
Please refer to Supplementary Tables XX and XX for ENA, NCBI, and Tree of Life sample IDs.

This repo in particular contains scripts to:

`00.` Clean accession data
`01.` Execute a phylogenomic reconstruction of the genus Linum based on nuclear (angiosperms353) and plastid genes via maximum likelihood, coalescence, and network methods
`02-03.` Execute a phylogeographic reconstruction of L. bienne across its distribution based on whole plastomes, as well as occurrence (GBIF) and climate data.
`04.` Downstream processing and plotting of outputs from above analyses (some of these scripts also make use of outputs from demographic pipeline at XXXX)


# DIRECTORY STRUCTURE & PACKAGES

Scripts are organised into four numbered folders under `scripts/`, each corresponding to a major analysis step making use of R or shell scripts, and Julia in a couple of instances. Within `scripts/` subdirectories, individual scripts are also numbered so that everything should be run in the order indicated by the numbers since intermediate outputs are produced which are used in downstream steps.

Steps `00`, `03`, and `04` are run entirely locally in R Studio. Steps `01` and `02` are run mostly on a server (SLURM), with some local post-processing. Julia is only used in steps `01` and `04`. For shell scripts, conda was used to install and manage most programs, with some exceptions.

Three symbols are used throughout:
- 📥 **load manually** — file or directory must be provided by the user before running
- 📦 **load from this repo** — file is provided in this repository
- ⚙️ **generated** — file or directory is created by the pipeline

---

## Overview of top-level structure

```
workDir/
├── scripts/                       # 📦 all analysis scripts (this repository)
│   ├── 00/                        # coordinate correction (R, local)
│   ├── 01/                        # phylogenomics: nuclear + plastid (shell + R, server/local)
│   ├── 02/                        # plastome mapping + BEAST dating (shell + R, server/local)
│   ├── 03/                        # niche modelling (R, local)
│   └── 04/                        # output processing and plotting (R + Julia, local)
│
├── [programs/]                    # 📥 third-party programs (see per-step package lists below)
│
└── data/                          # all input and output data (see below)
```

---

## scripts/00 — Sample coordinate correction

**Language:** R (local, R Studio)

**R packages:**
`tidyverse`, `ggplot2`, `leaflet`, `GGally`, `parzer`, `xlsx`

> Run **once before anything else**. The corrected Excel file is used as input throughout
> `scripts/02` and `scripts/04`.

```
workDir/
└── data/
    ├── populations_locations_1.xlsx            # 📥 raw sample metadata (uncorrected coordinates)
    └── populations_locations_1_coordFixed.xlsx # ⚙️ coordinate-corrected output
```

---

## scripts/01 — Phylogenomics (nuclear & plastid)

**Language:** shell (server, conda), R (local, R Studio), Julia (local, conda)

### Packages — QC and read trimming (scripts 00, server)
`fastqc`, `fastp`, `multiqc`

### Packages — HybPiper target capture assembly (scripts 01–03, 07–09, server)
`hybpiper`

### Packages — Alignment and tree building (scripts 04–06, 10, local)
`seqkit`, `mafft`, `trimal`, `iqtree`, `phyx`,
`ConcatFasta.py` — https://github.com/santiagosnchez/ConcatFasta

### Packages — Coalescence, dating and conflict analyses (scripts 05–06, 10–12, local)
`Astral-III` — https://github.com/smirarab/ASTRAL  
`phyparts` — https://bitbucket.org/blackrim/phyparts  
`phypartspiecharts.py` — https://github.com/mossmatters/MJPythonNotebooks  
`treePL`  
`phytop`  
`PhyloNet` — https://github.com/NakhlehLab/PhyloNet  
`QuIBL.py` — https://github.com/miriammiyagi/QuIBL  
`PhyloNetworks` *(Julia)*, `PhyloPlots` *(Julia)*, `RCall` *(Julia)*

### Server-side directories and files (scripts 00–03, 07–09)

> Create the directories marked below before submitting SLURM jobs.
> Raw reads should be placed in the appropriate `*_fqs/` directories before running.
> `gene_sequences_nuclear/` and `gene_sequences_plastid/` must be downloaded locally
> after steps 03 and 09 respectively, before continuing with local steps.

```
workDir/
│
├── logs/                               # 📥 create in advance (receives SLURM logs)
│
├── ang353_NewTargets/                  # 📥 NewTargets scripts and files
│   ├── filter_mega353.py               #    https://github.com/chrisjackson-pellicle/NewTargets
│   ├── mega353.fasta                   #    https://github.com/chrisjackson-pellicle/NewTargets
│   └── select_file.txt                 # 📦 taxon selection file for filtering
│
├── angiosperm353_paftol_fqs/           # 📥 raw reads tp download from Tree of Life
│   └── SraRunTable_downloadLinks.txt   # 📦 ENA download links (see Supplementary Table XX)
│
├── angiosperm353_juan_fqs/             # 📥 raw reads generated for this project (download from ENA)
│
├── trimmed_fqs/                        # ⚙️ created by pipeline (receives trimmed reads)
│
├── nuclear_targets.fasta               # ⚙️ filtered Angiosperm353 target file (script 01)
├── plastid_targets.faa                 # 📥 protein target file for plastid assembly
├── namelist.txt                        # ⚙️ sample accession names (script 01)
├── namelist_nuclear.txt                # 📦 as above, with "_nuclear" suffix
├── namelist_plastid.txt                # 📦 as above, with "_plastid" suffix
│
├── <SAMPLE>_nuclear/                   # ⚙️ HybPiper output per sample, nuclear
├── <SAMPLE>_plastid/                   # ⚙️ HybPiper output per sample, plastid
├── stats_nuclear/                      # ⚙️ HybPiper QC stats, nuclear (script 02)
├── stats_plastid/                      # ⚙️ HybPiper QC stats, plastid (script 08)
├── gene_sequences_nuclear/             # ⚙️ retrieved nuclear sequences (script 03)
└── gene_sequences_plastid/             # ⚙️ retrieved plastid sequences (script 09)
```

### Local directories and files (scripts 04–06, 10–12)

> `gene_sequences_nuclear/` and `gene_sequences_plastid/` must be downloaded from the server
> into `data/Angiosperm353/hybpiper/` before running these steps.

```
workDir/
├── ConcatFasta.py                      # 📥 https://github.com/santiagosnchez/ConcatFasta
├── phypartspiecharts.py                # 📥 https://github.com/mossmatters/MJPythonNotebooks
├── PhyloNetv3_8_2.jar                  # 📥 https://github.com/NakhlehLab/PhyloNet
├── Astral/
│   └── astral.5.7.8.jar                # 📥 https://github.com/smirarab/ASTRAL
├── phyparts/
│   └── target/
│       └── phyparts-*.jar              # 📥 https://bitbucket.org/blackrim/phyparts
│
├── QuIBL/
│   ├── QuIBL.py                        # 📥 https://github.com/miriammiyagi/QuIBL
│   ├── myparameters.txt                # 📦 QuIBL parameter file (genus-level analysis)
│   ├── myparameters_lbNode1.txt        # 📦 QuIBL parameter file (L. bienne node, all taxa)
│   ├── myparameters_lbNode2.txt        # 📦 QuIBL parameter file (L. bienne node, subset 1)
│   └── myparameters_lbNode3.txt        # 📦 QuIBL parameter file (L. bienne node, subset 2)
│
└── data/
    └── Angiosperm353/
        ├── hybpiper/
        │   ├── gene_sequences_nuclear/ # 📥 download from server after script 03
        │   ├── gene_sequences_plastid/ # 📥 download from server after script 09
        │   ├── genes2remove_nuclear.txt # 📦 manually curated gene exclusion list based on hybpiper stats output
        │   └── genes2remove_plastid.txt # 📦 manually curated gene exclusion list based on hybpiper stats output
        │
        ├── aln_nuclear/                # ⚙️ nuclear alignments (script 04)
        ├── aln_plastid/                # ⚙️ plastid alignments (script 10)
        │
        └── iqtrees/
            ├── nuclear_gtr_g/          # ⚙️ nuclear ML trees and ASTRAL results (scripts 05–06, 11–12)
            │   ├── supermx/
            │   │   └── config_supermx_s*_linumCN*.txt  # 📦 treePL config files
            │   └── Astral/
            │       ├── config_astral_s*_linumCN*.txt   # 📦 treePL config files
            │       └── nodeLB/         # ⚙️ PhyloNet conflict analyses (script 11)
            └── plastid/                # ⚙️ plastid ML trees (script 10)
                ├── config_supermx_s*_linumCN*.txt      # 📦 treePL config files
                ├── partitions2test.txt                 # 📦 topology test partition file
                └── topologies2test_lb.nwk              # 📦 alternative topologies for AU test
```

---

## scripts/02 — L. bienne whole-plastome mapping & BEAST dating

**Language:** shell (server, conda), R (local, R Studio)

### Packages — Read mapping, variant calling and consensus (scripts 01, 03–04, server)
`fastqc`, `fastp`, `multiqc`, `bowtie2`, `samtools`, `qualimap`, `picard`, `bcftools`,
`vcftools`, `seqkit`

### Packages — Alignment and tree building (script 01, server)
`mafft`, `trimal`, `iqtree`, `phyx`

### Packages — BEAST divergence time estimation (scripts 03–04, server)
`beast`

### Packages — VCF quality control and visualisation (scripts 02–03, local)
`tidyverse` *(R)*, `vcfR` *(R)*

### Server-side directories and files (scripts 01, 03–04)

> Create the directories marked below before submitting SLURM jobs.

```
workDir/   [= Lus_Plastome_paper2024/]
│
├── logs/                               # 📥 create in advance
│
├── refs/                               # 📥 create in advance; load reference files here
│   ├── NC_036356_1_Lus_chloroplast_ref.fasta  # 📥 L. usitatissimum chloroplast reference
│   └── annotation.gff3                        # 📥 chloroplast genome annotation
│
├── raw_fasq_files_plastome_linum/      # 📥 raw L. bienne plastome reads
│
├── beast_analysis/                     # 📦 BEAST XML config files
│   ├── run_empty_plastid.xml           #    prior-only run, plastid data
│   ├── run_plastid.xml                 #    full BEAST run, plastid data
│   ├── run_empty_superMX5.xml          #    prior-only run, nuclear supermatrix
│   └── run_superMX5.xml               #    full BEAST run, nuclear supermatrix
│
└── mapping/                            # ⚙️ created by pipeline
    └── all/                            # ⚙️ merged VCF, consensus fastas, alignments, trees
```

### Local directories and files (scripts 02–03)

> vcftools stats files, the filtered VCF, and BEAST outputs must be downloaded from
> the server before running local steps.

```
workDir/
└── data/
    └── Lb_plastomes_raw/
        ├── refs/                       # 📥 local copies of server reference files
        └── Consensus/
            ├── stats_vcf_unfiltered/   # 📥 vcftools stats files downloaded from server
            ├── aln/all/
            │   ├── all_plastomes.vcf.gz        # 📥 filtered VCF downloaded from server
            │   └── newRes2026/                 # 📥 BEAST output files downloaded from server
            ├── genetic_distances_tree_nuclear.nwk  # 📥 collaborator input (nuclear distance tree)
            └── plastid_genes_extractedVillarianum/ # ⚙️ manually aligned villarianum sequences (script 02)
```

---

## scripts/03 — GBIF occurrence data, climate data & niche modelling

**Language:** R (local, R Studio)

**R packages:**
`tidyverse`, `ggplot2`, `ggfortify`, `patchwork`, `rgbif`, `CoordinateCleaner`, `maps`,
`countrycode`, `sf`, `raster`, `fuzzySim`, `modEvA`, `dismo`, `maxnet`, `gam`,
`randomForest`, `gbm`, `corrplot`, `ecospat`, `blockCV`, `plotmo`, `sdmpredictors`, `pastclim`

> All scripts run locally in R Studio.
> WorldClim rasters and the countries shapefile must be downloaded manually before running.
> Barreto2023 palaeoclimate data are downloaded automatically by `pastclim` on first run
> (requires internet connection and several GB of disk space).

```
workDir/
└── data/
    └── niche_model/
        ├── species_occurrence/         # ⚙️ GBIF occurrence data (scripts 01–02)
        ├── climate/
        │   ├── wc2.1_10m_bio/          # 📥 WorldClim 2.1 bioclimatic variables (10 min)
        │   │                           #    https://www.worldclim.org/data/bioclim.html
        │   └── past/                   # ⚙️ palaeoclimate data (downloaded by pastclim, script 08)
        └── countries/
            └── world_countries.shp     # 📥 world countries shapefile (+ sidecar files)
```

> All niche model outputs (`.rds`, `.csv`, `.RData`, `.tif` files) are written directly
> into `data/niche_model/` and generated by scripts 04–08.

---

## scripts/04 — Output processing and plotting

**Language:** R (local, R Studio), Julia (local, conda `snaq_env`)

**R packages:**
`tidyverse`, `ggplot2`, `ggpubr`, `patchwork`, `gridExtra`, `ggnewscale`, `ggtree`, `quiblR`,
`ape`, `hash`, `tanggle`, `geiger`, `phytools`, `treeio`, `deeptime`, `rwty`, `raster`,
`geosphere`, `xlsx`

**Julia packages:**
`PhyloNetworks`, `PhyloPlots`, `RCall`

> All scripts run locally in R Studio, except `02b` which runs in Julia.
> All input files are generated by earlier pipeline steps or provided in the repo.
> BEAST output files must be downloaded from the server (see `scripts/02` above).

```
workDir/
└── data/
    ├── populations_locations_1_coordFixed.xlsx  # 📦/⚙️ from scripts/00
    │
    ├── Angiosperm353/iqtrees/          # ⚙️ all trees and QuIBL outputs from scripts/01
    │
    ├── Lb_plastomes_raw/Consensus/aln/all/newRes2026/  # 📥 BEAST outputs from server
    │
    └── Lb_genomes/
        └── yann_diversity_analysis/    # 📥 collaborator input files (diversity analysis)
            ├── Distance_matrix_new.txt
            └── list_individuals_tokeep_selection.txt
```
