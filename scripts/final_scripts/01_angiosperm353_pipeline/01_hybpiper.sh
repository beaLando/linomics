# ENVIRONMENTS 
## readNmap
conda create -n readNmap fastqc multiqc trimmomatic fastp bowtie2 samtools bcftools star qualimap conda-forge::parallel

## hybpiper
conda create hybpiper hybpiper 
conda activate hybpiper
conda install bioconda::rename conda-forge::biopython anaconda::pandas
conda deactivate


# GET DATA from PAFTOL
## NB: SraRunTable_downloadLinks.txt can be created from excel file SraRunTable.xlsx in same scripts folder (Table 2)

conda activate readNmap

cd ./angiosperm353_paftol_fqs/ #I prefer to be inside directory to avoid problems with wget and parallel command and paths

wget -i SraRunTable_downloadLinks.txt

nohup parallel --joblog parallel_jobs.log --nice 10 'bzcat {} | gzip -c > {.}.gz' ::: *bz2 & 
rm -r *.bz2

cd ..


# CHECK QUALITY of JUAN & PAFTOL's READS & CLEAN
## CHECK before trimming
nohup fastqc ./angiosperm353_juan_fqs/*fastq.gz -o ./angiosperm353_juan_fqs -t 3 -noextract & #jobID: 25459 #there seem to be adapters + instances of polyAs or polyGs
nohup fastqc ./angiosperm353_paftol_fqs/*fq.gz -o ./angiosperm353_paftol_fqs -t 3 -noextract & #jobID: 25384 #maybe polyAs or polyGs

multiqc -ip ./angiosperm353_juan_fqs/*_fastqc.zip ./angiosperm353_paftol_fqs/*_fastqc.zip -o ./multiqc_angiosperm353_juan_paftol

## CLEANING with FASTp (instead of trimmomatic) 
### Need to trim polyA(X) and polyG from all files, and adapters from Juan's only. Duplicate reads do not need to be removed.
### fastp filters by quality and length automatically (to disable -Q -L) > for now I leave it on and see what happens, it seems ok

for first in `ls ./angiosperm353_juan_fqs/*_1.fastq.gz`;
do
    second=`echo ${first} | sed 's/_1/_2/g'`
    first_out=${first%".fastq.gz"}
    second_out=${second%".fastq.gz"}

    fastp --in1 ${first} --in2 ${second} --out1 ${first_out}_trimmed.fastq.gz --out2 ${second_out}_trimmed.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --html fastp_juan.html --json fastp_juan.json
done


mv fastp_juan* ./angiosperm353_juan_fqs

for first in `ls ./angiosperm353_paftol_fqs/*_1.fq.gz`;
do
    second=`echo ${first} | sed 's/_1/_2/g'`
    first_out=${first%".fq.gz"}
    second_out=${second%".fq.gz"}

    fastp --in1 ${first} --in2 ${second} --out1 ${first_out}_trimmed.fastq.gz --out2 ${second_out}_trimmed.fastq.gz --detect_adapter_for_pe --trim_poly_g --trim_poly_x --html fastp_paftol.html --json fastp_paftol.json
done

mv fastp_paftol* ./angiosperm353_paftol_fqs

## CHECK after timming
nohup fastqc ./angiosperm353_juan_fqs/*_trimmed.fastq.gz -o ./angiosperm353_juan_fqs -t 3 -noextract & #jobID: 74531
nohup fastqc ./angiosperm353_paftol_fqs/*_trimmed.fastq.gz -o ./angiosperm353_paftol_fqs -t 3 -noextract & #jobID: 74703

multiqc -ip ./angiosperm353_juan_fqs/*trimmed_fastqc.zip ./angiosperm353_paftol_fqs/*trimmed_fastqc.zip -o ./multiqc_angiosperm353_juan_paftol

conda deactivate


# HYBPIPER - RETRIEVE ANGIOSPERM353 GENES
conda activate hybpiper

## organize directories
mkdir ./trimmed_fqs
mv ./angiosperm353_juan_fqs/*_trimmed.fastq.gz ./trimmed_fqs
mv ./angiosperm353_paftol_fqs/*_trimmed.fastq.gz ./trimmed_fqs

rename 's/-read_/_/g' ./trimmed_fqs/*

## create list of accession names
for fq in `ls ./trimmed_fqs/*_1_trimmed.fastq.gz`;
do
    fl=$(basename -- "$fq")
    accession=${fl%"_1_trimmed.fastq.gz"}
    echo $accession | cat >> namelist.txt
done

## create target genes file by filtering angiosperm353 whole catalog
### see here https://hackmd.io/@mossmatters/Sk82Jujhc#15-filtering-and-correcting-your-target-file-with-%60hybpiper-fix_targetfile%60
### and see here: https://github.com/chrisjackson-pellicle/NewTargets
### for my purpose I only used reference sequences from linaceae, phyllantaceae, and euphorbiaceae

python ./ang353_NewTargets/filter_mega353.py ./ang353_NewTargets/mega353.fasta ./ang353_NewTargets/select_file.txt
mv ./ang353_NewTargets/filtered_target_file.fasta ./filtered_target_file.fasta

dos2unix filtered_target_file.fasta
dos2unix namelist.txt #in scripts folder

## run hybpiper
### NB1: both blastx and bwa have problem when file path contains space, so hybpiper does not work if data are accessed from googleDrive "Shared drives" for example
### NB2: while loops are a bit tricky to run in background, if at some point you realize you want to stop everything.
###      to stop everything, do "ps xjf" to see all parallel processess and their parents and then kill, possibly a few times, most "parent" processes (kill -9 PID xxx)

#line below is if you have problems with paths to sample fqs (you can do them one by one)
#hybpiper assemble -t_dna filtered_target_file.fasta -r ./trimmed_fqs/C95_1_trimmed.fastq.gz ./trimmed_fqs/C95_2_trimmed.fastq.gz --prefix C95 --cpu 10

while read name; 
do hybpiper assemble -t_dna filtered_target_file.fasta -r ./trimmed_fqs/${name}_1_trimmed.fastq.gz ./trimmed_fqs/${name}_2_trimmed.fastq.gz --prefix $name --cpu 10; 
done < namelist.txt &

hybpiper stats -t_dna filtered_target_file.fasta gene namelist.txt
hybpiper recovery_heatmap seq_lengths.tsv
hybpiper paralog_retriever namelist.txt -t_dna filtered_target_file.fasta
hybpiper retrieve_sequences dna -t_dna filtered_target_file.fasta --sample_names namelist.txt --stats_file ./stats/hybpiper_stats.tsv --filter_by GenesAt50pct greater 49 --fasta_dir gene_sequences

conda deactivate