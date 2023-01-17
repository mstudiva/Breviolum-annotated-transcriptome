# Transcriptome Annotation, version January 11, 2023
# Created by Misha Matz (matz@utexas.edu), modified by Michael Studivan (studivanms@gmail.com) for use on FAU's HPC (KoKo)
# for use in generating transcriptome annotations for Breviolum


#------------------------------
# BEFORE STARTING, replace, in this whole file:
#	- studivanms@gmail.com by your actual email;
#	- mstudiva with your KoKo user name.

# The idea is to copy the chunks separated by empty lines below and paste them into your cluster
# terminal window consecutively.

# The lines beginning with hash marks (#) are explanations and additional instructions -
# please make sure to read them before copy-pasting.


#------------------------------
# setup

# To install Bioperl as a conda environment
conda create -y -n bioperl perl-bioperl

# getting scripts
cd ~/bin
git clone https://github.com/z0on/annotatingTranscriptomes.git
mv annotatingTranscriptomes/* .
rm -rf annotatingTranscriptomes
rm launcher_creator.py

git clone https://github.com/z0on/emapper_to_GOMWU_KOGMWU.git
mv emapper_to_GOMWU_KOGMWU/* .
rm -rf emapper_to_GOMWU_KOGMWU

git clone https://github.com/mstudiva/Breviolum-annotated-transcriptome.git
mv Breviolum-annotated-transcriptome/* .
rm -rf Breviolum-annotated-transcriptome

# creating backup directory
mkdir backup

# creating annotation directory
cd
mkdir annotate
cd annotate


#------------------------------
# getting transcriptomes

# Rivera (April 2020)
wget https://sites.bu.edu/davieslab/files/2020/04/B_psygmophilum_transcriptome.zip
unzip B_psygmophilum_transcriptome.zip

# uses fast-x toolkit to wrap each line to 60 characters
module load fastx-toolkit-0.0.14-gcc-8.3.0-ombppo2
srun fasta_formatter -i B_psygmophilum_transcriptome/B_psygmophilum_transcriptome.fasta -w 60 -o Breviolum.fasta

# use the stream editor to find and replace the first instance of component designations with the species name in each line
sed -i 's/Sym_/Breviolum/' Breviolum.fasta

# Camp (June 2021)
# this transcriptome has already been converted to protein translations, so skip to line 143
# download from https://osf.io/gsn8p/ and scp to your annotate directory
mv B1.annotated.fa Breviolum.fasta
sed -i 's/TRINITY_DN/Breviolum/' Breviolum.fasta

# Chen (2019)
# get from https://espace.library.uq.edu.au/view/UQ:8279c9a
gunzip Breviolum_minutum.CDS.fna.gz
mv Breviolum_minutum.CDS.fna Breviolum.fasta
sed -i 's/Bmin.gene/Breviolum/' Breviolum.fasta

# Avila-Magana (August 2021)
# from https://datadryad.org/stash/dataset/doi:10.5061/dryad.k3j9kd57b
gunzip Symbio_Sider.fna.gz
mv gunzip Symbio_Sider.fna Breviolum.fasta
sed -i 's/TRINITY_DN/Breviolum/' Breviolum.fasta

# transcriptome statistics
conda activate bioperl
echo "seq_stats.pl Breviolum.fasta > seqstats_Breviolum.txt" > seq_stats
launcher_creator.py -j seq_stats -n seq_stats -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch seq_stats.slurm

# transcriptome statistics from https://sites.bu.edu/davieslab/data-code/
Breviolum.fasta (Davies)
-------------------------
31970 sequences.
1291 average length.
8436 maximum length.
500 minimum length.
N50 = 1458
41.3 Mb altogether (41260850 bp).
0 ambiguous Mb. (0 bp, 0%)
0 Mb of Ns. (0 bp, 0%)
-------------------------

Breviolum.fasta (Camp)
-------------------------
98555 sequences.
335 average length.
8046 maximum length.
99 minimum length.
N50 = 490
33 Mb altogether (33040441 bp).
25.4 ambiguous Mb. (25354444 bp, 76.7%)
0.9 Mb of Ns. (945134 bp, 2.9%)
-------------------------


#------------------------------
# uniprot annotations with blast

# getting uniprot_swissprot KB database
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz

# unzipping
gunzip uniprot_sprot.fasta.gz &

# indexing the fasta database
echo "makeblastdb -in uniprot_sprot.fasta -dbtype prot" >mdb
launcher_creator.py -j mdb -n mdb -q shortq7 -t 6:00:00 -e studivanms@gmail.com
sbatch mdb.slurm

# splitting the transcriptome into 100 chunks
splitFasta.pl Breviolum.fasta 100

# blasting all 200 chunks to uniprot in parallel, 4 cores per chunk
ls subset* | perl -pe 's/^(\S+)$/blastx -query $1 -db uniprot_sprot\.fasta -evalue 0\.0001 -num_threads 4 -num_descriptions 5 -num_alignments 5 -out $1.br/'>bl
launcher_creator.py -j bl -n blast -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch blast.slurm

# watching progress:
grep "Query= " subset*.br | wc -l
# you should end up with the same number of queries as sequences from the seq_stats script

# combining all blast results
cat subset*br > myblast.br
mv subset* ~/annotate/backup/

# for trinity-assembled transcriptomes: annotating with isogroups
grep ">" Breviolum.fasta | perl -pe 's/>Breviolum(\d+)(\S+)\s.+/Breviolum$1$2\tBreviolum$1/'>Breviolum_seq2iso.tab
cat Breviolum.fasta | perl -pe 's/>Breviolum(\d+)(\S+).+/>Breviolum$1$2 gene=Breviolum$1/'>Breviolum_iso.fasta


#-------------------------
# extracting coding sequences and corresponding protein translations:
echo "perl ~/bin/CDS_extractor_v2.pl Breviolum_iso.fasta myblast.br allhits bridgegaps" >cds
launcher_creator.py -j cds -n cds -l cddd -t 6:00:00 -q shortq7 -e studivanms@gmail.com
sbatch cddd


#------------------------------
# GO annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

# the Camp transcriptome already has protein translations, so run this first before proceeding
# cat Breviolum.fasta | perl -pe 's/>Breviolum(\d+)(\S+).+/>Breviolum$1$2 gene=Breviolum$1/'>Breviolum_iso_PRO.fas

# selecting the longest contig per isogroup
fasta2SBH.pl Breviolum_iso_PRO.fas >Breviolum_out_PRO.fas

# scp your *_out_PRO.fas file to laptop, submit it to
http://eggnog-mapper.embl.de
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*_out_PRO.fas .

# copy link to job ID status and output file, paste it below instead of current link:
# symB (Davies) status: http://eggnog-mapper.embl.de/job_status?jobname=MM_7hc47__h
# symB (Camp) status: http://eggnog-mapper.embl.de/job_status?jobname=MM_0a0rsfk1

# once it is done, download results to HPC:
wget http://eggnog-mapper.embl.de/MM_7hc47__h/out.emapper.annotations # Davies
wget http://eggnog-mapper.embl.de/MM_0a0rsfk1/out.emapper.annotations # Camp

# GO:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$10 }' out.emapper.annotations | grep GO | perl -pe 's/,/;/g' >Breviolum_iso2go.tab
# gene names:
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$8 }' out.emapper.annotations | grep -Ev "\tNA" >Breviolum_iso2geneName.tab


#------------------------------
# KOG annotation
# updated based on Misha Matz's new GO and KOG annotation steps on github: https://github.com/z0on/emapper_to_GOMWU_KOGMWU

cp ~/bin/kog_classes.txt .

#  KOG classes (single-letter):
awk -F "\t" 'BEGIN {OFS="\t" }{print $1,$7 }' out.emapper.annotations | grep -Ev "[,#S]" >Breviolum_iso2kogClass1.tab
# converting single-letter KOG classes to text understood by KOGMWU package (must have kog_classes.txt file in the same dir):
awk 'BEGIN {FS=OFS="\t"} NR==FNR {a[$1] = $2;next} {print $1,a[$2]}' kog_classes.txt Breviolum_iso2kogClass1.tab > Breviolum_iso2kogClass.tab


#------------------------------
# KEGG annotations:

# selecting the longest contig per isogroup:
srun fasta2SBH.pl Breviolum_iso.fasta >Breviolum_4kegg.fasta

# scp *4kegg.fasta to your laptop
cd /path/to/local/directory
scp mstudiva@koko-login.hpc.fau.edu:~/path/to/HPC/directory/*4kegg.fasta .
# use web browser to submit _4kegg.fasta file to KEGG's KAAS server http://www.genome.jp/kegg/kaas/
# select SBH method, upload nucleotide query
# symB (Davies) status: https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1672971915&key=JaGOMN7D
# symB (Camp) status: https://www.genome.jp/kaas-bin/kaas_main?mode=user&id=1673036132&key=n9SFDxIG

# Once it is done, download to HPC - it is named query.ko by default
wget https://www.genome.jp/tools/kaas/files/dl/1672971915/query.ko # Davies
wget https://www.genome.jp/tools/kaas/files/dl/1673035930/query.ko # Camp

# selecting only the lines with non-missing annotation:
cat query.ko | awk '{if ($2!="") print }' > Breviolum_iso2kegg.tab

# the KEGG mapping result can be explored for completeness of transcriptome in terms of genes found,
# use 'html' output link from KAAS result page, see how many proteins you have for conserved complexes and pathways, such as ribosome, spliceosome, proteasome etc


#------------------------------
# file transfer

# copy all files to local machine
cd /path/to/local/directory
scp mstudiva@koko-login.fau.edu:~/path/to/HPC/directory/* .
