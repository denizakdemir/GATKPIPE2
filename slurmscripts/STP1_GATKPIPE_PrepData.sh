# This code unzips the fastq files and creates a list of the fastq files to be used in the next steps of the pipeline

# This is a bash script to be run in the terminal
# to run this script type the following command in the terminal, from the directory containing the script:
# bash slurmscripts/STP1_GATKPIPE_PrepData.sh

# directory structure:
# .
# |____0_index
# | |____README.md
# |____4_processing
# | |____README.md
# |____README.md
# |____1_data
# | |____reference3D7
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.1.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.9.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.9.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.toplevel.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.12.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.9.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.2.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.12.fa.gz
# | | |____CHECKSUMS
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.toplevel.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.2.fa.gz
# | | |____README
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.1.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.12.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.4.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.2.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.1.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_rm.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna_sm.chromosome.4.fa.gz
# | | |____Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.chromosome.4.fa.gz
# | |____README.md
# | |____referenceIPO323
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.12.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.4.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.toplevel.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.18.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.1.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.15.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.toplevel.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.14.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.21.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.8.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.15.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.10.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.1.fa.gz
# | | |____CHECKSUMS
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.14.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.9.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.4.fa.gz
# | | |____README
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.2.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.14.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.20.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.18.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.4.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.16.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.11.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.2.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.6.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.15.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.2.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.3.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.12.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.21.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.19.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.5.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.1.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.17.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.18.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.13.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.9.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.7.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_sm.chromosome.21.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna.chromosome.12.fa.gz
# | | |____Zymoseptoria_tritici.MG2.dna_rm.chromosome.9.fa.gz
# | |____X204SC22100973-Z01-F001.csv
# | |____X204SC22100973-Z01-F001.tar
# | |____X204SC22100973-Z01-F002.tar
# |____2_fastqc
# | |____README.md
# |____.gitignore
# |____slurmscripts
# | |____README.md
# |____documentation
# | |____README.md
# |____.git
# | |____COMMIT_EDITMSG
# | |____config
# | |____hooks
# | | |____pre-commit.sample
# | | |____update.sample
# | | |____applypatch-msg.sample
# | | |____pre-push.sample
# | | |____post-update.sample
# | | |____pre-rebase.sample
# | | |____commit-msg.sample
# | | |____pre-applypatch.sample
# | | |____prepare-commit-msg.sample
# | |____index
# | |____logs
# | | |____refs
# | | | |____heads
# | | | | |____main
# | | | |____remotes
# | | | | |____origin
# | | | | | |____main
# | | | | | |____HEAD
# | | |____HEAD
# | |____packed-refs
# | |____info
# | | |____exclude
# | |____description
# | |____objects
# | | |____pack
# | | |____a1
# | | | |____5522a188f55f4e1d41cf20c5427e24c231d47f
# | | |____5a
# | | | |____411577c16ae502c6d7a24d0a97328634b7a86e
# | | |____cd
# | | | |____1e40c7763a3835c6820698ae9f8dac4ee4285b
# | | |____c7
# | | | |____8f4ab9e413c565ae454e27b8d776254ee327b7
# | | |____info
# | | |____36
# | | | |____b99649ae797a08ebb4d91887aa1e754f36c7ac
# | | |____f7
# | | | |____e1b8f9b84b84c6a9f24556f39ea1d0b1b47ab8
# | | |____6b
# | | | |____055fd64400cc05e4ed3814d4ab8b2c399ce4c4
# | | |____1e
# | | | |____b3de567d1100ae145d01fdef12ae914788ebea
# | | |____e6
# | | | |____9de29bb2d1d6434b8b29ae775ad8c2e48c5391
# | | |____b1
# | | | |____ba6cbf55738e1a3bfcc2a5387bbb15b2de1035
# | | | |____f1f2430d4cca10296634af3e36e028b1f415a1
# | | |____f9
# | | | |____3e3a1a1525fb5b91020da86e44810c87a2d7bc
# | | |____d7
# | | | |____f928f764beb69f1409a3c76cf9f7bf112daa0f
# | | |____17
# | | | |____4d740c1d28d3c46b4a2f0d2147d484c74bff91
# | | |____9f
# | | | |____c95617f6b8e9a8941601dbdf3a0e658ddca57e
# | | |____a6
# | | | |____39db314498fb027f84b576abc6329d967d10e7
# | | |____e7
# | | | |____5435c1b1aa229e6b2e849b5da20696779bf734
# | | |____97
# | | | |____5a314adec1fdd93de1da4d3ea22fba215ecfc8
# | | |____d4
# | | | |____f0afdd1c0c1595eb0eb9ea6045a2a5d20acaf0
# | | |____90
# | | | |____9cf411df84978e346a6998a467c7710590aab6
# | | |____5c
# | | | |____2c0df841097638b2c75b59e4eb0a8619321bd9
# | | |____29
# | | | |____4a0868ba4f933e4d5691c9d75ae7834de5bfb1
# | | |____ca
# | | | |____255b34726d37b67e9df857de73761e834f93d1
# | |____refs
# | | |____heads
# | | | |____main
# | | |____remotes
# | | | |____origin
# | | | | |____main
# | | | | |____HEAD
# | | |____tags
# | |____branches
# | |____HEAD
# |____5_annotation
# | |____README.md
# |____3_mapping
# | |____README.md

# 1- unzip the fastq file :X204SC22100973-Z01-F001.tar under the directory 1_data/fastq_set1

mkdir -p 1_data/fastq_set1
tar -xf 1_data/X204SC22100973-Z01-F001.tar -C 1_data/fastq_set1

# 2- unzip the fastq file :X204SC22100973-Z01-F002.tar under the directory 1_data/fastq_set2

mkdir -p 1_data/fastq_set2
tar -xf 1_data/X204SC22100973-Z01-F002.tar -C 1_data/fastq_set2

# 3- Unzip the reference genome Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz under the directory 0_index/referenceIPO323

mkdir -p 0_index/referenceIPO323
gunzip -c 1_data/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa.gz > 0_index/referenceIPO323/Zymoseptoria_tritici.MG2.dna.toplevel.fa

# 4- Unzip the reference genome Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz under the directory 0_index/reference3D7

mkdir -p 0_index/reference3D7
gunzip -c 1_data/reference3D7/Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa.gz > 0_index/reference3D7/Zymoseptoria_tritici_st99ch_3d7_gca_900091695.Zt_ST99CH_3D7.dna.toplevel.fa
