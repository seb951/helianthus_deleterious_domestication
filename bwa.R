#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
from = as.numeric(args[1]) # Specify which sequences in "list_ind" file you want to align, directlty from the shell. Alternatively, you can do this from the "alignments" function itself.
to = as.numeric(args[2])

###plan to attack the domestication load question
#Align all the annuus sequences (wild, elite, improved)
#Call SNP
#Create snp file for each wild, elite and improved groups
#Check if they are syn, non-syn, non coding
#run provean to see which are deleterious
#confirm (or test at least) that there are more deleterious alleles in the domesticates.

#setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/domestication_load/")
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/")

#
###
###Step 1 set up working directory and individuals of interest##
###
alignments = function(list_ind= "reference/annuus.txt", from = 1, to = 1) {

###
individuals = read.delim(list_ind, header = T, stringsAsFactors = F)
individuals = cbind(individuals,0)
individuals[,1] = gsub("/home/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/texanus/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/weedy/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/wild/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("~/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/",individuals[,1])


for(i in 1:nrow(individuals)) {individuals[i,5] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}

###
### Step 2. Index reference
###
if(file.exists("reference/HA412_trinity_noAltSplice_400bpmin.fa.pac") == F) system("bwa index reference/HA412_trinity_noAltSplice_400bpmin.fa")
if(file.exists("reference/HA412_trinity_noAltSplice_400bpmin.fa.fai") == F) system("samtools faidx reference/HA412_trinity_noAltSplice_400bpmin.fa")
###create sequence dictionnary with picard tools
if(file.exists("reference/HA412_trinity_noAltSplice_400bpmin.dict") == F) system("/Library/Internet\\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -jar /Users/costanza/applications_seb/picard/picard.jar CreateSequenceDictionary  R=reference/HA412_trinity_noAltSplice_400bpmin.fa O= reference/HA412_trinity_noAltSplice_400bpmin.dict")

###
###Step 3. BWA ALIGNMENTS create fai files.
###

#for(i in 1: nrow(individuals)
for(i in c(from:to))
{
bwa_aln1 = bwa_aln2 = bwa_sampe =  "ls"

if((individuals[i,3] == "ill") & (individuals[i,4] == "sanger")) # new illumina quality format
{bwa_aln1 = paste("bwa aln -t 3 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 3 -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

	if((individuals[i,3] == "ill") & (individuals[i,4] == "Ill1.3")) # old illumina quality format
{bwa_aln1 = paste("bwa aln -t 3 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_1.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai", sep = "");
bwa_aln2 = paste("bwa aln -t 3 -I -q 20 reference/HA412_trinity_noAltSplice_400bpmin.fa ",sub(".fq","_2.fq",individuals[i,1], fixed = T),"  >alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai", sep = "")}

bwa_sampe = paste("bwa sampe -r '@RG\tID:bwa\tSM:",individuals[i,5],"' reference/HA412_trinity_noAltSplice_400bpmin.fa alignments/",sub(".fq","_1.fq",individuals[i,5], fixed = T), ".sai"," alignments/",sub(".fq","_2.fq",individuals[i,5], fixed = T), ".sai ",sub(".fq","_1.fq",individuals[i,1], fixed = T)," ", sub(".fq","_2.fq",individuals[i,1], fixed = T)," >alignments/",individuals[i,5], ".sam", sep = "")

sam_view = paste("/usr/local/bin/samtools view -bSh -o alignments/",individuals[i,5], ".bam", " alignments/",individuals[i,5], ".sam", sep = "")
sam_sort = paste("/usr/local/bin/samtools sort alignments/",individuals[i,5], ".bam ","alignments/",individuals[i,5], ".sorted",sep = "")
sam_index = paste("/usr/local/bin/samtools index alignments/",individuals[i,5],".sorted.bam",sep = "")
sam_rm = paste("rm alignments/",individuals[i,5], ".bam ", "alignments/",individuals[i,5], ".sam alignments/",sub(".fq","_1.fq.sai",individuals[i,5])," alignments/",sub(".fq","_2.fq.sai",individuals[i,5]),sep = "")
#####
gatk_RTC = paste("/Library/Internet\\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -Xmx4g -jar /Users/costanza/applications_seb/GenomeAnalysisTK.jar -T RealignerTargetCreator -I alignments/",individuals[i,5],".sorted.bam -R reference/HA412_trinity_noAltSplice_400bpmin.fa -o alignments/",individuals[i,5],".intervals",sep = "")
gatk_realign = paste("/Library/Internet\\ Plug-Ins/JavaAppletPlugin.plugin/Contents/Home/bin/java -Xmx4g -jar /Users/costanza/applications_seb/GenomeAnalysisTK.jar -T IndelRealigner -I alignments/",individuals[i,5],".sorted.bam -R reference/HA412_trinity_noAltSplice_400bpmin.fa -targetIntervals alignments/",individuals[i,5],".intervals"," -o alignments/",individuals[i,5],".realign.sorted.bam 2>alignments/",individuals[i,5],".log",sep = "")


#java -Xmx4g -jar GenomeAnalysisTK.jar -T IndelRealigner -R ref.fasta -I input.bam -targetIntervals intervalListFromRTC.intervals -o realignedBam.bam
#run commands

aaa = Sys.time()
system(bwa_aln1)
system(bwa_aln2) 
system(bwa_sampe)
system(sam_view)
system(sam_sort)
system(sam_index)
system(gatk_RTC)
system(gatk_realign)
system(sam_rm)
Sys.time() - aaa
}

}

alignments(from =  from, to =  to) ### align 
