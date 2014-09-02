#!/usr/bin/Rscript --verbose
args = commandArgs(TRUE)
genes = as.numeric(args[1]) # how many genes you want to check
cpu = as.numeric(args[2]) #how many cpus you have


###
###call SNPs
###
###########################
### defining the SNP calling function ###
###########################
snp = function(list_ind= "reference/annuus.txt",user = "seb",genes = 1000,cpu = 16,wd = "~/Documents/domestication_load") {

setwd(wd) ### setup the working directory in the right format with the proper subdirectories.

###create "split" object for mpileup if it doesn't exist###
 split = read.table("reference/split", stringsAsFactors = F)

###load individuals' names
individuals = read.delim(list_ind, header = T, stringsAsFactors = F);individuals = cbind(individuals,0) # reference file
for(p in 1:nrow(individuals))
{individuals[p,5] = strsplit(individuals[p,1], split = "/")[[1]][length(strsplit(individuals[p,1], split = "/")[[1]])];
if(length(grep("fq",individuals[p,1])) == 1) {individuals[p,5] = paste(individuals[p,5],".sorted.bam",sep = ""); individuals[p,1] = paste(individuals[p,1],".sorted.bam",sep = "")}
if(length(grep("fq",individuals[p,1])) != 1) {individuals[p,5] = paste(individuals[p,5],".fq.sorted.bam",sep = ""); individuals[p,1] = paste(individuals[p,1],".fq.sorted.bam",sep = "")}
} # proper format to further process with idxstats

all_final1 = paste("alignments/", individuals[1:21,5], collapse = " ", sep = "") # remove the fake individuals (the last 6)
all_final2 = paste("/home/transfer/Documents/new/bam_files_jan2012/", individuals[22:41,5], collapse = " ", sep = "") # remove the fake individuals (the last 6)
all_final = paste(all_final1,all_final2, sep = " ")
#all_final = paste(individuals[,1], collapse = " ", sep = "")

###Pre-cleaning step###
	system("wc -l mpileup/0_51k_variants.clean.vcf >mpileup/wordcount1")
	wordcount = read.table("mpileup/wordcount1")
	system("awk 'NR == 1' mpileup/0_51k_variants.raw.vcf  >mpileup/header")
	#system("awk 'NR == 1' mpileup/variants2.raw.vcf  >mpileup/header")
	
	header = read.delim("mpileup/header", header = F, stringsAsFactors = F)
	header = gsub("/home/seb/Documents/repeatability/","",header, fixed = T)
	header = gsub("/home/seb/Documents/repeatability/alignments/","",header, fixed = T)
	header = gsub("/media/seb_1TB/hybrid_species_MER_GRNnote_transposon/alignments/","",header, fixed = T)
	header = gsub("/home/transfer/Documents/new/bam_files_jan2012/","",header, fixed = T)
	header = gsub("alignments/","",header[10:length(header)], fixed = T)
	header = gsub(".fq.sorted.bam","",header, fixed = T)
	header = gsub(".sorted.bam","",header, fixed = T)
	header = gsub("_R1.fq.bz2","",header, fixed = T)
	header = gsub("/SciBorg/array0/renaut/speciation_islands_individuals/","", header, fixed = T)
	write.table(t(as.matrix(c("reference", "position","reference_annuus_allele",header))),"results/snp_table_all_new13",row.names = F, col.names = F, quote = F, sep = " ")

	temp = NULL
	mis = NULL
	temp_geno = c("RR","LR","LL")
	con = file("mpileup/0_51k_variants.clean.vcf")
	open(con)
	for(i in 1:as.numeric(wordcount[1]))
	#for(i in 1:100000)
		{
		x = readLines(con,1) #you will now read the mega big 0-51k_variants.clean.vcf file line by line in R. 
		xx = strsplit(x, split = "\t")[[1]]
		xxx = xx[(length(xx)-(length(header) - 1)):length(xx)]
		
		for(q in 1:length(xxx))	# this is to replace low quality calls (maximum Phred-scaled genotype likelihoods below 20). ie. essentially, you need at least 2 reads to call a SNP!
		{
		temp = strsplit(xxx[q], split = ":|,")[[1]]
		if(max(as.numeric(temp[2:length(temp)])) < 13) xxx[q] = "XX" else {xxx[q] = temp[1]; xxx[q] = gsub("/","",xxx[q],fixed = T); xxx[q] = gsub("0","R",xxx[q],fixed = T); xxx[q] = gsub("1","L",xxx[q],fixed = T)} #new way of doing things... phred score of genotype of at least 7 (90% sure this is the correct genotype,usually at least 3 reads). Essentially genotypes are coded as such: "GT:PL:GQ        0/1:0,0,0:3			0 refers to the ref allele, 1: the alternate
		#PL: you want the lowest number (probabilities of each genotypes), but you'd also want a large difference between the lowest and highest number, otherwise it means both are equally possible. #GQ the highest numbers.

		}
		xxx = gsub("R",xx[4],xxx)
		xxx = gsub("L",xx[5],xxx)
		xxx = c(xx[c(1,2,4)],xxx)
		if(nchar(xx[5]) > 1) xxx[4:length(xxx)] = rep("XX",length(xxx)-3)
		
		if((length(xxx[xxx == "XX"]) / length(xxx[4:length(xxx)])) < 0.5)	cat(t(as.matrix(xxx)),file = "results/snp_table_all_new13",append = T, fill = F, "\n") else mis = c(mis, i) #only cat loci with less than 50% missing data. 
	#	if(regexpr(xxx[3],paste(xxx[4:21],collapse = "")) < 0) temp = rbind(temp,xxx)
		if(i %% 100000 == 0) print(paste(i,"of", wordcount[1], Sys.time()))
		}
close(con)

system("cat results/snp_table_all_new13 | sed 's/[ \t]*$//' >results/snp_table_all2_new13")
#system("rm mpileup/0_51* mpileup/header mpileup/wordcount1 log prog ") ### rm undesirable outputs
}

############
### running ###
############

snp(list_ind= "reference/annuus.txt",user = "seb",genes = genes,cpu = cpu,wd = "~/Documents/domestication_load")
