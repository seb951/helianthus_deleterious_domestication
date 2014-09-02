#!/usr/bin/Rscript --verbose

#args = commandArgs(TRUE)
#gene = as.numeric(args[1])


###
###load snp file0
###
setwd("~/Documents/domestication_load")#set up working directory 
library("Kmisc")

results_deleterious = as.matrix(read.table("results/new_results_deleterious", header = T))  #file already exists
results_deleterious[results_deleterious[,11] == "dodgerblue3",11] = "dodgerblue4"


cutoff = c(0,-2.5,-4,-6)
#snp_table = read.table("results/snp_table_noX_new13",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...
snp_table = read.table("results/snp_table_all3_new13",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...

colnames(snp_table) = gsub("X.home.transfer.Documents.new.bam_files_jan2012.","",colnames(snp_table)) #shorten names of the columns.
type = colnames(snp_table)[4:ncol(snp_table)]

individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!!
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])
individuals = individuals[-c(20),]  #remove MEN

provean_results_all = read.delim("provean/new_provean_results_all", sep = " ", stringsAsFactors = F)#provean results
 
unique_orf =  as.matrix(read.table("reference/unique_orf",header = T, stringsAsFactors = F))

ns_s = cbind(rep(0,nrow(snp_table)),0,0) #add a column to see what kind of SNP it is and the SNP EFFECT FROM PROVEAN...
names = paste(c(colnames(snp_table),"ns_site","score","freq_wild","freq_land","freq_eli","freq_weed","freq_all"),collapse = " ") #column names for the result file.

snp_table_effect = read.table("results/new_snp_table_effect", header = T, stringsAsFactors = F);snp_table_effect = snp_table_effect[,-c(23)] #remove ME


snp_table_greg = snp_table; colnames(snp_table_greg)[4:43] = paste(colnames(snp_table_greg)[4:43],individuals[,5], sep = "_")
snp_table_effect_greg = snp_table_effect; colnames(snp_table_effect_greg)[4:43] = paste(colnames(snp_table_effect_greg)[4:43],individuals[,5], sep = "_")

write.table(snp_table_greg,"results/snp_table_greg",col.names = T,row.names = F,sep = "\t", quote = F)#snp table for greg.
write.table(snp_table_effect_greg,"results/snp_table_effect_greg",col.names = T,row.names = F,sep = "\t", quote = F)#snp table for greg.


#look only at ha/ rha individuals
snp_table_comp = snp_table[1:10000,c(1:5,8,10:11,16:18)]
snp_table_effect_comp = snp_table_effect[1:10000,]

c_index = rep(0,nrow(snp_table_comp))
#complementation index
#1 == AA AA	#condition 1
#2 == AC AC	#condition 2
#3 == AC AA	#condition 3
#4 == CC AA	#condition 4
ncomp = matrix(0,nrow = 8, ncol = 8,dimnames  = list(colnames(snp_table_comp[4:11]),colnames(snp_table_comp[4:11]))) #number of comparisons performed
ncomp2 = matrix(0,nrow = 8, ncol = 8,dimnames  = list(colnames(snp_table_comp[4:11]),colnames(snp_table_comp[4:11]))) #number of comparisons performed
ncomp3 = matrix(0,nrow = 8, ncol = 8,dimnames  = list(colnames(snp_table_comp[4:11]),colnames(snp_table_comp[4:11]))) #number of comparisons performed

# are fixed mutations between two lines more deleterious when btwn RHA an HA, than 2 HAs? NO
# are deleteriou mutations more fixed when looking btwen RHA and HA, than 2 HAs? NO
for(line1 in 1:8)
	{
		for(line2 in 1:8)
			{
				for(i in 1: nrow(snp_table_comp))
					{
						read1 = snp_table_comp[i,line1+3]
						read2 = snp_table_comp[i,line2+3]
						if((read1 == read2) & (strsplit(read1,"")[[1]][1] == strsplit(read1,"")[[1]][2])) c_index[i]= 1 #condition 1
						if((read1 == read2) & (strsplit(read1,"")[[1]][1] != strsplit(read1,"")[[1]][2])) c_index[i] = 2	#condition 2
						if((read1 != read2) & (strsplit(read2,"")[[1]][1] != strsplit(read2,"")[[1]][2] | strsplit(read1,"")[[1]][1] != strsplit(read1,"")[[1]][2])) c_index[i] = 3	#condition 3
						if((read1 != read2) & strsplit(read2,"")[[1]][1] == strsplit(read2,"")[[1]][2] & strsplit(read1,"")[[1]][1] == strsplit(read1,"")[[1]][2]) c_index[i] = 4	#condition 3
					}
			ncomp[line1,line2] = mean(snp_table_effect_comp[c_index == 4,45]) #mean deleriousness of mutations which are fixed between 2 lines 
			ncomp2[line1,line2] = rle(sort(c_index[snp_table_effect_comp[,45] > -2]))[[1]][4] / length(c_index) #percentage of mutations non-deleterious that are fixed
			ncomp3[line1,line2] = rle(sort(c_index[snp_table_effect_comp[,45] < -2]))[[1]][4] / length(c_index) #percentage of non-deleterious mutations that are fixed.
			}
			print(paste(line1,"of",line2,"done, Time is:",Sys.time()))	
	}
	
	
###are deleterious mutations more heterozygous? No. but not even sure it makes sense as there is definitevely a correlation btwn allele frequency and level of deleteriousness.
ho = rep(0,nrow(snp_table)) 
#for (i in 1:nrow(snp_table))
for (i in 1:100000)

	{
		temp = unlist(apply(snp_table[i,4:43],2,strsplit,""))
		temp2 = (temp[seq(1,length(temp),by = 2)] == temp[seq(2,length(temp),by = 2)])
		ho[i] =	length(temp2[temp2 == F]) / length(snp_table[i,4:43][snp_table[i,4:43] != "XX"])
		if(i %% 10000 == 0) print(paste(i,"of",nrow(snp_table)))
	}

plot(ho[1:100000],snp_table_effect[1:100000,45])

