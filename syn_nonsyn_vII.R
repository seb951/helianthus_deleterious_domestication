#!/usr/bin/Rscript --verbose

#args = commandArgs(TRUE)
#gene = as.numeric(args[1])


###
###load snp file0
###
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/")#set up working directory 
library("Kmisc")


if(file.exists("results/snp_table_noX") == F) {
snp_table = read.table("results/snp_table_all3",stringsAsFactors = F, header = T)
colnames(snp_table) = gsub("X.home.transfer.Documents.new.bam_files_jan2012.","",colnames(snp_table)) #shorten names of the columns.
snp_table2 = as.matrix(snp_table) # ,nrow = nrow(snp_table), ncol = ncol(snp_table))

for(i in 1:nrow(snp_table))
#for(i in 1:1000)
	{
	majority_temp = rle(sort(strsplit(paste(snp_table[i,c(4:ncol(snp_table))],collapse = ""),"")[[1]])) #replace the consensus read with a majority rule consensus nucleotide.
	majority_rule = cbind(majority_temp[[1]][majority_temp[[2]] != "X"],majority_temp[[2]][majority_temp[[2]] != "X"])

	snp_table2[i,3] = majority_rule[order(as.numeric(majority_rule[,1]),decreasing = T),2][1] #replace the consensus read with a majority rule consensus nucleotide.
	snp_table2[i,] = gsub("X",snp_table2[i,3],snp_table2[i,]) #replace missing data with the reference
	if(i %% 10000 == 0) print(paste(i,Sys.time()))
	}
		snp_table[1:nrow(snp_table),1:ncol(snp_table)] = snp_table2
	write.table(snp_table,"results/snp_table_noX", row.names = T, col.names =T,quote = F)
																	} else snp_table = read.table("results/snp_table_noX",stringsAsFactors = F, header = T)		
					
reference_matrix =  as.matrix(read.table("reference/reference_matrix",header = T, stringsAsFactors = F))
individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!!
individuals = individuals[-20,]
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])

provean_results_all = read.delim("provean/new_provean_results_all", sep = " ", stringsAsFactors = F)#provean results
 
unique_orf =  as.matrix(read.table("reference/unique_orf",header = T, stringsAsFactors = F))

ns_s = cbind(rep(0,nrow(snp_table)),0,0) #add a column to see what kind of SNP it is and the SNP EFFECT FROM PROVEAN...
names = paste(c(colnames(snp_table),"ns_site","score","freq_wild","freq_land","freq_eli","freq_weed","freq_all"),collapse = " ") #column names for the result file.

###
### annotate reference matrix with the consensus object.
###
if(file.exists("results/new_snp_table_effectXXX") == F) {
###system(paste("echo",names,">results/new_snp_table_effect")) #result file with column names...
for(nss in 473673:nrow(snp_table))
#for(nss in 1:25000)
	{
		snp_table_effect = matrix(0,nrow = 1, ncol = 50) ### vector with info on what kind of snp it is... EMPTY IT EVERYTIME JUST TO BE SURE...
		reference_gene = reference_matrix[reference_matrix[,1] == snp_table[nss,1],c(1,2,2,2)] #First colum is the name of the gene, second is the reference, 3rd is for allele 1, 4rth for allele 2.
		y =  snp_table[nss, ]
				
		y_alleles = rle(sort(strsplit(paste(y[,-c(1,2,3)],collapse = ""),"")[[1]]))$values # look at all the alleles present...
		#annotate SNP in all haplotypes. two alleles per snps. Also annotate the consensus...
		substring(reference_gene[2], as.numeric(y[2]), as.numeric(y[2])) = substring(y[3],1,1) #annotate the consensus
		substring(reference_gene[3], as.numeric(y[2]), as.numeric(y[2])) = y_alleles[1]
		substring(reference_gene[4], as.numeric(y[2]), as.numeric(y[2])) = y_alleles[length(y_alleles)]		
						
	

###
#Interrogate each SNP, find what codon it belongs to and see whether it is syn or nonsyn. add this information in the consensus table. 
###
reference_orf = reference_gene; reference_aa = reference_gene ; reference_aa[3:length(reference_aa)] = ""
aa = as.matrix(read.delim("reference/aa_code.txt", header = T))
aa_temp = c(0,0,0,0)
symbols =      c("A","C","G","T","M","R","W","S","Y","K","V","H","D","B")
replacements = c("t","g","c","a","k","y","w","s","r","m","b","d","h","v")
codon = NULL

	#x_seq = reference_gene[,1] %in% snp_table[i,1]
	x_orf = unique_orf[,6] %in% reference_gene[1]
	#y_seq = reference_matrix[x_seq == T, ]
	y_orf = unique_orf[x_orf == T, ]
	if(length(y_orf) > 0) { # make sure there is an ORF		
				if(as.numeric(y_orf[2]) < as.numeric(y_orf[4]) & (length(y_orf) > 0))  # the easy case
				{
					orf1 = substring(reference_gene[3],as.numeric(y_orf[2]),as.numeric(y_orf[4])) 
					orf2 = substring(reference_gene[4],as.numeric(y_orf[2]),as.numeric(y_orf[4])) 
					reference_orf[3] = orf1; reference_orf[4] = orf2 #brings sequences to orf object
				}

		#						
				if(as.numeric(y_orf[2]) > as.numeric(y_orf[4]) & (length(y_orf) > 0))  #reverse complement case
				{
				 	orf1 = substring(reference_gene[3],as.numeric(y_orf[4]),as.numeric(y_orf[2])) ; orf1_r = str_rev(orf1)
					orf2 = substring(reference_gene[4],as.numeric(y_orf[4]),as.numeric(y_orf[2])) ; orf2_r = str_rev(orf2)
					for(s in 1:length(symbols)) {orf1_r  = gsub(symbols[s], replacements[s], orf1_r);orf2_r  = gsub(symbols[s], replacements[s], orf2_r)}#complement sequence.
					orf1_r = gsub("(\\w)", "\\U\\1", orf1_r, perl=TRUE);  	orf2_r = gsub("(\\w)", "\\U\\1",orf2_r, perl=TRUE) #gsub for small caps to capital letters. 
					reference_orf[3] = orf1_r; reference_orf[4] = orf2_r #brings sequences to orf object
				}
			
		#add information of provean and also were the snp is...		
		provean_results_pergene = provean_results_all[regexpr(y[,1],provean_results_all[,3]) >0,] #snp effect per gene...
				if(nrow(provean_results_pergene) > 0) 
				{
					rename = gsub("[1234567890]","",provean_results_pergene[,1]) #only do this if provean found something for that gene.
					rename2 = t(matrix(unlist(strsplit(rename,split = "")), ncol = length(rename), nrow = 2))
					rename3 = t(apply(rename2,1,sort))
					rename4 = apply(rename3,1,paste,collapse = "")
					provean_results_pergene = cbind(provean_results_pergene,rename4,0)			
					for(v in 1:nrow(provean_results_pergene))
						{
							if(as.numeric(y_orf[2]) > as.numeric(y_orf[4]) & (length(y_orf) > 0)) provean_results_pergene[v,5] = as.numeric(y_orf[4]) + as.numeric(y_orf[3]) - (as.numeric(gsub("^.|.$","",provean_results_pergene[v,1]))*3) + 1 #reverse complement
							if(as.numeric(y_orf[2]) < as.numeric(y_orf[4]) & (length(y_orf) > 0)) provean_results_pergene[v,5] = as.numeric(y_orf[2]) + (as.numeric(gsub("^.|.$","",provean_results_pergene[v,1]))*3) - 1 #the easy case
						}
				}
																														
			###get AA sequence for everyone.
			codon_positions = seq(1, nchar(y_orf[5]), 3)
			temp_ns = NULL
						for(a in 1:length(codon_positions)) 		#populate AAs	
				
						{
							temp_codon1 = substring(reference_orf[3],codon_positions[a],(codon_positions[a]+2))
							temp_codon2 = substring(reference_orf[4],codon_positions[a],(codon_positions[a]+2))
						#	if(temp_codon1 != temp_codon) print(c(temp_codon1,temp_codon2))
							temp_aa1 = aa[temp_codon1 == aa[,2],3]
							temp_aa2 = aa[temp_codon2 == aa[,2],3]
							if((temp_aa1 != temp_aa2) & (nrow(provean_results_pergene) > 0)) 
								{ 
									temp_ns = paste(sort(c(temp_aa1,temp_aa2)),  collapse = ""); #find out the snp effect...
									temp_effect = provean_results_pergene[provean_results_pergene[,4] == temp_ns,]
									if(nrow(temp_effect) == 1) ns_s[nss,2:3] = as.character(provean_results_pergene[provean_results_pergene[,4] == temp_ns,1:2]) #there is only a single snp which corresponds.
									if(nrow(temp_effect) > 1)  ns_s[nss,2:3] =  as.character(temp_effect[which.min(abs(y[,2] - temp_effect[,5])),1:2]) #there are several snps which corresponds, but now you must sort them according to position too...
								}
													
							if(length(temp_aa1) > 0)
								{
									reference_aa[3] = paste(reference_aa[3],temp_aa1,sep = "");
									reference_aa[4] = paste(reference_aa[4],temp_aa2,sep = "")
								} else
								{
									reference_aa[3] = paste(reference_aa[3],0,sep = ""); #append all AAs or a 0 is the AA is not there...
									reference_aa[4] = paste(reference_aa[4],0,sep = "")
								}
						}
				if(reference_aa[3] != reference_aa[4])   ns_s[nss,1] = "ns" # non-syn
				if((reference_aa[3] != reference_aa[4]) & ((gregexpr("X",reference_aa[3])[[1]] > 0) |  (gregexpr("X",reference_aa[4])[[1]] > 0)))   ns_s[nss,1] = "STOP"  #alternate stop codon was introduced...
				#
				if(length(gregexpr("X",reference_aa[3])[[1]]) > 1) print(paste(nss,reference_aa[1],reference_aa[3],sep = "###"))
		
				if((reference_aa[3] == reference_aa[4]) & (orf1 != orf2)) ns_s[nss,1] = "s" # syn
				if((reference_aa[3] == reference_aa[4]) & (orf1 == orf2)) ns_s[nss,1] = "nc" # non-coding
										
				}
				
																

snp_table_effect[1,(ncol(snp_table_effect)-6):(ncol(snp_table_effect)-5)] = ns_s[nss,2:3] ###annotate the SNP effect info into a new matrix, per individual.###

		for(j in 4:ncol(snp_table))
			{
			snp_temp = snp_table[nss,j]
			if(substring(snp_temp,1,1) != substring(snp_temp,2,2)) snp_table_effect[1,j] = ns_s[nss,1] #there is an hetero snp...
			if((substring(snp_temp,1,1) == substring(snp_temp,2,2)) &  (substring(snp_temp,1,1) == snp_table[nss,3]) ) snp_table_effect[1,j] = 0 #there is no snp so it is missing data...
			if((substring(snp_temp,1,1) == substring(snp_temp,2,2)) &  (substring(snp_temp,1,1) != snp_table[nss,3]) ) snp_table_effect[1,j] = ns_s[nss,1] #there is a fixed snp...
			}
		
###frequency of the rare allele for everyone and per group...
reads = y[,-c(1:3)]
counts_all = rle(sort(strsplit(paste(reads,collapse = ""), split = "")[[1]]))
counts_wild = rle(sort(strsplit(paste( reads[individuals[,5] == "wild"],collapse = ""), split = "")[[1]]))
counts_land = rle(sort(strsplit(paste(reads[individuals[,5] == "landrace"],collapse = ""), split = "")[[1]]))
counts_eli = rle(sort(strsplit(paste(reads[individuals[,5] == "elite"],collapse = ""), split = "")[[1]]))
counts_weed = rle(sort(strsplit(paste(reads[individuals[,5] == "weed"],collapse = ""), split = "")[[1]]))

temp1 = signif(counts_wild$length[counts_wild$values != snp_table[nss,3]] / sum(counts_wild$lengths),4)	 #frequency of the derived allele		
temp2 = signif(counts_land$length[counts_land$values != snp_table[nss,3]] / sum(counts_land$lengths),4)	
temp3 = signif(counts_eli$length[counts_eli$values != snp_table[nss,3]] / sum(counts_eli$lengths),4)	
temp4 = signif(counts_weed$length[counts_weed$values != snp_table[nss,3]] / sum(counts_weed$lengths),4)
temp5 = signif(counts_all$length[counts_all$values != snp_table[nss,3]] / sum(counts_all$lengths),4)	

if(length(temp1) != 0) snp_table_effect[1,length(snp_table_effect)-4] = temp1
if(length(temp2) != 0) snp_table_effect[1,length(snp_table_effect)-3] = temp2
if(length(temp3) != 0) snp_table_effect[1,length(snp_table_effect)-2] = temp3
if(length(temp4) != 0) snp_table_effect[1,length(snp_table_effect)-1] = temp4
if(length(temp5) != 0) snp_table_effect[1,length(snp_table_effect)] = temp5

snp_table_effect[1,is.na(snp_table_effect[1,])] = 0 # remove potential NAs and change to 0...
snp_table_effect[1,1:3] = as.character(snp_table[nss,1:3])
###output to file... CAT CAT CAT
write.table(snp_table_effect,"new_snp_table_effect_temp",row.names = F, col.names = F, quote = F)

system("cat new_snp_table_effect_temp >>results/new_snp_table_effect")

if(nss %% 1000 == 0) print(paste(nss,"of",nrow(snp_table),"Time is:",Sys.time()))
			}
system("rm new_snp_table_effect_temp")
} else snp_table_effect = read.table("results/new_snp_table_effect", header = T, stringsAsFactors = F)

#;snp_table_effect = snp_table_effect[,-c(23)] #remove MEN


