#!/usr/bin/Rscript --verbose

args = commandArgs(TRUE)
gene = args[1]
#gene = 100, or specify "all" if you want all genes.


	#April 30th : NOW PROVEAN RUNS WITH ONLY THE VIRIDIPLANTAE DATASET,  THE NEW NR DATABASE AND VERSION 1.1.4 (WHICH MEANS THAT YOU CAN SPECIFY THE NUM_THREADS ARGUMENT IN PSIBLAST)
	#April 30th: provean is only ran once per gene, then the results are divided per indidual... 
###
###load snp file
###
setwd("~/Documents/domestication_load")#set up working directory 
library("Kmisc")


###you must install PROVEAN.
###FIRST THING, to the blast command, you need to ADD THE -gilist /location/of/gi/restriction to the Sequence.cpp file in the src directory of provean. (needs to be added in 2 places)
###./configure PSIBLAST=/usr/bin/psiblast BLASTDBCMD=/usr/bin/blastdbcmd CDHIT=/usr/bin/cdhit
###make
###sudo make install
###modify the /usr/local/bin/provean.sh script to add the location of the blast database.


if(file.exists("results/snp_table_noX_new13") == F) {
snp_table = read.table("results/snp_table_all3_new13",stringsAsFactors = F, header = T)
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
	write.table(snp_table,"results/snp_table_noX_new13", row.names = T, col.names =T,quote = F)
																	} else snp_table = read.table("results/snp_table_noX_new13",stringsAsFactors = F, header = T)		
					
reference_matrix =  as.matrix(read.table("reference/reference_matrix",header = T, stringsAsFactors = F))
individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!! 
unique_orf =  as.matrix(read.table("reference/unique_orf",header = T, stringsAsFactors = F))

### prepare the files for each category ###
type = colnames(snp_table)[4:44]
#type = c("wild","weed","land","wild_tex","eli")
for(ec in 1: length(type)) 
	{
	echo_seq = paste("echo all_sequences_",type[ec]," >provean/new_all_sequences_",type[ec],sep = "") #per individuals
	echo_snp = paste("echo aa_change score gene >provean/new_provean_results_",type[ec],sep = "") #per individuals
	echo_snp_all = paste("echo aa_change score gene >provean/new_provean_results_all",sep = "") #for all genes and all individuals
	
	system(echo_seq)
	system(echo_snp)
	system(echo_snp_all)
	}


### annotate reference matrix with the consensus object.
###
#######
##################################################################################################################
##################################################################################################################
##################################################################################################################
if(gene == "all") gene = nrow(reference_matrix) else gene = as.numeric(gene)

#for(g in 1:24)
for(g in 1:gene) #run it for each gene
	{
		reference_gene = reference_matrix[g,c(1,2,rep(2,2*(ncol(snp_table) - 3)))] #First colum is the name of the gene, second is the reference each individual gets two haplotypes.
			{
				x = snp_table[,1] %in% reference_matrix[g,1]
				y =  snp_table[x == T, ]
			
				if(nrow(y) == 1) #annotate SNP in all haplotypes. two alleles per snps. Also annotate the consensus...
					{
						for(ind in 1:(ncol(snp_table) - 3)) 
							{
								if(ind == 1) substring(reference_gene[2], as.numeric(y[2]), as.numeric(y[2])) = substring(y[3],1,1) #annotate the consensus
								substring(reference_gene[(ind*2)+2], as.numeric(y[2]), as.numeric(y[2])) = substring(y[ind+3],1,1); 
								substring(reference_gene[(ind*2)+3], as.numeric(y[2]), as.numeric(y[2])) = substring(y[ind+3],2,2)
							}
					}
	
				if(nrow(y) > 1) 
				for(j in 1:nrow(y)) #annotate all SNPs...
					{
						for(ind in 1:(ncol(snp_table) - 3))  #...in all haplotypes. Also annotate the consensus...
							{
								if(ind == 1) substring(reference_gene[2], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,3],1,1) #annotate the consensus
								substring(reference_gene[(ind*2)+1], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,ind+3],1,1); 
								substring(reference_gene[(ind*2)+2], as.numeric(y[j,2]), as.numeric(y[j,2])) = substring(y[j,ind+3],2,2)
							}
					}
			}

	### now you have a single gene done...
if(nrow(y) > 0) {
###
#Interrogate each SNP, find what codon it belongs to and see whether it is syn or nonsyn. add this information in the consensus table. 
###
reference_orf = reference_gene; reference_aa = reference_gene ; reference_aa[2:length(reference_aa)] = ""
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
		for(ind in 1:(ncol(snp_table) - 3)) # get orf for each haplotype
			{
				if(as.numeric(y_orf[2]) < as.numeric(y_orf[4]) & (length(y_orf) > 0))  # the easy case
				{
					if(ind == 1) orf_cons = substring(reference_gene[2],as.numeric(y_orf[2]),as.numeric(y_orf[4])) #consensus sequence
					orf1 = substring(reference_gene[(ind*2)+1],as.numeric(y_orf[2]),as.numeric(y_orf[4]))
					orf2 = substring(reference_gene[(ind*2)+2],as.numeric(y_orf[2]),as.numeric(y_orf[4])) 
					if(ind == 1) reference_orf[2] = orf_cons #bring consensus to orf object
					reference_orf[(ind*2)+1] = orf1; reference_orf[(ind*2)+2] = orf2 #bring sequences to orf object
				}

		#						
				if(as.numeric(y_orf[2]) > as.numeric(y_orf[4]) & (length(y_orf) > 0))  #reverse complement case
				{
				 	if(ind == 1) {orf_cons = substring(reference_gene[2],as.numeric(y_orf[4]),as.numeric(y_orf[2])) ; orf_cons_r = str_rev(orf_cons)}
				 	orf1 = substring(reference_gene[(ind*2)+1],as.numeric(y_orf[4]),as.numeric(y_orf[2])) ; orf1_r = str_rev(orf1)
					orf2 = substring(reference_gene[(ind*2)+2],as.numeric(y_orf[4]),as.numeric(y_orf[2])) ; orf2_r = str_rev(orf2)
					for(s in 1:length(symbols)) {orf1_r  = gsub(symbols[s], replacements[s], orf1_r);orf2_r  = gsub(symbols[s], replacements[s], orf2_r); if(ind == 1) {orf_cons_r  = gsub(symbols[s], replacements[s], orf_cons_r)}}#complement sequence.
					orf1_r = gsub("(\\w)", "\\U\\1", orf1_r, perl=TRUE); orf2_r = gsub("(\\w)", "\\U\\1",orf2_r, perl=TRUE); if(ind == 1) {orf_cons_r = gsub("(\\w)", "\\U\\1",orf_cons_r, perl=TRUE)}#gsub for small caps to capital letters. 
					reference_orf[(ind*2)+1] = orf1_r; reference_orf[(ind*2)+2] = orf2_r; if(ind == 1) {reference_orf[2] = orf_cons_r} #brings sequences to orf object
				}
			}
																														
			###get AA sequence for everyone.
			codon_positions = seq(1, nchar(y_orf[5]), 3)
			
			for(ind in 1:(ncol(snp_table) - 3))    # get orf for each haplotype
				{
					for(a in 1:length(codon_positions)) 		#populate AAs
						{
							if(ind == 1) temp_codon_cons = substring(reference_orf[2],codon_positions[a],(codon_positions[a]+2))
							temp_codon1 = substring(reference_orf[(ind)*2 +1],codon_positions[a],(codon_positions[a]+2))
							temp_codon2 = substring(reference_orf[(ind)*2 +2],codon_positions[a],(codon_positions[a]+2))
							if(ind == 1) temp_aa_cons = aa[temp_codon_cons == aa[,2],3];if(length(aa[temp_codon_cons == aa[,2],3]) ==0)  temp_aa_cons = "B"  #in case there is a ambiguity, it will be replaced by a "B", which doesn't exist...
							temp_aa1 = aa[temp_codon1 == aa[,2],3];if(length(temp_aa1) ==0) temp_aa1 = "B"  #in case there is a ambiguity, it will be replaced by a "B", which doesn't exist...
							temp_aa2 = aa[temp_codon2 == aa[,2],3];if(length(temp_aa2) ==0) temp_aa2 = "B"  #in case there is a ambiguity, it will be replaced by a "B", which doesn't exist...
							if(length(temp_aa1) > 0) {
							if(ind == 1) reference_aa[2] = paste(reference_aa[2],temp_aa_cons,sep = "");
							reference_aa[(ind)*2+1] = paste(reference_aa[(ind)*2+1],temp_aa1,sep = "");
							reference_aa[(ind)*2+2] = paste(reference_aa[(ind)*2+2],temp_aa2,sep = "")
							} else {
							if(ind == 1) reference_aa[2] = paste(reference_aa[2],0,sep = "");
							reference_aa[(ind)*2+1] = paste(reference_aa[(ind)*2+1],0,sep = ""); #append all AAs or a 0 is the AA is not there...
							reference_aa[(ind)*2+2] = paste(reference_aa[(ind)*2+2],0,sep = "")}
						}
				}
						
###prepare a file for provean###						
###prepare a file for provean###						
###prepare a file for provean###	



type_ind = rep(0,ind*2)
#type_ind[seq(1,ind*2,by =2)] = individuals[,5]; type_ind[seq(2,ind*2,by =2)] = individuals[,5]
type_ind[seq(1,ind*2,by =2)] = type; type_ind[seq(2,ind*2,by =2)] = type
####
seq_all.var = NULL #file which will contain all mutations concatenated.
seq_ind.var = as.list(rep(0,length(type))) #list which will contain all mutations per individuals.

###for all four type FILE

for(t in 1:length(type))
	{
		seq = reference_aa[3:length(reference_aa)][type_ind == type[t]]
		seq.fasta = c(paste(">",reference_aa[1],"###gene",g,sep = ""),reference_aa[2])
		seq.var.temp = NULL
		
		seq.var.temp1 = c(1:nchar(seq[1]))[(strsplit(seq[1],"")$consensus == strsplit(reference_aa[2],"")$consensus) == F] 
		seq.var.temp2 = c(1:nchar(seq[1]))[(strsplit(seq[2],"")$consensus == strsplit(reference_aa[2],"")$consensus) == F]
			if(length(seq.var.temp1) > 1) for(ii in 1:length(seq.var.temp1))
				{
				seq.var.temp = c(seq.var.temp,paste(strsplit(reference_aa[2],"")$consensus[seq.var.temp1[ii]],seq.var.temp1[ii],strsplit(seq[1],"")$consensus[seq.var.temp1[ii]],sep = "")) #add all variation into a single vector.
				} else {seq.var.temp = c(seq.var.temp,paste(strsplit(reference_aa[2],"")$consensus[seq.var.temp1],seq.var.temp1,strsplit(seq[1],"")$consensus[seq.var.temp1],sep = ""))}
			if(length(seq.var.temp2) > 1) for(ii in 1:length(seq.var.temp2)) #second haplotype
				{
				seq.var.temp = c(seq.var.temp,paste(strsplit(reference_aa[2],"")$consensus[seq.var.temp2[ii]],seq.var.temp2[ii],strsplit(seq[2],"")$consensus[seq.var.temp2[ii]],sep = "")) #add all variation into a single vector.
				} else {seq.var.temp = c(seq.var.temp,paste(strsplit(reference_aa[2],"")$consensus[seq.var.temp2],seq.var.temp2,strsplit(seq[2],"")$consensus[seq.var.temp2],sep = ""))}
					
					
			seq.var.temp = unique(seq.var.temp)
				
			seq.var = unique(seq.var.temp[seq.var.temp != "NA"])			
			
			if(length(seq.var)>0) { #make sure there is variation in each individual
			seq_ind.var[[t]] = sort(seq.var)	#individual variation
			seq_all.var = c(seq_all.var,seq.var) 	#all mutations variation
											}		
		}
	if(length(seq_all.var) > 0) {seq_all.var = unique(sort(seq_all.var));write.table(seq_all.var,paste("provean/seq_all.var"),row.names = F, col.names = F, quote = F)}	#write all variation, if there is any
	write.table(seq.fasta,"provean/seq.fasta",row.names = F, col.names = F, quote = F)  #write sequence a single time for all individuals.

#################
#################run provean once for all individuals, then sort out the results per individuals
#################
		
		if(length(seq_all.var) > 0) { #only run if variation is present...

		if(file.exists(paste("provean/blast/set",g, sep = ""))) {system(paste("provean.sh --num_threads 7 -q provean/seq.fasta -v provean/seq_all.var --supporting_set provean/blast/set",g," >provean/provtest.log",sep = ""))} else {system(paste("provean.sh --num_threads 7 -q provean/seq.fasta -v provean/seq_all.var --save_supporting_set provean/blast/set",g," >provean/provtest.log", sep = ""))} #reuse the blast results if it already exists for that gene... otherwise create it for the other individuals to use it... 
		#April 30th : NOW PROVEAN RUNS WITH ONLY THE VIRIDIPLANTAE DATASET,  THE NEW NR DATABASE and version 1.1.4
			
		variants = read.csv("provean/provtest.log", header = F,sep = "\t", stringsAsFactors = F) #load variants
		
		if(length(variants[regexpr("VARIATION",variants[,1]) >0,1]) > 0 ) { #set up another filter because even though we already have if(length(seq_all.var) > 0) above, provean can crash and give no results
		variants = variants[(grep("VARIATION",variants[,1])+1):nrow(variants),];
		variants = cbind(variants,paste(">",reference_aa[1],"###gene",g,sep = "")); #add the name of the sequence were the variant is from 
		write.table(variants,"provean/variants",row.names = F, col.names = F, quote = F);
		system("cat provean/variants >>provean/new_provean_results_all");  #write sequence a single time for all individuals.system(#cat all the variants in a single file
		
		for(t in 1:length(type))
			{
				temp = variants[,1] %in% seq_ind.var[[t]]
				seq_ind.var[[t]] = variants[temp == T,] 
				write.table(seq_ind.var[[t]],"provean/variants",row.names = F, col.names = F, quote = F); #write it back
				if(nrow(seq_ind.var[[t]]) > 0) { 
				system(paste("cat provean/variants >>provean/new_provean_results_",type[t],sep = "")); # keep track of all MUTATIONS only if variants are present (seq_ind.var[[t]])
				system(paste("cat provean/seq.fasta >>provean/new_all_sequences_",type[t],sep = ""))} # keep track of all sequences only if variants are present (seq_ind.var[[t]])
			} ###result file
																											}
											}
										
	}	#only run if you have an ORF								
	
	} #only run if you have SNPs.
	if(g %% 10 == 0) print(paste("done", "gene",g,"out of",gene,", AA length =",nchar(seq.fasta[2])," , Time is:",Sys.time()))
}#run this per gene...


