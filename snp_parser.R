#!/usr/bin/Rscript --verbose


setwd("~/Documents/domestication_load")#set up working directory 
###########################
###curating snp file ######
###########################
system(paste("wc -l results/snp_table_all2_new13 >wordcount1", sep = ""))
wordcount = read.table("wordcount1")
mis = 0.50 #missing data threshold
con <- file(paste("results/snp_table_all2_new13", sep = ""))

open(con) ; file = readLines(con,1); n_ind = length(strsplit(file," ")[[1]]) -3;close(con); write.table(paste(file),"results/snp_table_all3_new13", row.names = F, col.names = F, quote = F)
 
con <- file(paste("results/snp_table_all2_new13", sep = "")) 
open(con) ; file = readLines(con,1)

for(i in 1:(as.numeric(wordcount[1])-1) )
#for(i in 1:100000)
		{
		file = readLines(con,1)
		file2 = strsplit(file," ")[[1]]
		missing = length(c(1:n_ind)[grepl("XX",file2[4:(n_ind+3)]) == T]) / n_ind

##############################
### trim based on 2pq (Ht) ###
##############################
#counting the alleles.
	counter_ACGTX = rep(0,5)
	names(counter_ACGTX) = c("A","C","G","T","X")
	ht = 0

		x = paste(file2[4:(n_ind+3)], collapse = "")
		xx = strsplit(x, split = "")
		
		counter_ACGTX[1] = length(c(1:((n_ind)*2))[grepl("A", xx[[1]])])
		counter_ACGTX[2] = length(c(1:((n_ind)*2))[grepl("C", xx[[1]])])
		counter_ACGTX[3] = length(c(1:((n_ind)*2))[grepl("G", xx[[1]])])
		counter_ACGTX[4] = length(c(1:((n_ind)*2))[grepl("T", xx[[1]])])
		counter_ACGTX[5] = length(c(1:((n_ind)*2))[grepl("X", xx[[1]])])
		p = sort(counter_ACGTX[1:4])[4] 
		q = sort(counter_ACGTX[1:4])[3] 
		al = sum(counter_ACGTX[1:4])
		ht = 2 * (p / al) * (q / al)

###################################
### trim based on Ho (paralogs) ###
###################################	
		ho  = 0
		a1 = substring(file2[4:(n_ind+3)],1,1)
		a2 = substring(file2[4:(n_ind+3)],2,2)		
		for(h in 1:length(a1))
		{
		if((a1[h] != a2[h]) & (a1[h] != "X") & (a2[h] != "X") & !is.na(a1[h]) & !is.na(a2[h])) ho = (ho + 1) #count the heterozygotes
		}
	ho = ho / n_ind # observed heterozygosity


######## 

if((missing < mis) & (ht > 0.095) & (ho < 0.6)) cat(file,file = paste("results/snp_table_all3_new13"),append = T, fill = 1)

if(i %% 2000 ==0) print(paste(i, Sys.time()))
		}
close(con)


#####################
### FIND BEST ORF ###
#####################
if(file.exists("reference/reference_matrix") == F) {		
### RUN GETORF PACKAGE FROM THE EMBOSS PIPELINE ###
reference_transcriptome = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F,sep = "\t"))
reference = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F,sep = "\t"))
reference_transcriptome_unique_ID = reference_transcriptome
reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1] = paste(">",c(1:length(reference_transcriptome_unique_ID[(regexpr(">", reference_transcriptome_unique_ID[,1],fixed = T) > 0),1])), sep = "") 
write.table(reference_transcriptome_unique_ID,"reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa", row.names = F, col.names = F, quote = F)

### run getorf ###
system("/usr/local/bin/getorf -sequence reference/HA412_trinity_noAltSplice_400bpmin_unique_ID.fa -minsize 300 -find 2 -outseq reference/HA412_trinity_noAltSplice_orf.fasta")

### PARSE THE OUTPUT OF GETORF###
out = as.matrix(read.delim("reference/HA412_trinity_noAltSplice_orf.fasta", header = F, sep = " "))
out = out[,1:5]
out[,2] = gsub("^.","", out[,2])
out[,4] = gsub(".$","", out[,4])
out = gsub("SENSE)","", out, fixed = T)
out = cbind(gsub("(REVERSE","", out, fixed = T), 0)
out = rbind(out,">")

x = c(1:nrow(out))[(regexpr(">",out[,1],fixed = T) > 0)]

for(i in 1: (length(x)-1)) {out[x[i],5]  = paste(out[(x[i]+1): (x[(i+1)]-1),1], collapse = "")}

out_seq = out[(regexpr(">",out[,1],fixed = T) > 0),]
out_seq[,3] =  abs(as.numeric(out_seq[,2]) - as.numeric(out_seq[,4]))
out_seq[,6] = apply(out_seq,2,substring,(regexpr(">", out_seq[,1],fixed = T)+1),(regexpr("_", out_seq[,1],fixed = T)-1))[,1]
unique_gene = unique(out_seq[,6])
unique_orf = NULL

for(i in 1:length(unique_gene))
	{
	x = out_seq[,6] %in% unique_gene[i]
	temp = out_seq[x == T, ]
	if(length(temp) == 6) longest_temp = temp[as.numeric(temp[3]) == max(as.numeric(temp[3]))] else longest_temp = temp[as.numeric(temp[,3]) == max(as.numeric(temp[,3])),]
	if(length(longest_temp) == 6) unique_orf = rbind(unique_orf, longest_temp) else unique_orf = rbind(unique_orf,longest_temp[1,])
	}

	unique_orf = unique_orf[1:(nrow(unique_orf)-1),] #These are the longest ORF.

### get the proper names back
	names = reference_transcriptome[(regexpr(">", reference_transcriptome[,1],fixed = T) > 0),1] 
	for(i in 1:nrow(unique_orf)) {unique_orf[i,1] = names[as.numeric(unique_orf[i,6])]}
	unique_orf[,6] = gsub(">","",unique_orf[,1], fixed = T)

###	
#create a objet reference_matrix which contains the annotated consensus sequences. 
###
reference_vector = c(1:nrow(reference))[(regexpr(">", reference[,1],fixed = T) > 0)]
reference_vector = cbind(reference[(regexpr(">", reference[,1],fixed = T) > 0),1], reference_vector,0)
reference_vector = rbind(reference_vector,c("null",nrow(reference),0))

for(i in 1:(nrow(reference_vector)-1))
{reference_vector[i,3] = paste(reference[(as.numeric(reference_vector[i,2])+1): (as.numeric(reference_vector[(i+1),2])-1),1],collapse = "")}

reference_matrix = cbind(gsub(">","",reference_vector[1:nrow(reference_vector),1]), reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3],reference_vector[1:nrow(reference_vector),3])
colnames(reference_matrix) = c("name","consensus","ONE_haplo1","ONE_haplo2","TWO_haplo1","TWO_haplo2")

consensus_final_per_gene_template = cbind(reference_matrix[,1],0,0,0,nchar(reference_matrix[,2]),0)
colnames(consensus_final_per_gene_template) = c("gene","s","ns","nc","length_total","length_ORF")

write.table(consensus_final_per_gene_template,"reference/consensus_final_per_gene_template.txt", row.names =F, col.names = T)
write.table(reference_matrix,"reference/reference_matrix", row.names =F, col.names = T)
write.table(unique_orf,"reference/unique_orf", row.names =F, col.names = T)
} else {
reference_matrix =  as.matrix(read.table("reference/reference_matrix",header = T, stringsAsFactors = F));
unique_orf =  as.matrix(read.table("reference/unique_orf",header = T, stringsAsFactors = F))
}






