#!/usr/bin/Rscript --verbose

###Step 0.1 set up working directory ###
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/SITES_macosx/deleterious_load/")

######################
###get the length of the genes###
###################### 
ref = read.delim("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/reference/HA412_trinity_noAltSplice_400bpmin.fa", header = F, stringsAsFactors = F)
x = c(1:nrow(ref))[gregexpr(">",ref[,1]) > 0]
len = matrix(0, nrow = length(x), ncol = 22); len[,1] =  ref[x,1]; x = c(x,nrow(ref))

colnames(len)  = c("name","len", "pi_all_eli", "pi_all_lan","pi_all_wee","pi_all_wil", 
"pi_s_eli", "pi_s_lan","pi_s_wee","pi_s_wil", 
"pi_ns_eli", "pi_ns_lan","pi_ns_wee","pi_ns_wil", 
"pi_STOP_eli", "pi_STOP_lan","pi_STOP_wee","pi_STOP_wil", 
"pi_nsdel_eli", "pi_nsdel_lan","pi_nsdel_wee","pi_nsdel_wil")

for(i in 1:nrow(len)) {len[i,2] = sum(nchar(ref[(x[i]+1):(x[i+1]-1),1]))}
len[,1] = gsub(">","", len[,1],fixed = T)
###

comp4 = read.delim("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/results_old/snp_table_all3_new13", header = T,stringsAsFactors = F, sep =  " ")
comp4 = comp4[,-23]
individuals = read.delim("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/reference/annuus.txt", header = T, stringsAsFactors = F)
individuals = cbind(individuals,0)
individuals[,1] = gsub("/home/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/texanus/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/weedy/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("/Linux/Loren/Seq/EST/annuus/wild/illumina/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/sequences/",individuals[,1])
individuals[,1] = gsub("~/","/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/",individuals[,1])
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])
individuals = individuals[-20,]#remove 20.

for(i in 1:nrow(individuals)) {individuals[i,6] = strsplit(individuals[i,1], split = "/")[[1]][length(strsplit(individuals[i,1], split = "/")[[1]])]}

#############################
#####categorise the SNPs#####
#############################
snp_table_effect = read.delim("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/results_old/new_snp_table_effect", header = T,stringsAsFactors = F, sep =  " ")
snp_table_effect = snp_table_effect[,-23] #remove MEN.
snp_table_effect_v = rep(0,nrow(snp_table_effect))
###single vector with the SNP table effect.
for(i in 1:nrow(snp_table_effect))
	{
		temp = snp_table_effect[i,4:43]
		if(length(temp[temp != 0][1]) != 0) snp_table_effect_v[i] = temp[temp != 0][1]
		if(is.na(snp_table_effect_v[i])) snp_table_effect_v[i] = 0
		if(snp_table_effect_v[i] == "ns") snp_table_effect_v[i] = ifelse(snp_table_effect[i,45] < -2.5,"nsdel","ns")
		if(i %% 10000 == 0) print(paste(i,"of",nrow(snp_table_effect),Sys.time()))
	}

all = individuals[1:40,]

###
aaa = Sys.time()
un_comp4 = unique(comp4[,1]) # unique set of genes, this is just faster than looking through the whole matrix everytime.

for(categ in 1:5) #do the loop for each category of mutation
	{
		if(categ == 1) comp5 = comp4###comp4 only for ALL mutations
		if(categ == 2) comp5 = comp4[snp_table_effect_v == "s",] ###comp4 only for SYN mutations
		if(categ == 3) comp5 = comp4[snp_table_effect_v == "ns",]
		if(categ == 4) comp5 = comp4[snp_table_effect_v == "STOP",]
		if(categ == 5) comp5 = comp4[snp_table_effect_v == "nsdel",]
	
	for(i in 1:nrow(len)) 
#		for(i in 1:100)
	{
	
	
	x1 = comp5[comp5[,1] == len[i,1],4:43];
	
	if(nrow(x1) > 0)  {
	
	if(is.na(x1[1,1]))	{x1 = matrix("AA", nrow = 2, ncol = 43)}; # just make sure the data is in a matrix format, even if nothing is there. 
	if(nrow(x1) == 1)	x1 = rbind(x1,"AA"); # if only one row, add another one to treat as matrix.

	x1 =  gsub("0","XX",as.matrix(x1));	
	x1 =  gsub("X","N",as.matrix(x1));	
	
	x1_one = x1; #x1[,regexpr(sps[o], colnames(x1))> 0]
	one_presites = matrix(0,nrow = (ncol(x1_one)*2), ncol = 2);
	seq = NULL; for(r in 1:40) {seq = c(seq,rep(r,2))}; one_presites[,1] = substring(all[seq,5],1,3); #
	for(n in 1:ncol(x1_one))
	{one_presites[(n*2)-1,2] = paste(substring(x1_one[,n],1,1),collapse = "");one_presites[(n*2),2] = paste(substring(x1_one[,n],2,2),collapse = "")};
	one_presites[,2] =  paste(one_presites[,1],paste(one_presites[,2],paste(rep("A",as.numeric(len[i,2])-nrow(x1_one)),collapse = ""),sep = ""),sep = "       ");
  
  one_presites = one_presites[order(one_presites[,1]),];
	
	#write tex_deb file 
	line1 = "SITES sample input" ;
	line2 = paste(ncol(x1_one)*2, nchar(one_presites[1,2])-10); ###16 is to remove the FUCKING f1s...
	line3 = paste(1,1);
	line4 = paste(1,nchar(one_presites[1,2]));
 	line5	= 4; #$four different groups...
  ###pi for each of the 4 groups...
	line6 = paste(rle(one_presites[,1])$values[1],rle(one_presites[,1])$lengths[1]);
	line7 = paste(rle(one_presites[,1])$values[2],rle(one_presites[,1])$lengths[2]);
	line8 = paste(rle(one_presites[,1])$values[3],rle(one_presites[,1])$lengths[3]);
	line9 = paste(rle(one_presites[,1])$values[4],rle(one_presites[,1])$lengths[4]);
	
	sites = c(line1,line2,line3,line4,line5,line6,line7,line8,line9,one_presites[,2]);

	write.table(sites, "sites_sample", row.names = F, col.names = F, quote = F);

	system("./sites -isites_sample -rsites_results -mmy_data -sa -asp -ogc -cf 1>log");

###read output###
#system("grep 'Fst' sites_results.SIT -A 5 | tail -1 >fst") 
	system("grep 'Theta' sites_results.SIT -A 6 | tail -5 >pi") ;

#fst = as.numeric(read.table("fst")[3])
	pi = read.table("pi")[,8];

###update the table 
	if(categ == 1) len[i,3:6] = pi
	if(categ == 2) len[i,7:10] = pi
	if(categ == 3) len[i,11:14] = pi
	if(categ == 4) len[i,15:18] = pi
	if(categ == 5) len[i,19:22] = pi
	
	}
	
	
	
		if(i %% 10 == 0) print(paste(i,Sys.time()))
	}
	}
Sys.time() - aaa

#write.table(len,"/home/seb/Documents/helianthus_analysis/92samples_new_alignments_sept11/results/fst_pergene", row.names = F, col.names =T)
write.table(len,"SITES_pi", row.names = F, col.names =T)


###
###PLOTTING
###
sites = read.table("SITES_pi",header = T, stringsAsFactors = F)
dim(sites)
head(sites)

results_pi = matrix(0,nrow = 4, ncol = 4)

colnames(results_pi) = c("pi_all","pi_ns/pi_s","pi_stop/pi_ns","pi_nsdel/pi_ns")
rownames(results_pi) = c("wild","weed","landrace","elite")

for(i in c(1:4))
{
	i_f = c(4:1)
	results_pi[i_f[i],1] = mean(sites[sites[,i+2] != 0,i+2])
	results_pi[i_f[i],2] = mean(sites[sites[,i+2] != 0,i+10]) / mean(sites[sites[,i+2] != 0,i+6])
	results_pi[i_f[i],3] = mean(sites[sites[,i+2] != 0,i+14]) / mean(sites[sites[,i+2] != 0,i+10])
	results_pi[i_f[i],4] = mean(sites[sites[,i+2] != 0,i+18]) / mean(sites[sites[,i+2] != 0,i+10])	
}
print(results_pi)

par(mfrow = c(2,2))
par(mar = c(5,5,4,2))
barplot(results_pi[order(results_pi[,1]),1], col = c("dodgerblue4","gold","firebrick1","burlywood4")[order(results_pi[,1])],xaxt = "n",main = expression(Nucleotide~diversity~(bolditalic(pi))))

barplot(results_pi[order(results_pi[,2]),2],col = c("dodgerblue4","gold","firebrick1","burlywood4")[order(results_pi[,2])],xaxt = "n",main = expression(over(bolditalic(pi)[bold(n)],bolditalic(pi)[bold(s)])))

barplot(results_pi[order(results_pi[,3]),3], col = c("dodgerblue4","gold","firebrick1","burlywood4")[order(results_pi[,3])],xaxt = "n",main = expression(over(bolditalic(pi)[bold(n~(nonsense))],~bolditalic(pi)[bold(n~(sense))])))
par(xpd = T)
barplot(results_pi[order(results_pi[,4]),4], col = c("dodgerblue4","gold","firebrick1","burlywood4")[order(results_pi[,4])],xaxt = "n",main = expression(over(bolditalic(pi)[bold(n~(deleterious))],~bolditalic(pi)[bold(n~(neutral))])))
###legend down right center...
par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), mar=c(0, 0, 0, 0), new=TRUE)
plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
legend(x = -0.1, 0.2, cex = 1.4, legend = c("wild","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")


dev.print(device=pdf, "/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/results/Figure_pi.pdf", onefile=FALSE)
dev.off()


### ### ### ###
### SANDBOX ###
### ### ### ###
if(1 ==2)
{
pi_mean = rep(0,20)
for(i in 1:20)
{
	pi_mean[i] = mean(as.numeric(len[len[,(i+2)]!= 0,(i+2)]))
	
}
}
