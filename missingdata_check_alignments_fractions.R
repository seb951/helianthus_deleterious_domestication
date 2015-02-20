#!/usr/bin/Rscript --verbose

#args = commandArgs(TRUE)
#gene = as.numeric(args[1])


### ###
###TEST 4 DIFFERENT MISSING DATA THRESHOLDS!!!!###
### ### 
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/")#set up working directory 
snp_table = read.table("results/snp_table_all3",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...
snp_table_effect = read.table("results/new_snp_table_effect",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...

type = colnames(snp_table)[4:ncol(snp_table)]

individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!!
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])
individuals = individuals[-c(20),]
provean_results_all = read.delim("provean/new_provean_results_all", sep = " ", stringsAsFactors = F)#provean results


individuals = cbind(individuals,0)
colnames(individuals)[6] = "missing"
 
 for(i in 4:ncol(snp_table))
{
	individuals[i-3,6] = length(snp_table[snp_table[,i] == "XX",i])/ length(snp_table[,i]) 
}


boxplot(individuals[,6]~individuals[,5])
t.test(individuals[individuals[,5] == "wild",6],individuals[individuals[,5] == "elite",6])
anova(lm(individuals[,6]~individuals[,5]))

results1 = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
results5 = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
results10 = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
results20 = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

###redo dn(del) /dn(neutral) for sites missing less than 1%/5%/10%/20% of the data.
for(i in 1:nrow(snp_table))
#for(i in 1:1000)
{
	missing_filter = snp_table[i,c(3:ncol(snp_table))]
	missing_filter = length(missing_filter[missing_filter == "XX"]) / length(missing_filter)
	for(t in 1:40)	
		{
			if((snp_table_effect[i,t+3] == "ns") & (missing_filter < 0.01)) results1[[t]] = c(results1[[t]], snp_table_effect[i,45])
			if((snp_table_effect[i,t+3] == "ns") & (missing_filter < 0.05)) results5[[t]] = c(results5[[t]], snp_table_effect[i,45])
			if((snp_table_effect[i,t+3] == "ns") & (missing_filter < 0.10)) results10[[t]] = c(results10[[t]], snp_table_effect[i,45])
			if((snp_table_effect[i,t+3] == "ns") & (missing_filter < 0.20)) results20[[t]] = c(results20[[t]], snp_table_effect[i,45])
		}
	if(i %% 10000 == 0) print(paste(i,Sys.time()))
}

results_deleterious1 = matrix(0,nrow = nrow(individuals), ncol = 11)
results_deleterious5 = matrix(0,nrow = nrow(individuals), ncol = 11)
results_deleterious10 = matrix(0,nrow = nrow(individuals), ncol = 11)
results_deleterious20 = matrix(0,nrow = nrow(individuals), ncol = 11)
for(t in 1:40)
{
	results_deleterious1[t,3:2] = c(length(results1[[t]][results1[[t]] < 0]),length(results1[[t]]))
	results_deleterious1[t,4] = length(results1[[t]][results1[[t]] < -2.5])
	results_deleterious1[t,5] = length(results1[[t]][results1[[t]] < -4])
	results_deleterious1[t,6] = length(results1[[t]][results1[[t]] < -6])
	results_deleterious5[t,3:2] = c(length(results5[[t]][results5[[t]] < 0]),length(results5[[t]]))
	results_deleterious5[t,4] = length(results5[[t]][results5[[t]] < -2.5])
	results_deleterious5[t,5] = length(results5[[t]][results5[[t]] < -4])
	results_deleterious5[t,6] = length(results5[[t]][results5[[t]] < -6])
	results_deleterious10[t,3:2] = c(length(results10[[t]][results10[[t]] < 0]),length(results10[[t]]))
	results_deleterious10[t,4] = length(results10[[t]][results10[[t]] < -2.5])
	results_deleterious10[t,5] = length(results10[[t]][results10[[t]] < -4])
	results_deleterious10[t,6] = length(results10[[t]][results10[[t]] < -6])
	results_deleterious20[t,3:2] = c(length(results20[[t]][results20[[t]] < 0]),length(results20[[t]]))
	results_deleterious20[t,4] = length(results20[[t]][results20[[t]] < -2.5])
	results_deleterious20[t,5] = length(results20[[t]][results20[[t]] < -4])
	results_deleterious20[t,6] = length(results20[[t]][results20[[t]] < -6])
}
	results_deleterious1[,7] = round(c(results_deleterious1[,3] / results_deleterious1[,2]),5)
	results_deleterious1[,8] = round(c(results_deleterious1[,4] / results_deleterious1[,2]),5)
	results_deleterious1[,9] = round(c(results_deleterious1[,5] / results_deleterious1[,2]),5)
	results_deleterious1[,10] = round(c(results_deleterious1[,6] / results_deleterious1[,2]),5)
	results_deleterious1[,1] =  individuals[,5]
	results_deleterious1[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
	colnames(results_deleterious1) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")
	results_deleterious5[,7] = round(c(results_deleterious5[,3] / results_deleterious5[,2]),5)
	results_deleterious5[,8] = round(c(results_deleterious5[,4] / results_deleterious5[,2]),5)
	results_deleterious5[,9] = round(c(results_deleterious5[,5] / results_deleterious5[,2]),5)
	results_deleterious5[,10] = round(c(results_deleterious5[,6] / results_deleterious5[,2]),5)
	results_deleterious5[,1] =  individuals[,5]
	results_deleterious5[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
	colnames(results_deleterious5) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")
	results_deleterious10[,7] = round(c(results_deleterious10[,3] / results_deleterious10[,2]),5)
	results_deleterious10[,8] = round(c(results_deleterious10[,4] / results_deleterious10[,2]),5)
	results_deleterious10[,9] = round(c(results_deleterious10[,5] / results_deleterious10[,2]),5)
	results_deleterious10[,10] = round(c(results_deleterious10[,6] / results_deleterious10[,2]),5)
	results_deleterious10[,1] =  individuals[,5]
	results_deleterious10[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
	colnames(results_deleterious10) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")
	results_deleterious20[,7] = round(c(results_deleterious20[,3] / results_deleterious20[,2]),5)
	results_deleterious20[,8] = round(c(results_deleterious20[,4] / results_deleterious20[,2]),5)
	results_deleterious20[,9] = round(c(results_deleterious20[,5] / results_deleterious20[,2]),5)
	results_deleterious20[,10] = round(c(results_deleterious20[,6] / results_deleterious20[,2]),5)
	results_deleterious20[,1] =  individuals[,5]
	results_deleterious20[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
	colnames(results_deleterious20) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")



write.table(results_deleterious1,"results/results_deleterious1", col.names = T)
write.table(results_deleterious5,"results/results_deleterious5", col.names = T)
write.table(results_deleterious10,"results/results_deleterious10", col.names = T)
write.table(results_deleterious20,"results/results_deleterious20", col.names = T)


###plots DIFFERENT MISSING DATA THRESHOLDS 1, 5, 10, 20
results_deleterious1 = read.table("results/results_deleterious1", header = T, stringsAsFactors = F)
results_deleterious5 = read.table("results/results_deleterious5", header = T, stringsAsFactors = F)
results_deleterious10 = read.table("results/results_deleterious10", header = T, stringsAsFactors = F)
results_deleterious20 = read.table("results/results_deleterious20", header = T, stringsAsFactors = F)


par(mfrow = c(2,2))
par(mar = c(5,5,4,2))
barplot(as.numeric(results_deleterious1[order(as.numeric(results_deleterious1[,8])),8]), col = results_deleterious1[order(as.numeric(results_deleterious1[,8])),11], font = 2,main =   expression(over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious1[,8]))],font = 3, las = 3, cex = 0.3)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(crt = 90, side = 1, at = -3,  text = c("SNPs with <","1% missing data"), font = 1, line = c(-16,-15), cex = 1, crt= 90)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
legend(x = 0,0.2, cex = 1, legend = c("wild","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")
#
par(mar = c(5,5,4,2))
barplot(as.numeric(results_deleterious5[order(as.numeric(results_deleterious5[,8])),8]), col = results_deleterious5[order(as.numeric(results_deleterious5[,8])),11], font = 2,main =   expression(over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious5[,8]))],font = 3, las = 3, cex = 0.3)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(crt = 90, side = 1, at = -3,  text = c("SNPs with <","5% missing data"), font = 1, line = c(-16,-15), cex = 1, crt= 90)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
#
par(mar = c(5,5,4,2))
barplot(as.numeric(results_deleterious10[order(as.numeric(results_deleterious10[,8])),8]), col = results_deleterious10[order(as.numeric(results_deleterious10[,8])),11], font = 2,main =   expression(over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious10[,8]))],font = 3, las = 3, cex = 0.3)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
mtext(crt = 90, side = 1, at = -3,  text = c("SNPs with <","10% missing data"), font = 1, line = c(-16,-15), cex = 1, crt= 90)
#
par(mar = c(5,5,4,2))
barplot(as.numeric(results_deleterious20[order(as.numeric(results_deleterious20[,8])),8]), col = results_deleterious20[order(as.numeric(results_deleterious20[,8])),11], font = 2,main =   expression(over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious20[,8]))],font = 3, las = 3, cex = 0.3)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
mtext(crt = 90, side = 1, at = -3,  text = c("SNPs with <","20% missing data"), font = 1, line = c(-16,-15), cex = 1, crt= 90)


###
dev.print(device=pdf, "results/FigureS_missing.pdf", onefile=FALSE)
dev.off()






###
###alignments checks
###ARE THERE MORE READS ALIGNED FOR WILD, WEED, ELITE, LAND
system("ls -1 alignments/*realign.sorted.bam >list")
list = read.table("list", header = F, stringsAsFactors = F)

for(i in 1:40)
{
	system(paste("/usr/local/bin/samtools idxstats ",list[i,1]," >idxstats_temp", paste = ""))
	idxstats_temp = read.table("idxstats_temp", header = F, stringsAsFactors = F,sep = "\t")
	idxstats_temp[51469,3] = idxstats_temp[51469,4]
	if(i == 1) {idxstats = idxstats_temp[,c(2,3)]; colnames(idxstats) = c("sequence_length",list[i,1]);rownames(idxstats) = idxstats_temp[,1]}
	if(i != 1) {idxstats = cbind(idxstats, idxstats_temp[,3]); colnames(idxstats)[ncol(idxstats)] = list[i,1]}
}

fraction_aligned = (colSums(idxstats) - idxstats[51469,]) / colSums(idxstats)
barplot(as.matrix(fraction_aligned)[2:41], las = 3, names.arg = individuals[,5]) 

lm_total = lm(as.matrix(fraction_aligned)[2:41]~individuals[,5])
anova(lm_total)


}