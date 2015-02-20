#!/usr/bin/Rscript --verbose


library(qvalue)
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/")#set up working directory 
snp_table_effect = read.table("results/new_snp_table_effect", header = T, stringsAsFactors = F);snp_table_effect = snp_table_effect[,-c(23)] #remove MEN

unique_map_transcript = read.delim("../speciation_islands_dec2011/reference/recombination_rates", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###

transcript = read.delim("../speciation_islands_dec2011/reference/recombination_rates_all_sites", header = T, stringsAsFactors = F,sep = " ") ###recombination rates###

snp_table_effect_r = cbind(snp_table_effect,0); colnames(snp_table_effect_r)[ncol(snp_table_effect_r)] = "recomb_rate"

snp_table_effect_r = snp_table_effect_r[rowSums(snp_table_effect_r[,45:48]) != 0, ] #keep only the genes polymorphic in the domesticates
#snp_table_effect_r = snp_table_effect_r[rowSums(snp_table_effect_r[,c(45,45)]) != 0, ] #keep only the genes polymorphic in the wild
snp_table_effect_r = snp_table_effect_r[snp_table_effect_r[,43] != 0, ] #keep only the ns SNPs

###for each gene, how many deleterious and total
transcript = cbind(transcript,0,0); colnames(transcript)[(ncol(transcript)-1):ncol(transcript)]  = c("del2_5","total")
	
for(r in 1:nrow(transcript))
	{
		temp = snp_table_effect_r[,1] %in% transcript[r,1] #what are the mutations for transcript r
		if(length(temp[temp == T]) > 0) 
			{temp2 = snp_table_effect_r[temp == T,]; 
			temp2 = temp2[rowSums(temp2[,45:48]) > 0,] #make sure the mutation is found in the domesticates
			#temp2 = temp2[rowSums(temp2[,c(45,45)]) > 0,] #make sure the mutation is found in the 	WILD
			 transcript[r,8:9] = c(length(temp2[temp2[,44] < -2.5,1]), length(temp2[temp2[,44] != 0,1]))} # fraction of deleterious and neutrals per gene
	
		if(r %% 1000 == 0) print(paste(r,"of",nrow(transcript),"Time is:",Sys.time()))
	}

###then for each position, how many deleterious and total
unique_map_transcript = cbind(unique_map_transcript,0,0,0); colnames(unique_map_transcript)[(ncol(unique_map_transcript)-2):ncol(unique_map_transcript )]  = c("del2_5","total","prop")

for(r in 1:nrow(unique_map_transcript))
	{
upperbound = transcript[(unique_map_transcript[r,4] +0.5)  >= transcript[,4] ,] #0.5 centimorgan window on each side
lowerbound = upperbound[upperbound[,4] > (unique_map_transcript[r,4] - 0.5) , ] #0.5 centimorgan window on each side

unique_map_transcript[r,8:10] = c(sum(lowerbound[,8]), sum(lowerbound[,9]), (sum(lowerbound[,8]) / sum(lowerbound[,9]) )    ) 
	
		if(r %% 1000 == 0) print(paste(r,"of",nrow(unique_map_transcript),"Time is:",Sys.time()))
	}

###resampling test to see if the regions with more deleterious SNP have a lower recombination rate...
resample_vector = rep(0,nrow(unique_map_transcript))
for(r in 1:nrow(unique_map_transcript))
	{
	for(s in 1:1000000)
		{
			ss = sample(snp_table_effect_r[,44],unique_map_transcript[r,9]) #random sample
			if((length(ss) != 0) & (unique_map_transcript[r,10] > (length(ss[ss < -2.5]) / length(ss))  ) ) resample_vector[r] = resample_vector[r] + 1
		}
	if(r %% 1000000 == 0) print(paste(r,"of",nrow(unique_map_transcript),"Time is:",Sys.time()))
	resample_vector[r] = signif((1000000 - resample_vector[r]) / 1000000,4) # pvalue
	if(resample_vector[r] == 0) resample_vector[r] = 1e-06
  }

unique_map_transcript_delmutations = cbind(unique_map_transcript,resample_vector)
resample_vector = unique_map_transcript_delmutations[,11]
unique_map_transcript_delmutations = cbind(unique_map_transcript_delmutations,qvalue(unique_map_transcript_delmutations[,11])$qvalues)
colnames(unique_map_transcript_delmutations)[11:12] = c("resampled_pvalue1000","qvalue")

write.table(unique_map_transcript_delmutations , "results/new_unique_map_transcript_delmutations", row.names = F)

			unique_map_transcript_delmutations = read.table("results/new_unique_map_transcript_delmutations", header = T)
	#		unique_map_transcript_delmutations = unique_map_transcript_delmutations[!is.na(unique_map_transcript_delmutations[,10]),]
			results  = rep(0,7);names(results)  = c("number_regions","ttest","wilcox","mean_noexcess","median_noexcess","mean_excess","median_excess");pval = 0.05
			temp_region = rle(unique_map_transcript_delmutations[,11] < pval)$values ; results[1] = length(temp_region[temp_region == T])
			results[2] = 	t.test(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] > pval,6],unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] < pval,6])$p.value # significance
			results[3] = 	wilcox.test(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] > pval,6],unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] < pval,6])$p.value # significance
			results[4] = 	signif(mean(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] > pval,6],na.rm = T),4) #mean recomb. rate for regions not containing excess of del. mutations
			results[5] = 	signif(median(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] > pval,6],na.rm = T),4) #mean recomb. rate for regions  containing excess of del. mutations
			results[6] = 	signif(mean(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] < pval,6],na.rm = T),4) #mean recomb. rate for regions  containing excess of del. mutations
			results[7] = 	signif(median(unique_map_transcript_delmutations[unique_map_transcript_delmutations[,11] < pval,6],na.rm = T),4) #mean recomb. rate for regions  containing excess of del. mutations
			print(results)

write.table(results,"results/new_table1_recombrate_deleterious")

###
###chromosome 10
###


if(1 == 2)
{
unique_map_transcript_delmutations = read.table("results/new_unique_map_transcript_delmutations", header = T)

c = 10
lg_x = unique_map_transcript_delmutations[unique_map_transcript_delmutations[,2] == c,]
#lg_x = unique_map_transcript_delmutations
par(mar = c(5,5,4,6))

plot(lg_x[lg_x[,12] > 0.05,3],log(lg_x[lg_x[,12] > 0.05,11],10),type = "p", pch = 19, col = "#00990070",ylim = c(-5,2.5),xlim = c(0,90),xaxt = "n",yaxt = "n",xlab = "",main = "Linkage Group 10",ylab = "")
#plot(lg_x[,3],(lg_x[,10]),type = "p", pch = 19, col = "#00990099",ylim = c(-5,2.5),xlim = c(0,90),xaxt = "n",yaxt = "n",xlab = "",ylab = "")
lines(lg_x[,3],log(lg_x[,6],10),type = "l",col = "black", lwd = 2,ylim = c(-5,2.5),xlim = c(0,100))
points(lg_x[lg_x[,12]<=0.05,3],ifelse(log(lg_x[lg_x[,12]<=0.05,11],10) == -Inf,-5,log(lg_x[lg_x[,12]<=0.05,11],10)),type = "p", pch = 17, col = 
c("#00990070","#00990070","#00990070","#009900FF","#00990070","#009900FF"), cex=1.5)
axis(side = 2, at = c(-3,-2,-1,0,1,2), labels = c(0.001,0.01,0.1,1,10,100),las = 1,col = "black") #left 
mtext(side = 2, "Recombination rate (cM / MB)",col = "black",line = 3.3,lwd = 2)

axis(side = 1,  at = c(0,10,20,30,40,50,60,70,80,90), labels = c(0,10,20,30,40,50,60,70,80,90),col = "black") #bottom 
mtext(side = 1, "Distance (cM)",col = "black",line = 3)

axis(side = 4,at = c(-5,-4,-3,-2,-1,0), labels = c(0.00001,0.0001,0.001,0.01,0.1,1),las = 1,col.lab = "#009900",col.axis = "#009900") #right 
mtext(side = 4,expression(bolditalic(p)~bold(-value)~(~over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))]))),col = "#009900",line = 4.9,at = -3)

dev.print(device=pdf, "results/Figure4_recomb_rate_LG10_domesticated.pdf", onefile=FALSE)
dev.off()


#plot(unique_map_transcript[,10],unique_map_transcript[,6])

map = unique_map_transcript[!is.na(unique_map_transcript[,10]),]
map = map[map[,10] != 0,]
map = map[!is.na(map[,6]),]
map = map[map[,6] > 0.000001,]


map = map[map[,6] > 0.000001,]

plot(unique_map_transcript[,10],log(unique_map_transcript[,6]))
plot(log(map[,6],10),map[,10])

}


