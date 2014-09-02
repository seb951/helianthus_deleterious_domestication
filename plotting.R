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
snp_table = read.table("results/snp_table_all3_new13",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...

colnames(snp_table) = gsub("X.home.transfer.Documents.new.bam_files_jan2012.","",colnames(snp_table)) #shorten names of the columns.
type = colnames(snp_table)[4:ncol(snp_table)]

individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!!
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])

provean_results_all = read.delim("provean/new_provean_results_all", sep = " ", stringsAsFactors = F)#provean results
 
unique_orf =  as.matrix(read.table("reference/unique_orf",header = T, stringsAsFactors = F))

ns_s = cbind(rep(0,nrow(snp_table)),0,0) #add a column to see what kind of SNP it is and the SNP EFFECT FROM PROVEAN...
names = paste(c(colnames(snp_table),"ns_site","score","freq_wild","freq_land","freq_eli","freq_weed","freq_all"),collapse = " ") #column names for the result file.

snp_table_effect = read.table("results/new_snp_table_effect", header = T, stringsAsFactors = F);snp_table_effect = snp_table_effect[,-c(23)] #remove MEN

###
###summarize results###
###
individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### ndividuals = individuals[1:6,] #for the test run

individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])
individuals = individuals[-c(20),]  #remove MEN

#snp_table_effect = snp_table_effect[1:10000,]
type = colnames(snp_table_effect)[4:(nrow(individuals)+3)]
results_ns_s = matrix(0,nrow = length(type), ncol = 8)
results_ns_s[,8] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
colnames(results_ns_s) = c("accessions","total","non-syn","syn","nc","ns/s","STOP","color")

for(t in 1:length(type))
{
valid_snp = snp_table_effect[snp_table_effect[,(t+3)] != 0,(t+3)]

results_ns_s[t,2] =  length(valid_snp) #total length
results_ns_s[t,3] = length(valid_snp[valid_snp == "ns"])  #ns 
results_ns_s[t,4] = length(valid_snp[valid_snp == "s"])  #s 
results_ns_s[t,6] = as.numeric(results_ns_s[t,3]) / (as.numeric(results_ns_s[t,4])) # ns /  s 
results_ns_s[t,5] = length(valid_snp[valid_snp == "nc"])  #nc
results_ns_s[t,7] = length(valid_snp[valid_snp == "STOP"])  #nc
}
results_ns_s[,1] = individuals[,5]



###
###plotting total and ns total
###
par(mfrow = c(2,2))
par(mar = c(5,5,4,2))
stacks = results_ns_s[order(as.numeric(results_ns_s[,2])),]
stacks = cbind(as.numeric(stacks[,3]),as.numeric(stacks[,4]),as.numeric(stacks[,5]),as.numeric(stacks[,7]))
barplot(as.numeric(results_ns_s[order(as.numeric(results_ns_s[,2])),2]), col = results_ns_s[order(as.numeric(results_ns_s[,2])),8], font = 2,main = "Polymorphic sites",  space = 0,border = 0, yaxt = "n",ylim = c(0,250000))	#all
barplot(t(stacks[,1:2]), space = 0, add = T, col = c("#FFFFFFFF","#FFFFFF99"), yaxt = "n")

mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_ns_s[,2]))], font = 3, las = 3,cex = 0.3)
mtext(side = 2, at = max(as.numeric(results_ns_s[,2])) / 2,  text = "Number of sites (x 1,000)", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20.5, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
text(x = -12, y = 300000,  "A", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)
axis(side = 2, at = seq(0,250000,by = 50000), labels = c(0,"50","100","150","200","250"), las = 1,font = 2,crt = 90)
text(x = 41.5, y = 25000,  "Non-syno\nnymous", font = 1,  cex = 0.8, adj = 0,srt = 45,col = "black",xpd = T)
text(x = 41.5, y = 100000,   "Synony-\nmous", font = 1, cex = 0.8, adj = 0,srt = 45,col = "black",xpd = T)
text(x = 41.5, y = 200000,  "Non\ncoding ", font = 1, cex = 0.8, adj = 0,srt = 45,col = "black",xpd = T)
legend(x = 0, 250000, cex = 1, legend = c("wild","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")

####
lm_total = lm(as.numeric(results_ns_s[,2])~as.factor(results_ns_s[,1])) #classes effect on number of sites
anova(lm_total)
t.test(as.numeric(results_ns_s[results_ns_s[,1] == "landrace",2]),as.numeric(results_ns_s[results_ns_s[,1] == "elite",2])) #landrace versus elites
t.test(as.numeric(results_ns_s[results_ns_s[,1] == "landrace" | results_ns_s[,1] == "elite",2]),as.numeric(results_ns_s[results_ns_s[,1] == "wild",2])) #domesticated versus wild

sum(as.numeric(results_ns_s[,3])) / sum(as.numeric(results_ns_s[,2])) #ns proportions
sum(as.numeric(results_ns_s[,4])) / sum(as.numeric(results_ns_s[,2])) #s proportions
sum(as.numeric(results_ns_s[,5])) / sum(as.numeric(results_ns_s[,2])) #nc proportions
sum(as.numeric(results_ns_s[,7])) / sum(as.numeric(results_ns_s[,2])) #STOP proportions


###
###plotting ns / s
###
par(mar = c(5,5,4,2))
barplot(as.numeric(results_ns_s[order(as.numeric(results_ns_s[,6])),6]), col = results_ns_s[order(as.numeric(results_ns_s[,6])),8], font = 2,main = expression(over(bolditalic(D)[bold(n)],bolditalic(D)[bold(s)])), ylim = c(0,0.51),space = 0,border = 0)
mtext(side = 2, at = max(as.numeric(results_ns_s[,6])) / 2,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20.5, text =  expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
mtext(side = 1, at =seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_ns_s[,6]))],font = 3, las = 3, cex = 0.3) #0.3
text(x = -12, y = 0.64,  "B", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)
lm_ns_s = lm(as.numeric(results_ns_s[,6])~as.factor(results_ns_s[,1])) #classes effect on number of sites
anova(lm_ns_s)
t.test(as.numeric(results_ns_s[results_ns_s[,1] == "landrace",6]),as.numeric(results_ns_s[results_ns_s[,1] == "elite",6])) #landrace versus elites
t.test(as.numeric(results_ns_s[results_ns_s[,1] == "landrace" | results_ns_s[,1] == "elite",6]),as.numeric(results_ns_s[results_ns_s[,1] == "wild",6])) #domesticated versus wild
###
###plotting STOP / ns
###
par(mar = c(5,5,4,2))
stop = cbind(results_ns_s[,1],100*as.numeric(results_ns_s[,7])/as.numeric(results_ns_s[,3]))
barplot(as.numeric(stop[order(as.numeric(stop[,2])),2]), col = results_ns_s[order(as.numeric(stop[,2])),8], font = 2,main = expression(over(bolditalic(D)[bold(n~(missense))],~bolditalic(D)[bold(n)])),space = 0,border = 0, ylim = c(0,1))
mtext(side = 2, at = max(as.numeric(stop[,2])) / 2,  text = "Proportions (X 100)", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20.5, text =  expression(italic(H.~annuus)~~accessions), font = 2, line = 3, cex = 1)
mtext(side = 1, at =seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(stop[,2]))],font = 3, las = 3, cex = 0.3)
text(x = -12, y = 1.26,  "C", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)
lm_stop = lm(as.numeric(stop[,2])~as.factor(stop[,1])) #classes effect on number of sites
anova(lm_stop)
t.test(as.numeric(stop[stop[,1] == "landrace",2]),as.numeric(stop[stop[,1] == "elite",2])) #landrace versus elites
t.test(as.numeric(stop[stop[,1] == "landrace" | stop[,1] == "elite",2]),as.numeric(stop[stop[,1] == "wild",2])) #domesticated versus wild
###
###
###	###	###	###
###cut off -2.5###
###	###	###	###
par(mar = c(5,5,4,2))
barplot(as.numeric(results_deleterious[order(as.numeric(results_deleterious[,8])),8]), col = results_deleterious[order(as.numeric(results_deleterious[,8])),11], font = 2,main =   expression(over(bolditalic(D)[bold(n~(deleterious))],~bolditalic(D)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious[,8]))],font = 3, las = 3, cex = 0.3)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 3, cex = 1)
text(x = -12, y = 0.25,  "D", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)
lm_del = lm(as.numeric(results_deleterious[,8])~as.factor(results_deleterious[,1])) #classes effect on number of sites
anova(lm_del)
t.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace",8]),as.numeric(results_deleterious[results_deleterious[,1] == "elite",8])) #landrace versus elites
t.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace" | results_deleterious[,1] == "elite",8]),as.numeric(results_deleterious[results_deleterious[,1] == "wild",8])) #domesticated versus wild
###


dev.print(device=pdf, "results/Figure1.pdf", onefile=FALSE)
dev.off()





###
###frequency distribution of the deleterious and neutral mutations
###

sn = snp_table_effect[snp_table_effect[,44] != 0, 44:50]
par(mfrow = c(2,2))
sn_sub = sn[sn[,3] != 0,] #make sure they are present in the population you are looking at...
x = hist(sn_sub[sn_sub[,2] > -2.5,3],breaks = seq(0,1,by = 0.01), plot = F)
x1 = hist(sn_sub[sn_sub[,2] < -2.5 ,3],breaks = seq(0,1,by = 0.01), plot = F)
plot(x1$mids+0.00,x1$density, col = "red", lwd = 4, type = "h",xlab = "Allele frequency", ylim = c(0,20),ylab = "Density",main = "Wild") # ns, non deleterious, freq in wild
points(x$mids+0.01,x$density,  lwd = 4, type = "h",col = "#000000FF")
legend(x = 0.4, 20, cex = 1, legend = c(expression(italic(D)[n~(neutral)]), expression(italic(D)[n~(deleterious)])), fill = c("black","red"), bty = "o", bg = "#FFFFFF80")
ks.test(sn_sub[sn_sub[,2] > -2.5,3],sn_sub[sn_sub[,2] < -2.5,3])
text(x = -0.1, y = 23.2,  "A", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)

sn_sub = sn[sn[,6] != 0,] #make sure they are present in the population you are looking at...
x = hist(sn_sub[sn_sub[,2] < -2.5,6],breaks = seq(0,1,by = 0.01), plot = F)
y = hist(sn_sub[sn_sub[,2] > -2.5,6],breaks = seq(0,1,by = 0.01), plot = F)
plot(x$mids,x$density, col = "red", lwd = 4, type = "h",ylim = c(0,40), xlab = "Allele frequency",ylab = "Density",main = "Weed") # ns, non deleterious, freq in wild
points(y$mids+0.02,y$density, col = "#000000FF", lwd = 4, type = "h")
ks.test(sn_sub[sn_sub[,2] > -2.5,6],sn_sub[sn_sub[,2] < -2.5,6])
text(x = -0.1, y = 46,  "B", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)

sn_sub = sn[sn[,4] != 0,] #make sure they are present in the population you are looking at...
x = hist(sn_sub[sn_sub[,2] < -2.5,4],breaks = seq(0,1,by = 0.01), plot = F)
y = hist(sn_sub[sn_sub[,2] > -2.5,4],breaks = seq(0,1,by = 0.01), plot = F)
plot(x$mids,x$density, col = "red", lwd = 4, type = "h",ylim = c(0,15), xlab = "Allele frequency", ylab = "Density",main = "Landrace") # ns, non deleterious, freq in wild
points(y$mids+0.02,y$density, col = "#000000FF", lwd = 4, type = "h")
ks.test(sn_sub[sn_sub[,2] > -2.5,4],sn_sub[sn_sub[,2] < -2.5,4])
text(x = -0.1, y = 18,  "C", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)

sn_sub = sn[sn[,5] != 0,] #make sure they are present in the population you are looking at...
x = hist(sn_sub[sn_sub[,2] < -2.5,5],breaks = seq(0,1,by = 0.01), plot = F)
y = hist(sn_sub[sn_sub[,2] > -2.5,5],breaks = seq(0,1,by = 0.01), plot = F)
plot(x$mids,x$density, col = "red", lwd = 4, type = "h",ylim = c(0,35),xlab = "Allele frequency", ylab = "Density",main = "Elite") # ns, non deleterious, freq in wild
points(y$mids+0.02,y$density, col = "#000000FF", lwd = 4, type = "h")
ks.test(sn_sub[sn_sub[,2] > -2.5,5],sn_sub[sn_sub[,2] < -2.5,5])
text(x = -0.1, y = 41,  "D", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)

dev.print(device=pdf, "results/Figure3_allele frequency NS neutral and deleterious.pdf", onefile=FALSE)
dev.off()

###correlation between effect of del mutation and frequency
x = sn[,2]
y = sn[,3]
plot(x,y,ylim=c(0,1), xlim = c(-15,15), xlab = "Provean score", ylab = "Frequency (wild)")
mylm<-lm(y~x)
abline(mylm,col="red", lwd = 4)
newx<- seq(-16,16)
prd<-predict(mylm,newdata=data.frame(x=newx),interval = c("confidence"), level = 0.95,type="response")
lines(newx,prd[,2],col="blue",lty=2)
lines(newx,prd[,3],col="blue",lty=2)
text(x = 9, col = "red",y = 0.95,label = paste("r = ",signif(cor.test(sn[,2],sn[,3])$estimate,2),", p-value < 2e-16",sep = ""),cex = 1.2)
dev.print(device=pdf, "results/Figure_S2_delerious_frequency_correlation.pdf", onefile=FALSE)
dev.off()


###

###
###PLOT THE DIFFERENCE FOR THE PRIVATE ALLELES (exclude the weedy guys)
###
results_private_deleterious = as.matrix(read.table("results/new_results_private_deleterious", header = T))  #file already exists
results_deleterious_w = results_deleterious[results_deleterious[,1] != "weed", ]
results_private_deleterious_w = results_private_deleterious[results_private_deleterious[,1] != "weed", ]
type_w = type[results_deleterious[,1] != "weed"]

for(t in 1:length(type_w))
{
if((results_private_deleterious_w[t,3] != "0") & (chisq.test(cbind(as.numeric(results_deleterious_w[t,c(3:2)]),as.numeric(results_private_deleterious_w[t,c(3:2)])))$p.value < (0.05/length(type_w)))) results_private_deleterious_w[t,12] = 1 #no FDR correction, otherwise, just divided 0.05 / nb samples (36)
}


results_deleterious[results_deleterious[,11] == "dodgerblue3",11] = "dodgerblue4"
results_private_deleterious[results_private_deleterious[,11] == "dodgerblue3",11] = "dodgerblue4"

par(mar = c(7,5,4,2))
barplot(as.numeric(results_deleterious_w[order(as.numeric(results_deleterious_w[,8])),8]), col = results_deleterious_w[order(as.numeric(results_deleterious_w[,8])),11], font = 2,main =  "Private alleles", ylim = c(0,0.3),space = 0,border = 0, cex.main = 1.8)
change = results_private_deleterious_w[order(as.numeric(results_deleterious_w[,8])),c(7,12)]   # change in frequency for private mutations
barplot(as.numeric(change[,1]), col = "#99999980", font = 2,add = T,border = 0,space = 0,axes = F)
mtext(side = 1, at = seq(0.5,35,by = 1),line = 1, text = type_w[order(as.numeric(results_deleterious_w[,8]))],font = 3, las = 3,cex = 0.8)

mtext(side = 2, at = max(as.numeric(results_private_deleterious_w[,7])) / 2,  text = "Proportions", font = 1, line = 2.8, cex = 1.4)
mtext(side = 1, at = 18.5, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 5.6, cex = 1.4)

legend(x = 0, y = 0.3, cex = 1, legend = c("wild","landrace","elite",expression(italic(D)[n~(deleterious~bold(private))]~~"/"~~italic(D)[n~(neutral~bold(private))])), fill = c("dodgerblue4","firebrick1","burlywood4", "#99999980"), bty = "o", bg = "#FFFFFF80")


text(x = 14.7, y = 0.275, labels = expression(italic(D)[n~(deleterious)]~~"/"~~italic(D)[n~(neutral)]))
text(x = 8, y = 0.275, "}", cex = 3.5)

dev.print(device=pdf, "results/Figure2_mutation_private_frequencies.pdf", onefile=FALSE)
dev.off()


###
###
###
results_deleterious[results_deleterious[,11] == "dodgerblue3",11] = "dodgerblue4"

par(mfrow = c(2,2), mar = c(2,4,4,2))
for(i in 7:10)
{
#lim = c(0.2,0.1,0.04,0.015)
cutoff = c(0,-2.5,-4,-6)
lim = c(0.8,0.2,0.1,0.04)
barplot(as.numeric(results_deleterious[order(as.numeric(results_deleterious[,i])),i]), col = results_deleterious[order(as.numeric(results_deleterious[,i])),11], font = 2,main = paste("Cut-off = ",cutoff[i-6],sep = ""), font.main = 2, cex.main = 2,space = 0,border = 0, ylim = c(0,lim[i-6]))

mtext(side = 1, at = 20.5, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 0.6, cex = 1.2)
mtext(side = 2, at = max(as.numeric(results_deleterious[,i])) / 2,  text = "proportions (deleterious)", font = 1, line = 2.2, cex = 1.2)

if(i == 7) {
#legend(x = 0, y = 0.8, cex = 1, legend = c("  wild / wild (Texas)","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")
legend(x = 0, y = 0.8, cex = 1, legend = c("wild","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")
#legend(x = 2.5, y = 0.2, cex = 1, legend = c(""), fill = c("dodgerblue3"), bty = "n")
}


}

dev.print(device=pdf, "results/Figure_S1_mutation_frequencies_cutoff.pdf", onefile=FALSE)
dev.off()


###
###
###
###SANDBOX
###
###
###



par(mar = c(5,5,2,2))
for(t in 1:length(type))
	{	
		if(t == 1) {plot(0,0,xlim = c(-12,8),ylim = c(0,0.55), xlab = "effect of non-synonymous mutations", ylab = "Frequency", col = "white", cex.lab = 1.5, font.lab = 2);rect(-15,0,0,0.55,col = "#99000070", border = NA)}
		lines(density(results[[t]][ ,2], n = 500),lwd =6, col = results_deleterious[t,11])
	}

legend(x = 1, y = 0.5, cex = 1, legend = c("  wild / wild (Texas)","weed","landrace","elite"), fill = c("dodgerblue4","gold","firebrick1","burlywood4"), bty = "o", bg = "#FFFFFF80")
legend(x = 1.6, y = 0.5, cex = 1, legend = c(""), fill = c("dodgerblue3"), bty = "n")

text(-7,0.5, "deleterious alleles", cex = 1.5, font = 2, col = "black")

dev.print(device=pdf, "results/mutation_distribution.pdf", onefile=FALSE)
dev.off()


####correlations
plot(as.numeric(results_deleterious[,2]),as.numeric(results_deleterious[,7]), lwd = 4, xlab = "", ylab = "")
mtext(side = 2, at = 0.14,  text = "proportions (ns delet.)", font = 1, line = 2.5, cex = 1.8)
mtext(side = 1, at = max(as.numeric(results_deleterious[,2])) / 2, text =  "Numbers (ns total)", font = 1, line = 3, cex = 1.8)
dev.print(device=pdf, "results/correlations.pdf", onefile=FALSE)
dev.off()



results_frequency = cbind(results_frequency,0)


for(u in 1:nrow(results_frequency))
	{
		temp = effect[u,effect[u,] != 0]
		if(length(temp) > 0) results_frequency[u,5] = min(temp)

		if(u %% 10000 == 0) print(u)
	}
	
	wild_h = hist(results_frequency[results_frequency[,3] != 0,3], breaks = seq(0,1,by = 0.05), plot = F)
	elite_h = hist(results_frequency[results_frequency[,2] != 0,2], breaks = seq(0,1,by = 0.05), plot = F)
	land_h = hist(results_frequency[results_frequency[,1] != 0,1], breaks = seq(0,1,by = 0.05), plot = F)
	
	plot(wild_h$mids[wild_h$density != 0],wild_h$density[wild_h$density != 0], type = "l", lwd = 14, ylim = c(0,9), ylab = "Frequency", xlab = "Population frequency of deleterious mutations")
	points(elite_h$mids[elite_h$density != 0]+ 0.01,elite_h$density[elite_h$density != 0], type = "l", lwd = 14, col = "red")
	points(land_h$mids[land_h$density != 0] + 0.02,land_h$density[land_h$density != 0], type = "l", lwd = 14, col = "blue")








