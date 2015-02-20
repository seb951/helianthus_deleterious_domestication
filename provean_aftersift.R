##########
setwd("/Volumes/backup_ubuntu_vancouver/backup_vancouver/seb/Documents/domestication_load/")#set up working directory 
individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]; individuals = individuals[-c(20),] #remove the test run individuals ### individuals = individuals[1:6,] #for the test run AND MEN LANDRACE!!!!!!!!!!!!!!!!!!11
results = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])

snp_table = read.table("results/snp_table_all3",stringsAsFactors = F, header = T)### note FEBRUARY 6, 2015;MEN LANDRACE HAS ALREADY BEEN DISCARDED....snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...
type = colnames(snp_table)[4:ncol(snp_table)]

if(file.exists("results/new_results_deleterious_sift") == F) {
results_deleterious = matrix(0,nrow = length(type), ncol = 11)
cut_offs = c(0,-2.5,-4,-6)
for(t in 1:length(type))
{
results_temp = read.table(paste("provean/new_sift_results_",type[t],sep = ""), stringsAsFactors = F, sep = " ", skip = 1, header = F)
###PATCH
###this is because there may be instance were the mutations was called twice...
results_temp_sorter = apply(results_temp,1,paste, collapse = "@@@")
results_temp = rle(sort(results_temp_sorter))$values
results_temp = t(matrix(unlist(strsplit(results_temp,"@@@")),nrow = 7))
results_temp = as.data.frame(results_temp, stringsAsFactors = F)
#results_temp[,2] = as.numeric(results_temp[,2])
#results_temp = results_temp[!is.na(results_temp[,2]),]

#results_temp =  results_temp[1:4000,]
results[[t]] = results_temp
results_deleterious[t,3:2] = c(length(results[[t]][results[[t]][,2] ==  "TOLERATED",2]),length(results[[t]][,2]))
results_deleterious[t,4] = length(results[[t]][results[[t]][,2] == "DELETERIOUS",2])
results_deleterious[t,5] = length(results[[t]][results[[t]][,2] == "DELETERIOUS",2])
results_deleterious[t,6] = length(results[[t]][results[[t]][,2] == "DELETERIOUS",2])

#print(paste(mean(results_temp[,2]),type[t],individuals[t,5])) #same number of mutations for everyone...
}
results_deleterious[,7] = round(c(results_deleterious[,3] / results_deleterious[,2]),5)
results_deleterious[,8] = round(c(results_deleterious[,4] / results_deleterious[,2]),5)
results_deleterious[,9] = round(c(results_deleterious[,5] / results_deleterious[,2]),5)
results_deleterious[,10] = round(c(results_deleterious[,6] / results_deleterious[,2]),5)

results_deleterious[,1] =  individuals[,5]

results_deleterious[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue4",9),rep("gold",4),rep("dodgerblue4",7))
colnames(results_deleterious) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")
write.table(results_deleterious,"results/new_results_deleterious_sift",row.names = F, col.names = T, quote = F)} else
results_deleterious = as.matrix(read.table("results/new_results_deleterious_sift", header = T))  #file already exists

####
####PLOTS
####
library(betareg)
par(mar = c(8,5,4,2))
barplot(as.numeric(results_deleterious[order(as.numeric(results_deleterious[,8])),8]), col = results_deleterious[order(as.numeric(results_deleterious[,8])),11], font = 2,main =   expression(over(bolditalic(P)[bold(n~(deleterious))],~bolditalic(P)[bold(n~(neutral))])), ylim = c(0,0.2),space = 0,border = 0)
mtext(side = 1, at = seq(0.5,40,by = 1),line = 1, text = type[order(as.numeric(results_deleterious[,8]))],font = 3, las = 3, cex = 0.7)
mtext(side = 2, at = 0.1,  text = "Proportions", font = 1, line = 2.8, cex = 1.5)
mtext(side = 1, at = 20, text = expression(italic(H.~annuus)~~accessions), font = 1, line = 6, cex = 1.5)
text(x = -12, y = 0.25,  "D", font = 1, cex = 1.8, adj = 0,srt = 0,col = "black",xpd = T)
lm_del = lm(as.numeric(results_deleterious[,8])~as.factor(results_deleterious[,1])) #classes effect on number of sites
anova(lm_del)
TukeyHSD(aov(as.numeric(results_deleterious[,8])~as.factor(results_deleterious[,1])))
beta = betareg(as.numeric(results_deleterious[,8])~as.integer(as.factor(results_deleterious[,1])),link = c("probit"));summary(beta)#beta regression
wilcox.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace",8]),as.numeric(results_deleterious[results_deleterious[,1] == "elite",8])) #landrace versus elites
wilcox.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace",8]),as.numeric(results_deleterious[results_deleterious[,1] == "wild",8])) #domesticated versus wild
wilcox.test(as.numeric(results_deleterious[results_deleterious[,1] == "elite",8]),as.numeric(results_deleterious[results_deleterious[,1] == "wild",8])) 
t.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace",8]),as.numeric(results_deleterious[results_deleterious[,1] == "elite",8])) #landrace versus elites
t.test(as.numeric(results_deleterious[results_deleterious[,1] == "landrace" | results_deleterious[,1] == "elite",8]),as.numeric(results_deleterious[results_deleterious[,1] == "wild",8])) #domesticated versus wild
###

dev.print(device=pdf, "results/Figure1d_sift.pdf", onefile=FALSE)
dev.off()


####correlations between PROVEAN AND SIFT
t = 1
results_temp = read.table(paste("provean/new_sift_results_",type[t],sep = ""), stringsAsFactors = F, sep = " ", skip = 1, header = F)
results_temp_sorter = apply(results_temp,1,paste, collapse = "@@@")
results_temp_matcher = apply(results_temp[,c(1,7)],1,paste, collapse = "@@@")
results_temp = rle(sort(results_temp_sorter))$values
results_temp_matcher_sift = rle(sort(results_temp_matcher))$values
results_temp = t(matrix(unlist(strsplit(results_temp,"@@@")),nrow = 7))
results_temp_sift = as.data.frame(results_temp, stringsAsFactors = F)
results_temp_sift = results_temp[!is.na(results_temp[,2]),]
results_temp_matcher_sift = results_temp_matcher_sift[!is.na(results_temp[,2])]


results_temp = read.table(paste("provean/new_provean_results_",type[t],sep = ""), stringsAsFactors = F, sep = " ", skip = 1, header = F)
results_temp_sorter = apply(results_temp,1,paste, collapse = "@@@")
results_temp_matcher = apply(results_temp[,c(1,3)],1,paste, collapse = "@@@")
results_temp = rle(sort(results_temp_sorter))$values
results_temp_matcher_provean = rle(sort(results_temp_matcher))$values
results_temp = t(matrix(unlist(strsplit(results_temp,"@@@")),nrow = 3))
results_temp = as.data.frame(results_temp, stringsAsFactors = F)
results_temp[,2] = as.numeric(results_temp[,2])
results_temp_provean = results_temp[!is.na(results_temp[,2]),]
results_temp_matcher_provean = results_temp_matcher_provean[!is.na(results_temp[,2])]



sift_provean_match = cbind(results_temp_matcher_provean,0,0,0)

for(i in 1: nrow(results_temp_provean))
{
  temp = results_temp_matcher_sift %in%  results_temp_matcher_provean[i]
  sift_provean_match[i,2] = results_temp_provean[i,2]
  if(length(temp[temp == T]) == 1) sift_provean_match[i,c(3,4)] = results_temp_sift[temp == T,c(2,3)]
}

###boxplot###
sift_provean_match = sift_provean_match[sift_provean_match[,3] != "0",]
boxplot(as.numeric(sift_provean_match[,2])~sift_provean_match[,3], ylab = "PROVEAN score", xlab = "SIFT category")
t.test(as.numeric(sift_provean_match[sift_provean_match[,3] == "TOLERATED",2]),as.numeric(sift_provean_match[sift_provean_match[,3] == "DELETERIOUS",2]))
text(1,11, label = "t-test, p-value < 2e-16")
dev.print(device=pdf, "results/Figure1d_sift_provean_correlation.pdf", onefile=FALSE)
dev.off()





