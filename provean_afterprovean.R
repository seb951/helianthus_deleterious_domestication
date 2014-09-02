##########
setwd("~/Documents/domestication_load")#set up working directory 
individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]; individuals = individuals[-c(20),] #remove the test run individuals ### individuals = individuals[1:6,] #for the test run AND MEN LANDRACE!!!!!!!!!!!!!!!!!!11
results = list(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)

individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])

snp_table = read.table("results/snp_table_all3_new13",stringsAsFactors = F, header = T);snp_table = snp_table[,-c(23)] #remove the MEN landrace. It shouldn't be there...
colnames(snp_table) = gsub("X.home.transfer.Documents.new.bam_files_jan2012.","",colnames(snp_table)) #shorten names of the columns.
type = colnames(snp_table)[4:ncol(snp_table)]

if(file.exists("results/new_results_deleterious") == F) {
results_deleterious = matrix(0,nrow = length(type), ncol = 11)
cut_offs = c(0,-2.5,-4,-6)
for(t in 1:length(type))
{
results_temp = read.table(paste("provean/new_provean_results_",type[t],sep = ""), stringsAsFactors = F, sep = " ", skip = 1, header = F)
#results_temp =  results_temp[1:4000,]
results[[t]] = results_temp
results_deleterious[t,3:2] = c(length(results[[t]][results[[t]][,2] < 0,2]),length(results[[t]][,2]))
results_deleterious[t,4] = length(results[[t]][results[[t]][,2] < -2.5,2])
results_deleterious[t,5] = length(results[[t]][results[[t]][,2] < -4,2])
results_deleterious[t,6] = length(results[[t]][results[[t]][,2] < -6,2])

print(paste(mean(results_temp[,2]),type[t],individuals[t,5])) #same number of mutations for everyone...
}
results_deleterious[,7] = round(c(results_deleterious[,3] / results_deleterious[,2]),5)
results_deleterious[,8] = round(c(results_deleterious[,4] / results_deleterious[,2]),5)
results_deleterious[,9] = round(c(results_deleterious[,5] / results_deleterious[,2]),5)
results_deleterious[,10] = round(c(results_deleterious[,6] / results_deleterious[,2]),5)

results_deleterious[,1] =  individuals[,5]

results_deleterious[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue3",9),rep("gold",4),rep("dodgerblue4",7))
colnames(results_deleterious) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-0","prop_del-2.5","prop_del-4","prop_del-6","color")
write.table(results_deleterious,"results/new_results_deleterious",row.names = F, col.names = T, quote = F)} else
results_deleterious = as.matrix(read.table("results/new_results_deleterious", header = T))  #file already exists

#################################
#################################
#################################FIND THE PRIVATE MUTATIONS IN THE SET.
#################################
#################################
if(file.exists("results/new_results_private_deleterious") == F) {

results_matrix = NULL
for(t in 1:nrow(individuals))
	{	
		results_matrix = rbind(results_matrix,cbind(results[[t]],paste(type[t],individuals[t,5],sep = "_"),paste(results[[t]][,1],results[[t]][,3], sep = "_"),t))
	}
colnames(results_matrix) = c("mutations","effect","gene","line","unique_name","individual")
results_matrix[,5] = as.character(results_matrix[,5])

unique_mutations = unique(sort(results_matrix[,5]))
#unique_mutations = unique_mutations[1:10000]
effect = matrix(0, nrow = length(unique_mutations),ncol = length(unique(sort(results_matrix[,4])))) #effect of each mutation
colnames(effect) = unique(sort(results_matrix[,4]))
for(u in 1:length(unique_mutations))
	{	
		temp = results_matrix[results_matrix[,5] == unique_mutations[u],]
		
		effect[u,temp[,6]] = temp[,2]
		if(u %% 1000 == 0) print(paste(u,"of",length(unique_mutations),"Time is:", Sys.time()))
	}

#unique_mutations = unique_mutations[1:10000]
### get the frequency of each mutations in each of the 4 categories. See if they are private to that category
results_frequency = matrix(0,nrow = length(unique_mutations), ncol = 4) 
results_private = matrix(0,nrow = nrow(effect), ncol = ncol(effect)); colnames(results_private) = colnames(effect) # matrix telling you if the mutation is private in the land / eli or wild individuals.

v_landeli = c(1:length(type))[regexpr("land|elite", unique(sort(results_matrix[,4]))) > 0]
v_wild = c(1:length(type))[regexpr("wild", unique(sort(results_matrix[,4]))) > 0]

for(u in 1:length(unique_mutations))
	{
		landeli = effect[u,regexpr("landrace|elite",colnames(effect)) > 0]
		land = effect[u,regexpr("landrace",colnames(effect)) > 0]
		elite = effect[u,regexpr("elite",colnames(effect)) > 0]
		wild = effect[u,regexpr("wild",colnames(effect)) > 0]
		weed = effect[u,regexpr("weed",colnames(effect)) > 0]

		results_frequency[u,1] = length(land[land < -2.5]) / length(land) #proportion of individuals which have a deleterious mutations....
		results_frequency[u,2] = length(elite[elite < -2.5]) / length(elite)
		results_frequency[u,3] = length(wild[wild < -2.5]) / length(wild)
		results_frequency[u,4] = length(weed[weed < -2.5]) / length(weed)

		if( (length(landeli[landeli != 0]) > 0) &   (length(wild[wild != 0]) == 0) )  results_private[u,v_landeli] =  1  # private in the landrace/elite
		if( (length(wild[wild != 0]) > 0) & (length(landeli[landeli != 0]) == 0) )  results_private[u,v_wild] =  1  # private in the wild
		
		if(u %% 10000 == 0) print(paste(u,"of",length(unique_mutations),"Time is:", Sys.time()))
	}
########################
### SUMMARIZE THE RESULTS ###
########################

results_private_deleterious = matrix(0,nrow = length(type), ncol = 12)
for(t in 1:length(type))
{
effect_no0 = effect[effect[,t] != 0,t]
private_no0 = results_private[effect[,t] != 0,t]
private_effect = effect_no0[private_no0 ==1]

results_private_deleterious[t,2] = length(private_effect)  #total number of privates
results_private_deleterious[t,3] = length(private_effect[private_effect < -2.5])  #total number of privates deleterious
results_private_deleterious[t,4] = length(private_effect[private_effect < -4])  #total number of privates deleterious
results_private_deleterious[t,5] = length(private_effect[private_effect < -6])  #total number of privates deleterious
results_private_deleterious[t,6] = length(private_effect[private_effect < -8])  #total number of privates deleterious
}

results_private_deleterious[,7] = round(c(results_private_deleterious[,3] / results_private_deleterious[,2]),5)
results_private_deleterious[,8] = round(c(results_private_deleterious[,4] / results_private_deleterious[,2]),5)
results_private_deleterious[,9] = round(c(results_private_deleterious[,5] / results_private_deleterious[,2]),5)
results_private_deleterious[,10] = round(c(results_private_deleterious[,6] / results_private_deleterious[,2]),5)

results_private_deleterious[,1] =  individuals[,5]

results_private_deleterious[,11] = c(rep("burlywood4",10),rep("gold",1),rep("firebrick1",9),rep("dodgerblue3",9),rep("gold",4),rep("dodgerblue4",7))
colnames(results_private_deleterious) = c("accessions","non-synonymous_total","ns_del-2.5","ns_del-4","ns-del-6","ns-del-8","prop_del-2.5","prop_del-4","prop_del-6","prop_del-8","color","chisqtest")

write.table(results_private_deleterious,"results/new_results_private_deleterious",row.names = F, col.names = T, quote = F)
} else results_private_deleterious = as.matrix(read.table("results/new_results_private_deleterious", header = T))  #file already exists











