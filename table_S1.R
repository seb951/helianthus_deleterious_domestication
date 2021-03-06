setwd("~/Documents/domestication_load")#set up working directory 

system("ls -1 sequences/*_1.fq >list")

list = read.table("list"); list = cbind(list,0)

for(i in 1:nrow(list))
	{
		system(paste("wc", list[i,1],">count"))
		count = read.table("count")
		list[i,2] = count[1,1] / 4
	}
list[,1] = gsub("sequences/","",list[,1])
list[,1] = gsub("_1","",list[,1])
list = list[-c(1:6),];list = list[-c(22),]

individuals = read.table("reference/annuus.txt", header = T, stringsAsFactors = F)
individuals =  individuals[c(1:(nrow(individuals)-6)),]  #remove the test run individuals ### individuals = individuals[1:6,] were for the test run. Nuke 'em !!!!
individuals[,5] = gsub("eli","elite",individuals[,5])
individuals[,5] = gsub("land","landrace",individuals[,5])
individuals[,5] = gsub("wild_tex","wild",individuals[,5])
individuals[,5] = gsub("weed","weed",individuals[,5])


system("rm count list")

write.table(list,"results/table_S1",row.names = F, col.names = F, quote = F)


###
setwd("~/Documents/domestication_load/sequences") ### setup the working directory in the right format

system("ls -1 *_1.fq >list_tta")
list_tta = read.table("list_tta", stringsAsFactors = F)

#####################
###stats alignment function###
####################
wc.out = NULL
        for(i in 1:nrow(list_tta)) {

        seq_1 = list_tta[i,1]
        system(paste("wc -l ",seq_1," >wc_seq",sep = "")) #how many raw sequences

wc = read.table("wc_seq", header = F, stringsAsFactors = F)[1,1]
wc.out = c(wc.out, wc/4)

                                                        }
system("rm list_tta")
wc.out = signif(wc.out / 1000000,4)

write.table(cbind(list_tta,wc.out),"../results/number_of_raw_sequences", row.names = F, col.names = T)


