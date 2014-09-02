


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



