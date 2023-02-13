# study the contradiction
#  create the truncated table  ------
G = c(1,0,-1)
E = c(1,0,-1)
G. = c(1,0,-1)
E. =c(1,0,-1)
a = expand.grid(G, E, G., E.)
colnames(a) = c("E.", "G.", "E", "G")
a = a[,c("G", "E", "G.", "E.")]
Truth = a[-which(a$G==1 & a$E==1 | a$G==0 & a$E==0 | a$G==-1 & a$E==-1),]

logic.table = data.frame(gps.era = Truth$G- Truth$E, 
                         gps.gps = Truth$G- Truth$G.,
                         gps.era1 = Truth$G- Truth$E.,
                         era.era = Truth$E- Truth$E.,
                         gps1.era1 = Truth$G.- Truth$E.,
                         gps1.era = Truth$G.- Truth$E)
truncate.table = logic.table
truncate.table.m = sapply(c(1:6), function(x) {
  a = truncate.table[,x]
  a[which(a==2)]=1
  a[which(a==-2)]=-1
  return(a)
})
truncate.table.mo = unique(truncate.table.m)
truncate.table.5 = as.data.frame(unique(truncate.table.mo[,-1]))
colnames(truncate.table.5) = colnames(logic.table)[2:6]
save(truncate.table.5, file = paste0(path_results, "attribution/truncated.table.RData"))

# check the contradiction function ---------------
check_contradict <- function(y, table.selected){
  names(y) = NULL
  colnames(table.selected) = NULL
  res = sapply(c(1:36), function(x) identical(unlist(table.selected[x,]), y))
  ind.o = which(res==TRUE)
  out = ifelse(length(ind.o)>0, ind.o, 0)
  return(out)
}

# convert to coded table 
# read test result 
Total.res = get(load(file = paste0(path_results,"attribution0/stats_test_real_data.RData")))
# result for all G-E
list.name = substr(Total.res$station, 1, 15)
Total.res$brp = list.name
for (i in c(1:length(unique(list.name)))) {
  ind.i = which(Total.res$brp == unique(list.name)[i])
  t.i = Total.res$`tG-E`[ind.i]
  Total.res$`tG-E`[ind.i]= rep(t.i, length(ind.i))
}



convert_coded <- function(x, significance.level, length.x){
  sapply(c(1:length(x)), function(i) ifelse( (2*pnorm(-abs(as.numeric(x[i])))) < significance.level, 1*sign(x[i]), 0)) 
}



#check contradiction for different level of significance 
trunc.table = get(load(file = paste0(path_results, "attribution0/truncated.table.RData")))
contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table))
table(unlist(contra))

sig.list = seq(0.00, 0.1,0.002)[-1]
nb.contradicted = rep(NA, length(sig.list))
for (j in 1:length(sig.list)) {
  Total.coded = data.frame(matrix(NA, ncol = 5, nrow = nrow(Total.res)))
  for (i in c(1:nrow(Total.res))) {
    case.i = Total.res[i, c(paste0("t", list.name.test[2:6]))]
    Total.coded[i,] = convert_coded(case.i, significance.level = sig.list[j])
  }
  contra = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded[x,]), trunc.table))
  nb.contradicted[j] = length(which(contra == 0))
  a = data.frame(table(contra))
  png(file = paste0(path_results,"attribution/pie", sig.list[j],".png"))
  pie(a$Freq,a$contra)
  dev.off()
}


plot(sig.list, nb.contradicted)


