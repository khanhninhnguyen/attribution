# correct the t values according to the E-E' 
convert_coded <- function(x){
  sapply(c(1:length(x)), function(i) ifelse(abs(x[i])>1.96, 1*sign(x[i]), 0)) 
}
check_contradict <- function(y, table.selected){
  names(y) = NULL
  colnames(table.selected) = NULL
  res = sapply(c(1:38), function(x) identical(unlist(table.selected[x,]), y))
  ind.o = which(res==TRUE)
  out = ifelse(length(ind.o)>0, ind.o, 0)
  return(out)
}

# Correct the horizontal distance based on E-E' ---------------------------

Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))

## Only consider the non-collocated pairs
list.non.collocated = list.name.test[c(-1,-5)]
jump.old = Total.res[,paste0("jump", list.non.collocated )]
sd.old = Total.res[,paste0("jump", list.non.collocated)]/Total.res[,paste0("t", list.non.collocated )]
jump.new = jump.old - Total.res$`jumpE-E'`
jump.new$`jumpG'-E` = jump.new$`jumpG'-E`+2*Total.res$`jumpE-E'`
t.new = jump.new/sd.old
Total.res1 = Total.res
Total.res1[,paste0("jump", list.non.collocated )] = jump.new
Total.res1[,paste0("t", list.non.collocated )] = t.new 
save(Total.res1, file = paste0(path_results, "attribution0/stats_test_real_data_corrected_dist.RData"))

## plot histogram to compare
dat.p = reshape2:: melt(data.frame(old = Total.res$`tG'-E`, new = t.new$`jumpG'-E`))
dat.p %>%
  ggplot( aes(x=value, fill=variable)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 
## check the contradicted cases 
trunc.table = get(load(file = paste0(path_results, "attribution0/truncated.table.RData")))

Total.coded.old = as.data.frame(sapply(c(1:6), function(x) convert_coded(Total.res[,paste0("t", list.name.test[x])])))
colnames(Total.coded.old) = paste0("old", list.name.test)
contra.old = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded.old[x,c(1:6)]), trunc.table))

Total.coded.new = as.data.frame(sapply(c(1:6), function(x) convert_coded(Total.res1[,paste0("t", list.name.test[x])])))
colnames(Total.coded.new) = paste0("new", list.name.test)
contra.new = sapply(c(1:nrow(Total.coded)), function(x) check_contradict(unlist(Total.coded.new[x,c(1:6)]), trunc.table))
Total.res1$config = contra.new
Total.res1[,list.name.test] = Total.coded.new
all.t = cbind(Total.coded.old, Total.coded.new) %>% dplyr:: mutate(con1 = contra.old, new = contra.new)

dat.p = data.frame(table(contra.old)) %>% rename("config" = "contra.old") %>%
  mutate(name = rep("old",length(table(contra.old)))) %>%
  rbind( data.frame(table(contra.new)) %>% rename("config" = "contra.new") %>%
           mutate(name = rep("new",length(table(contra.new)))))

ggplot(data = dat.p, aes(x = config, y = Freq, fill = name))+
  theme_bw()+
  geom_bar(stat="identity", position=position_dodge(), width = 0.5)
 





