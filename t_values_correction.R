# correct the t values according to the E-E' 
convert_coded <- function(x, cri){
  sapply(c(1:length(x)), function(i) ifelse(abs(x[i])>cri, 1*sign(x[i]), 0)) 
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
# save file for prediction 
station = substr(Total.res1$station, 1, 4)
brp = substr(Total.res1$station, 6, 15)
nearby = substr(Total.res1$station, 17, 20)
data.save = data.frame(main = station, brp = brp, nearby = nearby)
data.save = cbind(data.save, Total.res1[,c(7:12)], Total.res1[,c(22,23,20)])
write.table(data.save, file = paste0(path_results, "attribution/predictive_rule/stats_test_real_data_corrected_dist.txt"))

## plot histogram to compare
dat.p = reshape2:: melt(data.frame(old = Total.res$`tG'-E`, new = t.new$`jumpG'-E`))
dat.p %>%
  ggplot( aes(x=value, fill=variable)) + theme_bw()+
  geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity', bins = 100) +
  scale_fill_manual(values=c("#69b3a2", "#404080")) 
## check the contradicted cases 
trunc.table = get(load(file = paste0(path_results, "attribution0/truncated.table.RData")))

Total.coded.old = as.data.frame(sapply(c(1:6), function(x) convert_coded(Total.res[,paste0("t", list.name.test[x])], cri = 1.96)))
colnames(Total.coded.old) = paste0("old", list.name.test)
contra.old = sapply(c(1:nrow(Total.coded.old)), function(x) check_contradict(unlist(Total.coded.old[x,c(1:6)]), trunc.table))

Total.coded.new = as.data.frame(sapply(c(1:6), function(x) convert_coded(Total.res1[,paste0("t", list.name.test[x])], cri = 1.96)))
colnames(Total.coded.new) = paste0("new", list.name.test)
contra.new = sapply(c(1:nrow(Total.coded.new)), function(x) check_contradict(unlist(Total.coded.new[x,c(1:6)]), trunc.table))
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
 
table(all.t$con1)
table(all.t$new)


# fdr correction ----------------------------------------------------------

## without bias correction -------------------------------------------------
Keep.pval = as.data.frame(
  sapply(c(1:6), function(x) round(pnorm(-abs(unlist(Total.res[paste0("t", list.name.test[x])])), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 5)))
pvalA <- purrr::map(1:dim(Keep.pval)[1],~p.adjust(Keep.pval[.x,], method = "fdr")) %>% Reduce(c,.) %>% matrix(ncol=6,nrow=dim(Keep.pval)[1],byrow=TRUE)
colnames(pvalA) <- list.name.test
signif <- ifelse(pvalA<0.05,1,0)*sign(as.matrix(Total.res[paste0("t", list.name.test)])) %>% as.data.frame()
contra.1 = sapply(c(1:nrow(signif)), function(x) check_contradict(unlist(signif [x,c(1:6)]), trunc.table))
### bar plot of configuration
dat.p = data.frame(table(contra.old)) %>% rename("config" = "contra.old") %>%
  mutate(name = rep("raw",length(table(contra.old)))) %>%
  rbind( data.frame(table(contra.1)) %>% rename("config" = "contra.1") %>%
           mutate(name = rep("fdr",length(table(contra.1)))))

ggplot(data = dat.p, aes(x = config, y = Freq, fill = name))+
  theme_bw()+
  geom_bar(stat="identity", position=position_dodge(), width = 0.5)
### investigate why different 
two.config = data.frame(raw = contra.old, fdr = contra.1, name = Total.res$station)

Keep.pval = as.data.frame(
  sapply(c(1:6), function(x) round(pnorm(-abs(unlist(Total.res1[paste0("t", list.name.test[x])])), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 4)))
pvalA <- purrr::map(1:dim(Keep.pval)[1],~p.adjust(Keep.pval[.x,], method = "fdr")) %>% Reduce(c,.) %>% matrix(ncol=6,nrow=dim(Keep.pval)[1],byrow=TRUE)
colnames(pvalA) <- list.name.test
signif <- ifelse(pvalA<0.05,1,0)*sign(as.matrix(Total.res1[paste0("t", list.name.test)])) %>% as.data.frame()
contra.1 = sapply(c(1:nrow(signif)), function(x) check_contradict(unlist(signif [x,c(1:6)]), trunc.table))

dat.p = data.frame(table(contra.new)) %>% rename("config" = "contra.new") %>%
  mutate(name = rep("bias.cor",length(table(contra.new)))) %>%
  rbind( data.frame(table(contra.1)) %>% rename("config" = "contra.1") %>%
           mutate(name = rep("fdr",length(table(contra.1)))))

ggplot(data = dat.p, aes(x = config, y = Freq, fill = name))+
  theme_bw()+
  geom_bar(stat="identity", position=position_dodge(), width = 0.5)


# limit at 1% -------------------------------------------------------------

Total.coded.new1 = as.data.frame(sapply(c(1:6), function(x) convert_coded(Total.res[,paste0("t", list.name.test[x])], cri = 2.58)))
colnames(Total.coded.new1) = paste0("new", list.name.test)
contra.new1 = sapply(c(1:nrow(Total.coded.new1)), function(x) check_contradict(unlist(Total.coded.new1[x,c(1:6)]), trunc.table))

dat.p = data.frame(table(contra.old)) %>% rename("config" = "contra.old") %>%
  mutate(name = rep("5%",length(table(contra.old)))) %>%
  rbind( data.frame(table(contra.new1)) %>% rename("config" = "contra.new1") %>%
           mutate(name = rep("1%",length(table(contra.new1)))))

ggplot(data = dat.p, aes(x = config, y = Freq, fill = name))+
  theme_bw()+
  geom_bar(stat="identity", position=position_dodge(), width = 0.5)

## with FDR+bias correction
Keep.pval = as.data.frame(
  sapply(c(1:6), function(x) round(pnorm(-abs(unlist(Total.res1[paste0("t", list.name.test[x])])), mean = 0, sd = 1, lower.tail = TRUE)*2, digits = 5)))
pvalA <- purrr::map(1:dim(Keep.pval)[1],~p.adjust(Keep.pval[.x,], method = "fdr")) %>% Reduce(c,.) %>% matrix(ncol=6,nrow=dim(Keep.pval)[1],byrow=TRUE)
colnames(pvalA) <- list.name.test

### convert to t-value : keep the large value from the original 
t.fdr = Total.res1[,c(7:12)]
for(i in c(1:nrow(pvalA))){
  p.i = pvalA[i,]
  for(j in c(1:6)){
    if(p.i[j] != 0){
      sign.ij = sign(t.fdr[i,j])
      t.c = sign.ij * abs(qt(p.i[j]/2, df = 1000))
      t.fdr[i,j] = t.c
    }
  }
}

# save file for prediction 
Total.res1[,c(7:12)] <- t.fdr
station = substr(Total.res1$station, 1, 4)
brp = substr(Total.res1$station, 6, 15)
nearby = substr(Total.res1$station, 17, 20)
data.save = data.frame(main = station, brp = brp, nearby = nearby)
data.save = cbind(data.save, Total.res1[,c(7:12)], Total.res1[,c(22,23,20)])
write.table(data.save, file = paste0(path_results, "attribution/predictive_rule/stats_test_real_data_corrected_dist_fdr.txt"))

### 
signif <- ifelse(pvalA<0.05,1,0)*sign(as.matrix(Total.res[paste0("t", list.name.test)])) %>% as.data.frame()
contra.1 = sapply(c(1:nrow(signif)), function(x) check_contradict(unlist(signif [x,c(1:6)]), trunc.table))


