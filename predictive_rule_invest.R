# to investigate result of predictive rule
final.t = get(load(file = paste0(path_results, "attribution/Final.Table.RData")))
prob.t = get(load(file = paste0(path_results, "attribution/Post.Prob.List.RData")))
prob.tn = get(load(file = paste0(path_results, "attribution/Post.Prob.Listn.RData")))
tot = get(load(file = paste0(path_results,"attribution/stats_test_real_data.RData")))

# plot t-values of each configuration 
config1 = final.t[which(final.t$pred.y==1),]
dat.p = config1[,c(5:9)]
res1 = reshape2::melt(dat.p)
res1$value = abs(res1$value)
ggplot(res1, aes(x = variable, y = value))+ geom_boxplot()+geom_hline(yintercept = 1.96)

# All 0
d = rep(0,5)

a =final.t[which(apply(final.t[,c(13:17)], 1, function(x) all(x==d))),c(5:9,19)]
res1 = reshape2::melt(a, id = "pred.y")
res1=res1[which(res1$pred.y %in% c(1,8)),]
res1$pred.y = as.factor(res1$pred.y)
# res1$value = abs(res1$value)
ggplot(res1, aes(x = variable, y = value, fill = pred.y))+ geom_boxplot()+geom_hline(yintercept = 1.96)


