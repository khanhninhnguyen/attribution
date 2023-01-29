
final.t = get(load(file = paste0(path_results, "attribution/Final.Table.RData")))
prob.t = get(load(file = paste0(path_results, "attribution/Post.Prob.List.RData")))
prob.tn = get(load(file = paste0(path_results, "attribution/Post.Prob.Listn.RData")))

valid = get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
last.name = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[1])
last.brp = sapply(c(1:114), function(x) prob.t[[x]]$MainBreak[2])
last.pre = sapply(c(1:114), function(x) prob.t[[x]]$Config.Pred.Post)
valid.list = sapply(c(1:114), function(x) valid$valid[which(valid$name == last.name[x] & valid$detected == last.brp[x])])

res = data.frame(name = last.name, brp = last.brp, pred = last.pre, valid = valid.list)