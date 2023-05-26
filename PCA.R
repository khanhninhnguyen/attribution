# to investigate result of predictive rule
offset = 0
GE = 0
number.pop =3
one = get(load(file = paste0(path_results, "attribution/predictive_rule/details/Final.Table", significance.level = 0.01, offset, GE, number.pop, ".RData")))
five = get(load(file = "/home/knguyen/Documents/PhD/paper/attribution/result/attribution/Final.Table.RData"))

all = left_join(one, five, by = c("main", "brp", "nearby")) %>%
  dplyr::select(c("main", "brp", "nearby",  "tGGp.x", "tGEp.x", "tEEp.x" ,"tGpEp.x", "tGpE.x", "pred.y.x", "pred.y.y"))


all$pred.y.x=as.factor(all$pred.y.x)
all$pred.y.y=as.factor(all$pred.y.y)

all.PCA=all[,4:10]
library(FactoMineR)
pca.test=PCA(all.PCA,scale.unit = TRUE,quali.sup = c(6,7),graph = FALSE)

library(factoextra)

fviz_pca_var(pca.test,axes = c(1,2))
fviz_pca_ind(pca.test,axes = c(1,2))


pca.test.8.22=
fviz_pca_biplot(pca.test,axes = c(1,2),habillage = "pred.y.x")
fviz_pca_biplot(pca.test,axes = c(1,2),habillage = "pred.y.y")


