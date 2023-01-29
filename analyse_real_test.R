
tot = get(load(file = paste0(path_results,"attribution/stats_test_real_data.RData")))
t.table = cbind(reduced.list[,c(1:3)], tot[,c(7:12, 21,22,20)])

write.table(t.table,
            file = paste0(path_results,"attribution/FGLS_on_real_data_t.txt"),
            row.names = FALSE, sep = "\t", quote = FALSE)
# to analyse results of test ---n

tot$n = tot$n1+tot$n2

dat.a = data.frame(matrix(NA, ncol = 3, nrow = 0))
for (i in c(2,3,4,6)) {
  jump = tot[,c(paste0("jump", list.name.test[i]))]
  sig = tot[,c(list.name.test[i])]
  jump.sig = jump[which(sig != 0)]
  n1 = tot$n[which(sig != 0)]
  namei = tot$station[which(sig != 0)]
  dat.i = data.frame( jump.s = jump.sig, n = n1, name = namei)
  dat.a = rbind(dat.a, dat.i)
}

a = dat.a[which(dat.a$n<600),]
l1 = c(seq(500,1000,100))
sapply(c(1:length(l1)), function(x) {
  a = dat.a[which(dat.a$n<l1[x]),]
  mean(abs(a$jump.s))
})

# plot CDF jump -----------------------------------------------------------

jumps = tot[,c(1:6)]
names(jumps) = list.name.test
# jump1 = jumps
jump1 = jumps[which(tot$distance>50),]

dat.p = reshape2::melt(jump1)
dat.p$value = abs(dat.p$value)
dat.p$variable = factor(dat.p$variable,  levels = reoder.list.name1)

library(RColorBrewer)
p <- ggplot(dat.p, aes(x = value, col = variable ))+ theme_bw()+
  stat_ecdf(lwd = 0.3)+
  scale_x_continuous(breaks = seq(0, 1.5, 0.3), limits = c(0,1.5))+
  # geom_hline(yintercept = 0.5, size = 0.3) +
  scale_color_manual(values = brewer.pal(n = 6, name = 'Dark2'))+
  labs(y = "CDF", x = "Jump", linetype = "")+
  theme(axis.text = element_text(size = 5),legend.text=element_text(size=4.5),
        axis.title = element_text(size=5), legend.key.size = unit(0.2, "cm"),
        legend.box.spacing = unit(0, "pt"), legend.title= element_blank())
# print(p)
# dev.off()
ggsave(paste0(path_results,"attribution/jumps1.jpg" ), plot = p, width = 8.8, height = 6, units = "cm", dpi = 1200)

j1 = abs(jumps[which(tot$distance>50),])
j2 = abs(jumps[which(tot$distance<50),])
d1 = reshape2::melt(j1)
d2 = reshape2::melt(j2)

d = rbind(d1,d2)
d$distance = c(rep(">50",nrow(d1)),rep("<50", nrow(d2)))
d$distance[which(d$variable=="G-E")] = "<50"
d$variable = factor(d$variable,  levels = reoder.list.name1)

p <- ggplot(data = d, aes(x = variable, y = value, col = distance))+theme_bw()+
  geom_boxplot(lwd=0.3,outlier.size = 0.2) + 
  labs(y = "Jump", x = "", linetype = "")+
  theme(axis.text = element_text(size =6),legend.text=element_text(size=5.5),
        axis.title = element_text(size=6), legend.key.size = unit(0.2, "cm"),
        legend.box.spacing = unit(0, "pt"), legend.title=element_text(size=6))
# print(p)
# dev.off()
ggsave(paste0(path_results,"attribution/jumps_boxplot.jpg" ), plot = p, width = 8, height = 6, units = "cm", dpi = 600)


