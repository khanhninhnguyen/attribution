# Plot for the paper functions 
# convert to coded table 
convert_coded <- function(x, thres){
  sapply(c(1:length(x)), function(i) ifelse((2*pnorm(-abs(x[i])))<thres, 1*sign(x[i]), 0)) 
}

Total.coded = data.frame(matrix(NA, ncol = 5, nrow = nrow(Total.res)))
for (i in c(1:nrow(Total.res))) {
  case.i = unlist(Total.res[i, c(paste0("t", list.name.test[2:6]))])
  Total.coded[i,] = convert_coded(case.i, 0.01)
}


# Test result for 6 series  -----------------------------------------------
Total.res <- read.table(paste0(path_main,"paper/attribution/result/attribution/", "FGLS_on_real_data_t.txt"), header = TRUE)
colnames(Total.res)[4:9] <- paste0("t", list.name.test)
text1 = "Distance < 50 km"
text2 = "Distance > 50 km"
colnames(Total.coded) = list.name.test[2:6]
data1 = Total.coded[which(Total.res$distance<50),]
# data1 = Total.coded
data.p = reshape2::melt(data1) 
data.p = rbind(data.p, data.frame(variable = rep(list.name.test[1], length(t.value.GE)), value = convert_coded(t.value.GE)))
data.p$c = 1
data.plot = aggregate(c~., data = data.p, sum)
data.plot$S= nrow(data1)
data.plot$S[which(data.plot$variable == list.name.test[1])] = length(t.value.GE)
data.plot$fre = data.plot$c/data.plot$S
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name1)
library(gtable)    
library(grid)
library(gridExtra) 
unicode_minus = function(x) sub('^-', '\U2212', format(x))
p1 <- ggplot(data.plot, aes(fill=value, y=c, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = c),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(a)", subtitle = text1) + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

data1 = Total.coded[which(reduced.list$distances>50),]
# data1 = Total.coded
data.p = reshape2::melt(data1) 
data.p = rbind(data.p, data.frame(variable = rep(list.name.test[1], length(t.value.GE)), value = convert_coded(t.value.GE)))
data.p$c = 1
data.plot = aggregate(c~., data = data.p, sum)
data.plot$S= nrow(data1)
data.plot$S[which(data.plot$variable == list.name.test[1])] = length(t.value.GE)
data.plot$fre = data.plot$c/data.plot$S
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name1)

p2 <- ggplot(data.plot, aes(fill=value, y=c, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = c),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(b)", subtitle = text2) + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)

gB$widths <- gA$widths
# Arrange the two charts.

grid.newpage()
p = (grid.arrange(gA, gB,nrow = 1))

ggsave(paste0(path_results,"attribution/pop_significance_level1.jpg" ), plot = p, width = 14.4, height = 5, units = "cm", dpi = 1200)

