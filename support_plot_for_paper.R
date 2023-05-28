# Plot for the paper functions 
# convert to coded table 
convert_coded <- function(x, thres){
  sapply(c(1:length(x)), function(i) ifelse((2*pnorm(-abs(x[i])))<thres, 1*sign(x[i]), 0)) 
}

# Test result for 6 series  -----------------------------------------------
Total.res <- read.table(paste0(path_main,"paper/attribution/result/attribution/", "FGLS_on_real_data_t.txt"), header = TRUE)
colnames(Total.res)[4:9] <- paste0("t", list.name.test)

Total.coded = data.frame(matrix(NA, ncol = 6, nrow = nrow(Total.res)))
for (i in c(1:nrow(Total.res))) {
  case.i = unlist(Total.res[i, c(paste0("t", list.name.test[1:6]))])
  Total.coded[i,] = convert_coded(case.i, 0.01)
}

Total.coded = Total.coded %>% 
  `colnames<-`(list.name.test) %>% 
  mutate(distance = sapply(c(1:nrow(Total.res)), function(x) ifelse(Total.res$distance[x]<50, 1,2))) 

count.total =  Total.coded[,c(1:6)] %>% 
  reshape2::melt() %>% 
  mutate(count = 1) %>% 
  aggregate(count~.,data = ., sum)

count.GE = count.total[which(count.total$variable == "G-E"),]
# Fig 1:
text1 = "Distance < 50 km"

data1 = Total.coded[which(Total.res$distance<50),] 

data.plot = data1[,c(1:6)] %>% 
  reshape2::melt() %>% 
  mutate(count = 1) %>% 
  aggregate(count~.,data = ., sum)

data.plot[which(data.plot$variable == "G-E"),] <- count.GE
data.plot$value = as.factor(data.plot$value)
data.plot$variable = factor(data.plot$variable,  levels = reoder.list.name1)

unicode_minus = function(x) sub('^-', '\U2212', format(x))
p1 <- ggplot(data.plot, aes(fill=value, y=count, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = count),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(c)", subtitle = text1) + 
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 6),plot.subtitle = element_text(size = 6),
        legend.title=element_blank(), legend.box.spacing = unit(0, "pt"), plot.margin = rep(unit(0,"null"),4))

# Fig 2:
text2 = "Distance > 50 km"

data2 = Total.coded[which(Total.res$distance>50),] 

data.plot1 = data2[,c(1:6)] %>% 
  reshape2::melt() %>% 
  mutate(count = 1) %>% 
  aggregate(count~.,data = ., sum)

data.plot1[which(data.plot1$variable == "G-E"),] <- count.GE
data.plot1$value = as.factor(data.plot1$value)
data.plot1$variable = factor(data.plot1$variable,  levels = reoder.list.name1)

unicode_minus = function(x) sub('^-', '\U2212', format(x))
p2 <- ggplot(data.plot1, aes(fill=value, y=count, x=variable)) + 
  geom_bar(position="stack", stat="identity")+theme_bw()+
  geom_text(aes(label = count),
            colour = "black",  size=1.8,
            position = position_stack(vjust = 0.5)) +
  labs(x = NULL, y ="Count", tag = "(d)", subtitle = text2) + 
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

ggsave(paste0(path_results,"attribution/pop_significance_level2.jpg" ), plot = p, width = 14.4, height = 5, units = "cm", dpi = 1200)



