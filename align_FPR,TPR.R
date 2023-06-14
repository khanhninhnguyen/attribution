# plot the FPR,TPR aligned
library(ggplot2)   
library(gtable)    
library(grid)
library(gridExtra) 
unicode_minus = function(x) sub('^-', '\U2212', format(x))
rate.ext <- function(Res, list.sample, off.set){
  tot.tpr = data.frame(matrix(NA, ncol = 4, nrow = length(list.sample)))
  for (j in c(1:length(list.sample))) {
    #FGLS
    fgls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$Estimate[1])
    mu.fgls = fgls.beta + off.set
    fgls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fgls$t.table$`Std. Error`[1])
    SD.jump.fgls = fgls.sd.beta 
    t.values.fgls = mu.fgls/SD.jump.fgls
    tpr.fgls = length(which(abs(t.values.fgls) > 1.96))
    
    # GLS 
    gls.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$Estimate[1])
    mu.gls = gls.beta + off.set
    gls.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$gls$t.table$`Std. Error`[1])
    SD.jump.gls = gls.sd.beta 
    t.values.gls = mu.gls/SD.jump.gls
    tpr.gls = length(which(abs(t.values.gls) > 1.96))
    
    #OLS
    ols.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.ols$Estimate[1])
    mu.ols = ols.beta + off.set
    ols.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.ols$`Std. Error`[1])
    SD.jump.ols = ols.sd.beta 
    t.values.ols = mu.ols/SD.jump.ols
    tpr.ols = length(which(abs(t.values.ols) > 1.96))
    
    #HAC
    hac.sd.beta = sapply(c(1:nb.sim), function(x) Res$total[[j]][[x]]$fit.hac$`Std. Error`[1])
    SD.jump.hac = hac.sd.beta 
    t.values.hac = mu.ols/SD.jump.hac
    tpr.hac = length(which(abs(t.values.hac) > 1.96))
    
    tot.tpr[j,] = c(tpr.fgls, tpr.gls, tpr.ols, tpr.hac)
  }
  colnames(tot.tpr) = c("FGLS", "GLS", "OLS", "OLS-HAC")
  out = tot.tpr[,c("GLS","FGLS", "OLS-HAC", "OLS")]
  return(out)
  
}

list.param.ar = seq(0,0.9,0.15)
list.param.sig = seq(0, 0.8, 0.2)

N = 400
nb.sim=1000
unit1 = "cm"
unit2 = "pt"
unit3 = "null"
tag.list = c("a)","b)","c)","d)","e)","f)")
# plot 1------------
hetero = 0
autocor = 1
y.axis = "FPR"
x.axis = "rho"
extrem = ""
# x.axis = "sig"
list.sample = list.param.ar
# list.sample = list.param.sig

noise.name = ifelse(autocor==0, "white", "ar1")
Res = get(load(file = paste0(path_results,"attribution/N1000n/",hetero,"auto",autocor, x.axis, y.axis,noise.name,extrem,"R.Data")))
rate.ext(Res, list.sample = list.sample, off.set = 0)

FPR = rate.ext(Res, list.sample = list.sample, off.set = 0)/nb.sim
TPR = rate.ext(Res, list.sample = list.sample, off.set = 0.356)/nb.sim
if(x.axis=="rho"){
  FPR$phi = list.param.ar
  TPR$phi = list.param.ar
  dat.fpr = reshape2::melt(FPR, id ="phi")
  dat.tpr = reshape2::melt(TPR, id ="phi")
  
}else{
  FPR$range = list.param.sig*100
  TPR$range = list.param.sig*100
  dat.fpr = reshape2::melt(FPR, id ="range")
  dat.tpr = reshape2::melt(TPR, id ="range")
}
# plot FPR 
x.axis.name = ifelse(x.axis=="rho", "phi","range")

y.axis.name = "False Positive Rate"
thres = 0.05
if(x.axis.name =="phi"){
  break1 = list.sample
}else{
  break1 = list.sample*100
}
x.axis.name1 = ifelse(x.axis=="rho", expression(phi) ,"range(%)")

legend.1 = "none"
tag.pos = "topleft"
p1 <- eval(parse(
  text=paste0("ggplot(dat.fpr, aes(x =", x.axis.name, ",y = value, col = variable))+
 geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[1])  +
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = 2, r = 2, l = 2),  plot.tag.position = tag.pos,
              legend.box.spacing = unit(0,unit2 ))")))

y.axis.name = "True Positive Rate"
thres = 0.95
legend2 = "none"
p2 <- eval(parse(
  text=paste0("ggplot(dat.tpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[2]) + 
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = 2, r = 2, l = 2), plot.tag.position = tag.pos,legend.box.spacing = unit(0,unit2 ))")))

#plot 2----------

hetero = 1
autocor = 1
y.axis = "FPR"
x.axis = "rho"
extrem = ""
x.axis = "sig"
# list.sample = list.param.ar
list.sample = list.param.sig

noise.name = ifelse(autocor==0, "white", "ar1")
Res = get(load(file = paste0(path_results,"attribution/N1000n/",hetero,"auto",autocor, x.axis, y.axis,noise.name,extrem,"R.Data")))
rate.ext(Res, list.sample = list.sample, off.set = 0)

FPR = rate.ext(Res, list.sample = list.sample, off.set = 0)/nb.sim
TPR = rate.ext(Res, list.sample = list.sample, off.set = 0.356)/nb.sim
if(x.axis=="rho"){
  FPR$phi = list.param.ar
  TPR$phi = list.param.ar
  dat.fpr = reshape2::melt(FPR, id ="phi")
  dat.tpr = reshape2::melt(TPR, id ="phi")
  
}else{
  FPR$range = list.param.sig*100
  TPR$range = list.param.sig*100
  dat.fpr = reshape2::melt(FPR, id ="range")
  dat.tpr = reshape2::melt(TPR, id ="range")
}
# plot FPR 

x.axis.name = ifelse(x.axis=="rho", "phi","range")

y.axis.name = "False Positive Rate"
thres = 0.05
if(x.axis.name =="phi"){
  break1 = list.sample
}else{
  break1 = list.sample*100
}
x.axis.name1 = ifelse(x.axis=="rho", expression(phi) ,"range(%)")

p3 <- eval(parse(
  text=paste0("ggplot(dat.fpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[3]) + 
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = 2, r = 2, l = 2), plot.tag.position = tag.pos,legend.box.spacing = unit(0,unit2 ))")))

y.axis.name = "True Positive Rate"
thres = 0.95
p4 <- eval(parse(
  text=paste0("ggplot(dat.tpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[4]) + 
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = 2, r = 2, l = 2),plot.tag.position = tag.pos, legend.box.spacing = unit(0,unit2 ))")))


# plot 3------
hetero = 1
autocor = 1
y.axis = "FPR"
x.axis = "rho"
extrem = 1
# x.axis = "sig"
list.sample = list.param.ar
# list.sample = list.param.sig

noise.name = ifelse(autocor==0, "white", "ar1")
Res = get(load(file = paste0(path_results,"attribution/N1000n/",hetero,"auto",autocor, x.axis, y.axis,noise.name,extrem,"R.Data")))
rate.ext(Res, list.sample = list.sample, off.set = 0)

FPR = rate.ext(Res, list.sample = list.sample, off.set = 0)/nb.sim
TPR = rate.ext(Res, list.sample = list.sample, off.set = 0.356)/nb.sim
if(x.axis=="rho"){
  FPR$phi = list.param.ar
  TPR$phi = list.param.ar
  dat.fpr = reshape2::melt(FPR, id ="phi")
  dat.tpr = reshape2::melt(TPR, id ="phi")
  
}else{
  FPR$range = list.param.sig*100
  TPR$range = list.param.sig*100
  dat.fpr = reshape2::melt(FPR, id ="range")
  dat.tpr = reshape2::melt(TPR, id ="range")
}
# plot FPR 

x.axis.name = ifelse(x.axis=="rho", "phi","range")

y.axis.name = "False Positive Rate"
thres = 0.05
if(x.axis.name =="phi"){
  break1 = list.sample
}else{
  break1 = list.sample*100
}
x.axis.name1 = ifelse(x.axis=="rho", expression(phi) ,"range(%)")

p5 <- eval(parse(
  text=paste0("ggplot(dat.fpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[5]) + 
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
   theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = 0, r = 2, l = 2), plot.tag.position = tag.pos,legend.box.spacing = unit(0,unit2 ))")))

y.axis.name = "True Positive Rate"
thres = 0.95
p6 <- eval(parse(
  text=paste0("ggplot(dat.tpr, aes(x =", x.axis.name, ",y = value, col = variable))+
  geom_point(size=0.3) + geom_line(lwd = 0.3) +theme_bw() +
  ylab(y.axis.name) + labs(tag = tag.list[6]) + 
  xlab(x.axis.name1) + geom_hline(yintercept =", thres,", lwd = 0.1)+
  scale_y_continuous(breaks=seq(0, 1, 0.1), limits =c(0,1))+
  scale_x_continuous(breaks=break1)+
  theme(axis.text.x = element_text(size = 6), axis.text.y = element_text(size = 6),legend.text=element_text(size=6),
        axis.title = element_text(size = 7), legend.key.size = unit(1, unit1), plot.tag = element_text(size = 7),
        legend.title=element_blank(),plot.margin = margin(t = 2, b = -3, r = 2, l = 2), plot.tag.position = tag.pos,legend.box.spacing = unit(0,unit2 ))")))

# 
# gA <- ggplotGrob(p1)
# gB <- ggplotGrob(p2)
# gC <- ggplotGrob(p3)
# gD <- ggplotGrob(p4)
# gE <- ggplotGrob(p5)
# gF <- ggplotGrob(p6)
# 
# gB$widths <- gA$widths
# gC$widths <- gA$widths
# gD$widths <- gA$widths
# gE$widths <- gA$widths
# gF$widths <- gA$widths
# 
# # Arrange the two charts.
# 
# grid.newpage()
# 
# arrangeGrob(p1, theme(legend.position="none"), ncol=2)
# p = (grid.arrange(gA, gB, gC, gD, gE, gF,nrow = 3))

# grid.arrange(grobs = list(p1, p2), ncol=2,top="Main Title", theme(legend.position="none"))


grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {

  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  # print(lheight)
  lwidth <- sum(legend$width)
  # print(lwidth)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)

  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  combined

}

p = grid_arrange_shared_legend(p1, p2, p3, p4,p5, p6, nrow = 3, ncol = 2)
print(p)
ggsave(paste0(path_results,"attribution/simulation_FGLS_HAC1.jpg" ), plot = p, width = 14.4, height = 12, units = "cm", dpi = 1200)






