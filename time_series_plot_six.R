# Plot time series oof the six test 
source(paste0(path_code_att,"FGLS.R"))
unicode_minus = function(x) sub('^-', '\U2212', format(x))

dat = get(load( file = paste0(path_results,"attribution/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
final.t = get(load(file = paste0(path_results, "attribution/Final.Table.RData")))
valid1 = get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))

valid2 = get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))
  
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)
name.i =  "fair.2017-10-04.clgo"
name.i1 =  "fair.2017-10-04.clgo"
ind.case = which(all.case %in% name.i  == TRUE)
brp = as.Date(substr(name.i, 6, 15))
station.name = substr(name.i, 1, 4)
# p1 GPS-ERA
# " check if there is restriction in data"
data.i = dat[[name.i]]

name.s = "gps.era"
datai = remove_na_2sides(data.i, name.s)
brp.ind = which(datai$date == brp)
begin.date = datai$date[1]
end.date = datai$date[length(datai$date)]
# read result FGLS 5 series 
station1 = get(load(file = paste0(path_results,"attribution/FGLS-GE/",name.i1, "fgls.RData")))

n1 = length(na.omit(datai[(1:brp.ind),name.s]))
n2 = length(na.omit(datai[-(1:brp.ind),name.s]))

datai$fit = station1[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station1[[name.s]]$t.table$Estimate[10] + station1[[name.s]]$t.table$Estimate[9]
text1 = paste0("Jump = ", round(station1[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station1[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station1[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR: 0.79, MA: 0.43",
               ", n1 = ", n1, ", n2 = ", n2)

# plot metadata 
meta.g = valid1[which(valid1$name == station.name),]
meta.g = meta.g[which(meta.g$known>begin.date & meta.g$known<end.date),]
list.meta = meta.g$known
if(nrow(meta.g)>0){
  datai$meta = NA
  datai$meta[which( datai$date %in% list.meta == TRUE)] = 4
}

p1 <- ggplot(data = datai, aes(x = date)) +
  theme_bw() + geom_line(aes(y = gps.era), col = "gray", lwd = 0.3)+ 
  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  geom_vline(xintercept = brp, lwd = 0.2)+ 
  labs(subtitle = text1)+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 E") + xlab("")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5), plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), legend.position = "none", plot.margin = unit(c(0, 0.5, 0, 0), "cm"))

if(nrow(meta.g)>0){
  p1<-p1+ geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
}
# GPS'-ERA' 

name.s = "gps1.era1"
datai = remove_na_2sides(data.i, name.s)
if(begin.date > datai$date[1]){
  begin.date =  datai$date[1]
}
if(end.date < datai$date[length(datai$date)]){
  end.date = datai$date[length(datai$date)]
}

brp.ind = which(datai$date == brp)
# read result FGLS 5 series 
station = get(load(file = paste0(path_results,"attribution/FGLS-full/",name.i, "fgls.RData")))

n1 = length(na.omit(datai[(1:brp.ind),name.s]))
n2 = length(na.omit(datai[-(1:brp.ind),name.s]))

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR : 0.84, MA: -0.58",
               ", n1 = ", n1, ", n2 = ", n2)

p2 <- ggplot(data = datai, aes(x = date, y = gps1.era1)) +
  theme_bw() + geom_line(col = "gray", lwd = 0.3)+
  geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G' \u2013 E'") +xlab("")+  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5), plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), legend.position = "none", plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
# four others

name.s = "gps.gps"
datai = remove_na_2sides(data.i, name.s)

# 
# delta = datai$date[length(datai$date)] - brp
# start.d = brp - delta
# datai = datai[which(datai$date>start.d),]

brp.ind = which(datai$date == brp)

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR: 0.61",
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])

meta.g = valid1[which(valid1$name == station.name),]
meta.g = meta.g[which(meta.g$known>begin.date & meta.g$known<end.date),]
list.meta = meta.g$known
if(nrow(meta.g)>0){
  datai$meta = NA
  datai$meta[which( datai$date %in% list.meta == TRUE)] = 4
}

p3 <- ggplot(data = datai, aes(x = date, y = gps.gps)) +
  theme_bw() + geom_line(col = "gray", lwd = 0.3)+
  geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 G'") +xlab("")+
  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5), plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), legend.position = "none", plot.margin = unit(c(0, 0.5, 0, 0), "cm"))

list.meta = meta.g$known
if(nrow(meta.g)>0){
  p3 <- p3 + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
}

name.s = "gps.era1"
datai = remove_na_2sides(data.i, name.s)
brp.ind = which(datai$date == brp)

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR: 0.45", 
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])

meta.g = valid1[which(valid1$name == station.name),]
meta.g = meta.g[which(meta.g$known>begin.date & meta.g$known<end.date),]
list.meta = meta.g$known
if(nrow(meta.g)>0){
  datai$meta = NA
  datai$meta[which( datai$date %in% list.meta == TRUE)] = 4
}

p4 <- ggplot(data = datai, aes(x = date, y = gps.era1)) +
  theme_bw() + geom_line(col = "gray", lwd = 0.3)+
  geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 E'") +xlab("")+
  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5),plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), legend.position = "none",plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
max.y <- ggplot_build(p4)$layout$panel_params[[1]]$y.range[2]
list.meta = meta.g$known
if(nrow(meta.g)>0){
  p4 <- p4 + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
}

name.s = "era.era"
datai = remove_na_2sides(data.i, name.s)
brp.ind = which(datai$date == brp)

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR: 0.31",
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])

p5 <- ggplot(data = datai, aes(x = date, y = era.era)) +
  theme_bw() + geom_line(col = "gray", lwd = 0.3)+
  geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("E \u2013 E'") +xlab("")+
  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5),plot.subtitle= element_text(size = 5),
        legend.title=element_blank(), legend.position = "none", plot.margin = unit(c(0, 0.5, 0, 0), "cm"))

name.s = "gps1.era"
datai = remove_na_2sides(data.i, name.s)
brp.ind = which(datai$date == brp)

datai$fit = station[[name.s]]$t.table$Estimate[10]
datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
               ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
               ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
               ", AR: 0.48",
               ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])

p6 <- ggplot(data = datai, aes(x = date, y = gps1.era)) +
  theme_bw() + geom_line(col = "gray", lwd = 0.3)+
  geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
  scale_y_continuous(breaks = seq(-7,5,1), limits = c(-4,4))+
  geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G' \u2013 E") +xlab("")+
  theme(axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),legend.text=element_text(size=4),
        axis.title = element_text(size = 5), legend.key.size = unit(0.3, "cm"), 
        plot.tag = element_text(size = 5),plot.subtitle = element_text(size = 5),
        legend.title=element_blank(), legend.position = "bottom", plot.margin = unit(c(0, 0, 0, 0.5), "cm"))

yearser = function(sd, ed){
   paste0(substr(seq(sd, ed, "years"), 1, 4)[-1], "-01-01")
}

list.seq = mapply(function(x, y) yearser(x, y), begin.date, end.date)

p1<- p1+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
p2<- p2+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
p3<- p3+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
p4<- p4+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
p5<- p5+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
p6<- p6+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                    breaks = function(x) as.Date(list.seq),
                    date_labels = "%Y" )
GE <- ggplotGrob(p1)
G.E. <- ggplotGrob(p2)
GG. <- ggplotGrob(p3)
GE. <- ggplotGrob(p4)
EE. <- ggplotGrob(p5)
G.E <- ggplotGrob(p6)
G.E.$widths <- GE$widths
GG.$widths <- GE$widths
GE.$widths <- GE$widths
EE.$widths <- GE$widths
G.E $widths <- GE$widths

# Arrange the two charts.

grid.newpage()
p = (grid.arrange(GE,  EE., GG.,  G.E.,  GE., G.E, nrow = 3))

ggsave(paste0(path_results,"attribution/time_sereis", name.i, ".jpg" ), plot = p, width = 14.4, height = 10, units = "cm", dpi = 1200)

