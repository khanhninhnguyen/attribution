# function to plot time series 
source(paste0(path_code_att,"FGLS.R"))

dat = get(load( file = paste0(path_results,"attribution0/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
final.t = get(load(file = paste0(path_results, "attribution0/Final.Table.RData")))
rownames(final.t) = NULL
lengthlist = get(load(file = paste0(path_results, "attribution0/lengthlist.RData")))
final.t$n1 = lengthlist$X1
final.t$n2 = lengthlist$X2
valid1 = get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
valid2 = get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))
all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)
yearser = function(sd, ed){
  paste0(substr(seq(sd, ed, "years"), 1, 4)[-1], "-01-01")
}

plot_six <- function(name.case){
  
  # name of main and 5 others
  name.GG =  name.case
  station.name = substr(name.GG, 1, 4)
  brp =  as.Date(substr(name.GG, 6, 15),format="%Y-%m-%d")
  station.nearby =  substr(name.GG, 17, 20)
  ind.case = which(all.case %in% name.GG == TRUE)
  name.GE = all.case[which(final.t$main == station.name & final.t$brp == brp & is.na(final.t$tGE)==FALSE)]
  
  # plot for G-E------------------------
  data.i = dat[[name.GG]]
  name.s = "gps.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  
  # read result FGLS for G-E
  station1 = get(load(file = paste0(path_results,"attribution0/FGLS-GE/",name.GE, "fgls.RData")))
  n1GE = length(na.omit(datai[(1:brp.ind),name.s]))
  n2GE = length(na.omit(datai[-(1:brp.ind),name.s]))
  
  if((n1GE+n2GE) >(length(station1$gps.era$fit)+50)){
    ind.wtna = which(is.na(datai[,name.s]) == FALSE)
    brp.ind = which(ind.wtna == brp.ind)
    min.ind = ind.wtna[max(1,(brp.ind-1000), na.rm = TRUE)]
    max.ind = ind.wtna[min((brp.ind+999), length(ind.wtna), na.rm = TRUE)]
    datai = datai[c(min.ind:max.ind),]
    brp.ind = which(datai$date == brp)
    n1GE = length(na.omit(datai[c(min.ind:brp.ind),name.s]))-1
    n2GE = length(na.omit(datai[c(brp.ind:max.ind),name.s]))
  }
  
  begin.date = datai$date[1]
  end.date = datai$date[length(datai$date)]
  
  set.margin = list(lower = floor(min(data.i[,c(1:6)], na.rm = TRUE)), 
                    upper = ceiling(max(data.i[,c(1:6)], na.rm = TRUE)))
  datai$fit = station1[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station1[[name.s]]$t.table$Estimate[10] + station1[[name.s]]$t.table$Estimate[9]
  text1 = paste0(toupper(station.name), ", Jump = ", round(station1[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station1[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station1[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station1$gps.era$coef.arma$phi, digits = 2), ", MA: ", round(station1$gps.era$coef.arma$theta, digits = 2),
                 ", n1 = ", n1GE, ", n2 = ", n2GE)
  
  # plot metadata 
  meta.g = valid1[which(valid1$name == station.name),]
  meta.g = meta.g[which(meta.g$known>begin.date & meta.g$known<end.date),]
  list.meta = meta.g$known
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  
  pGE <- ggplot(data = datai, aes(x = date)) +
    theme_bw() + 
    geom_line(aes(y = gps.era), col = "gray", lwd = 0.3)+ 
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+
    ylab("G \u2013 E") + xlab("")+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGE <- pGE +
      geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
 
  print("G-E")
  # PLOT for G'-E' -------------------------------------------------------
  # read result FGLS 5 series 
  n1 = final.t$n1[ind.case]
  n2 = final.t$n2[ind.case]
  
  station = get(load(file = paste0(path_results,"attribution0/FGLS-full/",name.GG, "fgls.RData")))
  
  name.s = "gps1.era1"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  if(begin.date > datai$date[1]){
    begin.date =  datai$date[1]
  }
  if(end.date < datai$date[length(datai$date)]){
    end.date = datai$date[length(datai$date)]
  }
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0(toupper(station.nearby), ", Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", n1, ", n2 = ", n2)
  # plot metadata 
  meta.g1 = valid2[which(valid2$name == station.nearby),]
  meta.g1 = meta.g1[which(meta.g1$known>begin.date & meta.g1$known<end.date),]
  list.meta1 = meta.g1$known
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pG1E1 <- ggplot(data = datai, aes(x = date, y = gps1.era1)) +
    theme_bw() + geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ 
    ylab("G' \u2013 E'") +
    xlab("")+  
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
  # meta for nearby
  if(nrow(meta.g1)>0){
    pG1E1 <- pG1E1 +
      geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  print("G-E")
  
  # PLOT G-G'
  
  name.s = "gps.gps"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case],
                 ", dist = ", round(final.t$distance[ind.case]) ,"(km)")
  
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pGG <- ggplot(data = datai, aes(x = date, y = gps.gps)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 G'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGG <- pGG + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
  if(nrow(meta.g1)>0){
    pGG <- pGG + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  
  # PLOT E-E'
  
  name.s = "era.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  pEE <- ggplot(data = datai, aes(x = date, y = era.era)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("E \u2013 E'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  # PLOT G-E'
  
  name.s = "gps.era1"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  pGE1 <- ggplot(data = datai, aes(x = date, y = gps.era1)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 E'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGE1 <- pGE1 + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
 
  # PLOT G'-E
  
  name.s = "gps1.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pG1E <- ggplot(data = datai, aes(x = date, y = gps1.era)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G' \u2013 E ") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g1)>0){
    pG1E <- pG1E + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  
  list.seq = mapply(function(x, y) yearser(x, y), begin.date, end.date)
  
  pGE <- pGE + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  pG1E1 <- pG1E1 + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  pGG <- pGG + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  pEE<- pEE+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  pGE1 <- pGE1+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  pG1E <- pG1E+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                        breaks = function(x) as.Date(list.seq),
                        date_labels = "%Y" )
  GE <- ggplotGrob(pGE)
  G.E. <- ggplotGrob(pG1E1)
  GG. <- ggplotGrob(pGG)
  GE. <- ggplotGrob(pGE1)
  EE. <- ggplotGrob(pEE)
  G.E <- ggplotGrob(pG1E)
  G.E.$widths <- GE$widths
  GG.$widths <- GE$widths
  GE.$widths <- GE$widths
  EE.$widths <- GE$widths
  G.E $widths <- GE$widths
  
  
  grid.newpage()
  p = (grid.arrange(GE,  EE., GG.,  G.E.,  GE., G.E, nrow = 3))
  
  ggsave(paste0(path_results,"attribution0/six_diff/", name.case,"config",final.t$pred.y[ind.case], ".jpg" ), plot = p, width = 15.5, height = 11.5, units = "cm", dpi = 1200)
  
}
plot_six_residual <- function(name.case){
  
  # name of main and 5 others
  name.GG =  name.case
  station.name = substr(name.GG, 1, 4)
  brp =  as.Date(substr(name.GG, 6, 15),format="%Y-%m-%d")
  station.nearby =  substr(name.GG, 17, 20)
  ind.case = which(all.case %in% name.GG == TRUE)
  name.GE = all.case[which(final.t$main == station.name & final.t$brp == brp & is.na(final.t$tGE)==FALSE)]
  
  # plot for G-E------------------------
  data.i = dat[[name.GG]]
  name.s = "gps.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  
  # read result FGLS for G-E
  station1 = get(load(file = paste0(path_results,"attribution0/FGLS-GE/",name.GE, "fgls.RData")))
  n1GE = length(na.omit(datai[(1:brp.ind),name.s]))
  n2GE = length(na.omit(datai[-(1:brp.ind),name.s]))
  
  if((n1GE+n2GE) >(length(station1$gps.era$fit)+50)){
    ind.wtna = which(is.na(datai[,name.s]) == FALSE)
    brp.ind = which(ind.wtna == brp.ind)
    min.ind = ind.wtna[max(1,(brp.ind-999), na.rm = TRUE)]
    max.ind = ind.wtna[min((brp.ind+999), length(ind.wtna), na.rm = TRUE)]
    datai = datai[c(min.ind:max.ind),]
    brp.ind = which(datai$date == brp)
    n1GE = length(na.omit(datai[c(min.ind:brp.ind),name.s]))-1
    n2GE = length(na.omit(datai[c(brp.ind:max.ind),name.s]))
  }
  
  begin.date = datai$date[1]
  end.date = datai$date[length(datai$date)]
  
  set.margin = list(lower = floor(min(data.i[,c(1:6)], na.rm = TRUE)), 
                    upper = ceiling(max(data.i[,c(1:6)], na.rm = TRUE)))
  datai$fit = station1[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station1[[name.s]]$t.table$Estimate[10] + station1[[name.s]]$t.table$Estimate[9]
  text1 = paste0(toupper(station.name), ", Jump = ", round(station1[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station1[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station1[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station1$gps.era$coef.arma$phi, digits = 2), ", MA: ", round(station1$gps.era$coef.arma$theta, digits = 2),
                 ", n1 = ", n1GE, ", n2 = ", n2GE)
  
  # plot metadata 
  meta.g = valid1[which(valid1$name == station.name),]
  meta.g = meta.g[which(meta.g$known>begin.date & meta.g$known<end.date),]
  list.meta = meta.g$known
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  
  pGE <- ggplot(data = datai, aes(x = date)) +
    theme_bw() + 
    geom_line(aes(y = gps.era), col = "gray", lwd = 0.3)+ 
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+
    ylab("G \u2013 E") + xlab("")+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGE <- pGE +
      geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
  
  print("G-E")
  # PLOT for G'-E' -------------------------------------------------------
  # read result FGLS 5 series 
  n1 = final.t$n1[ind.case]
  n2 = final.t$n2[ind.case]
  
  station = get(load(file = paste0(path_results,"attribution0/FGLS-full/",name.GG, "fgls.RData")))
  
  name.s = "gps1.era1"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  if(begin.date > datai$date[1]){
    begin.date =  datai$date[1]
  }
  if(end.date < datai$date[length(datai$date)]){
    end.date = datai$date[length(datai$date)]
  }
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0(toupper(station.nearby), ", Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", n1, ", n2 = ", n2)
  # plot metadata 
  meta.g1 = valid2[which(valid2$name == station.nearby),]
  meta.g1 = meta.g1[which(meta.g1$known>begin.date & meta.g1$known<end.date),]
  list.meta1 = meta.g1$known
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pG1E1 <- ggplot(data = datai, aes(x = date, y = gps1.era1)) +
    theme_bw() + geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ 
    ylab("G' \u2013 E'") +
    xlab("")+  
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none", 
          plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
  # meta for nearby
  if(nrow(meta.g1)>0){
    pG1E1 <- pG1E1 +
      geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  print("G-E")
  
  # PLOT G-G'
  
  name.s = "gps.gps"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case],
                 ", dist = ", round(final.t$distance[ind.case]) ,"(km)")
  
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pGG <- ggplot(data = datai, aes(x = date, y = gps.gps)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 G'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGG <- pGG + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
  if(nrow(meta.g1)>0){
    pGG <- pGG + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  
  # PLOT E-E'
  
  name.s = "era.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  pEE <- ggplot(data = datai, aes(x = date, y = era.era)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("E \u2013 E'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  # PLOT G-E'
  
  name.s = "gps.era1"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  pGE1 <- ggplot(data = datai, aes(x = date, y = gps.era1)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 E'") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g)>0){
    pGE1 <- pGE1 + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
  }
  
  # PLOT G'-E
  
  name.s = "gps1.era"
  datai = remove_na_2sides(data.i, name.s)
  brp.ind = which(datai$date == brp)
  datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
  brp.ind = which(datai$date == brp)
  
  datai$fit = station[[name.s]]$t.table$Estimate[10]
  datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
  text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
                 ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
                 ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
                 ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
  
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  
  pG1E <- ggplot(data = datai, aes(x = date, y = gps1.era)) +
    theme_bw() + 
    geom_line(col = "gray", lwd = 0.3)+
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G' \u2013 E ") +
    xlab("")+
    scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
                       limits = c(set.margin$lower, set.margin$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5),
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5),
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(),
          legend.position = "none", 
          plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  if(nrow(meta.g1)>0){
    pG1E <- pG1E + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
  }
  
  list.seq = mapply(function(x, y) yearser(x, y), begin.date, end.date)
  
  pGE <- pGE + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                            breaks = function(x) as.Date(list.seq),
                            date_labels = "%Y" )
  pG1E1 <- pG1E1 + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                                breaks = function(x) as.Date(list.seq),
                                date_labels = "%Y" )
  pGG <- pGG + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                            breaks = function(x) as.Date(list.seq),
                            date_labels = "%Y" )
  pEE<- pEE+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                          breaks = function(x) as.Date(list.seq),
                          date_labels = "%Y" )
  pGE1 <- pGE1+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                             breaks = function(x) as.Date(list.seq),
                             date_labels = "%Y" )
  pG1E <- pG1E+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
                             breaks = function(x) as.Date(list.seq),
                             date_labels = "%Y" )
  GE <- ggplotGrob(pGE)
  G.E. <- ggplotGrob(pG1E1)
  GG. <- ggplotGrob(pGG)
  GE. <- ggplotGrob(pGE1)
  EE. <- ggplotGrob(pEE)
  G.E <- ggplotGrob(pG1E)
  G.E.$widths <- GE$widths
  GG.$widths <- GE$widths
  GE.$widths <- GE$widths
  EE.$widths <- GE$widths
  G.E $widths <- GE$widths
  
  
  grid.newpage()
  p = (grid.arrange(GE,  EE., GG.,  G.E.,  GE., G.E, nrow = 3))
  
  ggsave(paste0(path_results,"attribution0/six_diff/", name.case,"config",final.t$pred.y[ind.case], ".jpg" ), plot = p, width = 15.5, height = 11.5, units = "cm", dpi = 1200)
  
}
