# function to plot time series 
source(paste0(path_code_att,"FGLS.R"))

Total.res = get(load(paste0(path_results,"attribution0/stats_test_real_data.RData")))
dat = get(load( file = paste0(path_results,"attribution0/data.all_", win.thres = 10,"years_", nearby_ver,"screened.RData")))
final.t = get(load(file = paste0(path_results, "attribution0/Final.Table.RData")))
rownames(final.t) = NULL
lengthlist = get(load(file = paste0(path_results, "attribution0/lengthlist.RData")))
final.t$n1 = lengthlist$X1
final.t$n2 = lengthlist$X2
# valid1 = get(load(file = paste0(path_results,"validation/",nb_test.ref,"-",criterion,"metacompa",screen.value="",".RData")))
valid2 = get(load(file = paste0(path_results,"validation/",nb_test.near,"-",criterion,"metacompa",screen.value="",".RData")))
meta1 = read.meta(nb_test = nb_test.ref, path_meta = path_data_support)
meta2 = read.meta(nb_test = nb_test.near, path_meta = path_data_support)

all.case = paste0(final.t$main,".",as.character(final.t$brp), ".", final.t$nearby)
yearser = function(sd, ed){
  paste0(substr(seq(sd, ed, "years"), 1, 4)[-1], "-01-01")
}
sliding.median <- function(Y, name.var, length.wind){
  Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
  x = unlist(Y1[name.var], use.names = FALSE)
  n = length(x)
  slide.med <- rep(NA, n)
  l = rep(NA, n)
  for (i in c(1:n)) {
    begin = max(i-(length.wind-1),1)
    end = min(n, i+length.wind)
    x.i = x[begin:end]
    thre = 30
    if(i < 30|i>(n-30)){thre = 16}
    if(length(na.omit(x.i))> thre){
      slide.med[i] <- median(x.i, na.rm = TRUE)
      l[i] <- length(na.omit(x.i))
    }
  }
  slide.med = slide.med[which(Y1$date %in% Y$date)]
  # l = l[which(Y1$date %in% Y$date)]
  return(slide.med)
  # return(l)
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
  meta.g = meta1[which(meta1$name == station.name),]
  meta.g = meta.g[which(meta.g$ymd > begin.date & meta.g$ymd < end.date),]
  list.meta = meta.g$ymd
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
  meta.g1 = meta2[which(meta2$name == station.nearby),]
  meta.g1 = meta.g1[which(meta.g1$ymd>begin.date & meta.g1$ymd<end.date),]
  list.meta1 = meta.g1$ymd
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  # mark if there is any break in the nearby
  nearby.brp = valid2[which(valid2$name == station.nearby),]
  nearby.brp = nearby.brp[which(nearby.brp$detected > begin.date & nearby.brp$detected < end.date),]
  list.nearby.brp = nearby.brp$detected
  if(nrow(nearby.brp)>0){
    datai$detected = NA
    datai$detected[which( datai$date %in% list.nearby.brp == TRUE)] = set.margin$upper
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
  if(nrow(nearby.brp)>0){
    pG1E1 <- pG1E1 +
      geom_point(aes(y = detected), colour="black",shape = 7, size = 0.5)
  }
  
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
                 ", h.dist = ", round(final.t$distance[ind.case]) ,"(km)",
                 ", v.dist = ", round(Total.res$ver.distance[ind.case]) ,"(m)")
  
  if(nrow(meta.g)>0){
    datai$meta = NA
    datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
  }
  if(nrow(meta.g1)>0){
    datai$meta1 = NA
    datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
  }
  if(nrow(nearby.brp)>0){
    datai$detected = NA
    datai$detected[which( datai$date %in% list.nearby.brp == TRUE)] = set.margin$upper
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
  if(nrow(nearby.brp)>0){
    pG1E1 <- pG1E1 +
      geom_point(aes(y = detected), colour="black",shape = 7, size = 0.5)
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
  if(nrow(nearby.brp)>0){
    datai$detected = NA
    datai$detected[which( datai$date %in% list.nearby.brp == TRUE)] = set.margin$upper
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
  if(nrow(nearby.brp)>0){
    pG1E1 <- pG1E1 +
      geom_point(aes(y = detected), colour="black",shape = 7, size = 0.5)
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

plot_FGLS_details <- function(name.case, name.series, lab.y){
  
  # name of main and 5 others
  brp =  as.Date(substr(name.case, 6, 15),format="%Y-%m-%d")
  station.name = substr(name.case, 1, 4)
  station.nearby =  substr(name.case, 17, 20)
  ind.case = which(all.case %in% name.case == TRUE)
  
  # plot for G-E------------------------
  station = get(load(file = paste0(path_results,"attribution0/FGLS/",name.case, "fgls.RData")))
  datai = station[[name.series]]$design.matrix
  brp.ind = which(datai$date == brp)
  
  # read result FGLS for G-E
  n1GE = length(na.omit(datai[(1:brp.ind),"signal"]))
  n2GE = length(na.omit(datai[-(1:brp.ind),"signal"]))
  
  begin.date = datai$date[1]
  end.date = datai$date[length(datai$date)]
  
  # # add the fit and Fourier series
  # datai$fit = station[[name.series]]$t.table$Estimate[10]
  # datai$fit[(brp.ind+1):nrow(datai)] = station[[name.series]]$t.table$Estimate[10] + station[[name.series]]$t.table$Estimate[9]
  # datai$fourier = as.matrix(datai[,c(2:9)]) %*% as.matrix(station[[name.series]]$t.table$Estimate[1:8])
  
  datai$fit = NA
  datai$fit[which(is.na(datai$signal) == FALSE)] = station[[name.series]]$fit
  # add infor of test 
  text1 = paste0(toupper(name.case), ", Jump = ", round(station[[name.series]]$t.table$Estimate[9], digits = 2), 
                 ", t = ", round(station[[name.series]]$t.table$`t value`[9], digits = 2), 
                 ", n1 = ", n1GE, ", n2 = ", n2GE)
  text2 = paste0(" mean.MSD = ", round(mean(sqrt(station[[name.series]]$var), na.rm = TRUE) , digits = 2), 
                 ", AR: ", round(station[[name.series]]$coef.arma$phi, digits = 2), ", MA: ", round(station[[name.series]]$coef.arma$theta, digits = 2))
  set.margin = list(lower = floor(min(datai$signal, na.rm = TRUE)),
                    upper = ceiling(max(datai$signal, na.rm = TRUE)))
  # datai$fourier = datai$fourier - abs(min(datai$signal, na.rm = TRUE)) -1
  # plot metadata 
  if ( name.series %in% list.test[c(1,3)] == TRUE){
    meta.data = meta1
  } else if (name.series %in% list.test[c(5,6)] == TRUE){
    meta.data = meta2
  } else if (name.series == list.test[2]){
    meta.data = rbind(meta1, meta2)
  } else{
    meta.data = data.frame()
  }
  print(nrow(meta.data))
  if (nrow(meta.data)>0){
    meta.g = meta1[which(meta1$name == station.name),]
    meta.g = meta.g[which(meta.g$ymd>begin.date & meta.g$ymd<end.date),]
    list.meta = meta.g$ymd
    if(nrow(meta.g)>0){
      datai$meta = NA
      datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
    }
  }
  p <- ggplot(data = datai, aes(x = date)) +
    theme_bw() + 
    geom_line(aes(y = signal), col = "gray", lwd = 0.3)+ 
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    labs(subtitle = text1)+
    geom_line(aes(y = fit), col = "red", lwd = 0.3)+
    # geom_line(aes(y = fourier), col = "blue", lwd = 0.3)+
    ylab(lab.y) + xlab("")+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none")
  
  
  
  if(nrow(meta.data)>0){
    if(nrow(meta.g)>0){
      p <- p +
        geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
    }
  }
  
  # plot the residual 
  # compute the moving median 
  datai$mm = sliding.median(datai, name.var = "residual", length.wind = 60)
  datai$sd = sqrt(station[[name.series]]$var) - abs(min(datai$residual, na.rm = TRUE)) -1 
  set.margin1 = list(lower = floor(min(datai$sd, na.rm = TRUE)),
                    upper = ceiling(max(datai$residual, na.rm = TRUE)))
  p1 <- ggplot(data = datai, aes(x = date)) +
    theme_bw() + 
    geom_line(aes(y = residual), col = "gray", lwd = 0.3)+ 
    geom_line(aes(y = sd), col = "blue", lwd = 0.3)+ 
    geom_line(aes(y = mm), col = "red", lwd = 0.3)+ 
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    ylab("Residual") + xlab("")+ labs(subtitle = text2)+
    scale_y_continuous(breaks = seq(set.margin1$lower, set.margin1$upper,1),
                       limits = c(set.margin1$lower, set.margin1$upper))+
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none",plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  
  # plot the normalized residual 
  order.aima = c(ifelse(station$gps.era$coef.arma$phi !=0, 1,0), 0, ifelse(station$gps.era$coef.arma$theta !=0, 1,0))
  names(order.aima) = NULL
  arima.fit = arima(datai$norm.res,order = order.aima)
  port.test = Box.test(arima.fit$residuals)
  p2 <- ggplot(data = datai, aes(x = date)) +
    theme_bw() + 
    geom_line(aes(y = norm.res), col = "gray", lwd = 0.3)+ 
    geom_vline(xintercept = brp, lwd = 0.2)+ 
    ylab("Normalized Residual") + xlab("")+
    labs(subtitle = paste("mean = ", round(mean(datai$norm.res, na.rm = TRUE), digits = 2), 
                          " sd = ", round(sd(datai$norm.res, na.rm = TRUE), digits = 2), 
                          " p.portmanteau =", round(port.test$p.value, digits = 2))) +
    theme(axis.text.x = element_text(size = 5), 
          axis.text.y = element_text(size = 5),
          legend.text=element_text(size=4),
          axis.title = element_text(size = 5), 
          legend.key.size = unit(0.3, "cm"), 
          plot.tag = element_text(size = 5), 
          plot.subtitle = element_text(size = 5),
          legend.title=element_blank(), 
          legend.position = "none",plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
  grid.newpage()
  G <- ggplotGrob(p)
  G1 <- ggplotGrob(p1)
  G2 <- ggplotGrob(p2)
  G$widths <- G2$widths
  G1$widths <- G2$widths
  p.all = (grid.arrange(p, p1, p2, nrow = 3))
  
  ggsave(paste0(path_results,"attribution0/six_diff/", name.case,"config",final.t$pred.y[ind.case], name.series,".jpg" ), plot = p.all, width = 8.8, height = 10, units = "cm", dpi = 1200)
  
}
# plot_FGLS_details(name.case = "pots.2016-06-05.d020", name.series = "gps.era", lab.y = "G-E")

plot_six_details <- function(name.case){
  for (name.s in c(1:6)) {
    name.series = list.test[name.s]
    lab.y = list.name.test[name.s]
    plot_FGLS_details(name.case, name.series, lab.y)
  }
}

# plot_six_residual <- function(name.case){
#   
#   # name of main and 5 others
#   name.GG =  name.case
#   station.name = substr(name.GG, 1, 4)
#   brp =  as.Date(substr(name.GG, 6, 15),format="%Y-%m-%d")
#   station.nearby =  substr(name.GG, 17, 20)
#   ind.case = which(all.case %in% name.GG == TRUE)
#   name.GE = all.case[which(final.t$main == station.name & final.t$brp == brp & is.na(final.t$tGE)==FALSE)]
#   
#   # plot for G-E------------------------
#   data.i = dat[[name.GG]]
#   name.s = "gps.era"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   
#   # read result FGLS for G-E
#   station1 = get(load(file = paste0(path_results,"attribution0/FGLS-GE/",name.GE, "fgls.RData")))
#   n1GE = length(na.omit(datai[(1:brp.ind),name.s]))
#   n2GE = length(na.omit(datai[-(1:brp.ind),name.s]))
#   
#   if((n1GE+n2GE) >(length(station1$gps.era$fit)+50)){
#     ind.wtna = which(is.na(datai[,name.s]) == FALSE)
#     brp.ind = which(ind.wtna == brp.ind)
#     min.ind = ind.wtna[max(1,(brp.ind-999), na.rm = TRUE)]
#     max.ind = ind.wtna[min((brp.ind+999), length(ind.wtna), na.rm = TRUE)]
#     datai = datai[c(min.ind:max.ind),]
#     brp.ind = which(datai$date == brp)
#     n1GE = length(na.omit(datai[c(min.ind:brp.ind),name.s]))-1
#     n2GE = length(na.omit(datai[c(brp.ind:max.ind),name.s]))
#   }
#   
#   begin.date = datai$date[1]
#   end.date = datai$date[length(datai$date)]
#   
#   set.margin = list(lower = floor(min(data.i[,c(1:6)], na.rm = TRUE)), 
#                     upper = ceiling(max(data.i[,c(1:6)], na.rm = TRUE)))
#   datai$fit = station1[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station1[[name.s]]$t.table$Estimate[10] + station1[[name.s]]$t.table$Estimate[9]
#   text1 = paste0(toupper(station.name), ", Jump = ", round(station1[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station1[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station1[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station1$gps.era$coef.arma$phi, digits = 2), ", MA: ", round(station1$gps.era$coef.arma$theta, digits = 2),
#                  ", n1 = ", n1GE, ", n2 = ", n2GE)
#   
#   # plot metadata 
#   meta.g = meta1[which(meta1$name == station.name),]
#   meta.g = meta.g[which(meta.g$ymd>begin.date & meta.g$ymd<end.date),]
#   list.meta = meta.g$ymd
#   if(nrow(meta.g)>0){
#     datai$meta = NA
#     datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
#   }
#   
#   pGE <- ggplot(data = datai, aes(x = date)) +
#     theme_bw() + 
#     geom_line(aes(y = gps.era), col = "gray", lwd = 0.3)+ 
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     geom_vline(xintercept = brp, lwd = 0.2)+ 
#     labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+
#     ylab("G \u2013 E") + xlab("")+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5), 
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5), 
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(), 
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
#   
#   if(nrow(meta.g)>0){
#     pGE <- pGE +
#       geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
#   }
#   
#   print("G-E")
#   # PLOT for G'-E' -------------------------------------------------------
#   # read result FGLS 5 series 
#   n1 = final.t$n1[ind.case]
#   n2 = final.t$n2[ind.case]
#   
#   station = get(load(file = paste0(path_results,"attribution0/FGLS-full/",name.GG, "fgls.RData")))
#   
#   name.s = "gps1.era1"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
#   if(begin.date > datai$date[1]){
#     begin.date =  datai$date[1]
#   }
#   if(end.date < datai$date[length(datai$date)]){
#     end.date = datai$date[length(datai$date)]
#   }
#   brp.ind = which(datai$date == brp)
#   
#   datai$fit = station[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
#   text1 = paste0(toupper(station.nearby), ", Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
#                  ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
#                  ", n1 = ", n1, ", n2 = ", n2)
#   # plot metadata 
#   meta.g1 = meta2[which(meta2$name == station.nearby),]
#   meta.g1 = meta.g1[which(meta.g1$ymd>begin.date & meta.g1$ymd<end.date),]
#   list.meta1 = meta.g1$ymd
#   if(nrow(meta.g1)>0){
#     datai$meta1 = NA
#     datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
#   }
#   
#   pG1E1 <- ggplot(data = datai, aes(x = date, y = gps1.era1)) +
#     theme_bw() + geom_line(col = "gray", lwd = 0.3)+
#     geom_vline(xintercept = brp, lwd = 0.2)+ labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+ 
#     ylab("G' \u2013 E'") +
#     xlab("")+  
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5), 
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5), 
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(), 
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0, 0, 0.5), "cm"))
#   # meta for nearby
#   if(nrow(meta.g1)>0){
#     pG1E1 <- pG1E1 +
#       geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
#   }
#   print("G-E")
#   
#   # PLOT G-G'
#   
#   name.s = "gps.gps"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
#   brp.ind = which(datai$date == brp)
#   
#   datai$fit = station[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
#   text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
#                  ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
#                  ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case],
#                  ", dist = ", round(final.t$distance[ind.case]) ,"(km)")
#   
#   if(nrow(meta.g)>0){
#     datai$meta = NA
#     datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
#   }
#   if(nrow(meta.g1)>0){
#     datai$meta1 = NA
#     datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
#   }
#   
#   pGG <- ggplot(data = datai, aes(x = date, y = gps.gps)) +
#     theme_bw() + 
#     geom_line(col = "gray", lwd = 0.3)+
#     geom_vline(xintercept = brp, lwd = 0.2)+ 
#     labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 G'") +
#     xlab("")+
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5),
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5),
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(),
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
#   
#   if(nrow(meta.g)>0){
#     pGG <- pGG + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
#   }
#   if(nrow(meta.g1)>0){
#     pGG <- pGG + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
#   }
#   
#   # PLOT E-E'
#   
#   name.s = "era.era"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
#   brp.ind = which(datai$date == brp)
#   
#   datai$fit = station[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
#   text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
#                  ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
#                  ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
#   
#   pEE <- ggplot(data = datai, aes(x = date, y = era.era)) +
#     theme_bw() + 
#     geom_line(col = "gray", lwd = 0.3)+
#     geom_vline(xintercept = brp, lwd = 0.2)+ 
#     labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("E \u2013 E'") +
#     xlab("")+
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5),
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5),
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(),
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
#   
#   # PLOT G-E'
#   
#   name.s = "gps.era1"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
#   brp.ind = which(datai$date == brp)
#   
#   datai$fit = station[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
#   text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
#                  ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
#                  ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
#   
#   if(nrow(meta.g)>0){
#     datai$meta = NA
#     datai$meta[which( datai$date %in% list.meta == TRUE)] = set.margin$upper
#   }
#   pGE1 <- ggplot(data = datai, aes(x = date, y = gps.era1)) +
#     theme_bw() + 
#     geom_line(col = "gray", lwd = 0.3)+
#     geom_vline(xintercept = brp, lwd = 0.2)+ 
#     labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G \u2013 E'") +
#     xlab("")+
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5),
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5),
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(),
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
#   
#   if(nrow(meta.g)>0){
#     pGE1 <- pGE1 + geom_point(aes(y = meta), colour="blue",shape = 2, size = 0.5)
#   }
#   
#   # PLOT G'-E
#   
#   name.s = "gps1.era"
#   datai = remove_na_2sides(data.i, name.s)
#   brp.ind = which(datai$date == brp)
#   datai = datai[c((brp.ind-n1+1):(brp.ind+n2)),]
#   brp.ind = which(datai$date == brp)
#   
#   datai$fit = station[[name.s]]$t.table$Estimate[10]
#   datai$fit[(brp.ind+1):nrow(datai)] = station[[name.s]]$t.table$Estimate[10] + station[[name.s]]$t.table$Estimate[9]
#   text1 = paste0("Jump = ", round(station[[name.s]]$t.table$Estimate[9], digits = 2), 
#                  ", t = ", round(station[[name.s]]$t.table$`t value`[9], digits = 2), 
#                  ", SD = ", round(mean(sqrt(station[[name.s]]$var), na.rm = TRUE) , digits = 2), 
#                  ", AR: ", round(station[[name.s]]$coef.arma$phi, digits = 2), 
#                  ", MA: ", round(station[[name.s]]$coef.arma$theta, digits = 2),
#                  ", n1 = ", final.t$n1[ind.case], ", n2 = ", final.t$n2[ind.case])
#   
#   if(nrow(meta.g1)>0){
#     datai$meta1 = NA
#     datai$meta1[which( datai$date %in% list.meta1 == TRUE)] = set.margin$upper
#   }
#   
#   pG1E <- ggplot(data = datai, aes(x = date, y = gps1.era)) +
#     theme_bw() + 
#     geom_line(col = "gray", lwd = 0.3)+
#     geom_vline(xintercept = brp, lwd = 0.2)+ 
#     labs(subtitle = text1)+
#     geom_line(aes(y = fit), col = "red", lwd = 0.3)+ ylab("G' \u2013 E ") +
#     xlab("")+
#     scale_y_continuous(breaks = seq(set.margin$lower, set.margin$upper,1), 
#                        limits = c(set.margin$lower, set.margin$upper))+
#     theme(axis.text.x = element_text(size = 5), 
#           axis.text.y = element_text(size = 5),
#           legend.text=element_text(size=4),
#           axis.title = element_text(size = 5),
#           legend.key.size = unit(0.3, "cm"), 
#           plot.tag = element_text(size = 5),
#           plot.subtitle = element_text(size = 5),
#           legend.title=element_blank(),
#           legend.position = "none", 
#           plot.margin = unit(c(0, 0.5, 0, 0), "cm"))
#   
#   if(nrow(meta.g1)>0){
#     pG1E <- pG1E + geom_point(aes(y = meta1), colour="green",shape = 3, size = 0.5)
#   }
#   
#   list.seq = mapply(function(x, y) yearser(x, y), begin.date, end.date)
#   
#   pGE <- pGE + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                             breaks = function(x) as.Date(list.seq),
#                             date_labels = "%Y" )
#   pG1E1 <- pG1E1 + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                                 breaks = function(x) as.Date(list.seq),
#                                 date_labels = "%Y" )
#   pGG <- pGG + scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                             breaks = function(x) as.Date(list.seq),
#                             date_labels = "%Y" )
#   pEE<- pEE+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                           breaks = function(x) as.Date(list.seq),
#                           date_labels = "%Y" )
#   pGE1 <- pGE1+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                              breaks = function(x) as.Date(list.seq),
#                              date_labels = "%Y" )
#   pG1E <- pG1E+ scale_x_date(limits = as.Date(c(begin.date, end.date)), 
#                              breaks = function(x) as.Date(list.seq),
#                              date_labels = "%Y" )
#   GE <- ggplotGrob(pGE)
#   G.E. <- ggplotGrob(pG1E1)
#   GG. <- ggplotGrob(pGG)
#   GE. <- ggplotGrob(pGE1)
#   EE. <- ggplotGrob(pEE)
#   G.E <- ggplotGrob(pG1E)
#   G.E.$widths <- GE$widths
#   GG.$widths <- GE$widths
#   GE.$widths <- GE$widths
#   EE.$widths <- GE$widths
#   G.E $widths <- GE$widths
#   
#   
#   grid.newpage()
#   p = (grid.arrange(GE,  EE., GG.,  G.E.,  GE., G.E, nrow = 3))
#   
#   ggsave(paste0(path_results,"attribution0/six_diff/", name.case,"config",final.t$pred.y[ind.case], ".jpg" ), plot = p, width = 15.5, height = 11.5, units = "cm", dpi = 1200)
#   
# }
# 
