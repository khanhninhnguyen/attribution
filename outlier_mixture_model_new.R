source(paste0(path_code_att,"sliding_variance.R"))
source(paste0(path_code_att,"support_screening.R"))

# # Normalize data first:  ------------------------------------------------
window.thres = 2
data.cr = get(load( file = paste0(path_results,"attribution/six_diff_series_rm_crenel_restricted_closed_brp_",
                                  window.thres,"year_", nearby_ver,".RData")))
name.var = list.test[4]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()

for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[i]]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 31){
    bef.norm = two.step.norm(Y = before, list.test[1])
    if(length(bef.norm)>30){
      bef.norm.all <- list()
      for (k in c(1:6)) {
        norm.dat = one.step.norm(Y = before, list.test[k])
        names(norm.dat) <- NULL
        bef.norm.all[[list.test[k]]] <- norm.dat
      }
      bef.norm.all[["date"]] = before$date
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef <- bef.norm.all
    }
  }
  if(nrow(na.omit(after)) > 31){
    aft.norm = two.step.norm(Y = after, list.test[1])
    if(length(aft.norm) > 30){
      aft.norm.all <- list()
      for (k in c(1:6)) {
        norm.dat <- one.step.norm(Y = after, list.test[k])
        names(norm.dat) <- NULL
        aft.norm.all[[list.test[k]]] <- norm.dat
      }
      aft.norm.all[["date"]] = after$date
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm.all
    }
  }
}
# save(dist.mean, file = paste0(path_results,"attribution/dist.mean_2years_", nearby_ver,".RData"))
save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,".RData"))

# test for station GOPE
case.name = "gope.2006-09-07.wtzt"
d <- data.frame(Group = rep(c("raw","norm"), each = length(data.all[[case.name]]$bef$gps.gps)), 
                Sample =c(rnorm(236,0,1), data.all[[case.name]]$bef$gps.gps))
d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group)) + 
  geom_density()+scale_x_continuous(breaks=seq(-6,7,1))+  theme_bw()

n1 = 1552
res <- data.frame(matrix(NA, nrow = 1552, ncol = 3))
length.rm = data.frame(matrix(NA, nrow = 1552, ncol = 2))
for (i in seq(1,n1,100)) {
  x = data.all[[i]]$bef$gps.gps
  date1 = data.all[[i]]$bef$date
  if(is.null(x)==FALSE & sum(is.na(x))!= length(x) ){
    a = screen.qn.o(x, thres = 2, sdt = 1/1.349)
    b = screen.diff.o(x, dif = 1, thres = 2, sdt = 1/1.349)
    s1 = rep(0, length(x))
    s1[a$point.rm] <- 1
    s2 = rep(0, length(x))
    s2[b$point.rm] <- 2
    x1 <- x
    x1[a$point.rm] <- NA
    x2 <- x
    x2[b$point.rm] <- NA
    res[i,] <- c(shapiro.test(x)$p.value, shapiro.test(x1)$p.value, shapiro.test(x2)$p.value)
    length.rm[i,] <- c(length(a$point.rm), length(b$point.rm))
    Y = data.frame( date = date1, val = x, s = as.factor(s1+s2))
    Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+geom_point( aes(col = s), size = 0.7)+
      scale_color_manual(values=c("black","green","red","blue"))+
      labs(subtitle = names(data.all)[i])+
      theme(axis.text = element_text(face="bold",size=14))
    jpeg(paste0(path_results,"attribution/test.scr/",names(data.all)[i], "b.jpeg"),
         width = 3000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
  }
}

plot(x, type = "l")
points(x = a$point.rm, y = x[a$point.rm], col = "red", cex=2)
points(x = b$point.rm, y = x[b$point.rm], col = "green")


d <- data.frame(Group = rep(c("norm", "s1","s2"), each = length(x)), 
                Sample = c(x, x1, x2))
                
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()

                
                
                