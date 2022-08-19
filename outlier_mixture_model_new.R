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
  if(nrow(na.omit(before)) > 30){
    bef.norm = one.step.norm(Y = before, list.test[1])
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
  if(nrow(na.omit(after)) > 30){
    aft.norm = one.step.norm(Y = after, list.test[1])
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
save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"onestep.RData"))
dat = get(load( file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"onestep.RData")))

# test for station GOPE --------------------------------------------
case.name = names(data.cr)[i]
d <- data.frame(Group = rep(c("raw","norm"), each = length(data.all[[case.name]]$bef$gps.era1)), 
                Sample =c(before$gps.era1, data.all[[case.name]]$bef$gps.era1))
d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group)) + 
  geom_density()+scale_x_continuous(breaks=seq(-6,7,1))+  theme_bw()


plot(x, type = "l")
points(x = a$point.rm, y = x[a$point.rm], col = "red", cex=2)
points(x = b$point.rm, y = x[b$point.rm], col = "green")


d <- data.frame(Group = rep(c("norm", "s1","s2"), each = length(x)), 
                Sample = c(x, x1, x2))
                
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()

                
                
                


# see the qqplot of all cases ---------------------------
name.var = list.test[2]

res.x = list()
res.y = list()
for (i in c(1:length(dat))) {
  x = dat[[i]]$bef[[name.var]]
  if(is.null(x)==FALSE & sum(is.na(x))!= length(x) ){
    a = qqnorm(x, plot.it = FALSE)
    res.x[[i]] <- a$x
    res.y[[i]] <- a$y
  }
}
all.x= unlist(res.x, use.names = FALSE)
lim = seq(-3.2,3.2,0.01)
limy = rep(NA, length(lim))
for (r in c(2:length(lim))) {
  mini = lim[r-1]
  maxi = lim[r]
  y <- c()
  for (l in c(1: length(res.y))) {
    ind = which(res.x[[i]] > mini & res.x[[i]] < maxi)
    if(length(ind)!=0){ y = c(y,res.y[[i]][ind])}
  }
  limy[r] <- mean(y, na.rm = TRUE)
}




# screening all data ------------------------------------------------------

n1 = length(dat)
res <- data.frame(matrix(NA, nrow = n1, ncol = 3))
length.rm = data.frame(matrix(NA, nrow = n1, ncol = 2))
for (i in seq(1,n1,1)) {
  x = dat[[i]]$bef$gps.gps
  date1 = dat[[i]]$bef$date
  if(is.null(x)==FALSE & sum(is.na(x))!= length(x) ){
    a = screen.O(x)
    b = screen.diff.o(x, dif = 1, thres1 = 3/(sqrt(2)), thres2 = 3,sdt = 1)
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
    Y$s <- sapply(Y$s, function(x) {
      if(x==0){y = "n"
      }else if (x ==1) {y="s1"
      }else if(x==2){y="s2"
      }else{y="both"}} )
    
    Y1 <- tidyr::complete(Y, date = seq(min(Y$date), max(Y$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+geom_point( aes(col = s), size = 0.7)+
      scale_color_manual(values=c("n"= "black", "s1"= " green", "s2" = "red", "both" ="blue"))+
      labs(subtitle = paste0(names(data.all)[i], "    P.val = ", paste(round( res[i,], digits = 3), collapse  = ",  ")))+
      theme(axis.text = element_text(face="bold",size=14))
    jpeg(paste0(path_results,"attribution/test.scr/",names(data.all)[i], "b.jpeg"),
         width = 4000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
  }
}

colnames(res) <- c("raw", "scr1", "scr2")
boxplot(res, main = "boxplot of p-value")
save(res, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,"p.val.RData"))

r = which(res$scr1<0.01)
which(names(dat) == "medi.2006-05-19.rovi")

d <- data.frame(Group= rep(c("raw","scr1","scr2"), each = length(bef.norm)), Sample=c(before$gps.era1, bef.norm))
qplot(sample=Sample, data=d, color=as.factor(Group))

# test screening method on simulation 

res <- data.frame(matrix(NA, ncol = 2, nrow = 1000))
list.rm <- c()
for (j in c(1:1000)) {
  set.seed(j)
  x = rnorm(1000,0,1)
  x[seq(100, 1000, 100)] <- rep(5, 10)
  r = screen.O(x)
  list.rm <- c(list.rm, r$point.rm)
  res[j,] <- c(shapiro.test(x)$p.value, shapiro.test(r$data)$p.value)
}
list.rm[! duplicated(list.rm)]
summary(res)

which(names(dat) == "auck.2005-11-07.hamt")



