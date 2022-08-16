# distribution of real data after restricted homogeneity
source(paste0(path_code_att,"sliding_variance.R"))

# # Normalize data first:  ------------------------------------------------
window.thres = 10
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
    bef.norm = two.step.norm(Y = before, name.var)
    if(length(bef.norm)>30){
      hist.b = hist(bef.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.b$counts)
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef<- bef.norm
    }
  }
  if(nrow(na.omit(after)) > 31){
    aft.norm = two.step.norm(Y = after, name.var)
    if(length(aft.norm) > 30){
      hist.a = hist(aft.norm, breaks = seq(-30,30,0.05), plot = FALSE)
      dist.mean <- rbind(dist.mean, hist.a$counts)
      data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm
    }
  }
}
save(dist.mean, file = paste0(path_results,"attribution/dist.mean_2years_", nearby_ver,".RData"))
save(data.all, file = paste0(path_results,"attribution/data.all_2years_", nearby_ver,".RData"))



res <- data.frame(t = seq(-29.95,30,0.05), count = colMeans(dist.mean))

d = res[c(400:800),]

ggplot(d, aes(x = t, y = count))+geom_line()
fitG = function(x,y,mu,sig,scale){
    
    f = function(p){
      d = p[3]*dnorm(x,mean=p[1],sd=p[2])
      sum((d-y)^2)
    }
    
    optim(c(mu,sig,scale),f)
  }
p1 <- fitG(x = d$t, y = d$count, mu = 0, sig = 1, scale =1)
p = p1$par
plot(d$t, d$count)
lines(d$t,p[3]*dnorm(d$t,p[1],p[2]))

# MEAN VALUE INDICATE THE PURE GAUSSIAN --> LOOK INTO A SMALLER REGION WHERE NB OF POINT >0 IN OUTLIER PERIOD
a = colSums(dist.mean)

res <- data.frame(t = seq(-29.95,30,0.05), count = colMeans(dist.mean))

d <- tibble(Group=rep(1:2694, each=1200), Sample=as.vector(t(dist.mean)))

e = unlist(sapply(c(21), function(x) lengths(data.all[[x]])))
# d <- tibble(Group=rep(1:(length(e)), times=e), Sample=unlist(data.all))
d <- data.frame(Group=rep((21), times=e[21]), Sample=unlist(data.all[c(21)]))
qplot(sample=Sample, data=d, color=as.factor(Group))

plot( before$gps.era1, type = "l")
lines(bef.norm, col="red")
abline(h = -2.5)

d <- data.frame(Group= rep(c("raw","norm"), each = length(bef.norm)), Sample=c(before$gps.era1, bef.norm))
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()

d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group, fill = Group)) + 
  geom_histogram(alpha = 0.5, position = "identity", bins = 200)+
  theme_bw()

# choose only europe stations 

list.pos = read.table(file = paste0(path_data_support,"gps_sta_CODE_REPRO_2015_OPER_combi.txt"),header = TRUE)
colnames(list.pos) <- c("name", "lat", "lon", "hei","alt")
list.f = list.pos[which(list.pos$lat>30 & list.pos$lat<50 & list.pos$lon >-10 & list.pos$lon <30),]
name.f = tolower(list.f$name)
list.name.eu = name_main[name_main%in%name.f]

name.var = list.test[3]
dist.mean <- data.frame(matrix(NA, nrow = 0, ncol = (length(seq(-30,30,0.05))-1)))
data.all <- list()
for (i in c(1:length(data.cr))) {
  case.name = names(data.cr)[i]
  station.ref = substr(case.name ,start = 1, stop = 4)
  if(station.ref %in% list.name.eu){
    station.near = substr(case.name ,start = 17, stop = 20)
    data.i = data.cr[[i]]
    breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
    before =  data.i[which(data.i$date <= breakpoint),]
    after =  data.i[which(data.i$date > breakpoint),]
    if(nrow(na.omit(before)) > 31){
      bef.norm = two.step.norm(Y = before, name.var)
      if(length(bef.norm)>30){
        # hist.b = hist(bef.norm, breaks = seq(-30,30,0.05), plot = FALSE)
        # dist.mean <- rbind(dist.mean, hist.b$counts)
        data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$bef<- bef.norm
      }
    }
    if(nrow(na.omit(after)) > 31){
      aft.norm = two.step.norm(Y = after, name.var)
      if(length(aft.norm) > 30){
        # hist.a = hist(aft.norm, breaks = seq(-30,30,0.05), plot = FALSE)
        # dist.mean <- rbind(dist.mean, hist.a$counts)
        data.all[[paste0(station.ref,".",as.character( breakpoint), ".", station.near)]]$aft <- aft.norm
      }
    }
  }
}

e = unlist(sapply(c(217:220), function(x) lengths(data.all[[x]])))
# d <- tibble(Group=rep(1:(length(e)), times=e), Sample=unlist(data.all))
d <- data.frame(Group=rep((1:length(e)), times=e), Sample=unlist(data.all[217:220]))
qplot(sample=Sample, data=d, color=as.factor(Group))+theme_bw()


d$Group <- as.factor(d$Group)
ggplot(d, aes(x = Sample, colour = Group)) + 
  geom_density()+scale_x_continuous(breaks=seq(-6,7,1))+
  # geom_histogram(alpha = 0.5, position = "identity", bins = 200)+
  theme_bw()+  theme(legend.position="none")


coun = 0
for (i in c(217)) {
  case.name = names(data.all)[i]
  l = which( names(data.cr) == case.name)
  station.ref = substr(case.name ,start = 1, stop = 4)
  station.near = substr(case.name ,start = 17, stop = 20)
  data.i = data.cr[[l]]
  breakpoint = as.Date(substr(case.name,start = 6, stop = 15) , format = "%Y-%m-%d")
  before =  data.i[which(data.i$date <= breakpoint),]
  after =  data.i[which(data.i$date > breakpoint),]
  if(nrow(na.omit(before)) > 31){
    coun = coun+1
    bef.norm = two.step.norm(Y = before, name.var)
    data.b = data.frame(date = before$date, val = bef.norm)
    Y1 <- tidyr::complete(data.b, date = seq(min(data.b$date), max(data.b$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+labs(subtitle = paste0("case", coun))
    jpeg(paste0(path_results,"attribution/",station.ref,".",as.character( breakpoint), ".", station.near, "b.jpeg"), width = 3000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
    
  }
  if(nrow(na.omit(after)) > 31){
    coun = coun+1
    aft.norm = two.step.norm(Y = after, name.var)
    data.a = data.frame(date = after$date, val = aft.norm)
    Y1 <- tidyr::complete(data.a, date = seq(min(data.a$date), max(data.a$date), by = "day"))
    p <- ggplot(data = Y1, aes(x = date, y= val))+
      geom_line()+theme_bw()+labs(subtitle = paste0("case", coun))
    jpeg(paste0(path_results,"attribution/",station.ref,".",as.character( breakpoint), ".", station.near, "a.jpeg"), width = 3000, height = 1800,res = 300) # change name
    print(p)
    dev.off()
  }
}