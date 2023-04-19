# add some tests according to the paper review 

# autocorrelation between GPS and ERA -------------------------------------


list.station = substr(list.files(path = path_CODE_aux_ERA5_v1b), 1,4)
correlation = rep(NA, length(list.station))
for (i in c(1:length(list.station))) {
  a = read.series(path_CODE_aux_ERA5_v1b, station = list.station[i], na.rm = 0, add.full = 0)
  correlation[i] = cor(a$ERAI, a$GPS, use='complete.obs')
}

