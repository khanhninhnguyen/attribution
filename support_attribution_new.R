diff.var.1 <- function(name.test){
  if(name.test == "gps.gps"){
    varname = c("GPS", "GPS1")
  }
  if(name.test == "gps.era"){
    varname = c("GPS", "ERA")
  }
  if(name.test == "gps1.era"){
    varname = c("GPS1", "ERA")
  }
  if(name.test == "gps.era1"){
    varname = c("GPS", "ERA1")
  }
  if(name.test == "gps1.era1"){
    varname = c("GPS1", "ERA")
  }
  if(name.test == "era.era"){
    varname = c("ERA", "ERA1")
  }
  return(varname)
}