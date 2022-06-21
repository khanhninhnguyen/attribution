rm(list = ls())
# Main program to read all data and its path 
library("gridExtra")
library("dplyr") 
#library("colorøspace")
library("ggplot2")
library("grid")
library("forecast")
library(lmtest)
# main paths  --------------------------------------------------------------

path_main="/home/knguyen/Documents/PhD/"
setwd(path_main)
path_code_seg = paste0(path_main,"Code/segment/segmentation/")
path_code_att = paste0(path_main,"Code/attribution/")
path_data = paste0(path_main,"Data/IWV-ERA/")
path_data_support = paste0(path_main,"Data/support/")
path_homo = paste0(path_main,"Results/homo/")
path_results=paste0(path_main,"Results/")
path.validation = paste0(path_results, "validation/")

# data paths --------------------------------------------------------------

path_data_era5 = paste0(path_data, "COMPARE_ERA5/") 
path_CODE_aux_ERA5_v1b = paste0(path_data_era5,"CODE_REPRO2015_2018.v1bR/" ) # 34 new version of CODE aux ERA5 - ERA5 changed sampling method
path_NGL_v2_R = paste0(path_data_era5,"NGL.v2R/") #43 is the old data( missing few years) # 44 is the new data ( full time series )

# source code -------------------------------------------------------------
# from main path
source(paste0(path_code_seg,"read.series.segment.R")) 
source(paste0(path_code_seg,"read.meta.R"))
# from Proj Attribution
source(paste0(path_code,"attribution-seg-support.R"))
source(paste0(path_code,"ar_coef_estimation.R"))


# parameters --------------------------------------------------------------

criterion = "BM_BJ"
nb_test.ref = 34
nb_test.near = 44

path_series_main = path_CODE_aux_ERA5_v1b
version_name = as.character(expression(path_CODE_aux_ERA5_v1b))
path_series_nearby = path_NGL_v2_R
nearby_ver = "NGL"
tolerance_noise = 10
tolerance_valid = 62
thres_hdist = 500
thres_vdist = 200
limit.type = "superposed" 
check_correlation = 5

thres_period = 365
no.day.per = 40/365
significant.level = 0.05 
screening.data.gps = 2
vertical.correction = 1
screen_value = 1
screening1 = 1

validation.file.ref = paste0(path.validation,nb_test = nb_test.ref,"-",screen_value,
                             criterion,tolerance_noise, tolerance_valid,".txt")
validation.file.near = paste0(path.validation,nb_test = nb_test.near,"-",screen_value,
                              criterion,tolerance_noise, tolerance_valid,".txt")

list.test = c( "gps.gps", "gps.era1", "gps1.era", "era.era", "gps.era", "gps1.era1")

# list of 81 stations -----------------------------------------------------

list.station = read.table(paste0(path_data_support,"list-cross-stations.txt"),header = TRUE)
name_main = as.character(list.station[,1])

