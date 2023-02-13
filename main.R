# main program for attribution 
rm(list = ls())
# in each step, it is necessary to create a new folder with the name of function
# initiate: read packages and functions -----------------------------------
path_code_att="/home/knguyen/Documents/PhD/Code/attribution/"
source(paste0(path_code_att,"init_info.R"))

# vertical correction -----------------------------------------------------
source(paste0(path_code_att,"vertical_correction.R"))
# vertical_correction(window.thres = 10,
#                     path_results,
#                     version_name,
#                     nearby_ver,
#                     path_series_nearby,
#                     path_series_main)
# check homogeneities  ----------------------------------------------------
source(paste0(path_code_att,"restrict_based_homogeneity.R"))

remove_cluster(path_results, nb_test.ref, criterion, window.thres, nearby_ver)
check_homo(path_results,window.thres, nearby_ver, criterion)

# screening ---------------------------------------------------------------
source(paste0(path_code_att,"screening_final.R"))
screening(path_code_att, path_results, window.thres, nearby_ver)

# FGLS --------------------------------------------------------------------


# predictive rule ---------------------------------------------------------
source(paste0(path_code_att,"predictive_rule.R"))
predictive_rule(path_results, significance.level = 0.05, offset = 0, GE = 0, number.pop =3, R=10)


