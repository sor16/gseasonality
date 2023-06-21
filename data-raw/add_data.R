diagnosis_data <- read.csv('data-raw/sim_dat.csv',header=T)
ICD_codes <- read.table('data-raw/ICD_codes.tsv',header=T,sep = '\t')

usethis::use_data(diagnosis_data, overwrite = TRUE)
usethis::use_data(ICD_codes, overwrite = TRUE)
