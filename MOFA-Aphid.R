## MOFA script for Aphids ##

input_file="Input/Aphids and honeydew from exposure to niv MOFA.xlsx"

OLD_NORMALISATION=TRUE
unlink(recursive = T, "./Results")
source("./MOFA-CORE.R");

outfile <- paste0("../",basename(dir),"-",format(Sys.Date(), "%d-%B-%Y"),".zip")

utils::zip(zipfile = outfile, files = ".", flags = "-r9Xq")

# Remove the copied file
unlink("./MOFA-CORE.R")