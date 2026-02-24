## MOFA script for Aphids ##

######################  Set working directory #######################
dir="/Users/svztg/Library/CloudStorage/OneDrive-TheUniversityofNottingham/Digital_Research/Source_Control/Rumiana_Ray_NIV_paper/Aphid"

input_file="Input/Aphids and honeydew from exposure to niv MOFA.xlsx"

setwd(dir)
OLD_NORMALISATION=TRUE
unlink(recursive = T, "./Results")
source("../Scripts/MOFA-CORE.R");

# Copy file from ../Scripts/ to current directory
file.copy("../Scripts/MOFA-CORE.R", "./MOFA-CORE.R", overwrite = TRUE)

outfile <- paste0("../",basename(dir),"-",format(Sys.Date(), "%d-%B-%Y"),".zip")

utils::zip(zipfile = outfile, files = ".", flags = "-r9Xq")

# Remove the copied file
unlink("./MOFA-CORE.R")

