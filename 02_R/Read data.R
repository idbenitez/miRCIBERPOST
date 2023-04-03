library(dplyr)
library(labelled)
# Read.data ---- 
dat <- haven::read_sav("07_Data/BBDD_05_09_22.sav",encoding ="UTF-8")
dat <- as.data.frame(unlabelled(dat))
# New variables and categorizations ----
is <- function(x){(is.character(x) | is.factor(x))}
# dat <- dat %>%  mutate_if(is,function(x){iconv(as.character(x),from = "UTF-8", to="ASCII//TRANSLIT")})

# Save data in Rdata 
save(dat,file = "07_Data/BBDD_05_09_22.RData")



dat <- foreign::read.spss("07_Data/BBDD_05_09_22.sav", to.data.frame = TRUE) 
dat <- 