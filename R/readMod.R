# #Example
# mod <- "
# Y ~ 1[s] + var(Y)[s] + res[s];
# random: 1, var(Y), res, s;
# between: var ~ x;
# x ~ res
# s==4;
# "

readMod <- function(mod){
  mod <- gsub(" ", "", mod, fixed = TRUE)
  modType <- c("ar", "var")[c(grepl("\\<ar\\>", mod),grepl("\\<var\\>", mod))]
  if(length(modType)>1){
    cat("please specify either VAR or AR effects (not both)")
  }
  markov <- grepl("\\<[s]\\>", mod)
  if(markov){
    m <- str_extract_all(str_extract_all(mod, "s==[0-9]+")[[1]][1], "[0-9]+")
    sDepPars <- str_extract_all(mod, "\\[([^]]+)\\].*")
    sDepPars <- str_split(sDepPars, "\\+")
  }
}






