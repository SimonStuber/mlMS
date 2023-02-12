#Example
mod <- "
Level 1:
c(y1, y2, y3) ~ 1[s] + var(c(y1,y2,y3))[s] + res[s];

Level 2:
1, var(c(y1, y2, y3)), res, s
x ~ res[s]
x2 ~ 1[s]
s ~ x2;

s==4;
"

readMod <- function(mod){
  mod <- gsub(" ", "", mod, fixed = TRUE)
  sel <- if(any(c(grepl("\\<ar\\>", mod),grepl("\\<var\\>", mod)))){
        which(c(grepl("\\<ar\\>", mod),
              grepl("\\<var\\>", mod)))}else{
         3}
  modType <- c("ar", "var", "hmm")[sel]
  if(length(modType)>1){
   cat("AR/VAR effects do not match between levels")
  }
  markov <- grepl("\\<[s]\\>", mod)
  if(markov){
    m <- str_extract_all(str_extract_all(mod, "s==[0-9]+")[[1]][1], "[0-9]+")
    if(is.na(m)){
      cat("Please specify the number of discrete latent states with s=='number'")
    }else if(m==0){
      cat("Number of latent states is 0. Please remove the latent states variable
          or specify s==>1")
    }else if(m==1){
      cat("Number of latent states is 1. Please remove the latent states variable
          or specify s==>1")
    }
    sDepPars <- strsplit(mod, c("[\\+\\;\\:]"))
    sDepPars <- str_extract_all(sDepPars, "\\s")
  }
}






