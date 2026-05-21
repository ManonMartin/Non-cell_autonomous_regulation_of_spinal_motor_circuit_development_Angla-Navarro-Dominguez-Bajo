rm(list = ls())

# Packages that need to be installed and loaded
library(tidyverse)
library(pander)
library(readxl)
library(npmv)


# Import the data

mutation <- "Ntf3"
Ntf3Files <- list.files("data/results_density2D/results_Ntf3", full.names = TRUE)
Ntf3Files <- Ntf3Files |>
  as.data.frame(Ntf3Files) |>
  dplyr::rename(Files = Ntf3Files) |>
  mutate(mutation = "Ntf3")


mutation <- "OC"
OCFiles <- list.files("data/results_density2D/results_OC", full.names = TRUE)
OCFiles <- OCFiles |>
  as.data.frame(OCFiles) |>
  dplyr::rename(Files = OCFiles) |>
  mutate(mutation = "OC")

allFiles <- rbind(Ntf3Files, OCFiles)


dat_list_all <- map(1:nrow(allFiles), function(i){
  df_Brach <- readxl::read_xlsx(path = allFiles$Files[i], sheet = "Brach")
  df_Lumb <- readxl::read_xlsx(path = allFiles$Files[i], sheet = "Lumb")
  df_Thor <- readxl::read_xlsx(path = allFiles$Files[i], sheet = "Thor")
  df_Brach$mutation <- df_Lumb$mutation <- df_Thor$mutation <-
    allFiles$mutation[i]
  return(list(Brach = df_Brach, Lumb = df_Lumb, Thor = df_Thor))
})


# function for multivariate tests

multivariate_tests <- function(dat){
  
  ## Non parametric approach ===========================================
  
  # Nonparametric Comparison of Multivariate Samples
  
  res <- npmv::nonpartest(ML|DV~cond,dat, 
                          plots=FALSE, permreps = 10000)
  pval_ANOVANP <- res$results["ANOVA type test p-value",
                              "Permutation Test p-value"]
  
  final_res <- c(
    ANOVA_NP = pval_ANOVANP)
  
  return(final_res)
}

# run multivariate_tests

res_all <- map(dat_list_all,
               function(x) {
                 res <- map(x, multivariate_tests)
                 return(as.data.frame(res))
               }
)

options(digits=22)


# reshape and export the results

nam <- sub("-distances\\.xlsx","",basename(allFiles$Files)) 
nam_mut <- paste(nam, allFiles$mutation, sep = "_")


res_multivariate_tests <- list_rbind(res_all)
rownames(res_multivariate_tests) <- nam_mut


res_multivariate_tests_adjusted <- res_multivariate_tests

for (i in 1:nrow(res_multivariate_tests)){
  res_multivariate_tests_adjusted[i,] <- p.adjust(res_multivariate_tests[i,], 
                                                  method = "holm")
  
}

res_multivariate_tests_adjusted$file <- nam
res_multivariate_tests_adjusted$mutation <- allFiles$mutation


for(i in c("Brach", "Lumb", "Thor")){
  res_multivariate_tests_adjusted[,i] <- format.pval(res_multivariate_tests_adjusted[,i])
}

pander(res_multivariate_tests_adjusted, 
       caption = "p-values non parametric ANOVA like test")

writexl::write_xlsx(x = res_multivariate_tests_adjusted, 
                    path =  paste0("results/S01/multivariate_tests.xlsx"))