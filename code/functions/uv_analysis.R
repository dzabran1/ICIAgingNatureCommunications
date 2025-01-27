uv_analysis <- function(data, cyt_imm_var, y_var, time_var="BL_cyt") {
  
  ### data should be in a long format stacked by cyt_imm_var ###
  y_col <- seq(1, dim(data)[2])[names(data) == y_var]
  cyt_imm_col <- seq(1, dim(data)[2])[names(data) == cyt_imm_var]
  
  names(data)[y_col] <- "y"
  names(data)[cyt_imm_col] <- "cyt_imm"
  
  cyt_imm_list <- unique(data[,cyt_imm_col])$cyt_imm

  # output for age coefficients
  age.uv_opt <- data.frame(
    time_text = rep(time_var, length(cyt_imm_list)),
    cyt_imm = c(cyt_imm_list),
   
    uv_coef_age = rep(NA, length(cyt_imm_list)),
    uv_low_age = rep(NA, length(cyt_imm_list)),
    uv_up_age = rep(NA, length(cyt_imm_list)),
    uv_p_age = rep(NA, length(cyt_imm_list)),
    
    uv_coef_age_cat = rep(NA, length(cyt_imm_list)),
    uv_low_age_cat = rep(NA, length(cyt_imm_list)),
    uv_up_age_cat = rep(NA, length(cyt_imm_list)),
    uv_p_age_cat = rep(NA, length(cyt_imm_list))
    
  )
  
  # Loop over each immune cell group
  for (i in 1:length(cyt_imm_list)) {
    temp_dat <- data %>% filter(cyt_imm == cyt_imm_list[i])
    
   
    #Uivariate Models
    uv_mod1 <- lm(y ~ age.per10, data = temp_dat)
    uv_mod2 <- lm(y ~ age65_status, data = temp_dat)
    
    # Extract Uivariate results for age
    age.uv_opt$uv_coef_age[i] <- round(summary(uv_mod1)$coef["age.per10", "Estimate"], 2)
    age.uv_opt$uv_low_age[i] <- round(confint(uv_mod1)["age.per10", "2.5 %"], 2)
    age.uv_opt$uv_up_age[i] <- round(confint(uv_mod1)["age.per10", "97.5 %"], 2)
    age.uv_opt$uv_p_age[i] <- round(summary(uv_mod1)$coef["age.per10", "Pr(>|t|)"], 4)
    
    
    age.uv_opt$uv_coef_age_cat[i] <- round(summary(uv_mod2)$coef["age65_statusyes", "Estimate"], 2)
    age.uv_opt$uv_low_age_cat[i] <- round(confint(uv_mod2)["age65_statusyes", "2.5 %"], 2)
    age.uv_opt$uv_up_age_cat[i] <- round(confint(uv_mod2)["age65_statusyes", "97.5 %"], 2)
    age.uv_opt$uv_p_age_cat[i] <- round(summary(uv_mod2)$coef["age65_statusyes", "Pr(>|t|)"], 4)
    
    } ### end for-loop i
  
  return(age.uv_opt)
}
