library(multcomp)

response_inter_age_mva_analysis <- function(data, cyt_imm_var, y_var, time_var="Y_BL_cyt") {
  
  ### data should be in a long format stacked by cyt_imm_var ###
  y_col <- seq(1, dim(data)[2])[names(data) == y_var]
  cyt_imm_col <- seq(1, dim(data)[2])[names(data) == cyt_imm_var]
  
  names(data)[y_col] <- "y"
  names(data)[cyt_imm_col] <- "cyt_imm"
  
  cyt_imm_list <- unique(data[,cyt_imm_col])$cyt_imm

  
  # output for response coefficients
  resp.mva_opt <- data.frame(
    time_text = rep(time_var, length(cyt_imm_list)),
    cyt_imm = c(cyt_imm_list),
    
    p.inter = rep(NA, length(cyt_imm_list)), ### interaction term p-value
    
    mv_coef_resp_y = rep(NA, length(cyt_imm_list)),
    mv_low_resp_y = rep(NA, length(cyt_imm_list)),
    mv_up_resp_y = rep(NA, length(cyt_imm_list)),
    mv_p_resp_y = rep(NA, length(cyt_imm_list)),
    
    mv_coef_resp_o = rep(NA, length(cyt_imm_list)),
    mv_low_resp_o = rep(NA, length(cyt_imm_list)),
    mv_up_resp_o = rep(NA, length(cyt_imm_list)),
    mv_p_resp_o = rep(NA, length(cyt_imm_list))
    
  )
  
  # Loop over each immune cell group
  for (i in 1:length(cyt_imm_list)) {
    temp_dat <- data %>% filter(cyt_imm == cyt_imm_list[i])
    
    
    #Multivariate Models
    mv_mod <- lm(y ~ orr*age65_status + pr_trt + I.gp_gi_gu_else, data = temp_dat)
    
    # Extract Multivariate results for response vs non-response
    resp.mva_opt$p.inter[i] <- round(summary(mv_mod)$coef["orr:age65_statusyes", "Pr(>|t|)"], 4)
    
    
    ### In young group, response vs. non-response
    resp.mva_opt$mv_coef_resp_y[i] <- round(summary(mv_mod)$coef["orr", "Estimate"], 2)
    resp.mva_opt$mv_low_resp_y[i] <- round(confint(mv_mod)["orr", "2.5 %"], 2)
    resp.mva_opt$mv_up_resp_y[i] <- round(confint(mv_mod)["orr", "97.5 %"], 2)
    resp.mva_opt$mv_p_resp_y[i] <- round(summary(mv_mod)$coef["orr", "Pr(>|t|)"], 4)
    
    ### In aged group, response vs. non-response
    K_uv <- matrix(c(0, 1, 0,0,0,0,1),1)
    t_uv <- glht(mv_mod, linfct = K_uv)
    
    resp.mva_opt$mv_coef_resp_o[i] <- round(tidy(confint(t_uv))$estimate[1],2)
    resp.mva_opt$mv_low_resp_o[i] <- round(tidy(confint(t_uv))$conf.low[1],2)
    resp.mva_opt$mv_up_resp_o[i] <- round(tidy(confint(t_uv))$conf.high[1],2)
    resp.mva_opt$mv_p_resp_o[i] <- round(tidy(t_uv)$adj.p.value[1],4) ### only 1 comp, adj dose not matter
    
    
  } ### end for-loop i
  
  return(resp.mva_opt)
}