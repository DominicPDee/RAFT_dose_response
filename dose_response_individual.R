dose_response_individual <- function(df1, LC_in = c(50, 99)){

  concentration_sim = seq(0,
                          sqrt(max(df1$Concentration)),
                          length.out = 100)

  N_1 <- length(concentration_sim)

  data_stan <- list(N=nrow(df1),
                    mortality=as.integer(df1$Test_mort_n),
                    tested=as.integer(df1$Test_n),
                    concentration=sqrt(as.vector(df1$Concentration)),
                    N_1=as.integer(N_1),
                    concentration_sim=as.vector(concentration_sim),
                    Y=n_distinct(df1$Covariate),
                    year_index=as.numeric(as.factor(df1$Covariate)))

  if(file.exists(paste0("outputs/fits/individual_",unique(df1$rlookup),".rds"))){
    fit <- readRDS(file=paste0("outputs/fits/individual_",unique(df1$rlookup),".rds"))

  } else {

    model <- rstan::stan_model("models/5PL_ftime_ind_NEW.stan")
    fit <- sampling(model, data=data_stan, iter=5000, chains=4) # rhat < 1.01

    # 0/ Checking convergence ----
    Rhat_df1 <- summary(fit)$summary[,10]
    modelcheck <- data.frame(not_converged=sum(Rhat_df1 > 1.01, na.rm = T),
                             is_NA=ifelse(length(names(table(is.na(Rhat_df1))))==1,"no","yes"),
                             divergent_iterations=rstan::get_num_divergent(fit),
                             mod_iter=5000)

    if(modelcheck$not_converged>0){
      fit <- sampling(model, data=data_stan, iter=10000, chains=4) # rhat < 1.01
      modelcheck <- data.frame(not_converged=sum(Rhat_df1 > 1.01, na.rm = T),
                               is_NA=ifelse(length(names(table(is.na(Rhat_df1))))==1,"no","yes"),
                               divergent_iterations=rstan::get_num_divergent(fit),
                               mod_iter=10000)
    }

    saveRDS(fit, file=paste0("outputs/fits/individual_",unique(df1$rlookup),".rds"))
    write.csv(modelcheck, file=paste0("outputs/csv/modelcheck_individual_",unique(df1$rlookup),".csv"))
  }

  # 0/ Extracting log likelihood & LOO for model comparison -----
  logLikelihood_fi <- extract_log_lik(fit, "LogLikelihood")
  LOO_fi <- loo(logLikelihood_fi)
  saveRDS(LOO_fi, file=paste0("outputs/LOO/individual_",unique(df1$rlookup),".rds"))

  # 1/ Generating simulated mortality curve ----
  mortality <- rstan::extract(fit, "mean_mortality_sim")[[1]]

  ylookup <- data.frame(Covariate=df1$Covariate,
                        ylu=as.numeric(as.factor(df1$Covariate))) %>%
    unique()

  mean_mort <- apply(mortality, c(2,3), mean) %>% # get the mean for each data point
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    reshape2::melt(id.vars="concentration") %>%
    mutate(variable=as.character(variable)) %>% # get the correct year
    mutate(variable=gsub("V","",variable)) %>%
    mutate(variable=as.numeric(variable)) %>%
    dplyr::rename(ylu=variable) %>%
    dplyr::rename(Test_mort_perc=value) %>%
    left_join(ylookup) %>%
    mutate(dat="sim")

  mean_lower <- apply(mortality, c(2,3), function(x) quantile(x, 0.025)) %>% # get the mean for each data point and year of data
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    reshape2::melt(id.vars="concentration") %>%
    mutate(variable=as.character(variable)) %>% # get the correct year
    mutate(variable=gsub("V","",variable)) %>%
    mutate(variable=as.numeric(variable)) %>%
    dplyr::rename(ylu=variable) %>%
    dplyr::rename(Test_mort_perc=value) %>%
    left_join(ylookup) %>%
    mutate(dat="sim")

  mean_upper <- apply(mortality, c(2,3), function(x) quantile(x, 0.975)) %>% # get the mean for each data point and year of data
    as.data.frame() %>%
    mutate(concentration=concentration_sim^2) %>%
    reshape2::melt(id.vars="concentration") %>%
    mutate(variable=as.character(variable)) %>% # get the correct year
    mutate(variable=gsub("V","",variable)) %>%
    mutate(variable=as.numeric(variable)) %>%
    dplyr::rename(ylu=variable) %>%
    dplyr::rename(Test_mort_perc=value) %>%
    left_join(ylookup) %>%
    mutate(dat="sim")

  mean_mort <- mean_mort %>% mutate(lower=mean_lower$Test_mort_perc,
                                    upper=mean_upper$Test_mort_perc)

  fit_mort <- tibble(Concentration =mean_mort$concentration,
                     Mortality_perc=mean_mort$Test_mort_perc*100,
                     lower=mean_mort$lower*100, upper=mean_mort$upper*100,
                     dat=mean_mort$dat,
                     Insecticide=unique(df1$Insecticide),
                     rlookup = unique(df1$rlookup),
                     Covariate=mean_mort$Covariate)

  # Combine and plot
  df1_s <- df1 %>%
    bind_rows(fit_mort)

  plot=ggplot(df1_s,
              aes(x=Concentration ,
                  y=Mortality_perc,
                  group=as.factor(Covariate))) +
    geom_point(data = filter(df1_s, is.na(dat)),
               aes(colour=as.factor(Covariate))) +
    geom_line(data = filter(df1_s, !is.na(dat)),
              aes(y=Mortality_perc,
                  colour=as.factor(Covariate))) +
    geom_ribbon(data = filter(df1_s, !is.na(dat)),
                aes(ymin=lower,
                    ymax=upper,
                    fill=as.factor(Covariate)),
                alpha=0.2) +
    theme_classic() +
    theme(panel.grid.major = element_line(colour = "grey93"),
          panel.grid.minor = element_line(colour = "grey97")) +
    ylab("Mortality (%)") +
    xlab("Insecticide concentration (%, ug/bottle, mg/m2)") +
    theme(legend.position="bottom") +
    labs(fill="Covariate", colour="Covariate") +
    scale_shape_manual(values=c(8,16)) +
    scale_x_sqrt() +
    guides(shape=guide_legend(nrow = 2, byrow = T))
  #
  # # 2/ LC ----
  lcx <- function(y,B,C,E){
    return(exp(((log((((-1)/(y-1))^(1/E))-1))/B)+C))
  }
  #
  Bfit <- rstan::extract(fit)[["B"]]
  Cfit <- rstan::extract(fit)[["C"]]
  Efit <- rstan::extract(fit)[["E"]]
  # #
  sapply(LC_in, function(a){

    message("Running LC summary ", a)

    temp = do.call(rbind,
                   sapply(1:data_stan$Y, function(x){
                     do.call(rbind,
                             sapply(1:length(Bfit), function(i){
                               data.frame("LC" = lcx(a/100, Bfit[i],Cfit[i,x],Efit[i]),
                                          "ylu" = x)
                             },simplify = FALSE))
                   }, simplify = FALSE))

    temp$LC_value <- a

    lc <- temp %>%
      left_join(ylookup) %>%
      mutate(LC=LC^2,
             Insecticide=unique(df1$Insecticide),
             rlookup=unique(df1$rlookup))

    lc_summ <- lc %>%
      group_by(Insecticide, rlookup, LC_value, Covariate) %>%
      summarise(LC_mean=mean(LC),
                LC_median=median(LC),
                LC_lower=quantile(LC,0.025),
                LC_upper=quantile(LC,0.975))

    write.csv(lc, file=paste0("outputs/csv/individual_LC", a, "_",unique(df1$rlookup),".csv"))
    write.csv(lc_summ, file=paste0("outputs/csv/individual_summ_LC", a, "_", unique(df1$rlookup),".csv"))

    NULL
  }, simplify = FALSE)

  # ## 2.2/ Concentrations at which mosquitoes die ----
  # modelfun <- function(x,B,C,E){
  #   return(1 -  (1 /((1+exp(B*(log(x)-C)))^E)))
  # }
  #
  # lower <- log10(1e-6)
  # upper <- log10(max(df1$Concentration))
  # x <- seq(lower,upper, length.out=1000)
  # x <- 10^x
  # x <- sqrt(x)
  #
  # getCDF = do.call(rbind, sapply(1:length(x), function(a){
  #   do.call(rbind,sapply(1:data_stan$Y, function(b){
  #     data.frame(sim_y = mean(sapply(1:length(Bfit), function(i){
  #       modelfun(x[a], Bfit[i],Cfit[i,b],Efit[i])
  #     })),
  #     ylu = b,
  #     sim_x=x[a])
  #   }, simplify = F))},simplify = F))
  #
  # getCDF <- getCDF %>%
  #   mutate(sim_x2=sim_x^2)
  #
  # # ggplot(getCDF, aes(x=sim_x2, y=sim_y,
  # #                    colour=as.factor(ylu)))+
  # #   geom_point() +
  # #   scale_x_sqrt()
  #
  # lgetCDF=sapply(1:data_stan$Y, function(i){
  #   getCDF %>% subset(ylu==i)
  # },simplify = F)
  #
  # mortx <- function(x,y){
  #   f <- approxfun(y,x)
  #   u <- runif(20000) # random deviates from uniform distribution
  #   lcdplot = data.frame("Mort_C"=f(u))
  # }
  #
  # getPDF=do.call(rbind,sapply(1:data_stan$Y, function(i){
  #   data.frame(Mortx=mortx(lgetCDF[[i]]$sim_x2,lgetCDF[[i]]$sim_y),
  #              ylu=i)
  # },simplify = F)) %>% left_join(ylookup)
  #
  # mortx_plot = ggplot(getPDF,aes(x=Mort_C)) +
  #   geom_density(aes(fill=as.factor(Covariate),
  #                    colour=as.factor(Covariate)),
  #                alpha=.5,position="identity") +
  #   labs(fill="Covariate",colour="Covariate",
  #        x=expression("Mortality at"~italic(x)~"concentration"), y="Density") +
  #   theme_classic() +
  #   ggtitle(paste0("Concentrations at which mosquitoes die (",unique(df1$rlookup),")"))
  #
  # # 2.3/ Background mortality estimate ----
  # BM = data.frame(mean_BM=mean(rstan::extract(fit)[["phi"]]),
  #                 median_BM=median(rstan::extract(fit)[["phi"]]),
  #                 lower_BM=quantile((rstan::extract(fit)[["phi"]]),0.025),
  #                 upper_BM=quantile((rstan::extract(fit)[["phi"]]), 0.975),
  #                 Village=unique(df1$Village),
  #                 Treatment="Deltamethrin")

  # 3/ Variability/predictive error ----
  ## 3.1/ Mortality variability (along y axis)

  ### STEPS:
  ### a) find all relevant x values (actual x values)
  ### b) extract actual y values for each of these x values
  ### c) generate simulated y value for each of these x values
  ### d) compute difference between actual y and simulated y (FOR EACH ITERATION)
  ### e) average these differences (FOR EACH ITERATION)
  ### f) compute mean, median and 95% CIs of variability estimate from all iterations

  # sim_y <- rstan::extract(fit,"mean_mortality_bm")[[1]] %>%
  #   as.data.frame()
  # act_x_year <- data.frame(act_x=(data_stan$concentration)^2,
  #                          ylu=data_stan$year_index)
  # act_y <- as.integer(df1$Mortality_perc)/100
  # n  <- nrow(df1)
  #
  # break_lookup <- act_x_year %>%
  #   rowid_to_column() %>%
  #   group_by(ylu) %>%
  #   mutate(breakp=max(rowid)) %>%
  #   pull(breakp) %>% unique()
  # begin=c(1,break_lookup+1)[1:length(break_lookup)]
  # end=break_lookup
  #
  # sim_y_list=sapply(seq_along(break_lookup), function(i){
  #   sim_y[,begin[i]:end[i]]
  # }, simplify = F)
  #
  # act_y_list <- sapply(seq_along(break_lookup), function(i){
  #   act_y[begin[i]:end[i]]
  # }, simplify = F)
  #
  # diff_it <- do.call(rbind,sapply(1:nrow(sim_y), function(i){
  #   do.call(rbind,sapply(1:length(sim_y_list), function(x){
  #     matrix(sum((1/(sqrt(sim_y_list[[x]][i,]*(1-sim_y_list[[x]][i,]))))*abs(act_y_list[[x]]-sim_y_list[[x]][i,]))/length(sim_y_list[[x]])) %>%
  #       cbind(it=i,
  #             ylu=x)
  #     #t(abs((sim_y[i,]*100) - act_y))
  #   }, simplify = F))
  # },simplify = F)) %>%
  #   as.data.frame() %>%
  #   left_join(ylookup)
  #
  # var_mort_new <- diff_it %>%
  #   group_by(Year) %>%
  #   summarise(mean_vmort=mean(V1)*100,
  #             median_vmort=median(V1)*100,
  #             lower_vmort=quantile(V1,0.025)*100,
  #             upper_vmort=quantile(V1,0.975)*100) %>%
  #   mutate(Village=unique(df1$Village),
  #          Treatment="Deltamethrin")

  # sim_y <- rstan::extract(fit,"mean_mortality_bm")[[1]] %>%
  #   as.data.frame()
  # act_x_year <- data.frame(act_x=(data_stan$concentration)^2,
  #                          ylu=data_stan$year_index)
  # act_y <- as.integer(df1$Mortality_perc)
  # diff_it=do.call(cbind,sapply(1:nrow(sim_y), function(i){
  #   t(abs((sim_y[i,]*100) - act_y))
  # },simplify = F)) %>%
  #   cbind(ylu=data_stan$year_index) %>%
  #   as.data.frame()
  #
  # varmort_it=do.call(rbind,sapply(1:data_stan$Y,
  #                   function(i){
  #                     data.frame(varmort=apply(diff_it[diff_it$ylu==i,c(1:length(Bfit))],2,sum)/nrow(act_x_year[act_x_year$ylu==i,]),
  #                                ylu=i)
  #                   },simplify = F)) %>%
  #   left_join(ylookup)
  #
  #   var_mort <- varmort_it %>% group_by(Year) %>%
  #     summarise(mean_vmort=mean(varmort),
  #               median_vmort=median(varmort),
  #               lower_vmort=quantile(varmort,0.025),
  #               upper_vmort=quantile(varmort,0.975)) %>%
  #     mutate(Village=unique(df1$Village),
  #            Treatment="Deltamethrin")

  # # 4/ Model fitting ----
  # ## 4.1/ Residuals plots ----
  # mean_mortality_bm <- rstan::extract(fit, "mean_mortality_bm")[[1]]
  #
  # df1_residuals <- apply(mean_mortality_bm,2,median) %>%
  #   melt(value.name = "sim_y") %>%
  #   bind_cols(Concentration=data_stan$concentration,
  #             act_y=df1$Mortality_perc/100,
  #             ylu=data_stan$year_index) %>%
  #   mutate(Residuals=act_y-sim_y) %>%
  #   left_join(ylookup)
  #
  # resplot = ggplot(df1_residuals,aes(x=Concentration, y=Residuals,
  #                                    colour=as.factor(Year))) +
  #   geom_point() +
  #   geom_hline(aes(yintercept=0), linetype="dashed") +
  #   ggtitle(paste0("Residuals for ",unique(df1$Village))) + labs(colour="Year") +
  #   geom_smooth(se=F)
  #
  # # ## 4.2/ Actual vs simulated ----
  # mort_sim <- df1 %>%
  #   bind_cols(Mortality_sim_perc=apply(mean_mortality_bm,2,median)*100)
  #
  # ### 4.2.1/ By year
  # mod_lm_y <- mort_sim %>%
  #   group_by(Village,Year) %>%
  #   do(mod = lm(Mortality_perc~Mortality_sim_perc, data=.))
  # df_coeff_y <- mod_lm_y %>%
  #   do(data.frame(
  #     Year=.$Year,
  #     Village=.$Village,
  #     var=names(coef(.$mod)),
  #     coef(summary(.$mod)),
  #     r2=summary(.$mod)$r.square,
  #     RMSE=sqrt(mean(.$mod$residuals^2))))
  #
  # ap_year <- ggplot(mort_sim,aes(x=Mortality_sim_perc, y=Mortality_perc, colour=as.factor(Year))) +
  #   geom_point(size=4) +
  #   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  #   ylab("Actual mortality (%)") +
  #   xlab("Predicted mortality (%)") +
  #   theme_bw() +
  #   ggtitle(paste0("Field model predictive accuracy (",unique(df1$Village),")")) +
  #   labs(colour="Year") +
  #   xlim(c(0,100)) + ylim(c(0,100))
  #
  # ### 4.2.2/ All years
  # mod_lm <- mort_sim %>%
  #   group_by(Village) %>%
  #   do(mod = lm(Mortality_perc~Mortality_sim_perc, data=.))
  # df_coeff <- mod_lm %>%
  #   do(data.frame(
  #     Village=.$Village,
  #     var=names(coef(.$mod)),
  #     coef(summary(.$mod)),
  #     r2=summary(.$mod)$r.square,
  #     RMSE=sqrt(mean(.$mod$residuals^2))))
  #
  # ap_all <- ggplot(mort_sim,aes(x=Mortality_sim_perc, y=Mortality_perc)) +
  #   geom_point(size=4) +
  #   geom_abline(intercept = 0, slope = 1, linetype="dashed") +
  #   ylab("Actual mortality (%)") +
  #   xlab("Predicted mortality (%)") +
  #   theme_bw() +
  #   ggtitle(paste0("Field model predictive accuracy (",unique(df1$Village),")")) +
  #   geom_text(data=df_coeff,
  #             aes(x=10,y=75,
  #                 label=paste0('R2=',round(r2,digits=2),"\nRMSE=",round(RMSE,digits=0),"%"),
  #                 fontface=3),
  #             size=5) +
  #   xlim(c(0,100)) + ylim(c(0,100))
  #
  # # 5/ Save objects----
  write.csv(df1_s, file=paste0("outputs/csv/df1_s_individual_",unique(df1$rlookup),".csv"))
  ggsave(plot,file=paste0("outputs/plots/individual_",unique(df1$rlookup),".png"),
         width = 7,height = 5)

  # write.csv(getPDF,file=paste0("outputs/csv/individual_PDF_",unique(df1$rlookup),"_CORRECTED.csv"))
  # ggsave(mortx_plot,file=paste0("outputs/plots/individual_mortxplot_",unique(df1$rlookup),"_CORRECTED.png"),
  #        width = 4.85,height = 4.39)
  # write.csv(BM,file=paste0("field/outputs/csv/paper/time_ind_BM_",unique(df1$Village),".csv"))
  # write.csv(var_mort,file=paste0("field/outputs/csv/paper/time_ind_var_mort_",unique(df1$Village),".csv"))
  # write.csv(var_mort_new,file=paste0("field/outputs/csv/paper/time_ind_var_mort_new_",unique(df1$Village),".csv"))
  #
  # ggsave(resplot,file=paste0("field/outputs/plots/paper/time_ind_residuals_",unique(df1$Village),".png"),
  #        width = 9.5,height = 4.39)
  # ggsave(ap_year,file=paste0("field/outputs/plots/paper/time_ind_AP_year_",unique(df1$Village),".png"),
  #        width = 5,height = 4.39)
  # write.csv(df_coeff_y,file=paste0("field/outputs/csv/paper/time_ind_APy-coeff_",unique(df1$Village),".csv"))
  # ggsave(ap_all,file=paste0("field/outputs/plots/paper/time_ind_AP_",unique(df1$Village),".png"),
  #        width = 4.85,height = 4.39)
  # write.csv(df_coeff,file=paste0("field/outputs/csv/paper/time_ind_AP-coeff_",unique(df1$Village),".csv"))
  # write.csv(df1_residuals,file=paste0("field/outputs/csv/paper/time_ind_residuals_",unique(df1$Village),".csv"))
  # write.csv(mort_sim,file=paste0("field/outputs/csv/paper/time_ind_mortsim_",unique(df1$Village),".csv"))

}
