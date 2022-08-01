
# 6. summarises the distribution for each biodiversity effect and compares it to the observed value
# - we assume that for an effect to be meaningful, the 95% HPDI interval of an effect must be
# within the interval of obs value +- 0.5*observed value. This means that the 95% HPDI interval
# contains the true value within a reasonable distance of the true value

# calculate accuracy metrics

# what percentage of the posterior falls between the 50% mean +-

# set the percentage threshold
thresh <- 0.5

df_unc_sum <- 
  
  full_join(
    
    df_unc %>% 
      group_by(Beff) %>%
      summarise(mu = mean(Value, na.rm = TRUE),
                HPDI_low = HPDI(Value, prob = 0.95)[1],
                HPDI_high = HPDI(Value, prob = 0.95)[2],),
    
    rename(BEF_obs, Value_obs = Value ), by = "Beff"
    
  ) %>%
  mutate(Value_obs_min = (1-thresh)*Value_obs,
         Value_obs_max = (1+thresh)*Value_obs ) %>%
  mutate(HPDI_true = if_else( (HPDI_low > Value_obs_min) & (HPDI_high < Value_obs_max), TRUE, FALSE ) ) %>%
  mutate(mu_deviation = abs((abs(mu - Value_obs)/Value_obs)*100) ) %>%
  mutate(mu_threshold = thresh)

# check if the different effects are consistent with the NBE
t1 <- sum(df_unc_sum[df_unc_sum$Beff %in% c("AS", "TC", "NO", "TI", "SI", "ST"), ][["mu"]])
t2 <- df_unc_sum[df_unc_sum$Beff %in% c("NBE"), ][["mu"]]

assertthat::assert_that(assertthat::see_if(near(t1, t2)),
                        msg = "biodiversity effects do not sum to NBD")

# add covariates regarding the parameters
df_unc_sum$start_abun <- paste(start_abun, collapse = "_")
df_unc_sum$inter_comp <- mean(comp_mat[lower.tri(comp_mat)])
df_unc_sum$optima <- paste(sp_att$optima, collapse = "_")
df_unc_sum$niche_breadth <- paste(sp_att$env_niche_breadth, collapse = "_")  
df_unc_sum$t_steps <- paste(t_sel, collapse = "_")
df_unc_sum$dispersal <- dispersal
df_unc_sum$mono_cor <- cor(mu_m1, MC.x.NA$M[is.na(MC.x.NA$M1)])

# definition of a bad monoculture
df_unc_sum$bad_mono_n <- sum( apply(PI_m1, 2, diff)/mean(mu_m1) > 1 )/length(mu_m1) 

# reorder the columns

# important: When HPDI_true is TRUE, it means that the 95% HPDI interval is within
# thresh*obs value +- obs value which means that 95% HPDI interval is narrow enough
# to be meaningful
df_unc_sum <- 
  df_unc_sum %>%
  select(t_steps, dispersal, start_abun, optima, niche_breadth, inter_comp, mono_cor, bad_mono_n,
         Beff, Value_obs, mu, mu_deviation,
         starts_with("HPDI"), mu_threshold,
         HPDI_true)

return(list(bind_rows(MC.x), df_unc_sum) )

} )
