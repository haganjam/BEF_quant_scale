
# 5. calculates Isbell et al.'s (2018) biodiversity effects using all the uncertainty present in the
# - modelled monocultures by drawing different samples from the posterior distribution
# - using starting relative abundances from a Dirichlet distribution

# which objects do I need in memory to run this?

# MC1_pred object

# MC1_NA object

# dr object



# fill in the monocultures with one sample from the posterior
# for that one sample from the posterior, calculate biodiversity effects with 100 samples from the Dirichlet distribution
Mono_reps <- 
  apply(MC1_pred[sample(x = 1:nrow(MC1_pred), 100), ], 1, function(x) {
    
    # fill in the missing monoculture data with one sample from the posterior
    df <- MC1_NA
    df[is.na(df$M1), ]$M1 <- x
    df <- df[,c("sample", "time", "place", "species", "M1", "Y")]
    df <- rename(df, M = M1)
    
    RYe_reps <- 
      apply(dr, 2, function(z) {
        
        a <- Isbell_2018_sampler(data = df, RYe = z, RYe_post = FALSE)
        return(bind_rows(a$Beff, rename(a$L.Beff, Beff = L.Beff)))
        
      } )
    
    return(RYe_reps)
    
  } )

df_unc <- bind_rows(Mono_reps, .id = "ID")



