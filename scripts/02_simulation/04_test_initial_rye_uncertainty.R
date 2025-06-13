#'
#' @title: Simulate metacommunities to test the analytical pipeline
#' 
#' @description: Calculates biodiversity effects when initial relative abundance data is missing
#' 
#' @details: This script uses a distribution of starting relative abundances from 
#' the Dirichlet distribution to generate distributions of biodiversity effects.
#' The code takes around 20 minutes to run on a regular desktop computer. 
#' 
#' @authors: James G. Hagan (james_hagan(at)outlook.com)
#' 

# load the relevant libraries
library(dplyr)
library(ggplot2)

# load the relevant functions
source("scripts/01_partition_functions/01_isbell_2018_partition.R")
source("scripts/Function_plotting_theme.R")

# read in the data
mc_sims <- readRDS(file = "results/mc_sim_list.rds")
start_RA <- readRDS(file = "results/mc_sim_list_start_RA.rds")

# check the length of the mc_sims
length(mc_sims)

# check an example dataset
# n <- sample(1:length(mc_sims), 1)
n <- 1437

# plot the mixtures
p1 <- 
  mc_sims[[n]]$mc_dat %>%
  mutate(Species = as.character(species)) %>%
  mutate(Patch = paste0("Patch ", place)) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = Y, colour = Species)) +
  geom_line() +
  scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", n = 3)) +
  facet_wrap(~Patch, scales = "free") +
  ylab("Mixture population size (N)") +
  xlab("Time point") +
  theme_meta() +
  theme(legend.position = "right",
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10),
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.key = element_rect(fill = NA))
plot(p1)

ggsave(filename = "manuscript/figures/app_3_fig_s9.png", p1, dpi = 600,
       unit = "cm", width = 20, height = 15)

# plot the monocultures
p2 <- 
  mc_sims[[n]]$mc_dat %>%
  mutate(Species = as.character(species)) %>%
  mutate(Patch = paste0("Patch ", place)) %>%
  ggplot(data = .,
         mapping = aes(x = time, y = M, colour = Species)) +
  geom_line() +
  scale_colour_manual(values = wesanderson::wes_palette("Darjeeling1", n = 3)) +
  facet_wrap(~Patch, scales = "free") +
  ylab("Monoculture population size (N)") +
  xlab("Time point") +
  theme_meta() +
  theme(legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(size = 11),
        legend.key = element_rect(fill = NA))
plot(p2)

ggsave(filename = "figures/app_3_fig_s10.png", p2, dpi = 600,
       unit = "cm", width = 20, height = 15)


# calculate the biodiversity incorporating uncertainty in RYe

# generate output list
output_list <- vector("list", length = length(mc_sims))

# loop over each each simulated dataset and over each sample from the Dirichlet
for(i in 1:length(mc_sims)) {

  # iterate over all possible initial RYE values
  RYE_rep <- 
    
    lapply(start_RA, function(RA) {
      
      RYE <- vector("list", length = nrow(RA))
      for(j in 1:nrow(RA)) { RYE[[j]] <- RA[j,] }
      
      # calculate the BEF effects
      BEF_post <- isbell_2018_part(data = mc_sims[[i]]$mc_dat, RYe = RYE)
      names(BEF_post[["L.Beff"]])[names(BEF_post[["L.Beff"]]) == "L.Beff"] <- "Beff"
      
      # combine the general biodiversity effects and local effects into one data.frame
      BEF_post <- rbind(BEF_post[["Beff"]], BEF_post[["L.Beff"]])
      
      # convert to a data.frame
      BEF_post <- as.data.frame(BEF_post, row.names = NULL)
      
      return(BEF_post)
      
    })
  
  
  # bind into a data.fram
  RYE_rep <- dplyr::bind_rows(RYE_rep, .id = "RYE")
  
  # summarise these data
  output <- 
    RYE_rep |>
    dplyr::group_by(Beff) |>
    dplyr::summarise(PI95_low = quantile(Value, 0.025),
                     PI95_high = quantile(Value, 0.975),
                     mean_BEF = mean(Value, na.rm = TRUE), .groups = "drop"
                     )
  
  # add the observed values
  names(mc_sims[[i]][["BEF_obs"]]) <- c("Beff", "BEF_obs")
  
  # join these datasets
  output <- dplyr::full_join(output, mc_sims[[i]][["BEF_obs"]], by = "Beff")
  
  # write output to the list
  output_list[[i]] <- output
  
}

# bind the output list into a data.frame
output_df <- bind_rows(output_list, .id = "sim_rep")

# save this as an RDS file
saveRDS(object = output_df, file = "results/mc_sim_test.rds")

### END
