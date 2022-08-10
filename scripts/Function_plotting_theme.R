
#' @title theme_meta
#' 
#' @description A customised plotting theme to equalise the formatting of all plots plotted with ggplot2
#' 
#' @details This script generates a customised ggplot theme and a function to output
#' the colours for the different types of biodiversity so that they remain consistent
#' throughout the publication
#' 
#' @author James G. Hagan (james_hagan(at)outlook.com)

# customised plotting theme
theme_meta <- function(base_size = 12, base_family = "") {
    
  theme(panel.background = element_rect(fill = "white"), 
          panel.border = element_rect(fill="NA", color="black", size=0.75, linetype="solid"),
          axis.line.x = element_line(color="black", size = 0.2),
          axis.line.y = element_line(color="black", size = 0.2),
          panel.grid.major =  element_blank(),
          panel.grid.minor =  element_blank(),
          axis.ticks.length = unit(-0.16, "cm"),
          axis.title.x = element_text(colour ="black", size = 12, face = "plain", margin=margin(5,0,0,0,"pt")),
          axis.title.y = element_text(colour = "black", size = 12, face = "plain", margin=margin(0,5,0,0,"pt")),
          axis.text.x = element_text(colour = "black", size=11, face = "plain",  margin=margin(10,0,0,0,"pt")),
          axis.text.y = element_text(colour ="black", size=11, face = "plain", margin=margin(0,10,0,0,"pt")),
          axis.ticks.x = element_line(colour = "black", size = 0.4),
          axis.ticks.y = element_line(colour = "black", size = 0.4),
          legend.text = element_text(colour = "black", size=10, face = "plain"),
          legend.title = element_text(colour = "black", size=10, face = "plain"))
  
}

# colour palette function

# this requires the viridis library
require(viridis)

v_col_BEF <- function(eff_in = c("LC", "LS", "TC", "TS", "NBE", "NO", "IT", "AS", "TI", "SI", "ST")) {
  
  # set-up the colour scale
  v.col <- c(c("black", "brown"), 
             viridis(option = "C", n = 4, alpha = 1, begin = 0, end = 0.4),
             viridis(option = "C", n = 5, alpha = 1, begin = 0.5, end = 1))
  
  # name the different colours based on the biodiversity effect
  names(v.col) <- c("NBE", "NO", "LC", "LS", "TC", "TS", "IT", "AS", "TI", "SI", "ST")
  
  # get the correct order
  v.col.sel <- v.col[ names(v.col) %in% eff_in ]
  v.col.sel <- v.col.sel[order(match(names(v.col.sel) , eff_in))]
  
  return(v.col.sel)
  
}

### END
