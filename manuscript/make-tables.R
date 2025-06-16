
# load relevant libraries
library(flextable)
library(ftExtra)
library(webshot)
library(here)

## make the tables functions (i.e. export to pdf)
save_flextable_as_pdf <- function(ft, pdf_path, vwidth = 1000, vheight = 300) {
  # Ensure required packages are loaded
  stopifnot(requireNamespace("flextable"),
            requireNamespace("webshot"),
            requireNamespace("here"))
  
  # Define intermediate HTML path (same dir as PDF, same basename)
  html_path <- sub("\\.pdf$", ".html", pdf_path)
  
  # Save flextable as HTML
  flextable::save_as_html(ft, path = html_path)
  
  # Render HTML to PDF via webshot
  invisible(webshot::webshot(
    url = html_path,
    file = pdf_path,
    vwidth = vwidth,
    vheight = vheight
  ))
}

## table s1
print("table-s1")

# create the table as a data.frame
tb <- data.frame(
  V1 = c(
    "Net biodiversity effect",
    "Total complementarity",
    "Total selection",
    "Non-random overyielding",
    "Total insurance",
    "Average selection",
    "Temporal insurance",
    "Spatial insurance",
    "Spatio-temporal insurance"
  ),
  V2 = c("NBE", "TC", "TS", "NO", "IT", "AS", "TI", "SI", "ST"),
  V3 = c(
    "The extent to which mixtures of species are higher functioning than expected based on their monoculture functioning and RYEs. \n\n**NBE = TS + TC**\n\nThe TC effect quantifies the extent to which multiple species at a given time and place are required to maximise ecosystem function. When the environment is homogeneous in space and time, species can still partition resources or facilitate one another, leading to positive biodiversity effects.",
    "(see above, part of NBE)",
    "**TS = NO + IT**",
    "Fraction of the TS that is due to the highest functioning species in monoculture overyielding the most in mixtures. \n\nThis effect is expected to be negative because species that are low functioning in monoculture can overyield in mixture more easily than high functioning species in monoculture.",
    "**IT = AS + TI + SI + ST**",
    "When there is a competitor that is highest functioning in monoculture and competitively dominant, mixtures are expected to be completely dominated by the superior competitor generating a positive NBE which is quantified by AS. High AS levels are expected when environmental conditions are similar in space and time.",
    "Quantifies the extent to which mixtures of species become dominated by species at times when their monoculture yields are highest. High TI levels are expected when there is considerable environmental variation through time and species differ in their responses to that variation.",
    "Quantifies the extent to which mixtures of species become dominated by species in places where their monoculture yields are highest. High SI levels are expected when there is considerable environmental variation in space and species differ in their responses to that variation.",
    "The combination of TI and SI i.e. where species become dominated by species at times and places where their monoculture yields are highest. High ST levels are expected when there is an interaction between spatial and temporal environmental variation that species respond to in different ways."
  ),
  stringsAsFactors = FALSE
)

# compose using flextable
ft <-
  flextable::flextable(tb) |>
  colformat_md(j = "V3") |>
  flextable::set_header_labels(
    V1 = "Effect",
    V2 = "Abbreviation",
    V3 = "Ecological interpretation") |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft, 
                      pdf_path = here::here("manuscript/figures/app_1_table_s1.pdf"),
                      vwidth = 1000, vheight = 300)


## table s2
print("table-s2")

# load the table data
tb <- readRDS(here::here("manuscript/figures/app_1_table_s2.rds"))

# rename the first column
names(tb)[1] <- "V1"

# create the flextable
ft <-
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  flextable::set_header_labels(
    V1 = "") |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_1_table_s2.pdf"),
                      vwidth = 1000, vheight = 300)

## table s3
print("table-s3")

# load the table
tb <- readRDS(here::here("manuscript/figures/app_1_table_s3.rds"))

# set the model labels
mod_labs <- paste0("Model ", 1:8)
mod_labs[8] <- paste0(mod_labs[8], " (Null)")
tb$model <- mod_labs

# round the values
tb <-
  tb |>
  dplyr::mutate(dplyr::across(.cols = -1, ~ round(.x, digits = 2)))

# add a parameters column
tb$pars <- c(40, 27, 18, 17, 14, 12, 11, 3)

# reorder the columns
tb <-
  tb |>
  dplyr::select(model, psis_loo, psis_loo_se, psis_loo_d, psis_loo_se_d,
                ploo, ploo_se, pars, kvals) 

# select and convert to table
ft <-
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  flextable::set_header_labels(
    model = "",
    psis_loo = "LOO",
    psis_loo_se = "LOO SE", 
    psis_loo_d = "ΔLOO",
    psis_loo_se_d = "ΔLOO SE",
    ploo = "p-LOO",
    ploo_se = "p-LOO SE",
    pars = "Pars.",
    kvals = "k > 0.5 (%)"
  ) |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# make the first row bold
ft <- flextable::bold(ft, i = 1, bold = TRUE, part = "body")

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_1_table_s3.pdf"),
                      vwidth = 1000, vheight = 300)

## table s4
print("table-s4")

# set-up the table as a data.frame
tb <- dplyr::tibble(
  "Place" = c(1, 1, 1, 1, 2, 2, 2, 2),
  "Time" = c(1, 1, 2, 2, 1, 1, 2, 2),
  "Species" = c(1, 2, 1, 2, 1, 2, 1, 2),
  "RY~E~" = rep(0.5, 8)
)

# make the flextable
ft <-
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  colformat_md(j = 4, part = "header") |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_3_table_s4.pdf"),
                      vwidth = 1000, vheight = 300)

## table s5
print("table-s5")

# set-up the table as a data.frame
tb <- dplyr::tibble(
  "Place" = c(1, 1, 1, 1, 2, 2, 2, 2),
  "Time" = c(1, 1, 2, 2, 1, 1, 2, 2),
  "Species" = c(1, 2, 1, 2, 1, 2, 1, 2),
  "RY~E~" = c(0.2, 0.8, 0.6, 0.4, 0.3, 0.7, 0.1, 0.9)
)

# make the flextable
ft <-
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  colformat_md(j = 4, part = "header") |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_3_table_s5.pdf"),
                      vwidth = 1000, vheight = 300)

## table s6
print("table-s6")

# load the table
tb <- readRDS(here::here("manuscript/figures/app_3_table_s6.rds"))

# round the values
tb <- 
  tb|>
  dplyr::mutate(dplyr::across(.cols = names(tb)[-c(1, 2)], ~ round(.x, digits = 1)))

# convert to a flextable
ft <-
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  flextable::fontsize(size = 14) |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_3_table_s6.pdf"),
                      vwidth = 1000, vheight = 300)

## table s7
print("table-s7")

# load the table
tb <- readRDS(here::here("manuscript/figures/app_3_table_s7.rds"))

# rename the columns
names(tb) <- paste0("V", 1:ncol(tb))

# round the values
tb <- 
  tb|>
  dplyr::mutate(dplyr::across(.cols = names(tb)[-c(1)], ~ round(.x, digits = 3)))

# combine the median and pi intervals
tb$Med_PI <- with(tb,
                  paste0(V3, " [", V4, " - ", V5, "]"))

# select the relevant columns
tb <-
  tb |>
  dplyr::select(V1, V2, Med_PI, V6)

# rename the columns
names(tb) <- c("V1", "Prop. within 95% PI", "Median abs. error, [PI95%]", "Pooled mean effect")

# convert to a flextable
ft <- 
  tb |>
  flextable::flextable(col_keys = names(tb)) |>
  flextable::set_header_labels(
    V1 = ""
  ) |>
  fontsize(size = 14, part = "all") |>
  flextable::autofit()

# save as a pdf
save_flextable_as_pdf(ft,
                      pdf_path = here::here("manuscript/figures/app_3_table_s7.pdf"),
                      vwidth = 1000, vheight = 300)

