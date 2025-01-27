# Master script

library(here)

here::i_am("code/01_run_render.R")

source(here("code/01.1_clinical.data.processing.R"))

source(here("code/01.1_cytokine.data.processing.R"))

source(here("code/01.2_cytof.data.processing.R"))
