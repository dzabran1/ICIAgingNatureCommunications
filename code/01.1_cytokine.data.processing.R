
library(here)
library(tidyverse)
library(openxlsx)
library(readxl)

# Set up

here::i_am("code/01.1_cytokine.data.processing.R")

output.path <- here("output")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

output.path.rds <- here("output/RDS files")
if (!dir.exists(here(output.path.rds))) dir.create(here(output.path.rds))

output.path.excel <- here("output/Excel files")
if (!dir.exists(here(output.path.excel))) dir.create(here(output.path.excel))

options(scipen =999)

# Read in clinical data from rds file
clinical.max.irae <- readRDS(here("output/RDS files/clinical.rds"))

### Import cytokine data

file_path_cytokine <- here("data/processed cytokine data/de.id_total_cytokine_aging_analysis.xlsx")

cytokine.import <- read_excel(file_path_cytokine)


# Cytokine Clinical tibble preparation

# Cytokine and clinical combine

total_cytokine <- cytokine.import %>%
  select(!pub_id:`pub_id paper`) %>%
  right_join(clinical.max.irae, by = "pt id")

# Creating baseline cytokine tibble

demo_baseline <- total_cytokine %>%
  filter(`blood timepoint` == "baseline")

demo_baseline_log2 <- demo_baseline %>%
  mutate(log2_conc = log(average_baseline_conc, 2)) %>%
  relocate(log2_conc, .after = "concentration") 

demo_baseline_log2_review <- demo_baseline_log2 %>%
  ungroup() %>%
  distinct(`pt id`, `blood timepoint`, .keep_all = TRUE) 

# Creating change in cytokine tibble

# select list of pts with baseline samples available

demo_baseline_run1 <- demo_baseline %>%
  filter(run == 1)

demo_baseline_run2 <- demo_baseline %>%
  filter(run == 2)

demo_baseline_run3 <- demo_baseline %>%
  filter(run == 3)

baseline_avail_pt_id_run1 <- demo_baseline_run1$`pt id`

baseline_avail_pt_id_run2 <- demo_baseline_run2$`pt id`

baseline_avail_pt_id_run3 <- demo_baseline_run3$`pt id`

#Create new tibble selecting for all samples available for pts with baseline samples available
run1_filtered_total_cytokine_pre.clean <- total_cytokine %>%
  filter(`pt id` %in% baseline_avail_pt_id_run1)

run2_filtered_total_cytokine_pre.clean <- total_cytokine %>%
  filter(`pt id` %in% baseline_avail_pt_id_run2)

run3_filtered_total_cytokine_pre.clean <- total_cytokine %>%
  filter(`pt id` %in% baseline_avail_pt_id_run3)

#Calculating changes for run 1 based on baseline conc with that run
run1_filtered_total_cytokine <- run1_filtered_total_cytokine_pre.clean %>%
  select(!average_baseline_conc) %>%
  filter(run == 1) %>%
  filter(!is.na(concentration)) %>%
  group_by(`pt id`,`cytokines`) %>%
  arrange(`pt id`,`blood timepoint order`) %>%
  mutate(log2_conc = log2(`concentration`),
         baseline_conc = concentration[row_number() ==1],
         delta = concentration - concentration[row_number() == 1],
         relative_change = (concentration - concentration[row_number() == 1])/concentration[row_number() == 1],
         fold_log2 = log2(concentration/concentration[row_number()==1]),
         fold = concentration/concentration[row_number() == 1]) %>%
  filter(!`blood timepoint` == "baseline") %>%
  relocate(log2_conc:fold, .after = concentration) %>%
  ungroup()

#Calculating changes for run 2 based on baseline conc with that run
run2_filtered_total_cytokine <- run2_filtered_total_cytokine_pre.clean %>%
  select(!average_baseline_conc) %>%
  filter(run == 2) %>%
  filter(!is.na(concentration)) %>%
  group_by(`pt id`,`cytokines`) %>%
  arrange(`pt id`,`blood timepoint order`) %>%
  mutate(log2_conc = log2(`concentration`),
         baseline_conc = concentration[row_number() ==1],
         delta = concentration - concentration[row_number() == 1],
         relative_change = (concentration - concentration[row_number() == 1])/concentration[row_number() == 1],
         fold_log2 = log2(concentration/concentration[row_number()==1]),
         fold = concentration/concentration[row_number() == 1]) %>%
  filter(!`blood timepoint` == "baseline") %>%
  relocate(log2_conc:fold, .after = concentration) %>%
  ungroup() 

#Calculating changes for run 3 based on baseline conc with that run
run3_filtered_total_cytokine <- run3_filtered_total_cytokine_pre.clean %>%
  select(!average_baseline_conc) %>%
  filter(run == 3) %>%
  filter(!is.na(concentration)) %>%
  group_by(`pt id`,`cytokines`) %>%
  arrange(`pt id`,`blood timepoint order`) %>%
  mutate(log2_conc = log2(`concentration`),
         baseline_conc = concentration[row_number() ==1],
         delta = concentration - concentration[row_number() == 1],
         relative_change = (concentration - concentration[row_number() == 1])/concentration[row_number() == 1],
         fold_log2 = log2(concentration/concentration[row_number()==1]),
         fold = concentration/concentration[row_number() == 1]) %>%
  filter(!`blood timepoint` == "baseline") %>%
  relocate(log2_conc:fold, .after = concentration) %>%
  ungroup() 

# Combining run 1 and run 2 for total change in cytokine tibble
demo_change <- rbind(run1_filtered_total_cytokine,run2_filtered_total_cytokine,run3_filtered_total_cytokine)

# Creating condensed tibble for unique pts, timepoints, and run (reviewing purposes)
demo_change_review <- demo_change %>%
  distinct(`pt id`, `blood timepoint`,`run`, .keep_all = TRUE) %>%
  relocate(plate, run, .after = `blood timepoint`)

# Creating RDS files

saveRDS(total_cytokine, file = here(output.path.rds, "total_cytokine.rds"))

saveRDS(demo_baseline_log2, file = here(output.path.rds, "demo_baseline_log2.rds"))

saveRDS(demo_change, file = here(output.path.rds, "demo_change.rds"))
