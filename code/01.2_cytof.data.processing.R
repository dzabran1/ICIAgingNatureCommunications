library(here)
library(tidyverse)
library(readxl)
library(openxlsx)

# Set up

here::i_am("code/01.2_cytof.data.processing.R")

output.path <- here("output")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

output.path.rds <- here("output/RDS files")
if (!dir.exists(here(output.path.rds))) dir.create(here(output.path.rds))

output.path.excel <- here("output/Excel files")
if (!dir.exists(here(output.path.excel))) dir.create(here(output.path.excel))

file_path_cytof_raw <- here("data/CYTOF/imCORE_CyTOF_batch1&2_post.normalization.xlsx")
cytof.raw.tb <- read_excel(file_path_cytof_raw) 

file_path_cytof_raw_v1 <- here("data/CYTOF/imCORE_CyTOF_batch1&2_post.normalization_v1.xlsx")
cytof.raw.tb_v1 <- read_excel(file_path_cytof_raw_v1) 

file_path_cytof_sample.labels <- here("data/CYTOF/CyTOF 1-2 runs cohort analysis designation for aging analysis.xlsx")
cytof.labels.tb <- read_excel(file_path_cytof_sample.labels) 

missing.cytof.sample <- cytof.raw.tb %>%
  anti_join(cytof.labels.tb, by =  "sample")

# detailed cytof clustering
# tibble cleaning
cytof.clean <- cytof.labels.tb %>%
  left_join(cytof.raw.tb, by = "sample") %>%
  filter(!(`pt id` == 15 & `blood timepoint number` == 0)) %>%
  filter(!(`pt id` == 18 & `blood timepoint number` == 3)) %>%
  filter(!(`pt id` == 43 & `blood timepoint number` == 2)) %>%
  mutate(`blood timepoint string` = ifelse(`pt id` == 15 & `blood timepoint string` == "baseline (2)", "baseline", `blood timepoint string`),
         `sample` = ifelse(`pt id` == 15 & `sample` == "15_1", "15_0", `sample`),
         `sample` = ifelse(`pt id` == 15 & `sample` == "15_2", "15_1", `sample`),
         `blood timepoint number` = ifelse(`pt id` == 15 & `blood timepoint number` == 1, 0, `blood timepoint number`),
         `blood timepoint number` = ifelse(`pt id` == 15 & `blood timepoint number` == 2, 1, `blood timepoint number`),
         `blood timepoint string` = ifelse(`pt id` == 18 & `blood timepoint string` == "irae (1)", "irae", `blood timepoint string`))

cytof.clean_baseline.sample_dates <- cytof.clean %>%
  filter(`blood timepoint number` == 0) %>%
  select(`pt id`,`blood timepoint string`, `sample date`)

# general cytof clustering (v1)
# tibble cleaning
cytof.clean_v1 <- cytof.labels.tb %>%
  left_join(cytof.raw.tb_v1, by = "sample") %>%
  filter(!(`pt id` == 15 & `blood timepoint number` == 0)) %>%
  filter(!(`pt id` == 18 & `blood timepoint number` == 3)) %>%
  filter(!(`pt id` == 43 & `blood timepoint number` == 2)) %>%
  mutate(`blood timepoint string` = ifelse(`pt id` == 15 & `blood timepoint string` == "baseline (2)", "baseline", `blood timepoint string`),
         `sample` = ifelse(`pt id` == 15 & `sample` == "15_1", "15_0", `sample`),
         `sample` = ifelse(`pt id` == 15 & `sample` == "15_2", "15_1", `sample`),
         `blood timepoint number` = ifelse(`pt id` == 15 & `blood timepoint number` == 1, 0, `blood timepoint number`),
         `blood timepoint number` = ifelse(`pt id` == 15 & `blood timepoint number` == 2, 1, `blood timepoint number`),
         `blood timepoint string` = ifelse(`pt id` == 18 & `blood timepoint string` == "irae (1)", "irae", `blood timepoint string`))


# importing absolute white count information
# cleaning and restructuring tibble

file_path_abs.counts <- here("data/clinical/final_absolute_cell_counts.xlsx")
wbc_counts_tb <- read_excel(file_path_abs.counts) %>%
  rename(`blood timepoint string` = "jci timepoint")

baseline_wbc_counts_for.cytof <- wbc_counts_tb %>%
  select(`pt id`, baseline_lymph_plus_monocytes) %>%
  mutate(`blood timepoint string` = "baseline") %>%
  relocate(`blood timepoint string`, .after = "pt id") %>%
  left_join(cytof.clean_baseline.sample_dates, by = c("pt id","blood timepoint string")) %>%
  relocate(`sample date`, .after = "blood timepoint string")  %>%
  rename(lymph.plus.monocytes = "baseline_lymph_plus_monocytes")

on.treatment_wbc_counts_for.cytof <- wbc_counts_tb %>%
  select(`pt id`, `blood timepoint string`, `sample date`, on.treatment_lymph_plus_monocytes) %>%
  rename(lymph.plus.monocytes = "on.treatment_lymph_plus_monocytes")

wbc_abs.count_with_sample.date <- rbind(baseline_wbc_counts_for.cytof, on.treatment_wbc_counts_for.cytof)

# Combine white count information with cytof

cytof.clean.long <- cytof.clean %>%
  left_join(wbc_abs.count_with_sample.date, by = c("pt id","blood timepoint string", "sample date")) %>%
  pivot_longer(`B_I`:UA, names_to = "immune_cell", values_to = "proportion") %>%
  mutate(cell_count = lymph.plus.monocytes * (proportion/100))

cytof.clean.long_v1 <- cytof.clean_v1 %>%
  left_join(wbc_abs.count_with_sample.date, by = c("pt id","blood timepoint string", "sample date")) %>%
  pivot_longer(`B`:UA, names_to = "immune_cell", values_to = "proportion") %>%
  mutate(cell_count = lymph.plus.monocytes * (proportion/100))

# import clinical rds file

clinical <- readRDS(here("output/RDS files/clinical.rds"))

# Combing cytof labels with clinical tibble

cytof_cohort_clinical <- cytof.labels.tb %>%
  left_join(clinical, by = "pt id")

cytof_cohort_clinical_review <- cytof_cohort_clinical %>%
  ungroup() %>%
  distinct(`pt id`, .keep_all = TRUE)

cytof.clinical.clean.long <- cytof.clean.long %>%
  left_join(clinical, by = "pt id") 

total_cohort_cytof_clinical <- cytof.clinical.clean.long %>%
  rename(`blood timepoint` = `blood timepoint string`)

cytof.clinical.clean.long_v1 <- cytof.clean.long_v1 %>%
  left_join(clinical, by = "pt id") 

total_cohort_cytof_clinical_v1 <- cytof.clinical.clean.long_v1 %>%
  rename(`blood timepoint` = `blood timepoint string`)

# For detailed CyTOF clustering

## CyTOF: Calculating change between on treatment timepoints compared to baseline with baseline samples removed

total_cytof_change <- total_cohort_cytof_clinical %>%
  group_by(`pt id`,`immune_cell`) %>%
  arrange(`pt id`,`blood timepoint`) %>%
  mutate(log2_proportion = log2(`proportion`),
         baseline_proportion = proportion[row_number() ==1],
         delta_prop = proportion - proportion[row_number() == 1],
         relative_change_prop = (proportion - proportion[row_number() == 1])/proportion[row_number() == 1],
         fold_log2_prop = log2(proportion/proportion[row_number()==1]),
         fold_prop = proportion/proportion[row_number() == 1]) %>%
  mutate(log2_cell.count = log2(`cell_count`),
         baseline_cell.count = cell_count[row_number() ==1],
         delta_cell.count = cell_count - cell_count[row_number() == 1],
         relative_change_cell.count = (cell_count - cell_count[row_number() == 1])/cell_count[row_number() == 1],
         fold_log2_cell.count = log2(cell_count/cell_count[row_number()==1]),
         fold_cell.count = cell_count/cell_count[row_number() == 1]) %>%
  filter(!`blood timepoint` == "baseline") %>%
  relocate(log2_proportion:fold_cell.count, .after = proportion) %>%
  ungroup() 

## CyTOF: Calculating change between on treatment timepoints compared to baseline but keeping baseline samples in tibble

total_cytof_change_with_baseline <- total_cohort_cytof_clinical %>%
  group_by(`pt id`,`immune_cell`) %>%
  arrange(`pt id`,`blood timepoint`) %>%
  mutate(log2_proportion = log2(`proportion`),
         baseline_proportion = proportion[row_number() ==1],
         delta_prop = proportion - proportion[row_number() == 1],
         relative_change_prop = (proportion - proportion[row_number() == 1])/proportion[row_number() == 1],
         fold_log2_prop = log2(proportion/proportion[row_number()==1]),
         fold_prop = proportion/proportion[row_number() == 1]) %>%
  mutate(log2_cell.count = log2(`cell_count`),
         baseline_cell.count = cell_count[row_number() ==1],
         delta_cell.count = cell_count - cell_count[row_number() == 1],
         relative_change_cell.count = (cell_count - cell_count[row_number() == 1])/cell_count[row_number() == 1],
         fold_log2_cell.count = log2(cell_count/cell_count[row_number()==1]),
         fold_cell.count = cell_count/cell_count[row_number() == 1]) %>%
  relocate(log2_proportion:fold_cell.count, .after = proportion) %>%
  ungroup() 

# For general clustering method with broader immune groups (naming convention in code is v1)

## CyTOF: change on treatment without baseline sample

total_cytof_change_v1 <- total_cohort_cytof_clinical_v1 %>%
  group_by(`pt id`,`immune_cell`) %>%
  arrange(`pt id`,`blood timepoint`) %>%
  mutate(log2_proportion = log2(`proportion`),
         baseline_proportion = proportion[row_number() ==1],
         delta_prop = proportion - proportion[row_number() == 1],
         relative_change_prop = (proportion - proportion[row_number() == 1])/proportion[row_number() == 1],
         fold_log2_prop = log2(proportion/proportion[row_number()==1]),
         fold_prop = proportion/proportion[row_number() == 1]) %>%
  mutate(log2_cell.count = log2(`cell_count`),
         baseline_cell.count = cell_count[row_number() ==1],
         delta_cell.count = cell_count - cell_count[row_number() == 1],
         relative_change_cell.count = (cell_count - cell_count[row_number() == 1])/cell_count[row_number() == 1],
         fold_log2_cell.count = log2(cell_count/cell_count[row_number()==1]),
         fold_cell.count = cell_count/cell_count[row_number() == 1]) %>%
  filter(!`blood timepoint` == "baseline") %>%
  relocate(log2_proportion:fold_cell.count, .after = proportion) %>%
  ungroup() 

## CyTOF: change on treatment with baseline sample

total_cytof_change_with_baseline_v1 <- total_cohort_cytof_clinical_v1 %>%
  group_by(`pt id`,`immune_cell`) %>%
  arrange(`pt id`,`blood timepoint`) %>%
  mutate(log2_proportion = log2(`proportion`),
         baseline_proportion = proportion[row_number() ==1],
         delta_prop = proportion - proportion[row_number() == 1],
         relative_change_prop = (proportion - proportion[row_number() == 1])/proportion[row_number() == 1],
         fold_log2_prop = log2(proportion/proportion[row_number()==1]),
         fold_prop = proportion/proportion[row_number() == 1]) %>%
  mutate(log2_cell.count = log2(`cell_count`),
         baseline_cell.count = cell_count[row_number() ==1],
         delta_cell.count = cell_count - cell_count[row_number() == 1],
         relative_change_cell.count = (cell_count - cell_count[row_number() == 1])/cell_count[row_number() == 1],
         fold_log2_cell.count = log2(cell_count/cell_count[row_number()==1]),
         fold_cell.count = cell_count/cell_count[row_number() == 1]) %>%
  relocate(log2_proportion:fold_cell.count, .after = proportion) %>%
  ungroup() 

# Creating rds and excel files

## detailed clustering

write.xlsx(cytof.clean, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_clinical.label_post.normalization.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(cytof_cohort_clinical_review, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof.cohorts_clinical.data.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(total_cohort_cytof_clinical, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_clinical.data_all.samples.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(total_cytof_change, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_clinical.data_calculated.change.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

saveRDS(total_cohort_cytof_clinical, file = here(output.path.rds, "total_cytof_clinical.rds"))

saveRDS(total_cytof_change, file = here(output.path.rds, "total_cytof_change.rds"))

saveRDS(total_cytof_change_with_baseline, file = here(output.path.rds, "total_cytof_change_with_baseline.rds"))

## general clustering (V1)

write.xlsx(cytof.clean_v1, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_v1_clinical.label_post.normalization.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(total_cohort_cytof_clinical_v1, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_v1_clinical.data_all.samples.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

write.xlsx(total_cytof_change_v1, 
           file = here(output.path.excel, paste0(Sys.Date(), "_", "cytof_v1_clinical.data_calculated.change.xlsx")), colNames = TRUE, rowNames = FALSE, append = FALSE)

saveRDS(total_cohort_cytof_clinical_v1, file = here(output.path.rds, "total_cytof_clinical_v1.rds"))

saveRDS(total_cytof_change_v1, file = here(output.path.rds, "total_cytof_change_v1.rds"))

saveRDS(total_cytof_change_with_baseline_v1, file = here(output.path.rds, "total_cytof_change_with_baseline_v1.rds"))
