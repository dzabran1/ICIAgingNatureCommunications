
library(here)
library(tidyverse)
library(openxlsx)
library(readxl)

# Set up

here::i_am("code/01.1_clinical.data.processing.R")

output.path <- here("output")
if (!dir.exists(here(output.path))) dir.create(here(output.path))

output.path.rds <- here("output/RDS files")
if (!dir.exists(here(output.path.rds))) dir.create(here(output.path.rds))

output.path.excel <- here("output/Excel files")
if (!dir.exists(here(output.path.excel))) dir.create(here(output.path.excel))

### Import Clinical Master File data

file_path <- here("data/clinical/Aging Master File data de.id cutoff 4.4.23.xlsx")
import.tb <- read_excel(file_path)

# Clinical data processing

clinical.max.irae <- import.tb %>%
  rename(`age65_status` = "age status (older than 65 years)",
         `autoimmune history` = "baseline autoimmune history",
         `regimen` = "immunotherapy regimen used for analysis",
         `cause_death` = "cause of death (cancer vs other)",
         `cause_death_non.cancer` = "cause of non.cancer death",
         `ae onset` = "ae onset for time to event analysis",
         `orr_character` = "objective response",
         `response` = "best response",
         `os` = "overall survival",
         `pfs` = "progression free survival") %>%
  mutate(race_w_b_o = case_when(
    str_detect(race, paste(c("white","black"), collapse = "|")) ~ race,
    is.na(race) ~ NA_character_,
    TRUE ~ "other")) %>%
  mutate(race_b = case_when(str_detect(race, "black") ~ "black",
                            is.na(race) ~ NA_character_,
                            TRUE ~ "other")) %>%
  mutate(race_w = case_when(str_detect(race, "white") ~ "white",
                            is.na(race) ~ NA_character_,
                            TRUE ~ "other")) %>%
  relocate(race_w_b_o, race_b, race_w, .after = race) %>%
  mutate(stage_status = case_when(stage == "Advanced/Metastatic" ~ "advanced/metastatic",
                                  is.na(stage) ~ NA_character_,
                                  TRUE ~ "neoadjuvant/adjuvant")) %>%
  mutate(`cancer group` = case_when(
    str_detect(`cancer type`, 
               paste(c("biliary", "pancreas","hepatocellular","colorectal", "^gastric$", "neuroendocrine", "^esophagogastric junction$"),
                     collapse = "|")) ~ "gi",
    str_detect(`cancer type`,
               paste(c("bladder","prostate","^renal cell carcinoma$"),
                     collapse = "|")) ~ "gu",
    str_detect(`cancer type`,
               paste(c("cervical","ovarian","endometrial","^vulvovaginal cancer$"),
                     collapse = "|")) ~ "gyn",
    str_detect(`cancer type`,
               paste(c("^squamous cell skin cancer$","^melanoma$","merkel"),
                     collapse = "|")) ~ "skin",
    str_detect(`cancer type`,
               paste(c("^head and neck cancer$","^lung cancer$"),
                     collapse = "|")) ~ "upper aerodigestive",
    str_detect(`cancer type`,
               paste(c("^adrenal carcinoma$","thyroid"),
                     collapse = "|")) ~ "endocrine",
    TRUE ~ "other")) %>%
  mutate(cancer_group = case_when(`cancer group` == "gi" ~ "gastrointestinal",
                                  `cancer group` == "gu" ~ "genitourinary",
                                  `cancer group` == "skin" ~ "skin",
                                  `cancer group` == "upper aerodigestive" ~ "upper aerodigestive",
                                  TRUE ~ "other")) %>%
  mutate(cancer_group_gi_other = case_when(`cancer group` == "gi" ~ "gi",
                                           TRUE ~ "other")) %>%
  mutate(cancer_group_gu_other = case_when(`cancer group` == "gu" ~ "gu",
                                           TRUE ~ "other")) %>%
  mutate(cancer_group_skin_other = case_when(`cancer group` == "skin" ~ "skin",
                                             TRUE ~ "other")) %>%
  mutate(cancer_group_UAD_other = case_when(`cancer group` == "upper aerodigestive" ~ "upper aerodigestive",
                                            TRUE ~ "other")) %>%
  mutate(irae_status = ifelse(irae == "yes", 1, 0)) %>%
  mutate(`autoimmune_status` = case_when(`autoimmune history` == "Yes" ~ 1,
                                         `autoimmune history` == "No" ~ 0,
                                         TRUE ~ NA_real_)) %>%
  mutate(`pfs_status` = case_when(`death` == "yes" ~ 1,
                                  `progression` == "yes" ~ 1,
                                  `progression` == "no" & !is.na(pfs) ~ 0,
                                  TRUE ~ NA_real_)) %>%
  mutate(`os_status` = case_when(`death` == "yes" ~ 1,
                                 `death` == "no" ~ 0,
                                 TRUE ~ NA_real_)) %>%
  mutate(`os_status_cancer` = case_when(`death` == "yes" & `cause_death` == "cancer" ~ 1,
                                        `death` == "no" | `cause_death` == "other" ~ 0,
                                        TRUE ~ NA_real_)) %>%
  mutate(regimen_type = case_when(str_detect(`regimen`, paste(c("^ipilimumab$","^relatlimab$"), collapse = "|")) ~ "Dual ICI",
                                  str_detect(`regimen`, paste(c("^pembrolizumab$",
                                                                "^nivolumab$", "^atezolizumab$",
                                                                "^avelumab$", "^cemiplimab$"
                                  ), collapse = "|")) ~ "ICI monotherapy",
                                  TRUE ~ "ICI with targeted or chemotherapy")) %>%
  relocate(regimen_type, .after = regimen) %>%
  mutate(orr = ifelse(orr_character == "yes", 1, 0)) %>%
  mutate(orr.age.status = case_when(age65_status == "yes" & orr_character == "yes" ~ "RO",
                                    age65_status == "yes" & orr_character == "no" ~ "NRO",
                                    age65_status == "no" & orr_character == "yes" ~ "RY",
                                    age65_status == "no" & orr_character == "no" ~ "NRY",
                                    TRUE ~ NA_character_))

# Save rds file

saveRDS(clinical.max.irae, file = here(output.path.rds, "clinical.rds"))
