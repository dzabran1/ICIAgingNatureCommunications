2025-01-27

Overview:

This file provides explanation for the post processed analysis of cytokine and CyTOF data associated with the publication “Age-related divergence of circulating immune responses in patients with solid tumors treated with Age-related divergence of circulating immune responses in patients with solid tumors treated with immune checkpoint inhibitors.”

Data availability:

The authors declare that the minimal data set for this study cannot be shared publicly due to ethical and legal restrictions on sharing de-identified data that aligns with the consent of research participants. Current JHU compliance policies require data with no direct consent for public open access sharing be under restricted access. We will provide access through Vivli, an established repository for clinical data that provides open access without a fee restricted to approved researchers under a Data Use Agreement. JHU compliance policy for Vivli requires additional anonymization of certain demographics, including use of age ranges and limiters to outlier values for weight, height, and certain rare diseases, while retaining sufficient value for reference and validation of results. Researchers can request more detailed data from the corresponding author shared though an approved collaboration arrangement.

Running code and analysis:

Once data has been obtained from Vivli, users wishing to re-analyze these data using scripts and R markdown files on GitHub should first deposit the data files into the “data” subfolder. The GitHub scripts and R markdown files should be correctly pathed to the data subfolder structure. If there are errors in running the scripts, please check if the pathing of the deposited data is in the correct sub folders.

In terms of running the scripts and markdown files, any output files (i.e. rds, excel, or figures) will be pathed to the “output” folder with corresponding subfolders.

In order to perform the analysis, please start off with data cleaning and processing by running the “01_run_render.R” script. Once this is performed, users can run the statistical analyses and figure generations by running the R markdown files 02.1_aging.total_cohort.table.Rmd to 02.8_aging_cytof.marker_analysis.Rmd which will roughly follow the flow of the paper.

