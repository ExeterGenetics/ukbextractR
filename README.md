# ukbextractR
`ukbextractR` — Reproducible Extraction of Linked Healthcare Records from UK Biobank

* Authors: Harry Green and Jiaqi Li
* Affiliation: University of Exeter
* Contact: h.d.green@exeter.ac.uk

`ukbextractR` is a lightweight suite of functions designed to reproducibly extract UK Biobank–linked healthcare records (HES inpatient, HES episodes, GP records, cancer registry, and OPCS procedure data) for downstream epidemiological and genetic analyses. This workflow is intended for use within the UK Biobank Research Analysis Platform (DNAnexus) and relies on Luke Pilling's `ukbrapR` package to export the raw tables into persistent project storage.

## Project Setup (create_tables.R)

Before using `ukbextractR` functions, you'll need to have extracted the healthcare records into your persistent project directory using `create_tables.R`. This will submit a number of `Table Exporter` jobs and create .tsv files in `/ukbrapr_data` in your project directory. This will pull from the latest version of the UK Biobank data. Make sure you only run this script **once** per project. More information on how this works can be found at https://lcpilling.github.io/ukbrapR/.

## Session Setup (session_setup.R)

Running `session_setup.R` at the start of an RStudio Workbench session will move the created .tsvs onto your instance and create a number of functions available for use:

## Available Functions

### read_ICD9 / read_ICD10

The functions `read_ICD9` and `read_ICD10` will pull ICD-10 linked data from the UK Biobank HES records. It uses a system level `grep` to identify diagnosis records that match your list of codes, then joins that with the main HES records table to link in the epistart and epiend dates. The grep is performed on the start of the code, so E11 will also pull out E11.0, E11.1, etc. If you do not want that, specify complete ICD10 codes. Multiple codes can be given in a list.

`frozen_shoulder=read_ICD10('M75')` will give a table of all ICD10 records matching M75*.

`diabetes=read_ICD10(c('E10','E11','E12'))` will give a table of all ICD10 records matching E10*, E11* etc.

### read_GP

The function `read_GP` works similarly to the above but can be pointed at either the clinical records (default) or the prescription data. The user can provide a mix of read 2 and read 3 codes (for clinical) or read 2, BNF and DMD codes (for prescriptions). The argument `table='scripts'` will point the function at the prescription records

`bmi=read_GP('22K..')` will give a table of all read 2 and read 3 records matching `22K..`. This is both the read 2 and read 3 code for BMI.

`metformin=read_GP('f41',table='scripts')` will give a table of all prescription records with read_2 code `f41` (metformin).

### read_cancer

The function `read_cancer` pulls data from the cancer registry data rather than the ICD10 records.

`read_cancer('C61')` will give a table of all cancer registry records for prostate cancer.

`read_cancer(c('C50','C61'))` will give a table of all cancer registry records for breast or prostate cancer.

### read_OPCS

The function `read_OPCS` reads from the operations and procedure codes in OPCS4 format. Note that the UK Biobank does NOT use decimal points, so Y11.0 is coded as Y110, and using Y11.0 will not return any results.

`ablation=read_OPCS('Y11')` will give a table of all OPCS4 records matching Y11* (destruction of organ).

## FAQs:

#### How do I quickly count the number of people with a particular code?

`read_cancer('C61')%>%select(eid)%>%distinct()%>%nrow()` will give the number of unique IDs matching your code.

#### How do I find all codes relating to an individual?

No individual eids should ever be present in any script. But this can be done by using `grep` on the tsv file you're interested in.

#### Do you have a quick way of running the setup scripts?

Not yet. I used to have a script that would directly run an r script from github, but this repo is currently private while it's still in development.

#### I used your old repo a lot, why the change?

The old version required me to produce csvs of all the healthcare records in my project, and everyone that wanted to use it needed viewership access (and also to be on project 103356). It also required the point and click interface, which is bad for reproducibility. Removing it means people from outside of Exeter can use this package in its full functionality without needing the point and click interface.

#### Do I need to do anything different?

You should be able to use `read_ICD9`, `read_ICD10`, `read_GP`, `read_cancer` and `read_OPCS` the same as before. The only difference is the columns that are outputted, notably it doesn't join with the baseline_table so it doesn't automatically attach the assessment centre date and age at record.

#### What happened to the baseline_table?

That has been removed from this version. It relied on a preformatted table which was held in my project. I haven't yet worked this into the new setup.

#### It triggered an error

That isn't a question. But please send me the exact command you ran, and what the error was, and I can try to fix it.
