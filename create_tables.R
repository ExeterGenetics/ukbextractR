# This short script needs to be run once and only once per project
# It sends off a few table exporter jobs and will add tsvs to your persistant storage space
# See documentation for the ukbrapR package for more details


if (!file.exists('/mnt/project/ukbrapr_data/hesin.tsv')){
  remotes::install_github("lcpilling/ukbrapR")
  library('ukbrapR')
  ukbrapR:::export_tables(submit=TRUE)
}