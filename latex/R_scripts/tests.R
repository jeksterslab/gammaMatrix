# find root directory
root <- rprojroot::is_rstudio_project
# set environment variables
Sys.setenv(R_LIBS_USER = ".library")
# source all scripts in R/
R_folder <- root$find_file(
  "R"
)
script_files <- list.files(
  R_folder,
  pattern = "\\.R$",
  full.names = TRUE,
  all.files = TRUE
)
for (i in script_files) {
  source(i)
}
rm(i, script_files, R_folder)
# load all data in data/
data_folder <- root$find_file(
  "data"
)
if (dir.exists(data_folder)) {
  data_files <- list.files(
    data_folder,
    pattern = "\\.rda$",
    full.names = TRUE,
    all.files = TRUE
  )
  for (i in data_files) {
    load(i)
  }
  rm(i, data_files)
}
rm(data_folder)
# run tests
test_folder <- root$find_file(
  "tests",
  "testthat"
)
test_files <- list.files(
  test_folder,
  pattern = "\\.R$",
  full.names = TRUE,
  all.files = TRUE
)
for (i in test_files) {
  source(i)
}
rm(i, test_files, test_folder)
