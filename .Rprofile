# Options
options(menu.graphics = FALSE)
# Set repos
local({
  r <- getOption("repos")
  r["CRAN"] <- "https://cran.rstudio.com/"
  options(repos = r)
})
# Set environment variables
Sys.setenv(R_LIBS_USER = ".library")
Sys.setenv(R_ENVIRON_USER = ".Renviron")
Sys.setenv(R_PROFILE_USER = ".Rprofile")
Sys.setenv(R_HISTFILE = ".Rhistory")
Sys.setenv(R_COMPLETION = TRUE)

.First <- function() {
  if (interactive()) {
    options(
      prompt = paste0(
        "[",
        basename(getwd()),
        "] R> "
      )
    )
    message(
      paste(
        R.version.string,
        Sys.time(),
        getwd(),
        sep = " | "
      )
    )
    if (Sys.getenv("TERM") == "xterm-256color") {
      if ("colorout" %in% rownames(utils::installed.packages())) {
        require(
          "colorout"
        )
      }
    }
  }
}
.Last <- function() {
  if (interactive()) {
    try(
      savehistory(
        file = file.path(getwd(), ".Rhistory")
      )
    )
  }
}
