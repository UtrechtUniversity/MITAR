#### Define required functions ####
# Code extracted from https://github.com/JesseAlderliesten/InstallPkgs.


##### Utility functions #####
# Utility function to write more succinct code in stopifnot().
all_characters <- function(x, allow_zero_char = FALSE) {
  is.character(x) &&
    (allow_zero_char || length(x) > 0L) &&
    all(nzchar(x, keepNA = FALSE)) && !anyNA(x)
}

is_logical <- function(x) {
  is.logical(x) && length(x) == 1L
}


##### check_BiocManager #####
# (re)install the BiocManager package from CRAN if it is not installed or not
# functional.
check_BiocManager <- function(lib_path) {
  if(!requireNamespace("BiocManager", lib.loc = lib_path, quietly = TRUE)) {
    message("Trying to install package 'BiocManager' in ", lib_path)
    install.packages("BiocManager", lib.loc = lib_path, lib = lib_path,
                     type = "binary")
    if(requireNamespace("BiocManager", lib.loc = lib_path, quietly = TRUE)) {
      stop("If no error message is printed below, the BiocManager package",
           " was\nsuccesfully installed. Restart R before proceeding.")
    } else {
      stop("Installation of the BiocManager package failed.\nIf a warning like",
           " 'lib = \"", lib_path[[1]][1], "\" is not writeable'\nwas issued,",
           " you most likely forgot to run R as administrator,\nor used an",
           " incorrect path (the warnings printed below might point to that).",
           "\nRestart R as administrator (e.g., right-click on the R or RStudio",
           " icon), select\n'Run as administrator', open the 'InstallPkgs'",
           " R-project file, and try again.")
    }
  }
}


##### check_OS_is_Windows #####
# Check if the operating system is windows.
# Input:
# - on_error: character string indicating the action if the operating system is
#   not Windows. Options are 'warn' (default), 'message' or 'quiet'.
# Return:
# - A logical value indicating if the operating system is Windows, returned
#   invisibly.
check_OS_is_Windows <- function(on_error = c("warn", "message", "quiet")) {
  on_error <- match.arg(on_error, several.ok = FALSE)
  # Using 'Sys.info()' which returns information about the platform R is running
  # on, not 'R.version()' which returns information about the platform R was
  # built on.
  system_name <- Sys.info()["sysname"]
  
  # Need tolower() because 'ignore_case' does not work if 'fixed' is TRUE.
  if(grepl("windows", tolower(system_name), fixed = TRUE)) {
    OS_is_Windows <- TRUE
  } else {
    OS_is_Windows <- FALSE
    text <- paste0("This script might fail because it is meant to be used on",
                   " Windows, whereas you\nare using '", system_name,
                   "' instead. See the platform-specific installation",
                   " instructions\nat https://cran.r-project.org/ and the",
                   " 'R installation and administration manual'\nat",
                   " https://cran.r-project.org/doc/manuals/r-release/R-admin.html",
                   " for help.")
    
    switch(on_error,
           message = message(text),
           quiet = NULL,
           warning(text))
  }
  invisible(OS_is_Windows)
}


##### get_paths #####
# Get paths where packages are (or should be) installed.
# Input:
# - path: character(0) (default) or a character vector indicating paths.
# - quietly: logical of length 1 (default FALSE) indicating if messages should
#   be suppressed.
# Notes:
# - A warning is issued if the working directory is returned as element
#   'first_path', because that implies that no paths are present in '.libPath'
#   and no path was provided in argument 'path'.
# Return:
# - A list with five character elements, with elements for which no path is
#   found set to character(0):
#   - 'first_path': the first non-empty path;
#   - 'argument_paths': the value of argument 'path';
#   - 'Rversion_paths': paths from .libPath that contain the same R version
#     number as the currently running R session;
#   - 'other_paths': paths from .libPath that do not contain the same R version
#     number as the currently running R session;
#   - 'wd_path': the working directory.
get_paths <- function(path = character(0), quietly = FALSE) {
  stopifnot(all_characters(path, allow_zero_char = TRUE), is_logical(quietly))
  
  paths_libPath <- .libPaths()
  bool_paths_Rversion <- grepl(pattern = paste0("R-", as.character(getRversion())),
                               x = paths_libPath, fixed = TRUE)
  
  # Using a list because elements other than 'wd_path' can have lengths larger
  # than one. The selection of paths returns character(0) if no path is present,
  # because 'paths_libPath' is character.
  paths_possible <- list(argument_paths = path,
                         Rversion_paths = paths_libPath[bool_paths_Rversion],
                         other_paths = paths_libPath[!bool_paths_Rversion],
                         wd_path = getwd())
  
  # With these settings all_characters() only returns 'TRUE' for non-empty,
  # non-NA_character_ character strings containing at least one character. The
  # first such element in paths_possible is selected, and the first element of
  # that string is included as element 'first_path' in the returned list.
  path_first_OK <- which(vapply(X = paths_possible, FUN = all_characters,
                                FUN.VALUE = logical(1), allow_zero_char = FALSE,
                                USE.NAMES = FALSE))[1]
  path_first_name <- names(paths_possible[path_first_OK])
  path_first <- paths_possible[path_first_OK][[1]][1]
  
  if(path_first_name == "wd_path") {
    warning("Returning the working directory ('", path_first,
            "')\nas first path because no path was supplied in argument 'path'",
            " and no path was\nfound in '.libPaths'. Installation of R",
            " packages might fail if the path points\nto a network drive.")
  } else {
    if(!quietly) {
      message("Path '", path_first,
              "' is element 'first_path' of the returned list.")
    }
  }
  c(list(first_path = path_first), paths_possible)
}


#### Code to run ####
# This code assumes it is run on Windows, so warn if it is not.
check_OS_is_Windows(on_error = "warn")

# Define a path for the folders were R packages have to be installed. The path
# is chosen from existing library paths (i.e., .libPaths()), unless argument
# 'path' is set: then that path will be used. Among the existing library paths,
# a path containing the R version number corresponding to the currently running
# version of R is preferred.
path_mitar <- get_paths(path = character(0), quietly = TRUE)$first_path

# Install BiocManager from CRAN if it is not installed or not functional
check_BiocManager(lib_path = path_mitar)

# Select and install missing packages
needed_pkgs <- c(
  "deSolve",   # https://CRAN.R-project.org/package=deSolve
  "rootSolve", # https://CRAN.R-project.org/package=rootSolve
  "TruncatedNormal", # https://CRAN.R-project.org/package=TruncatedNormal
  "dplyr",     # https://CRAN.R-project.org/package=dplyr
  "tidyr",     # https://CRAN.R-project.org/package=tidyr
  "ggplot2",   # https://CRAN.R-project.org/package=ggplot2
  "ggrepel",   # https://CRAN.R-project.org/package=ggrepel
  "scales",    # https://CRAN.R-project.org/package=scales
  "cowplot"    # https://CRAN.R-project.org/package=cowplot
)

# Install those packages from 'needed_pkgs' that are not installed or not
# functional.
missing_pkgs <- needed_pkgs[suppressWarnings(suppressPackageStartupMessages(
  which(!sapply(X = needed_pkgs, FUN = requireNamespace, lib.loc = path_mitar,
                quietly = TRUE, simplify = TRUE))
))]
if(length(missing_pkgs) > 0L) {
  BiocManager::install(pkgs = missing_pkgs, lib = path_mitar, type = "both",
                       update = FALSE, ask = FALSE, checkBuilt = TRUE)
} else {
  message("All packages are installed and functional")
}
