args <- commandArgs(TRUE)
if (!requireNamespace("yaml", quietly = TRUE)) stop("Missing package: yaml")
config <- yaml::read_yaml(args[1], eval.expr = FALSE)
for (section in names(config)) {
  for (key in names(config[[section]])) {
    if (!grepl("^[a-z][a-z0-9_]*$", section) || !grepl("^[a-z][a-z0-9_]*$", key)) stop("Invalid configuration key")
    value <- config[[section]][[key]]
    if (length(value) != 1L) stop("Configuration values must be scalar")
    cat(sprintf("%s=%s\n", paste(section, key, sep = "_"), shQuote(as.character(value))))
  }
}
