suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(ggalluvial))
suppressPackageStartupMessages(library(grid))
suppressPackageStartupMessages(library(maps))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ggplot2))
grDevices::pdfFonts(Arial = grDevices::Type1Font("Arial", rep("data/Arial.afm", 4)))
options(stringsAsFactors = FALSE, digits = 7, warn = 1)
args <- commandArgs(TRUE)
if (length(args)%%2) stop("Expected named argument pairs")
cli <- as.list(args[seq(2, length(args), 2)])
names(cli) <- sub("^--", "", args[seq(1, length(args), 2)])
for (key in names(cli)) {
    value <- cli[[key]]
    if (value %in% c("TRUE", "FALSE", "true", "false")) 
        value <- tolower(value) == "true"
    else if (grepl("^-?[0-9]+([.][0-9]+)?([eE][-+]?[0-9]+)?$", value)) 
        value <- as.numeric(value)
    assign(key, value)
}
HH_range <- as.integer(strsplit(as.character(levels), ",", fixed = TRUE)[[1]])
if (anyNA(HH_range) || any(HH_range < 1)) stop("Invalid levels")
if (jobs < 1 || jobs != as.integer(jobs)) stop("Invalid jobs")
set.seed(seed)
Sys.setenv(TMPDIR = temp)
dir.create(temp, recursive = TRUE, showWarnings = FALSE)
options(device = function(...) grDevices::pdf(file.path(temp,"graphics-device.pdf"),family="Arial",useDingbats=FALSE,...))
if (any(vapply(c(output, figure), function(folder) dir.exists(folder) && length(list.files(folder, recursive = TRUE)) > 
    0, logical(1))) && !overwrite) stop("Output exists; enable overwrite explicitly")
for (folder in c(output, figure)) dir.create(folder, recursive = TRUE, showWarnings = FALSE)
write_result <- function(x, file, row.names = FALSE, fileEncoding = "UTF-8", ...) {
    file <- sub("\\.csv$", ".tsv", file)
    file <- sub("\\.txt$", ".tsv", file)
    if (is.matrix(x)) 
        x <- data.frame(population = rownames(x), x, check.names = FALSE)
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    readr::write_tsv(as.data.frame(x), file, na = "NA")
}
read_table <- function(file, sheet = 1) {
    if (!file.exists(file)) 
        stop("Input does not exist: ", file)
    ext <- tolower(tools::file_ext(file))
    if (ext %in% c("xlsx", "xls")) 
        return(as.data.frame(readxl::read_excel(file, sheet = sheet)))
    if (ext == "csv") 
        return(read.csv(file, check.names = FALSE))
    read.delim(file, check.names = FALSE)
}
save_plot <- function(filename, plot = ggplot2::last_plot(), ...) {
    filename <- file.path(figure, basename(filename))
    if (exists("bootstrap_step") && bootstrap_step) 
        filename <- file.path(figure, paste0("HH", HH, "_", basename(filename)))
    dots <- list(...)
    dots$dpi <- NULL
    built <- ggplot2::ggplot_build(plot + ggplot2::theme(text = ggplot2::element_text(family = "Arial")))
    grob <- ggplot2::ggplot_gtable(built)
    do.call(ggplot2::ggsave, c(list(filename = filename, plot = grob,
        device = grDevices::pdf, family = "Arial", useDingbats = FALSE), dots))
    frames <- list()
    if (is.data.frame(plot$data) && nrow(plot$data)) 
        frames[[1]] <- mutate(plot$data, data_source = "source")
    for (i in seq_along(built$data)) {
        z <- built$data[[i]]
        z[] <- lapply(z, function(v) if (is.list(v)) 
            vapply(v, paste, character(1), collapse = ";")
        else if (is.numeric(v)) 
            sprintf("%.17g", v)
        else as.character(v))
        frames[[length(frames) + 1]] <- mutate(z, data_source = paste0("layer_", i))
    }
    frames <- lapply(frames, function(z) {
        z[] <- lapply(z, function(v) if (is.numeric(v)) 
            sprintf("%.17g", v)
        else if (is.numeric(v)) 
            sprintf("%.17g", v)
        else as.character(v))
        z
    })
    write_result(bind_rows(frames), sub("\\.pdf$", ".tsv", filename))
}
open_pdf <- function(file, ...) {
    target <- file.path(figure, basename(file))
    env <- parent.frame()
    for (key in c("D", "mat")) {
        if (exists(key, env, inherits = FALSE)) {
            write_result(get(key, env), sub("\\.pdf$", ".tsv", target))
            break
        }
    }
    grDevices::pdf(target, ..., family = "Arial", useDingbats = FALSE)
}
for (module in c("2-2-statistics", "3-2-distances", "4-2-plots", "5-2-resampling", "6-2-geography-library")) source(paste0("R/", 
    module, ".R"))
