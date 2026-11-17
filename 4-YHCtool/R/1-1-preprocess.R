source("R/1-2-runtime.R")
df1 <- read_table(input)
names(df1) <- tolower(names(df1))
aliases <- list(sample_id = c("sample_id", "sample", "id"), haplogroup_full = c("haplogroup_full", "haplogroup", 
    "hap"), population = c("population", "pop", "group"))
for (key in names(aliases)) {
    hit <- intersect(aliases[[key]], names(df1))
    if (length(hit) != 1) 
        stop("Missing or ambiguous column: ", key)
    names(df1)[match(hit, names(df1))] <- key
}
if (!nrow(df1)) stop("Empty sample input")
required <- names(aliases)
if (anyNA(df1[required]) || any(vapply(df1[required], function(x) any(!nzchar(trimws(x))), logical(1)))) stop("Missing sample values")
if (anyDuplicated(df1$sample_id)) stop("Duplicate sample identifiers")
df1$hap_clean <- sanitize_hap(df1$haplogroup_full)
if (any(!nzchar(df1$hap_clean))) stop("Invalid haplogroup")
if (use_equalize) {
    target <- min(table(df1$population))
    df1 <- bind_rows(lapply(split(df1, df1$population), function(z) z[sample.int(nrow(z), target), , drop = FALSE]))
}
write_result(df1, file.path(output, "samples.tsv"))
