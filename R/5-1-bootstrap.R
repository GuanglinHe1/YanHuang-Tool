source("R/1-2-runtime.R")
bootstrap_step <- TRUE
bootstrap_pop_levels <- NULL
freqplot_sig_ref_pop <- NULL
df1 <- read_table(input)
for (HH in HH_range[HH_range <= max(nchar(df1$hap_clean))]) {
    df <- mutate(df1, hap = hap_substr(hap_clean, HH))
    run_bootstrap_for_hh(df, HH, file.path(output, sprintf("HH%02d", HH)))
}
