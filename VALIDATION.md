# Validation

- Six-step synthetic pipeline completed successfully.
- Independently checked frequency tables and seven distance formulas at three haplogroup levels (tolerance 1e-12).
- Verified distance symmetry and diagonal values.
- Verified separate bootstrap results for each level and exact seeded reproducibility of bootstrap result TSV files.
- Verified 154 PDF/TSV pairs, Arial font resources and at most one text-show operation per PDF text object.
- Checked empty input, duplicate sample IDs, missing columns, overwrite protection, omitted optional inputs and missing configured metadata.
- Checked automatic graphics devices remain inside project temp.
- Checked R parsing and shell syntax.

Only synthetic input fixtures are retained. Generated outputs, logs, debugging files, original datasets, obsolete source versions and local Git metadata were removed after validation. Run the documented pipeline and validation command to regenerate and check the example outputs.

The installed R packages emit some version and deprecation warnings; the checked runs completed with successful exit codes.
