#' Build and Update Package Documentation
#'
#' This script automates the process of:
#' 1. Updating roxygen2 documentation (devtools::document)
#' 2. Generating the PDF Reference Manual
#' 3. Moving the PDF to inst/pdf/
#'
#' @description Run this script before tagging a new release to ensure the
#' PDF manual is up-to-date with the latest code changes.

# 1. Update Roxygen Documentation (Standard Help Files)
message("Updating standard documentation...")
devtools::document()

# 2. Define Paths
pkg_path <- getwd()
inst_pdf_dir <- file.path(pkg_path, "inst", "pdf")
manual_name <- "because_reference_manual.pdf"
output_path <- file.path(inst_pdf_dir, manual_name)

# Ensure directory exists
if (!dir.exists(inst_pdf_dir)) {
    dir.create(inst_pdf_dir, recursive = TRUE)
}

# 3. Generate PDF Manual
message("Generating PDF Reference Manual...")

# We use system call to R CMD Rd2pdf
# We intentionally ignore the output to avoid clutter, unless there's an error
cmd <- sprintf(
    "R CMD Rd2pdf . --output=%s --force --no-preview",
    shQuote(manual_name)
) # Generates in current dir first due to some Rd2pdf quirks or specify full path

# Actually Rd2pdf --output accepts full path, let's try direct
cmd_direct <- sprintf(
    "R CMD Rd2pdf %s --output=%s --force --no-preview",
    shQuote(pkg_path),
    shQuote(output_path)
)

status <- system(cmd_direct)

if (status == 0) {
    message(sprintf("SUCCESS: Manual generated at %s", output_path))
} else {
    stop(
        "FAILED: Could not generate PDF manual. Ensure you have LaTeX installed."
    )
}

message("Documentation update complete.")
