# Inspect installed function bodies to verify patch application
tryCatch(
    {
        library(because)

        cat("=== Inspecting because::because for latent filtering logic ===\n")
        # Look for the specific comment or logic added
        # "Remove occupancy variables from potential_latents"

        because_code <- as.character(body(because::because))
        found_filter <- any(grepl(
            "potential_latents <- setdiff\\(potential_latents, occ_vars\\)",
            because_code
        ))

        if (found_filter) {
            cat("✅ Latent filtering logic FOUND in because().\n")
        } else {
            cat("❌ Latent filtering logic NOT FOUND in because().\n")
            # Print snippet for debugging if missing
            cat("Snippet of potential_latents logic:\n")
            print(grep("potential_latents", because_code, value = TRUE))
        }

        cat(
            "\n=== Inspecting because:::dsep_with_latents for obs injection ===\n"
        )
        dsep_with_latents <- get(
            "dsep_with_latents",
            envir = asNamespace("because")
        )
        dsep_code <- as.character(body(dsep_with_latents))

        # Look for "ov_obs ~ ov + p_ov"
        # I commented it out, so it shouldn't be in the body?
        # Or comments might be stripped in installed package?
        # If it IS present (active code), grep will find it.

        found_injection <- any(grepl("ov_obs ~ ov", dsep_code))

        if (found_injection) {
            cat(
                "❌ Obs injection logic FOUND in dsep_with_latents (It should be gone!).\n"
            )
        } else {
            cat("✅ Obs injection logic NOT FOUND in dsep_with_latents.\n")
        }

        cat("\n=== Package Version ===\n")
        print(packageVersion("because"))
        print(system.file(package = "because"))
    },
    error = function(e) {
        cat("Error inspecting package: ", e$message, "\n")
    }
)
