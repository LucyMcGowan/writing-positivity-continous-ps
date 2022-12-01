#! /usr/bin/env Rscript

# This is a helper script to run the pipeline.
# Choose how to execute the pipeline below.
# See https://books.ropensci.org/targets/hpc.html
# to learn about your options.

# Use `cmd`/`ctrl` + `shift` + `B` to run in the build pane

cli::cli_alert_info("Starting targets pipeline")
targets::tar_make()
cli::cli_alert_success("targets pipeline completed")
