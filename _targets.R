library(targets)
library(tarchetypes)

# Set target options:
tar_option_set(
  # get packages from `DESCRIPTION` file
  packages = desc::desc_get_deps()$package
)

# to run parts of the pipeline in parallel
future::plan(
  future::multisession,
  workers = parallel::detectCores() - 1
)

# Run the R scripts in the R/ folder
tar_source()

tar_plan(
  parameters = params(),
  tar_file(
    data_path,
    simulate_arrow_dataset(parameters)
  ),
  tar_target(
    model_results,
    fit_simulated_models(data_path, parameters),
    format = "parquet"
  ),
  fig_a1 = plot_figure_a1(data_path),
  fig_a2 = plot_figure_a2(data_path),
  fig_weight = plot_figure_weight(data_path),
  fig_skew = plot_figure_skew(model_results),
  fig_sims = plot_figure_sims(model_results),
  fig_coverage = plot_figure_coverage(model_results),
  fig_bias = plot_figure_bias(model_results),
  fig_variance = plot_figure_variance(model_results),
  tar_quarto(readme, "README.qmd")
)
