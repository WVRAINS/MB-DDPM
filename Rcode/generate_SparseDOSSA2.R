if (!require("devtools")) install.packages("devtools")
if (!require("SparseDOSSA2")) devtools::install_github("biobakery/SparseDOSSA2")
if (!require("future.apply")) install.packages("future.apply")

library(SparseDOSSA2)
library(future.apply)

plan(multisession, workers = availableCores() - 1)

# setwd("./data")

load("NielsenHB_original_data.Rdata", verbose = TRUE)

Y.species <- t(NielsenHB_count)

Y.count.case <- Y.species[patient.info[["disease"]] == "IBD", ]
Y.count.ctrl <- Y.species[patient.info[["disease"]] == "healthy", ]

process_group_data <- function(data_matrix, group_name, n_samples = 1000) {
  data <- t(data_matrix)
  
  feature_counts <- rowSums(data > 0)
  valid_features <- feature_counts >= 2
  data_filtered <- data[valid_features, , drop = FALSE]
  
  set.seed(2174171)
  fitted_model <- tryCatch({
    fit_SparseDOSSA2(data = data_filtered)
  }, error = function(e) {
    cat("Model fitting error:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(fitted_model)) {
    stop("Model fitting failed, cannot proceed with data simulation")
  }
  
  sim_result <- SparseDOSSA2(
    template = fitted_model,
    n_sample = n_samples,
    new_features = FALSE,
    verbose = TRUE
  )
  
  simulated_data <- matrix(0, 
                           nrow = nrow(data), 
                           ncol = n_samples,
                           dimnames = list(rownames(data), NULL))
  
  simulated_features <- intersect(rownames(sim_result$simulated_data), rownames(data))
  if (length(simulated_features) > 0) {
    simulated_data[simulated_features, ] <- sim_result$simulated_data[simulated_features, ]
  }
  
  final_data <- t(simulated_data)
  
  if (any(grepl("s__", colnames(final_data)))) {
    final_data_spec <- final_data[, grep("s__", colnames(final_data)), drop = FALSE]
    final_compo <- t(apply(final_data_spec, 1, function(x) x / sum(x)))
  } else {
    final_compo <- t(apply(final_data, 1, function(x) x / sum(x)))
  }
  
  list(
    count_matrix = final_data,
    compo_matrix = final_compo,
    simulated_features = length(simulated_features)
  )
}

results <- future_lapply(list(
  list(data = Y.count.case, name = "case"),
  list(data = Y.count.ctrl, name = "ctrl")
), function(x) {
  process_group_data(x$data, x$name)
}, future.seed = TRUE)

case.result <- results[[1]]
ctrl.result <- results[[2]]

case.count.SparseDOSSA2 <- case.result$count_matrix
case.compo.SparseDOSSA2 <- case.result$compo_matrix
ctrl.count.SparseDOSSA2 <- ctrl.result$count_matrix
ctrl.compo.SparseDOSSA2 <- ctrl.result$compo_matrix

save(
  case.count.SparseDOSSA2,
  case.compo.SparseDOSSA2,
  ctrl.count.SparseDOSSA2,
  ctrl.compo.SparseDOSSA2,
  file = "Data_by_SparseDOSSA2_v1.0.Rdata"
)

plan(sequential)