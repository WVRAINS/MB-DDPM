clr <- function(x.f,
                base = exp(1),
                tol = .Machine$double.eps,
                ...) {
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG) / mean(nzero), 0.0)
}

vis.spearman <- function(..., taxa.name,
                         titles = c("Real", "MB-DDPM", "MB-GAN", 
                                    "MIDASim", "SparseDOSSA", "NorTA"),
                         color.palette = colorRampPalette(c("blue", "white", "red"))(200),
                         plot.type = "upper",
                         show.diag = FALSE,
                         method = "ellipse",
                         label = "(b)",  
                         label.cex = 2,  
                         main.cex = 1.8, 
                         number.cex = 3) {

  if (missing(taxa.name)) {
    stop("taxa.name must be specified")
  }

  mats <- list(...)
  if (length(mats) == 0) {
    stop("At least one input matrix must be provided")
  }

  if (length(titles) < length(mats)) {
    warning("Number of titles provided is less than number of matrices. Using default names.")
    titles <- c(titles, paste("Matrix", seq_len(length(mats) - length(titles))))
  }

  mats <- lapply(mats, function(mat) {
    mat <- as.matrix(mat)
    missing_taxa <- setdiff(taxa.name, colnames(mat))
    if (length(missing_taxa) > 0) {
      stop(paste("Taxa not found in matrix:", paste(missing_taxa, collapse=", ")))
    }
    mat[, taxa.name, drop = FALSE]
  })

  cor_list <- lapply(mats, function(mat) {
    cor(mat, method = "spearman", use = "pairwise.complete.obs")
  })

  cor_range <- range(unlist(cor_list), na.rm = TRUE)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(oma = c(2, 0, 0, 0),
      bg = "white",
      cex.axis = 1.5,
      cex.lab = 1.5)

  n <- length(cor_list)
  if (n == 3) {
    layout(matrix(1:3, nrow = 1))
  } else if (n == 6) {
    layout(matrix(1:6, nrow = 2, byrow = TRUE))
  } else {
    ncol <- min(3, ceiling(sqrt(n)))
    nrow <- ceiling(n / ncol)
    layout(matrix(1:(ncol * nrow), nrow = nrow, byrow = TRUE))
  }
  
  plot_params <- list(
    type = plot.type,
    diag = show.diag,
    method = method,
    is.corr = TRUE,
    col = color.palette,
    tl.pos = "n",
    cl.lim = cor_range,
    mar = c(0, 0, 3, 0),
    cex.main = main.cex,
    number.cex = number.cex
  )

  for (i in seq_along(cor_list)) {
    do.call(corrplot::corrplot, c(list(cor_list[[i]], main = titles[i]), plot_params))
  }
  
  if (!is.null(label)) {
    mtext(label, side = 3, outer = TRUE, line = -1.8, adj = 0, cex = label.cex)
  }

  layout(1)
}

scatter.spearman <- function(Xmat1, Xmat2, Xmat3 = NULL, Xmat4 = NULL, 
                             Xmat5 = NULL, Xmat6 = NULL, taxa.name,
                             method.names = c("MB-DDPM", "MB-GAN", "MIDASim",
                                              "SparseDOSSA", "NorTA"),
                             plot.hist = FALSE, nbreak = 10,
                             label = "(b)",  
                             label.cex = 2,        # Added parameter for label size
                             point.size = 1.2, text.size = 1.8,
                             plot.layout = NULL) {
  
  if (missing(Xmat1) || missing(Xmat2)) {
    stop("At least Xmat1 (real data) and Xmat2 (first simulation) must be provided")
  }
  if (missing(taxa.name)) {
    stop("taxa.name must be specified")
  }

  mats <- list(Xmat1, Xmat2)
  if (!is.null(Xmat3)) mats <- c(mats, list(Xmat3))
  if (!is.null(Xmat4)) mats <- c(mats, list(Xmat4))
  if (!is.null(Xmat5)) mats <- c(mats, list(Xmat5))
  if (!is.null(Xmat6)) mats <- c(mats, list(Xmat6))

  check_taxa <- function(mat) {
    all(taxa.name %in% colnames(mat))
  }
  
  for (i in seq_along(mats)) {
    if (!check_taxa(mats[[i]])) {
      stop(sprintf("Not all specified taxa found in matrix %d", i))
    }
  }

  extract_subset <- function(mat) {
    mat[, colnames(mat) %in% taxa.name, drop = FALSE]
  }
  
  mats_sub <- lapply(mats, extract_subset)

  calc_cor <- function(mat) {
    cor(mat, method = "spearman", use = "pairwise.complete.obs")
  }
  
  cor_mats <- lapply(mats_sub, calc_cor)

  extract_upper <- function(cor_mat) {
    as.vector(cor_mat[upper.tri(cor_mat)])
  }
  
  cor_vecs <- lapply(cor_mats, extract_upper)

  all_vals <- unlist(cor_vecs)
  axis_range <- range(all_vals, na.rm = TRUE)
  axis_range <- axis_range + c(-0.05, 0.05) * diff(axis_range)

  n_sim <- length(mats) - 1

  nrow <- 1
  ncol <- n_sim
  
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))

  par(mfrow = c(nrow, ncol),
      mar = c(5.5, 5, 4, 2),
      oma = c(1, 1, 3, 1),
      pty = "s",
      cex.lab = text.size,
      cex.main = text.size * 1.1,
      cex.axis = text.size * 0.9)

  create_scatter <- function(real_vec, sim_vec, title) {
    lm.fit <- lm(sim_vec ~ real_vec)
    lm.info <- summary(lm.fit)
    
    plot(real_vec, sim_vec,
         xlim = axis_range, ylim = axis_range,
         pch = 16, cex = point.size,
         col = adjustcolor("black", alpha.f = 0.6),
         xlab = "Real data correlation",
         ylab = paste(title, "correlation"),
         main = sprintf("%s vs Real\nR² = %.3f  MSE = %.3f",
                        title,
                        lm.info$adj.r.squared,
                        mean(lm.info$residuals^2)))
    
    abline(lm.fit, col = "red", lwd = 1.8)
    abline(0, 1, col = "blue", lty = 2, lwd = 1.5)
    
    legend("topleft", 
           legend = c("Regression", "y = x"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(1.8, 1.5),
           bty = "n",
           cex = text.size * 0.8)
  }

  real_vec <- cor_vecs[[1]]
  for (i in seq_len(n_sim)) {
    method_name <- if (i <= length(method.names)) method.names[i] else paste("Method", i)
    create_scatter(real_vec, cor_vecs[[i + 1]], method_name)
  }
  
  if (!is.null(label)) {
    mtext(label, side = 3, outer = TRUE, line = -1.8, adj = 0, cex = label.cex)
  }
}

vis.proportionality <- function(..., taxa.name, 
                                method = "color", 
                                col = colorRampPalette(c("red", "white"))(30),
                                titles = c("Real", "MB-DDPM", "MB-GAN", 
                                           "MIDASim", "SparseDOSSA", "NorTA"),
                                label = '(a)',          # Added parameter for flexible labeling
                                label.cex = 2,        # Added parameter for label size
                                main.cex = 1.8) {     # Added parameter for title size
  
  # Input validation
  if (missing(taxa.name)) {
    stop("taxa.name must be specified")
  }
  
  # Get all input matrices
  mats <- list(...)
  if (length(mats) == 0) {
    stop("At least one input matrix must be provided")
  }
  
  # Check titles length matches number of matrices
  if (length(titles) < length(mats)) {
    warning("Number of titles provided is less than number of matrices. Using default names.")
    titles <- c(titles, paste("Matrix", seq_len(length(mats) - length(titles))))
  }
  
  # Convert to matrices and check dimensions
  mats <- lapply(mats, function(x) {
    x <- as.matrix(x)
    if (nrow(x) == 0 || ncol(x) == 0) {
      stop("All input matrices must have positive dimensions")
    }
    x
  })
  
  # Zero imputation function
  impute_zeros <- function(mat) {
    non_zero_vals <- sort(unique(as.vector(mat[mat > 0])))
    impute_val <- ifelse(length(non_zero_vals) >= 2, non_zero_vals[2], min(non_zero_vals))
    mat[mat == 0] <- impute_val
    mat
  }
  
  # Apply zero imputation
  mats <- lapply(mats, impute_zeros)
  
  # Calculate proportionality matrices
  calc_phi <- function(mat) {
    if (!requireNamespace("compositions", quietly = TRUE)) {
      stop("Package 'compositions' required but not installed")
    }
    
    mat <- compositions::acomp(mat)
    vlr <- compositions::variation(mat, robust = "pearson")
    clr_trans <- t(apply(mat, 1, compositions::clr))
    clr_var <- apply(clr_trans, 2, var)
    phi <- sweep(vlr, 2, clr_var, FUN = "/")
    phi
  }
  
  phi_list <- lapply(mats, calc_phi)
  
  # Extract subset for target taxa
  phi_sub <- lapply(phi_list, function(phi) {
    keep_rows <- rownames(phi) %in% taxa.name
    keep_cols <- colnames(phi) %in% taxa.name
    
    if (sum(keep_rows) == 0 || sum(keep_cols) == 0) {
      stop("No matching taxa found in one or more matrices")
    }
    
    phi[keep_rows, keep_cols]
  })
  
  # Set common color range
  col_range <- range(unlist(phi_sub), na.rm = TRUE)
  
  # Set graphical parameters
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(oma = c(2, 0, 0, 0),  # Bottom margin for label
      bg = "white")         # Background color
  
  # Custom layout based on number of matrices
  n <- length(phi_sub)
  if (n == 3) {
    # 3 plots in 1 row
    layout(matrix(1:3, nrow = 1))
  } else if (n == 6) {
    # 6 plots in 2 rows
    layout(matrix(1:6, nrow = 2, byrow = TRUE))
  } else {
    # Automatic layout for other counts
    ncol <- min(3, ceiling(sqrt(n)))
    nrow <- ceiling(n / ncol)
    layout(matrix(1:(ncol * nrow), nrow = nrow, byrow = TRUE))
  }
  
  # Define consistent plot parameters
  plot_params <- list(
    type = "full",
    diag = FALSE,
    method = method,
    is.corr = FALSE,
    mar = c(0, 0, 1, 0),  # Margins for each plot
    tl.pos = "n",          # No text labels
    cl.lim = col_range,    # Common color limits
    col = col,
    cex.main = main.cex    # Title size
  )
  
  # Generate plots
  for (i in seq_along(phi_sub)) {
    do.call(corrplot::corrplot, c(list(phi_sub[[i]], main = titles[i]), plot_params))
  }
  
  if (!is.null(label)) {
    mtext(label, side = 3, outer = TRUE, line = -1.8, adj = 0, cex = label.cex)
  }
  
  # Reset layout to default
  layout(1)
}

scatter.proportionality <- function(Xmat1, Xmat2, Xmat3 = NULL, Xmat4 = NULL, 
                                    Xmat5 = NULL, Xmat6 = NULL, taxa.name,
                                    method.names = c("MB-DDPM", "MB-GAN", "MIDASim",
                                                     "SparseDOSSA", "NorTA"),
                                    point.size = 1.2, text.size = 1.8,
                                    label = NULL,
                                    label.cex = 2,        # Added parameter for label size
                                    plot.layout = NULL) {
  require(compositions, quietly = TRUE)

  if (missing(Xmat1) || missing(Xmat2)) {
    stop("At least Xmat1 (real data) and Xmat2 (first simulation) must be provided")
  }
  if (missing(taxa.name)) {
    stop("taxa.name must be specified")
  }

  mats <- list(Xmat1, Xmat2)
  if (!is.null(Xmat3)) mats <- c(mats, list(Xmat3))
  if (!is.null(Xmat4)) mats <- c(mats, list(Xmat4))
  if (!is.null(Xmat5)) mats <- c(mats, list(Xmat5))
  if (!is.null(Xmat6)) mats <- c(mats, list(Xmat6))

  mats <- lapply(mats, function(mat) {
    mat <- as.matrix(mat)
    impute.val <- unique(sort(as.vector(mat)))[2]
    mat[mat == 0] <- impute.val
    return(mat)
  })

  check_taxa <- function(mat) {
    all(taxa.name %in% colnames(mat))
  }
  
  for (i in seq_along(mats)) {
    if (!check_taxa(mats[[i]])) {
      stop(sprintf("Not all specified taxa found in matrix %d", i))
    }
  }

  compute_phi <- function(mat) {
    class(mat) <- 'acomp'
    vlr <- variation(mat)
    clr_mat <- t(apply(mat, 1, clr))
    var_clr <- apply(clr_mat, 2, var)
    sweep(vlr, 2, var_clr, FUN = "/")
  }
  
  phi.list <- lapply(mats, compute_phi)

  phi.sub.list <- lapply(phi.list, function(phi) {
    phi[rownames(phi) %in% taxa.name, colnames(phi) %in% taxa.name]
  })

  all_vals <- unlist(phi.sub.list)
  axis_range <- range(all_vals, na.rm = TRUE)
  axis_range <- axis_range + c(-0.05, 0.05) * diff(axis_range)

  n_sim <- length(mats) - 1
  
  nrow <- 1
  ncol <- n_sim
  
  old.par <- par(no.readonly = TRUE)
  on.exit(par(old.par))
  
  par(mfrow = c(nrow, ncol),
      mar = c(4.5, 5, 4, 2),
      oma = c(3, 1, 3, 1), 
      pty = "s",
      cex.lab = text.size,
      cex.main = text.size * 1.1,
      cex.axis = text.size * 0.9)

  create_scatter <- function(real_vec, sim_vec, title) {
    lm.fit <- lm(sim_vec ~ real_vec)
    lm.info <- summary(lm.fit)
    
    plot(real_vec, sim_vec,
         xlim = axis_range, ylim = axis_range,
         pch = 16, cex = point.size,
         col = adjustcolor("black", alpha.f = 0.6),
         xlab = "Real data proportionality",
         ylab = paste(title, "proportionality"),
         main = sprintf("%s vs Real\nR² = %.3f  MSE = %.3f",
                        title,
                        lm.info$adj.r.squared,
                        mean(lm.info$residuals^2)))
    
    abline(lm.fit, col = "red", lwd = 1.8)
    abline(0, 1, col = "blue", lty = 2, lwd = 1.5)
    
    legend("topleft", 
           legend = c("Regression", "y = x"),
           col = c("red", "blue"),
           lty = c(1, 2),
           lwd = c(1.8, 1.5),
           bty = "n",
           cex = text.size * 0.8)
  }

  real_vec <- as.vector(phi.sub.list[[1]])
  for (i in seq_len(n_sim)) {
    method_name <- if (i <= length(method.names)) method.names[i] else paste("Method", i)
    sim_vec <- as.vector(phi.sub.list[[i + 1]])
    create_scatter(real_vec, sim_vec, method_name)
  }
  if (!is.null(label)) {
    mtext(label, side = 3, outer = TRUE, line = -1.8, adj = 0, cex = label.cex)
  }
}