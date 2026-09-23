# Multiple imputations by predictive mean matching and regression trees

# The code contains the function to perform weighted multiple imputation of a single numeric variable.
# The methods implementable are predictive mean matching, single regression tree, bootstrap regression trees and Bayesian regression trees.
# In order to use the method, provide a dataframe. Text variable should be transformed to numeric (0-1) variable.
# The output consists in a matrix, with rows consisting in the different units and columns having the different imputations.

# The code uses the package rpart 4.1.24 available on CRAN, and has been implemented on R 4.5.1
# For questions please contact riccardo.desantis@unipd.it
# Last update 23/09/2026

# # Reproducible Example (please load the functions below before using)
# 
# set.seed(1234)
# y <- rnorm(25)
# x <- rnorm(25)
# z <- rnorm(25,0.5)
# v <- runif(25)
# x[sample(1:25,3)] <- NA
# dat.frame <- data.frame(y,x,z,v)
# 
# pmm_imputation(ampdf = dat.frame,
#                donors.pmm = 4,
#                no.imputs = 5,
#                seed = 123,
#                miss.ind = 2,
#                cov.ind = c(1,3),
#                weights.ind = 4)
# 
# cart_imputation_bayes_boot(ampdf = dat.frame,
#                            minbucket.cart = 3,
#                            no.imputs = 5,
#                            seed = 123,
#                            miss.ind = 2,
#                            cov.ind = c(1,3),
#                            weights.ind = 4)
# 
# cart_imputation_boot(ampdf = dat.frame,
#                      minbucket.cart = 3,
#                      no.imputs = 5,
#                      seed = 123,
#                      miss.ind = 2,
#                      cov.ind = c(1,3),
#                      weights.ind = 4)
# 
# cart_imputation(ampdf = dat.frame,
#                 minbucket.cart = 3,
#                 no.imputs = 5,
#                 seed = 123,
#                 miss.ind = 2,
#                 cov.ind = c(1,3),
#                 weights.ind = 4)


#functions

cart_imputation_bayes_boot <- function(
    ampdf,
    minbucket.cart,
    no.imputs,
    seed = NULL,
    miss.ind,
    cov.ind,
    weights.ind = NULL
) {
  
  if(!is.numeric(miss.ind) | !is.numeric(cov.ind)) stop("Indicators must be numbers")
  
  obs_idx <- !is.na(ampdf[,miss.ind])
  mis_idx <- is.na(ampdf[,miss.ind])
  
  Y_obs <- ampdf[obs_idx, miss.ind]
  
  df_obs <- ampdf[obs_idx, , drop = FALSE]
  
  xobs_data <- df_obs[
    , cov.ind,
    drop = FALSE
  ]
  
  xmis_data <- ampdf[
    mis_idx,
    cov.ind,
    drop = FALSE
  ]
  
  n_obs <- nrow(df_obs)
  n_mis <- sum(mis_idx)
  
  n_obs <- nrow(df_obs)
  n_mis <- sum(mis_idx)
  
  if (n_obs == 0L) {
    stop("None observed units are present.")
  }
  
  if (n_mis == 0L) {
    return(
      matrix(
        numeric(0),
        nrow = 0L,
        ncol = no.imputs
      )
    )
  }
  
  obs_weights <- if(!is.null(weights.ind)){df_obs[,weights.ind]}
  else {rep(1,n_obs)}
  
  if (
    length(obs_weights) != n_obs ||
    anyNA(obs_weights) ||
    any(!is.finite(obs_weights)) ||
    any(obs_weights <= 0)
  ) {
    stop(
      "weights_work must be positive and not missing."
    )
  }
  
  original_weight_sum <- sum(obs_weights)
  
  impute.df <- matrix(
    NA_real_,
    nrow = n_mis,
    ncol = no.imputs
  )
  
  for (m in seq_len(no.imputs)) {
    
    if (!is.null(seed)) {
      set.seed(seed + (m - 1L) * 100L)
    }
    
    # 1. Bayesian bootstrap:
    bb_weights <- rexp(
      n = n_obs,
      rate = 1
    )
    
    bb_weights <- bb_weights / sum(bb_weights)
    
    # 2. Moltiply with meta-analytics weights
    tree_weights <- obs_weights * bb_weights
    
    # 3. Make some of weights to K1
    tree_weights <- tree_weights *
      original_weight_sum / sum(tree_weights)
    
    # 4. One tree per imputation
    fit <- rpart::rpart(
      Y_obs ~ .,
      data = xobs_data,
      weights = tree_weights,
      method = "anova",
      control = rpart::rpart.control(
        minbucket = minbucket.cart,
        minsplit = 2L * minbucket.cart,
        cp = 1e-04,
        xval = 0L,
        maxsurrogate = 0L
      )
    )
    
    fit$frame$yval <- as.numeric(
      row.names(fit$frame)
    )
    
    
    leaf_obs <- as.numeric(
      predict(
        fit,
        newdata = xobs_data
      )
    )
    
    
    leaf_mis <- as.numeric(
      predict(
        fit,
        newdata = xmis_data
      )
    )
    
    # 5. donor pools from same leaf
    donor_pools <- lapply(
      leaf_mis,
      function(node) {
        
        in_leaf <- leaf_obs == node
        
        list(
          values = Y_obs[in_leaf]
        )
      }
    )
    
    # 6. extraction of the donor -- no weights
    impute.df[, m] <- vapply(
      donor_pools,
      function(pool) {
        
        n_donors <- length(pool$values)
        
        if (n_donors == 0L) {
          stop(
            "One terminal edge has no donors."
          )
        }
        
        if (n_donors == 1L) {
          return(
            as.numeric(pool$values[1L])
          )
        }
        
        donor_idx <- sample.int(
          n = n_donors,
          size = 1L
        )
        
        as.numeric(pool$values[donor_idx])
      },
      numeric(1)
    )
  }
  
  colnames(impute.df) <- paste0(
    "imp_",
    seq_len(no.imputs)
  )
  
  rownames(impute.df) <- row.names(ampdf)[mis_idx]
  
  impute.df
}

cart_imputation_boot <- function(
    ampdf,
    minbucket.cart,
    no.imputs,
    seed = NULL,
    miss.ind,
    cov.ind,
    weights.ind = NULL
) {
  if(!is.numeric(miss.ind) | !is.numeric(cov.ind)) stop("Indicators must be numbers")
  
  obs_idx <- !is.na(ampdf[,miss.ind])
  mis_idx <- is.na(ampdf[,miss.ind])
  
  Y_obs <- ampdf[obs_idx, miss.ind]
  
  df_obs <- ampdf[obs_idx, , drop = FALSE]
  
  xobs_data <- df_obs[
    , cov.ind,
    drop = FALSE
  ]
  
  xmis_data <- ampdf[
    mis_idx,
    cov.ind,
    drop = FALSE
  ]
  
  n_obs <- nrow(df_obs)
  n_mis <- sum(mis_idx)
  
  
  if (n_obs == 0L) {
    stop("None observed units are present.")
  }
  
  if (n_mis == 0L) {
    return(matrix(numeric(0), nrow = 0L, ncol = no.imputs))
  }
  
  obs_weights <- if(!is.null(weights.ind)){df_obs[,weights.ind]}
  else {rep(1,n_obs)}
  
  if (
    length(obs_weights) != n_obs ||
    anyNA(obs_weights) ||
    any(!is.finite(obs_weights)) ||
    any(obs_weights <= 0)
  ) {
    stop(
      "weights_work must be positive and not missing."
    )
  }
  
  original_weight_sum <- sum(obs_weights)
  
  impute.df <- matrix(
    NA_real_,
    nrow = n_mis,
    ncol = no.imputs
  )
  
  for (m in seq_len(no.imputs)) {
    
    if (!is.null(seed)) {
      set.seed(seed + (m - 1L) * 100L)
    }
    
    # 1. Uniform bootstrap on observed units
    boot_idx <- sample.int(
      n = n_obs,
      size = n_obs,
      replace = TRUE
    )
    
    df_boot <- df_obs[boot_idx, , drop = FALSE]
    
    Y_boot <- df_boot[, miss.ind]
    X_boot <- df_boot[, cov.ind, drop = FALSE]
    
    # 2. Use the correct weights
    boot_weights <- obs_weights[boot_idx]
    
    # 3. Weights normalization post-bootstrap:
    boot_weights <- boot_weights *
      original_weight_sum / sum(boot_weights)
    
    # 4. One tree per imputation
    fit <- rpart::rpart(
      Y_boot ~ .,
      data = X_boot,
      weights = boot_weights,
      method = "anova",
      control = rpart::rpart.control(
        minbucket = minbucket.cart,
        minsplit = 2L * minbucket.cart,
        cp = 1e-04,
        xval = 0L,
        maxsurrogate = 0L
      )
    )
    
    
    fit$frame$yval <- as.numeric(row.names(fit$frame))
    
    
    leaf_obs <- as.numeric(
      predict(fit, newdata = xobs_data)
    )
    
    leaf_mis <- as.numeric(
      predict(fit, newdata = xmis_data)
    )
    
    # 5. Creation donor pool 
    donor_pools <- lapply(
      leaf_mis,
      function(node) {
        
        in_leaf <- leaf_obs == node
        
        list(
          values = Y_obs[in_leaf]
        )
      }
    )
    
    # 6. Uniform estraction of the donor
    impute.df[, m] <- vapply(
      donor_pools,
      function(pool) {
        
        n_donors <- length(pool$values)
        
        if (n_donors == 0L) {
          stop(
            "One terminal edge has no donors."
          )
        }
        
        if (n_donors == 1L) {
          return(as.numeric(pool$values[1L]))
        }
        
        donor_idx <- sample.int(
          n = n_donors,
          size = 1L
        )
        
        as.numeric(pool$values[donor_idx])
      },
      numeric(1)
    )
  }
  
  colnames(impute.df) <- paste0("imp_", seq_len(no.imputs))
  rownames(impute.df) <- row.names(ampdf)[mis_idx]
  
  impute.df
}


cart_imputation<-function(ampdf,minbucket.cart,no.imputs,seed=NULL,
                          miss.ind,
                          cov.ind,
                          weights.ind = NULL){
  if(!is.numeric(miss.ind) | !is.numeric(cov.ind)) stop("Indicators must be numbers")
  
  obs_idx <- !is.na(ampdf[,miss.ind])
  mis_idx <- is.na(ampdf[,miss.ind])
  
  Y_obs <- ampdf[obs_idx, miss.ind]
  
  df_obs <- ampdf[obs_idx, , drop = FALSE]
  
  xobs_data <- df_obs[
    , cov.ind,
    drop = FALSE
  ]
  
  xmis_data <- ampdf[
    mis_idx,
    cov.ind,
    drop = FALSE
  ]
  
  n_obs <- nrow(df_obs)
  n_mis <- sum(mis_idx)
  
  obs_weights <- if(!is.null(weights.ind)){df_obs[,weights.ind]}
  else {rep(1,n_obs)}
  
  impute.df <- matrix(NA, nrow = sum(mis_idx), ncol = no.imputs)
  
  fit <- rpart::rpart(
    Y_obs ~ .,
    data = xobs_data,
    weights = obs_weights,
    control = rpart::rpart.control(
      minbucket = minbucket.cart,
      minsplit = 2 * minbucket.cart,
      cp = 1e-04
    )
  )
  
  fit$frame$yval <- as.numeric(row.names(fit$frame))
  
  leaf_obs <- fit$where
  
  leaf_mis <- match(
    as.character(predict(fit, newdata = xmis_data)),
    row.names(fit$frame)
  )
  
  donor_pools <- lapply(
    leaf_mis,
    function(s) {
      in_leaf <- (leaf_obs == s)
      list(
        values = Y_obs[in_leaf],
        weights = obs_weights[in_leaf]
      )
    }
  )
  
  impute.df <- matrix(
    NA,
    nrow = sum(mis_idx),
    ncol = no.imputs
  )
  
  for(m in seq_len(no.imputs)) {
    
    if(!is.null(seed))
      set.seed(seed + m*100)
    
    impute.df[, m] <- vapply(
      donor_pools,
      function(pool) {
        if (length(pool$values) == 1) {
          return(pool$values)
        }
        sample(pool$values, size = 1)
      },
      numeric(1)
    )
  }
  
  return(impute.df)
}



#funzioni
pmm_imputation<-function(ampdf,donors.pmm,no.imputs,seed=NULL,ridge = 1e-5,
                         miss.ind,
                         cov.ind,
                         weights.ind = NULL){
  if(!is.numeric(miss.ind) | !is.numeric(cov.ind)) stop("Indicators must be numbers")
  
  if(is.null(ampdf$weights_work)) 
  {ampdf$weights_work <- rep(1,nrow(ampdf))}
  
  obs_idx <- !is.na(ampdf[,miss.ind])
  mis_idx <- is.na(ampdf[,miss.ind])
  
  Y_obs <- ampdf[obs_idx, miss.ind]
  
  df_obs <- ampdf[obs_idx, , drop = FALSE]
  
  n_obs <- nrow(df_obs)
  
  xobs_data <- df_obs[
    , cov.ind,
    drop = FALSE
  ]
  
  xmis_data <- ampdf[
    mis_idx,
    cov.ind,
    drop = FALSE
  ]
  
  xobs <- cbind(as.matrix(xobs_data), intercept = 1)
  xmis <- cbind(as.matrix(xmis_data), intercept = 1)
  
  residual_df <- nrow(xobs) - ncol(xobs)
  
  obs_weights <- if(!is.null(weights.ind)){df_obs[,weights.ind]}
  else {rep(1,n_obs)}
  
  obs_weights <- nrow(xobs) * obs_weights / sum(obs_weights)
  
  impute.df <- matrix(NA, nrow = sum(mis_idx), ncol = no.imputs)
  
  sqrt_w <- sqrt(obs_weights)
  xobs_w <- (xobs * sqrt_w)
  yobs_w <- (sqrt_w * Y_obs)
  
  xtwx <- crossprod(as.matrix(xobs_w))
  penalized_xtwx <- xtwx + diag(ridge *  diag(xtwx))
  
  v_beta <- solve(penalized_xtwx)
  
  beta_hat <- v_beta %*% crossprod(xobs_w, yobs_w)
  
  residuals <- (Y_obs - xobs %*% beta_hat)
  weighted_sse <- sum(obs_weights * residuals^2)
  
  if (!is.null(seed)) set.seed(seed)
  
  for(j in seq_len(no.imputs)){
    
    chi_draw <- stats::rchisq(1, df = residual_df)
    sigma_draw <- sqrt(weighted_sse / chi_draw)
    
    eig <- eigen(v_beta)
    v_half <- eig$vectors %*%
      diag(sqrt(eig$values))
    
    beta_draw <- beta_hat +
      sigma_draw * v_half %*% rnorm(ncol(xobs))
    
    yhat_obs <- xobs %*% beta_hat
    yhat_target <- xmis %*% beta_draw
    
    imputed <- vapply(
      yhat_target,
      function(target_prediction) {
        distance <- abs(yhat_obs - target_prediction)
        
        ordered_donors <- order(
          distance,
          runif(length(distance))
        )
        donor_pool <- ordered_donors[seq_len(donors.pmm)]
        
        selected <- donor_pool[sample.int(length(donor_pool), 1)]
        Y_obs[selected]
      },
      numeric(1)
    )
    impute.df[,j] <- imputed
  }
  
  return(impute.df)
}