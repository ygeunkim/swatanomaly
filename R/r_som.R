#' Matrix Representation of Sliding Windows using PCA
#'
#' @description
#' Bind every window by row,
#' using principal components analysis in each window.
#' @param x mat multivariate series.
#' @param win int window size.
#' @param jump int jump size of sliding window.
#' @param parallel If TRUE, use parallel operation with \link[foreach]{foreach}.
#' To use this operation, parallel backend must be registered beforehand. By default, FALSE.
#' @return mat
#' @details
#' Each column represents an observation in each window.
#' Each row represents an distinct window.
#' Given the size of window win and the jump size jump, the number of the row becomes (nrow - win) / jump + 1.
#' The data set is multivariate series, in general, so we perform dimension reduction using PCA.
#' @import foreach
#' @export
reduce_pca <- function(x, win, jump, parallel = FALSE) {
  win_num <- as.integer((nrow(x) - win) / jump + 1)
  out <- matrix(nrow = win_num, ncol = win)
  win_mat <- matrix(nrow = win, ncol = ncol(x))
  i <- integer(win_num)

  if (!parallel) {
    out <-
      foreach(i = seq_len(win_num), .combine = rbind) %do% {
        win_mat <- x[(1 + (i - 1) * jump):((i - 1) * jump + win),]
        (win_mat %*% svd(win_mat)$v)[,1]
      }
  } else {
    out <-
      foreach(i = seq_len(win_num), .combine = rbind, .packages = c("swatanomaly")) %dopar% {
        win_mat <- x[(1 + (i - 1) * jump):((i - 1) * jump + win),]
        t((win_mat %*% svd(win_mat)$v)[,1])
      }
  }

  out
}

#' Matrix Representation of Sliding Windows using MDS
#'
#' @description
#' Bind every window by row,
#' using multidimensional scaling in each window.
#' @param x mat multivariate series.
#' @param win int window size.
#' @param jump int jump size of sliding window.
#' @param parallel If TRUE, use parallel operation with \link[foreach]{foreach}.
#' @return mat
#' @details
#' Each column represents an observation in each window.
#' Each row represents an distinct window.
#' Given the size of window win and the jump size jump, the number of the row becomes (nrow - win) / jump + 1.
#' The data set is multivariate series, in general, so we perform dimension reduction using MDS.
#' @import foreach
#' @export
reduce_mds <- function(x, win, jump, parallel = FALSE) {
  win_num <- as.integer((nrow(x) - win) / jump + 1)
  out <- matrix(nrow = win_num, ncol = win)
  win_mat <- matrix(nrow = win, ncol = ncol(x))
  i <- integer(win_num)

  b <- x %*% t(x) # subset i * jump x i * jump to i * jump + win - 1 x i * jump + win - 1 in the loop
  b_spec <- eigen(b[1:win, 1:win]) # pre-allocate

  if (!parallel) {
    out <-
      foreach(i = seq_len(win_num), .combine = rbind) %do% {
        b_spec <-
          eigen(
            x[(1 + (i - 1) * jump):((i - 1) * jump + win), (1 + (i - 1) * jump):((i - 1) * jump + win)]
          )
        b_spec$vectors[,1] %*% diag(sqrt(b_spec$values[1]))
      }
  } else {
    out <-
      foreach(i = seq_len(win_num), .combine = rbind, .packages = c("swatanomaly")) %dopar% {
        win_mat <- x[(1 + (i - 1) * jump):((i - 1) * jump + win),]
        t((win_mat %*% svd(win_mat)$v)[,1])
      }
  }

  out
}
