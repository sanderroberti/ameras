// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <Rmath.h>
#include <iterator>
#include <iostream>
#include <vector>
#include <algorithm>


// [[Rcpp::export]]
Eigen::RowVectorXd dldd_clogit(Eigen::Map<Eigen::MatrixXd> &designmat,
                                        Eigen::Map<Eigen::VectorXd> &RRs) {
  // Compute denominator vector: for each row i, sum over columns: designmat(i,j)*RRs[j]
  Eigen::VectorXd denom = designmat * RRs;  // (nrow x 1) vector
  
  
  // Create a matrix where each row i is designmat.row(i) / denom[i]
  Eigen::MatrixXd denom_mat = denom.replicate(1, designmat.cols()); // replicate denom by cols
  
  Eigen::MatrixXd frac = designmat.array() / denom_mat.array();
  
  // Sum over rows for each column -> tmp is row vector with length ncol
  Eigen::RowVectorXd tmp = frac.colwise().sum();
  
  return tmp;
}

// [[Rcpp::export]]
Eigen::RowVectorXd dldd_prophaz(Eigen::Map<Eigen::VectorXd> &entry_t,
                                         Eigen::Map<Eigen::VectorXd> &exit_t,
                                         Eigen::Map<Eigen::VectorXi> &status_ord,
                                         Eigen::Map<Eigen::VectorXd> &RRs) {
  int n = entry_t.size();  // Number of individuals
  Eigen::RowVectorXd result(n);  // Store the result for each individual
  result.setZero();
  
  // Loop through individuals who had an event (delta_i = 1)
  for (int i = 0; i < n; ++i) {
    if (status_ord(i) == 1) {  // Only consider individuals who had an event (delta_i = 1)
      double sum_RR = 0.0;
      // Compute the sum of relative risks in the risk set for individual i
      for (int j = 0; j < n; ++j) {
        if (exit_t(i) > entry_t(j) && exit_t(i) <= exit_t(j)) {
          sum_RR += RRs(j);
        }
      }
      
      // For each individual k, add the contribution of individual i to the result if k is in the risk set
      for (int k = 0; k < n; ++k) {
        if (exit_t(i) > entry_t(k) && exit_t(i) <= exit_t(k)) {
          // Only add if individual k is in the risk set of i
          result(k) += 1.0 / sum_RR;
        }
      }
    }
  }
  
  return result;
}


// [[Rcpp::export]]
Eigen::MatrixXd compute_ERCmatrix_prophaz(Eigen::Map<Eigen::VectorXd> &entry_t,
                                       Eigen::Map<Eigen::VectorXd> &exit_t,
                                       Eigen::Map<Eigen::VectorXi> &status_ord,
                                       Eigen::Map<Eigen::VectorXd> &RRs,
                                       Eigen::Map<Eigen::VectorXd> &drdd,
                                       Eigen::Map<Eigen::VectorXd> &drdd2) {
  int n = entry_t.size();
  
  // Collect all event indices
  std::vector<int> event_indices;
  event_indices.reserve(n);
  for (int i = 0; i < n; ++i) {
    if (status_ord(i) == 1) {
      event_indices.push_back(i);
    }
  }
  
  int n_events = event_indices.size();
  if (n_events == 0) return Eigen::MatrixXd::Zero(n, n);
  
  // Diagonal update terms
  Eigen::VectorXd status_d = status_ord.cast<double>();
  Eigen::VectorXd diag_update = (drdd2.array() / RRs.array()) * status_d.array() -
    (drdd.array().square() / RRs.array().square()) * status_d.array();
  
  Eigen::MatrixXd result = Eigen::MatrixXd::Zero(n, n);
  
  // Decide whether to pre-compute outer_drdd based on n
  // Threshold: if n*n*8 bytes > 500 MB, use blocking
  bool use_blocking = (n > 8000);  // ~500 MB threshold
  
  Eigen::MatrixXd outer_drdd;
  if (!use_blocking) {
    outer_drdd = drdd * drdd.transpose();
  }
  
  // Process events in batches
  int event_batch_size = 1000;
  
  for (int batch_start = 0; batch_start < n_events; batch_start += event_batch_size) {
    int batch_end = std::min(batch_start + event_batch_size, n_events);
    int current_batch_size = batch_end - batch_start;
    
    // Build atriskmat for this batch
    Eigen::MatrixXd atriskmat(current_batch_size, n);
    
    for (int row = 0; row < current_batch_size; ++row) {
      int event_idx = event_indices[batch_start + row];
      double exit_event = exit_t(event_idx);
      
      for (int col = 0; col < n; ++col) {
        atriskmat(row, col) = (exit_t(col) >= exit_event && entry_t(col) < exit_event) ? 1.0 : 0.0;
      }
    }
    
    // Compute atriskmat_prod_RRs
    Eigen::VectorXd atriskmat_prod_RRs = atriskmat * RRs;
    Eigen::VectorXd inv_sq = atriskmat_prod_RRs.array().square().inverse();
    
    // Scale rows
    Eigen::MatrixXd scaled_atrisk(current_batch_size, n);
    for (int row = 0; row < current_batch_size; ++row) {
      scaled_atrisk.row(row) = atriskmat.row(row) * inv_sq(row);
    }
    
    // Compute inner
    Eigen::MatrixXd inner = atriskmat.transpose() * scaled_atrisk;
    
    if (use_blocking) {
      // Use blocking for large n
      int block_size = 1000;
      for (int i_start = 0; i_start < n; i_start += block_size) {
        int i_end = std::min(i_start + block_size, n);
        for (int j_start = 0; j_start < n; j_start += block_size) {
          int j_end = std::min(j_start + block_size, n);
          
          Eigen::MatrixXd inner_block = inner.block(i_start, j_start, i_end - i_start, j_end - j_start);
          Eigen::VectorXd drdd_i = drdd.segment(i_start, i_end - i_start);
          Eigen::VectorXd drdd_j = drdd.segment(j_start, j_end - j_start);
          Eigen::MatrixXd outer_block = drdd_i * drdd_j.transpose();
          
          result.block(i_start, j_start, i_end - i_start, j_end - j_start).array() += 
            inner_block.array() * outer_block.array();
        }
      }
    } else {
      // Use pre-computed outer product for small n
      result.array() += inner.array() * outer_drdd.array();
    }
    
    // atriskmat scaled by drdd2
    Eigen::MatrixXd atrisk_drdd2(current_batch_size, n);
    for (int row = 0; row < current_batch_size; ++row) {
      for (int col = 0; col < n; ++col) {
        atrisk_drdd2(row, col) = atriskmat(row, col) * drdd2(col);
      }
    }
    
    // last_term
    Eigen::VectorXd inv_atriskmat_prod_RRs = atriskmat_prod_RRs.array().inverse();
    Eigen::RowVectorXd last_term = inv_atriskmat_prod_RRs.transpose() * atrisk_drdd2;
    
    result.diagonal().array() -= last_term.transpose().array();
  }
  
  result.diagonal() += diag_update;
  
  return result;
}




// [[Rcpp::export]]
Eigen::MatrixXd compute_ERCmatrix_clogit(Eigen::Map<Eigen::MatrixXd> &designmat,
                                      Eigen::Map<Eigen::VectorXd> &RRs,
                                      Eigen::Map<Eigen::VectorXd> &drdd,
                                      Eigen::Map<Eigen::VectorXd> &drdd2,
                                      Eigen::Map<Eigen::VectorXi> &status) {
  
  int ncol = designmat.cols();
  
  // Compute designmat_prod_RRs = designmat * RRs
  Eigen::VectorXd designmat_prod_RRs = designmat * RRs;
  
  // Compute inv_sq vector
  Eigen::VectorXd inv_sq = designmat_prod_RRs.array().square().inverse();
  
  // Scale rows of designmat by inv_sq
  Eigen::MatrixXd scaled_atrisk = designmat.array().colwise() * inv_sq.array();
  
  // Compute inner matrix
  Eigen::MatrixXd inner = designmat.transpose() * scaled_atrisk;
  
  // Outer product drdd * drdd.transpose()
  Eigen::MatrixXd outer_drdd = drdd * drdd.transpose();
  
  // Compute mymat before diag update
  Eigen::MatrixXd mymat = inner.array() * outer_drdd.array();
  
  // Diagonal update terms
  Eigen::VectorXd status_d = status.cast<double>();
  
  Eigen::VectorXd diag_update = (drdd2.array() / RRs.array()) * status_d.array() -
    (drdd.array().square() / RRs.array().square()) * status_d.array();
  
  // designmat scaled columns by drdd2
  Eigen::MatrixXd atrisk_drdd2 = designmat.array().rowwise() * drdd2.transpose().array();
  
  // last_term
  Eigen::RowVectorXd last_term = (designmat_prod_RRs.array().inverse().matrix().transpose()) * atrisk_drdd2;
  
  // Update diagonal elements
  for (int i = 0; i < ncol; i++) {
    mymat(i, i) += diag_update(i) - last_term(i);
  }

  return mymat;
}


// [[Rcpp::export]]
Eigen::VectorXd loglik_prophaz_rcpp(
    const Eigen::VectorXd& exit_t,
    const Eigen::VectorXd& entry_t,
    const Eigen::MatrixXd& RR_entry,
    const Eigen::MatrixXd& RR_exit,
    const Eigen::VectorXi& status_ord,
    double loglim) {
  
  const int n = exit_t.size();
  const int K = RR_exit.cols();
  
  Eigen::VectorXd loglik = Eigen::VectorXd::Zero(K);
  Eigen::VectorXd risk_sum = Eigen::VectorXd::Zero(K);
  
  int ptr_entry = 0;
  int i = 0;
  
  while (i < n) {
    
    const double ti = exit_t[i];
    
    // add entrants
    while (ptr_entry < n && entry_t[ptr_entry] < ti) {
      risk_sum.noalias() += RR_entry.row(ptr_entry).transpose();
      ptr_entry++;
    }
    
    // find tie group
    int idx_end = i;
    while (idx_end < n && exit_t[idx_end] == ti) idx_end++;
    
    const int tie_size = idx_end - i;
    if (tie_size > 0) {
      
      // Extract tie block
      Eigen::MatrixXd tie_RR = RR_exit.block(i, 0, tie_size, K);
      Eigen::VectorXd tie_mask = status_ord.segment(i, tie_size).cast<double>();
      
      // sum log(RR_exit) for events
      Eigen::VectorXd sum_log = (tie_RR.array().log().colwise() * tie_mask.array()).colwise().sum();
      
      int d = static_cast<int>(tie_mask.sum()); // number of events
      
      if (d > 0) {
        loglik += sum_log - d * risk_sum.array().max(loglim).log().matrix();
      }
      // remove all exiting at ti
      Eigen::VectorXd tie_exit_sum = tie_RR.colwise().sum();
      risk_sum.noalias() -= tie_exit_sum;
    }
    i = idx_end;
  }
  return loglik;
}








extern "C" {
  void C_dldd_clogit(double *R_atriskColVec, int *R_nr, int *R_nc, double *R_RR,
                              double *ret) {
    
    int nr = *R_nr;
    int nc = *R_nc;
    
    Eigen::Map<Eigen::VectorXd> RRs(R_RR, nc);
    Eigen::Map<Eigen::MatrixXd> designmat(R_atriskColVec, nr, nc);
    
    Eigen::RowVectorXd tmp = dldd_clogit(designmat, RRs);  //Eigen::RowVectorXd is returned
    
    for (int i=0; i<nc; i++) ret[i] = tmp(i);
    
    return;
  }
  
}

extern "C" {
  void C_dldd_prophaz(double *R_entry, double *R_exit, int *R_status, double *R_RR,
                               int *R_n, double *ret) {
    
    int n = *R_n;  // Number of individuals
    
    // Map the input data from R to Eigen vectors/matrices
    Eigen::Map<Eigen::VectorXd> entry_t(R_entry, n);
    Eigen::Map<Eigen::VectorXd> exit_t(R_exit, n);
    Eigen::Map<Eigen::VectorXi> status_ord(R_status, n);
    Eigen::Map<Eigen::VectorXd> RRs(R_RR, n);
    
    // Call the C++ function to compute the result
    Eigen::RowVectorXd result = dldd_prophaz(entry_t, exit_t, status_ord, RRs);
    
    // Copy the result back to the output pointer (ret)
    for (int i = 0; i < result.size(); ++i) {
      ret[i] = result(i);
    }
    
    return;
  }
}




extern "C" {
  void C_compute_ERCmatrix_clogit(double *R_atriskColVec, int *R_nr, int *R_nc, double *R_RR, 
                               double *R_drdd, double *R_drdd2, int *R_status, double *ret) {
    
    int nr = *R_nr;
    int nc = *R_nc;
    
    Eigen::Map<Eigen::VectorXd> RRs(R_RR, nc);
    Eigen::Map<Eigen::VectorXd> drdd(R_drdd, nc);
    Eigen::Map<Eigen::VectorXd> drdd2(R_drdd2, nc);
    Eigen::Map<Eigen::VectorXi> status(R_status, nc);
    Eigen::Map<Eigen::MatrixXd> designmat(R_atriskColVec, nr, nc);
    
    Eigen::MatrixXd mymat = compute_ERCmatrix_clogit(designmat, RRs, drdd, drdd2, status); //nc x nc
    
    // Return column vector
    double *ptr = ret;
    for (int j=0; j<nc; j++) {
      for (int i=0; i<nc; i++) *ptr++ = mymat(i, j);
    }
    
    return;
  }
}


extern "C" {
  void C_compute_ERCmatrix_prophaz(double *R_entry, double *R_exit, int *R_status, double *R_RR,
                                double *R_drdd, double*R_drdd2,  int *R_n, double *ret) {
    
    int n = *R_n;  // Number of individuals
    
    // Map the input data from R to Eigen vectors/matrices
    Eigen::Map<Eigen::VectorXd> entry_t(R_entry, n);
    Eigen::Map<Eigen::VectorXd> exit_t(R_exit, n);
    Eigen::Map<Eigen::VectorXi> status_ord(R_status, n);
    Eigen::Map<Eigen::VectorXd> RRs(R_RR, n);
    Eigen::Map<Eigen::VectorXd> drdd(R_drdd, n);
    Eigen::Map<Eigen::VectorXd> drdd2(R_drdd2, n);
    
    // Call the C++ function to compute the result
    Eigen::MatrixXd mymat = compute_ERCmatrix_prophaz(entry_t, exit_t, status_ord, RRs, drdd, drdd2);
    
    double *ptr = ret;
    for (int j=0; j<n; j++) {
      for (int i=0; i<n; i++) *ptr++ = mymat(i, j);
    }
    
    return;
  }
}



extern "C" {
  
  void C_loglik_prophaz_rcpp(
      double *exit_t,
      double *entry_t,
      double *RR_entry,
      double *RR_exit,
      int    *status_ord,
      int    *n,
      int    *K,
      double *loglim,
      double *out
  ) {
    const int N = *n;
    const int P = *K;
    
    // Map inputs
    Eigen::Map<Eigen::VectorXd> exit_t_e(exit_t, N);
    Eigen::Map<Eigen::VectorXd> entry_t_e(entry_t, N);
    Eigen::Map<Eigen::MatrixXd> RR_entry_e(RR_entry, N, P);
    Eigen::Map<Eigen::MatrixXd> RR_exit_e(RR_exit, N, P);
    Eigen::Map<Eigen::VectorXi> status_e(status_ord, N);
    
    // Call Rcpp function
    Eigen::VectorXd res =
      loglik_prophaz_rcpp(exit_t_e,
                          entry_t_e,
                          RR_entry_e,
                          RR_exit_e,
                          status_e,
                          *loglim);
    
    // Copy result back to C array
    for (int j = 0; j < P; ++j) {
      out[j] = res(j);
    }
  }
  
}




