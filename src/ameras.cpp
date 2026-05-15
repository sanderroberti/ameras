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
#include <numeric>


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

// [[Rcpp::export]]
double compute_ERCsum_clogit(
    Rcpp::List                  &set_members,
    Eigen::Map<Eigen::VectorXd> &RRs,
    Eigen::Map<Eigen::VectorXd> &drdd,
    Eigen::Map<Eigen::VectorXd> &drdd2,
    Eigen::Map<Eigen::VectorXi> &status,
    Eigen::Map<Eigen::MatrixXd> &Kmat) {
  
  int nsets = set_members.size();
  int ncol  = RRs.size();
  
  // Compute set_RR_sums inside C++
  Eigen::VectorXd set_RR_sums(nsets);
  for (int s = 0; s < nsets; s++) {
    Rcpp::IntegerVector members = set_members[s];
    double sum = 0.0;
    for (int i = 0; i < members.size(); i++) {
      sum += RRs(members[i]);
    }
    set_RR_sums(s) = sum;
  }
  
  Eigen::VectorXd inv_set_RR  = set_RR_sums.cwiseInverse();
  Eigen::VectorXd inv_set_RR2 = inv_set_RR.cwiseProduct(inv_set_RR);
  Eigen::VectorXd inv_RRs     = RRs.cwiseInverse();
  Eigen::VectorXd inv_RRs2    = inv_RRs.cwiseProduct(inv_RRs);
  Eigen::VectorXd status_d    = status.cast<double>();
  
  // Per-individual quantities
  Eigen::VectorXd ind_inv_RR  = Eigen::VectorXd::Zero(ncol);
  Eigen::VectorXd ind_inv_RR2 = Eigen::VectorXd::Zero(ncol);
  
  for (int s = 0; s < nsets; s++) {
    Rcpp::IntegerVector members = set_members[s];
    for (int i = 0; i < members.size(); i++) {
      ind_inv_RR(members[i])  = inv_set_RR(s);
      ind_inv_RR2(members[i]) = inv_set_RR2(s);
    }
  }
  
  // Compute dldd
  Eigen::VectorXd dldd = drdd.cwiseProduct(
    status_d.cwiseProduct(inv_RRs) - ind_inv_RR
  );
  
  // Accumulate sum(mymat * Kmat) over sets
  double sum_mymat_Kmat = 0.0;
  
  for (int s = 0; s < nsets; s++) {
    Rcpp::IntegerVector members = set_members[s];
    double inv2     = inv_set_RR2(s);
    int    set_size = members.size();
    
    for (int a = 0; a < set_size; a++) {
      int ia = members[a];
      for (int b = 0; b < set_size; b++) {
        int ib = members[b];
        sum_mymat_Kmat += inv2 * drdd(ia) * drdd(ib) * Kmat(ia, ib);
      }
    }
  }
  
  // Diagonal contribution
  Eigen::VectorXd diag_update = status_d.cwiseProduct(
    drdd2.cwiseProduct(inv_RRs) -
    drdd.cwiseProduct(drdd).cwiseProduct(inv_RRs2)
  );
  Eigen::VectorXd last_term  = drdd2.cwiseProduct(ind_inv_RR);
  Eigen::VectorXd diag_total = diag_update - last_term;
  
  for (int i = 0; i < ncol; i++) {
    sum_mymat_Kmat += diag_total(i) * Kmat(i, i);
  }
  
  // sum(tcrossprod(dldd) * Kmat) = dldd^T * Kmat * dldd
  double sum_dldd_Kmat = dldd.transpose() * Kmat * dldd;
  
  return sum_mymat_Kmat + sum_dldd_Kmat;
}



// Instead of passing the full Kmat, pass the raw dose matrix
// and compute Kmat entries on the fly as needed
// [[Rcpp::export]]
double compute_ERCsum_prophaz(
    Eigen::Map<Eigen::VectorXd> &entry_t,
    Eigen::Map<Eigen::VectorXd> &exit_t,
    Eigen::Map<Eigen::VectorXd> &status_ord,
    Eigen::Map<Eigen::VectorXd> &RRs,
    Eigen::Map<Eigen::VectorXd> &drdd,
    Eigen::Map<Eigen::VectorXd> &drdd2,
    Eigen::Map<Eigen::MatrixXd> &Xc,  // precomputed centered dose matrix
    Eigen::Map<Eigen::VectorXd> &Kmat_diag, // precomputed row variances
    Eigen::Map<Eigen::VectorXd> &dldd) {
  
  int n    = entry_t.size();
  int K    = Xc.cols();
  

  
  // Kmat(i,j) = (1/(K-1)) * sum_k Xc(i,k) * Xc(j,k)
  //           = (1/(K-1)) * Xc.row(i) * Xc.row(j)^T
  // Never form full Kmat, compute dot products on the fly
  
  // Diagonal update terms
  Eigen::VectorXd status_d = status_ord;
  Eigen::VectorXd inv_RRs  = RRs.cwiseInverse();
  Eigen::VectorXd inv_RRs2 = inv_RRs.cwiseProduct(inv_RRs);
  
  
  Eigen::VectorXd diag_total =
    status_d.cwiseProduct(
      drdd2.cwiseProduct(inv_RRs) -
      drdd.cwiseProduct(drdd).cwiseProduct(inv_RRs2)
    );
  
  // Diagonal contribution
  double sum_mymat_Kmat = 0.0;
  for (int i = 0; i < n; i++) {
    sum_mymat_Kmat += diag_total(i) * Kmat_diag(i);
  }
  
  // Sort entry indices
  std::vector<int> entry_order(n);
  std::iota(entry_order.begin(), entry_order.end(), 0);
  std::sort(entry_order.begin(), entry_order.end(),
            [&](int a, int b) { return entry_t(a) < entry_t(b); });
  
  // Risk set state
  std::vector<bool> at_risk(n, false);
  double            risk_RR_sum = 0.0;
  
  // Running sum: sum_{j in risk set} drdd(j) * Kmat(:,j)
  // Kmat(:,j) = Xc * Xc.row(j)^T / (K-1)
  // sum_{j in risk} drdd(j) * Kmat(:,j)
  // = Xc * (sum_{j in risk} drdd(j) * Xc.row(j)^T) / (K-1)
  // = Xc * risk_drdd_Xc^T / (K-1)
  // where risk_drdd_Xc = sum_{j in risk} drdd(j) * Xc.row(j)
  
  // So instead of maintaining an N-vector sum_drdd_Kmat_col
  // we maintain a K-vector risk_drdd_Xc
  // This reduces memory from O(N) to O(K) for this running sum
  Eigen::VectorXd risk_drdd_Xc = Eigen::VectorXd::Zero(K);
  
  int entry_ptr = 0;
  int exit_ptr  = 0;
  
  for (int i = 0; i < n; ++i) {
    double t = exit_t(i);
    
    // Remove exited individuals
    while (exit_ptr < i && exit_t(exit_ptr) < t) {
      int j = exit_ptr;
      if (at_risk[j]) {
        at_risk[j]      = false;
        risk_RR_sum    -= RRs(j);
        risk_drdd_Xc   -= drdd(j) * Xc.row(j).transpose();
      }
      exit_ptr++;
    }
    
    // Add new entrants
    while (entry_ptr < n && entry_t(entry_order[entry_ptr]) < t) {
      int j = entry_order[entry_ptr];
      if (!at_risk[j]) {
        at_risk[j]      = true;
        risk_RR_sum    += RRs(j);
        risk_drdd_Xc   += drdd(j) * Xc.row(j).transpose();
      }
      entry_ptr++;
    }
    
    if (status_ord(i) > 0.5 && risk_RR_sum > 0.0) {
      double inv2        = 1.0 / (risk_RR_sum * risk_RR_sum);
      double inv_risk_RR = 1.0 / risk_RR_sum;
      
      // sum_{a,b in risk} inv2 * drdd(a) * drdd(b) * Kmat(a,b)
      // = inv2 * sum_{a in risk} drdd(a) * sum_{b in risk} drdd(b) * Kmat(a,b)
      // = inv2 * sum_{a in risk} drdd(a) * (Xc.row(a) * risk_drdd_Xc) / (K-1)
      // = inv2 / (K-1) * sum_{a in risk} drdd(a) * Xc.row(a) * risk_drdd_Xc
      // = inv2 / (K-1) * risk_drdd_Xc^T * risk_drdd_Xc
      // = inv2 / (K-1) * ||risk_drdd_Xc||^2
      
      sum_mymat_Kmat += inv2 / (K - 1) * risk_drdd_Xc.squaredNorm();
      
      // last_term diagonal correction
      // sum_{a in risk} drdd2(a) * inv_risk_RR * Kmat(a,a)
      for (int j = 0; j < n; ++j) {
        if (at_risk[j]) {
          sum_mymat_Kmat -= drdd2(j) * inv_risk_RR * Kmat_diag(j);
        }
      }
    }
  }
  
  // sum(tcrossprod(dldd) * Kmat)
  // = dldd^T * Kmat * dldd
  // = dldd^T * Xc * Xc^T * dldd / (K-1)
  // = ||Xc^T * dldd||^2 / (K-1)
  // Never forms N x N Kmat
  Eigen::VectorXd Xc_dldd = Xc.transpose() * dldd;
  double sum_dldd_Kmat = Xc_dldd.squaredNorm() / (K - 1);
  
  return sum_mymat_Kmat + sum_dldd_Kmat;
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




