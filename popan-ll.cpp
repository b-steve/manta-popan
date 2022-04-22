#include <Rcpp.h>
#include <math.h>
#include <iostream>

using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
double popan_ll(int k, NumericVector phivec, NumericVector pvec, NumericVector ptrvec_full,
		       NumericVector pentvec, double Ns, double nhist, IntegerVector first_obs,
		       IntegerVector last_obs, IntegerMatrix det_dat){
  NumericVector chivec_trans(k);
  chivec_trans.fill(1);
  NumericVector chivec_resid(k);
  chivec_resid(k - 1) = 1;
  for (int i = k - 2; i >= 0; i--){
    chivec_resid(i) = 1 - phivec(i) + phivec(i)*(1 - pvec(i + 1))*chivec_resid(i + 1);
  }
  // Log probability of not being detected.
  double log_p_unseen = log(sum(((1 - pvec)*ptrvec_full + (1 - pvec)*chivec_resid*(1 - ptrvec_full))*pentvec));
  // Likelihood contribution due to nhist.
  double nllike = -lgamma(Ns + 1) + lgamma(Ns - nhist + 1) - (Ns - nhist) * log_p_unseen + lgamma(nhist + 1);
  LogicalVector first_f;
  // Considering all animals first detected in occasion f.
  for (int f = 0; f < k; f++){
    // Counting number of animals first detected in occasion f.
    int n_f = 0;
    for (int i = 0; i < nhist; i++){
      if (first_obs(i) == (f + 1)){
	n_f++;
      }
    }
    if (n_f > 0){
      // Grabbing indices for animals first detected in occasion f.
      IntegerVector f_indices(n_f);
      IntegerVector last_dat_f(n_f);
      int j = 0;
      for (int i = 0; i < nhist; i++){
	if (first_obs(i) == (f + 1)){
	  f_indices(j) = i;
	  last_dat_f(j) = last_obs(i);
	  j++;
	}
      }
      // Putting together capture histories for animals first detected in occasion f.
      IntegerMatrix dat_f(n_f, k);
      for (int i = 0; i < n_f; i++){
	dat_f(i, _) = det_dat.row(f_indices(i));
      }
      // Calculating the probability of not being detected on occasions i,
      // ..., f - 1, given alive at i.
      NumericVector psif_vec(f + 1);
      psif_vec(f) = 1; // Check this for f > 0.
      if (f > 0){
	for (int i = f - 1; i >= 0; i--){
	  psif_vec(i) = (1 - pvec(i))*phivec(i)*psif_vec(i + 1);
	}
      }
      // Log-likelihood contribution for each animal.
      for (int hst = 0; hst < n_f; hst++){
	int last_hst = last_dat_f(hst);
	// This is weird but C++ is annoying and wouldn't coerce a
	// single-row matrix into a vector in one line.
	IntegerMatrix datvec_hst_mat = dat_f(Range(hst, hst), Range(f, last_hst - 1));
	IntegerVector datvec_hst = datvec_hst_mat(0, _);
	// Getting history-specific parameter vectors.
	NumericVector p_hst(last_hst - f);
	j = 0;
	for (int i = f; i < last_hst; i++){
	  
	  p_hst(j) = pvec(i);
	  j++;
	}
	NumericVector phi_hst(last_hst - f - 1);
	j = 0;
	for (int i = f; i < last_hst - 1; i++){
	  phi_hst(j) = phivec(i);
	  j++;
	}
	NumericVector prob_resid_et_start(f + 1);
	for (int i = 0; i <= f; i++){
	  prob_resid_et_start(i) = pentvec(i)*(1 - ptrvec_full(i))*psif_vec(i);
	}
	double log_prob_resid_hst = 0;
	for (int i = 0; i < datvec_hst.size(); i++){
	  log_prob_resid_hst += datvec_hst(i)*log(p_hst(i)) + (1 - datvec_hst(i))*log(1 - p_hst(i));
	}
	if (sum(datvec_hst) > 1){
	  log_prob_resid_hst += sum(log(phi_hst));
	}
	double prob_resid_hst = exp(log_prob_resid_hst);
	double prob_resid_end = chivec_resid(last_hst - 1);
	NumericVector prob_resid_et(f + 1);
	for (int i = 0; i <= f; i++){
	  prob_resid_et(i) = prob_resid_et_start(i)*prob_resid_hst*prob_resid_end;
	}
	NumericVector prob_trans_et(f + 1);
	prob_trans_et.fill(0);
	if (sum(datvec_hst) == 1){
	  prob_trans_et(f) = pentvec(f)*ptrvec_full(f)*pvec(f);
	}
	nllike -= log(sum(prob_resid_et + prob_trans_et));
      }
    }
  }
  return nllike;
}
