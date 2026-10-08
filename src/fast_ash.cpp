// mr.ash.rss without per-coordinate containers (see R/susie_fast_ash.R).
//
// mr_ash_sufficient() in mr_ash_rss.h passes the results of each
// coordinate update through unordered_map<string, ...> containers (one
// map per mixture component and one per coordinate) and allocates two
// p-vectors per coordinate. At p = 1000 that bookkeeping takes most of the
// time of a sweep. The functions below are the same code with plain
// structs and buffers allocated once per call. Every arithmetic expression
// is unchanged and evaluated in the same order, so the results are
// bit-identical. R/susie_fast_ash.R checks this at run time against
// mr_ash_rss_cpp and uses mr_ash_rss_cpp if a platform ever disagrees.

#include <cpp11.hpp>
#include <cpp11armadillo.hpp>
#include "mr_ash_rss.h"

using namespace cpp11;
using namespace arma;
using namespace std;

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

// bayes_ridge_sufficient()
struct fast_ridge_out {
	double bhat, s2, mu1, sigma2_1, logbf;
};

inline void fast_bayes_ridge_sufficient(double xTx, double xTy, double sigma2_e, double sigma2_0, fast_ridge_out& r) {
	double bhat = xTy / xTx;
	double s2 = sigma2_e / xTx;
	double sigma2_1 = 1 / (1 / s2 + 1 / sigma2_0);
	double mu1 = sigma2_1 / s2 * bhat;
	double logbf = log(s2 / (sigma2_0 + s2)) / 2 + (pow(bhat, 2) / s2 - pow(bhat, 2) / (sigma2_0 + s2)) / 2;
	r.bhat = bhat;
	r.s2 = s2;
	r.mu1 = mu1;
	r.sigma2_1 = sigma2_1;
	r.logbf = logbf;
}

// bayes_mix_sufficient(). `out` is the K x 5 work matrix; log_w0 = log(w0),
// which changes only between sweeps.
struct fast_mix_out {
	vec w1, mu1_k, sigma2_1_k;
	double mu1, sigma2_1, logbf;
};

inline void fast_bayes_mix_sufficient(double xTx, double xTy, double sigma2_e, const vec& w0, const vec& log_w0,
                                      const vec& sigma2_0, mat& out, fast_mix_out& m) {
	int K = sigma2_0.n_elem;
	fast_ridge_out r;
	for (int i = 0; i < K; i++) {
		fast_bayes_ridge_sufficient(xTx, xTy, sigma2_e, sigma2_e * sigma2_0[i], r);
		out(i, 0) = r.bhat;
		out(i, 1) = r.s2;
		out(i, 2) = r.mu1;
		out(i, 3) = r.sigma2_1;
		out(i, 4) = r.logbf;
	}

	m.w1 = softmax_rss(out.col(4) + log_w0);

	m.mu1_k = out.col(2);
	m.sigma2_1_k = out.col(3);

	m.mu1 = sum(m.w1 % m.mu1_k);
	m.sigma2_1 = sum(m.w1 % (square(m.mu1_k) + m.sigma2_1_k)) - pow(m.mu1, 2);

	double u = max(out.col(4));
	m.logbf = u + log(sum(w0 % exp(out.col(4) - u)));
}

// mr_ash_sufficient(). The per-component posterior means and variances
// (mu1_k_t, sigma2_1_k_t upstream) are not returned, so they are not stored.
inline unordered_map<string, mat> fast_mr_ash_sufficient(const vec& XTy, const mat& XTX, double yTy, int n, double& sigma2_e,
                                                         const vec& sigma2_0, vec& w0, const vec& mu1_init, double tol,
                                                         int max_iter, bool update_w0, bool update_sigma,
                                                         bool compute_ELBO) {
	int p = XTX.n_cols;
	int K = sigma2_0.n_elem;
	vec mu1_t = mu1_init;
	vec sigma2_1_t(p, fill::zeros);
	mat w1_t(p, K, fill::zeros);
	int t = 0;
	double ELBO = 0;
	vec varobj_vec(max_iter, fill::zeros);
	bool converged = false;

	vec XTrbar_j(p);
	mat out(K, 5);
	fast_mix_out bfit;

	while (!converged) {
		double var_part_ERSS = 0;
		double neg_KL = 0;

		t++;

		if (t > max_iter) {
			t = max_iter;  // Clamp to valid index range
			cerr << "Max number of iterations reached. Try increasing max_iter." << endl;
			break;
		}

		vec mu1_tminus1 = mu1_t;

		vec XTrbar = XTy - XTX * mu1_t;
		vec log_w0 = log(w0);

		for (int j = 0; j < p; j++) {
			// Remove j-th effect from expected residuals
			XTrbar_j = XTrbar + XTX.col(j) * mu1_t[j];

			double xTrbar_j = XTrbar_j[j];
			double xTx = XTX(j, j);

			fast_bayes_mix_sufficient(xTx, xTrbar_j, sigma2_e, w0, log_w0, sigma2_0, out, bfit);

			mu1_t[j] = bfit.mu1;
			sigma2_1_t[j] = bfit.sigma2_1;
			w1_t.row(j) = bfit.w1.t();

			if (compute_ELBO) {
				var_part_ERSS += sigma2_1_t[j] * xTx;
				neg_KL += bfit.logbf + (1 / (2 * sigma2_e)) * (-2 * xTrbar_j * mu1_t[j] + (xTx * (sigma2_1_t[j] + pow(mu1_t[j], 2))));
			}

			// Update expected residuals
			XTrbar = XTrbar_j - XTX.col(j) * mu1_t[j];
		}

		// w0 aliases the caller's R vector, which upstream updates in place.
		if (update_w0) {
			w0 = sum(w1_t, 0).t() / p;
		}

		double beta_diff = norm(mu1_t - mu1_tminus1, 2);
		double beta_norm = norm(mu1_t, 2);

		double ERSS = yTy - 2 * dot(XTy, mu1_t) + as_scalar(mu1_t.t() * XTX * mu1_t) + var_part_ERSS;
		if (compute_ELBO) {
			ELBO = -0.5 * log(n) - 0.5 * n * log(2 * datum::pi * sigma2_e) - (1 / (2 * sigma2_e)) * ERSS + neg_KL;
		}
		varobj_vec[t - 1] = ELBO;

		if (update_sigma) {
			sigma2_e = (yTy - dot(XTy, mu1_t)) / n;
		}

		if (t >= 2 && beta_diff < tol * max(1.0, beta_norm)) {
			converged = true;
		}
	}

	return {{"mu1", mat(mu1_t)}, {"sigma2_1", mat(sigma2_1_t)}, {"w1", w1_t},
		{"sigma2_e", mat(1, 1, fill::value(sigma2_e))}, {"w0", mat(w0)},
		{"ELBO", mat(1, 1, fill::value(ELBO))},
		{"iter", mat(1, 1, fill::value((double)t))},
		{"varobj", mat(varobj_vec.subvec(0, t - 1))}};
}

// mr_ash_rss(), calling fast_mr_ash_sufficient().
inline unordered_map<string, mat> fast_mr_ash_rss([[maybe_unused]] const vec& bhat, const vec& shat, const vec& z, const mat& R,
                                                  double var_y, int n, double sigma2_e, const vec& s0, vec& w0,
                                                  const vec& mu1_init, double tol, int max_iter, bool update_w0,
                                                  bool update_sigma, bool compute_ELBO, bool standardize) {
	int p = z.n_elem;

	vec mu1_init_use = mu1_init;
	if (mu1_init.is_empty()) {
		mu1_init_use = vec(p, fill::zeros);
	}

	vec z_use = z;

	vec adj(p, fill::ones);
	if (std::isfinite(n)) {
		adj = (n - 1) / (square(z_use) + n - 2);
		z_use %= sqrt(adj);
	}
	mat XtX;
	vec Xty;
	if (std::isfinite(var_y) && !shat.is_empty()) {
		vec XtXdiag = var_y * adj / square(shat);
		XtX = diagmat(sqrt(XtXdiag)) * R * diagmat(sqrt(XtXdiag));
		XtX = 0.5 * (XtX + XtX.t());
		Xty = z_use % sqrt(adj) % (var_y / shat);
	} else {
		XtX = (n - 1) * R;
		Xty = z_use * sqrt(n - 1);
		var_y = 1.0;
	}

	vec sx(p, fill::ones);
	if (standardize) {
		vec dXtX = XtX.diag();
		sx = sqrt(dXtX / (n - 1));
		sx.replace(0, 1);
		XtX = diagmat(1 / sx) * XtX * diagmat(1 / sx);
		Xty /= sx;
		mu1_init_use %= sx;
	}

	unordered_map<string, mat> result = fast_mr_ash_sufficient(Xty, XtX, var_y * (n - 1), n, sigma2_e, s0, w0, mu1_init_use,
	                                                           tol, max_iter, update_w0, update_sigma, compute_ELBO);

	if (standardize) {
		unordered_map<string, mat> out_adj = rescale_post_mean_covar(vectorise(result["mu1"]), vectorise(result["sigma2_1"]), sx);
		result["mu1"] = out_adj["mu1_orig"];
		result["sigma2_1"] = out_adj["sigma2_1_orig"];
	}

	return {{"mu1", result["mu1"]}, {"sigma2_1", result["sigma2_1"]}, {"w1", result["w1"]},
		{"sigma2_e", result["sigma2_e"]}, {"w0", result["w0"]}, {"ELBO", result["ELBO"]},
		{"iter", result["iter"]}, {"varobj", result["varobj"]}};
}

// Same arguments and result as mr_ash_rss_cpp().
[[cpp11::register]]
writable::list fast_mr_ash_rss_cpp(const doubles& bhat, const doubles& shat, const doubles& z,
                                   const doubles_matrix<>& R, double var_y, int n, double sigma2_e,
                                   const doubles& s0, const doubles& w0, const doubles& mu1_init,
                                   double tol, int max_iter, bool update_w0, bool update_sigma,
                                   bool compute_ELBO, bool standardize) {
	vec bhat_vec = as_Col(bhat);
	vec shat_vec = as_Col(shat);
	vec z_vec = as_Col(z);
	mat R_mat = as_Mat(R);
	vec s0_vec = as_Col(s0);
	vec w0_vec = as_Col(w0);
	vec mu1_init_vec = as_Col(mu1_init);

	unordered_map<string, mat> result = fast_mr_ash_rss(bhat_vec, shat_vec, z_vec, R_mat, var_y, n, sigma2_e, s0_vec, w0_vec,
	                                                    mu1_init_vec, tol, max_iter, update_w0, update_sigma, compute_ELBO,
	                                                    standardize);

	writable::list ret;
	for (const auto& item : result) {
		cpp11::named_arg na(item.first.c_str());
		na = as_doubles_matrix(item.second);
		ret.push_back(na);
	}

	return ret;
}
