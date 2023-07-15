#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// Binary search find the first element that isn't smaller than x in a sorted array
int f_loc(arma::vec arry, double x, int low, int high)
{
  if (arry[0] >= x)
  {
    return 0;
  }
  int mid = 0;
  while (low <= high)
  {
    mid = (low + high) / 2;
    if (arry[mid] >= x)
    {
      if (arry[mid - 1] < x)
      {
        return mid;
      }
      high = mid - 1;
    }
    else
    {
      low = mid + 1;
    }
  }
  return mid;
}

// Function to center a design matrix with user-defined weight
// [[Rcpp::export]]
NumericMatrix center(NumericMatrix M, NumericVector w, int N)
{
  int p = M.ncol();
  int n = M.nrow();
  NumericVector cmean(p);
  for (int i = 0; i < n; i++)
  {
    double wi = w[i];
    for (int j = 0; j < p; j++)
    {
      cmean[j] += (M(i, j) / wi);
    }
  }
  cmean = cmean / n / N;
  NumericMatrix out(n, p);
  for (int i = 0; i < n; i++)
  {
    out(i, _) = M(i, _) - cmean;
  }
  return out;
}

// Calculate the KM estimator of the survival function
// [[Rcpp::export]]
List km(NumericVector e, NumericVector d,
        NumericVector p)
{
  int t = e.length();
  NumericVector s(t, 1.0);
  double den = 0;
  int i = t - 1;
  while (i >= 0)
  {
    double ei = e[i];
    int j = i;
    double nut = 0;
    while (j >= 0)
    {
      if (e[j] < ei)
      {
        break;
      }
      else
      {
        double pj = p[j];
        den += 1 / pj;
        if (d[j] == 1)
        {
          nut += 1 / pj;
        }
      }
      j--;
    }
    double tmp = 1 - nut / den;
    s[j + 1] = tmp;
    i = j;
  }
  NumericVector out(t);
  double start = 1;
  for (int i = 0; i < t; i++)
  {
    start = start * s[i];
    out[i] = start;
  }
  List out1 = List::create(e, out);
  return out1;
}

// Gehan type non-smooth estimating function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_ns(arma::mat x, arma::vec y, arma::uvec d,
                   arma::vec beta, arma::vec p, int n)
{
  arma::vec e = y - x * beta;
  int r = y.n_elem;
  int m = x.n_cols;
  arma::mat out(r, m);
  for (int i = 0; i < r; i++)
  {
    if (d[i] == 0)
    {
      continue;
    }
    arma::rowvec nut(m);
    for (int j = 0; j < r; j++)
    {
      if (e[i] <= e[j])
      {
        nut += (x.row(i) - x.row(j)) / p[j];
      }
    }
    out.row(i) = nut / p[i] / n / r / r;
  }
  return out;
}

// Gehan type non-smooth estimating function expressed by well-defined martingale
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_mtg(arma::mat x, arma::vec d, arma::vec e,
                    arma::uvec ind)
{
  int n = x.n_rows;
  int m = x.n_cols;
  int r0 = ind.n_elem;
  arma::mat out(n, m);
  // Define the pilot subsample
  arma::mat x_pt = x.rows(ind);
  arma::vec d_pt = d.elem(ind);
  arma::vec e_pt = e.elem(ind);
  // Calculate Zi and Ni in the pilot subsample
  arma::mat xsum(r0, m);
  arma::rowvec csum(r0);
  arma::rowvec xtmp(m);
  double ctmp = 0;
  int k = r0 - 1;
  while (k >= 0)
  {
    ctmp += 1;
    xtmp += x_pt.row(k);
    xsum.row(k) = xtmp;
    csum[k] = ctmp;
    int l = k + 1;
    // Dealing with ties
    while (e_pt[l] == e_pt[k])
    {
      xsum.row(l) = xsum.row(k);
      csum[l] = csum[k];
      l++;
    }
    k--;
  }
  // Calculate the contribution of the ith observation
  for (int i = 0; i < n; i++)
  {
    // Find the location of e[i] in e_pt (O{log(r_0)})
    // By binary search
    int loc = f_loc(e_pt, e[i], 0, r0 - 1);
    arma::rowvec cont_i(m);
    arma::rowvec xi = x.row(i);
    if (d[i] == 1)
    {
      cont_i = csum[loc] * xi - xsum.row(loc);
    }
    // Dealing with ties
    int nloc = loc;
    while (e[i] == e_pt[nloc])
    {
      nloc++;
    }
    if (e[i] > e_pt[nloc])
    {
      nloc++;
    }
    nloc--;
    for (int j = 0; j <= nloc; j++)
    {
      if (d_pt[j] == 1)
      {
        cont_i += xsum.row(j) / csum[j] - xi;
      }
    }
    out.row(i) = cont_i / r0;
  }
  return out;
}

// Gehan type smooth estimating function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_smth(const arma::mat &x, const arma::vec &y,
                     const arma::vec &d, const arma::vec &p,
                     const arma::vec &b, const int &n)
{
  int r = y.n_elem;
  int m = b.n_elem;
  arma::vec e = y - x * b;
  arma::mat out(r, m);
  for (int i = 0; i < r - 1; i++)
  {
    arma::mat xj = x.rows(i + 1, r - 1);
    arma::mat xdif = x.row(i) - xj.each_row();
    arma::vec edif = e.subvec(i + 1, r - 1) - e[i];
    arma::vec dij = d.subvec(i + 1, r - 1);
    arma::vec rij = sqrt(sum(xdif % xdif, 1));
    arma::vec pij = p.subvec(i + 1, r - 1);
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 1);
    out.row(i) = sum(xdif.each_col() % ((phij % (d[i] + dij) - dij) / pij) / p[i], 0) / n / n / r / r;
  }
  return out;
}

// Gehan type smooth estimating function expressed by well-defined martingale
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_s_mtg(const arma::mat &x, const arma::vec &y,
                      const arma::vec &d, const arma::vec &p,
                      const arma::vec &bt, const arma::uvec ind,
                      const int &n)
{
  int r = x.n_rows;
  int m = x.n_cols;
  arma::vec e = y - x * bt;
  int r0 = ind.n_elem;
  // Define the pilot subsample
  arma::mat x_pt = x.rows(ind);
  arma::vec d_pt = d.elem(ind);
  arma::vec e_pt = e.elem(ind);
  arma::vec p_pt = p.elem(ind);
  // Calculate Zi and Ni in the pilot subsample
  arma::mat xpsum(r, m);
  arma::vec psum(r);
  for (int i = 0; i < r; i++)
  {
    arma::mat xdif = x_pt.each_row() - x.row(i);
    arma::vec edif = e_pt - e[i];
    arma::vec rij = sqrt(arma::sum(xdif % xdif, 1));
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 0);
    psum[i] = arma::accu(phij / p_pt) / r0;
    xpsum.row(i) = arma::sum(x_pt.each_col() % (phij / p_pt), 0) / r0;
  }
  arma::mat xpsum_pt = xpsum.rows(ind);
  arma::vec psum_pt = psum.elem(ind);
  arma::mat tmp = xpsum_pt.each_col() / psum_pt;
  // Calculate the estimating function
  arma::mat out(r, m);
  arma::mat out1 = x.each_col() % psum - xpsum;
  arma::mat out2(r, m);
  for (int j = 0; j < r0; j++)
  {
    if (d_pt[j] == 0 || std::isnan(arma::accu(tmp.row(j))))
    {
      continue;
    }
    arma::mat xdif = x.each_row() - x_pt.row(j);
    arma::vec edif = e - e_pt[j];
    arma::vec rij = sqrt(arma::sum(xdif % xdif, 1));
    arma::vec phij = arma::normcdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 0);
    out2 += (x.each_row() - tmp.row(j)).each_col() % phij / p_pt[j] / r0;
  }
  out = (out1.each_col() % (d / p)  - out2.each_col() / p) / r0 / n / n;
  return out;
}

// The jacobian matrix of the gehan type smooth estimating function
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat gehan_s_jaco(const arma::mat &x, const arma::vec &y,
                       const arma::vec &d, const arma::vec &p,
                       const arma::vec &b, const int &n)
{
  int r = y.n_elem;
  int m = b.n_elem;
  arma::vec e = y - x * b;
  arma::mat out(m, m);
  for (int i = 0; i < r - 1; i++)
  {
    arma::mat xj = x.rows(i + 1, r - 1);
    arma::mat xdif = x.row(i) - xj.each_row();
    arma::vec edif = e.subvec(i + 1, r - 1) - e[i];
    arma::vec dij = d.subvec(i + 1, r - 1) + d[i];
    arma::vec rij = sqrt(sum(xdif % xdif, 1));
    arma::vec pij = p.subvec(i + 1, r - 1);
    arma::vec phij = arma::normpdf(sqrt(r) * edif / rij);
    phij.replace(arma::datum::nan, 1);
    arma::mat xdif2 = xdif.each_col() % (phij % dij / rij / pij);
    xdif2.replace(arma::datum::nan, 0);
    out += xdif.t() * xdif2 / p[i] / n / n / r / r * sqrt(r);
  }
  return out;
}
