//*****************************************************************
// Iterative template routine -- QMR
//
// QMR.h solves the unsymmetric linear system Ax = b using the
// Quasi-Minimal Residual method following the algorithm as described
// on p. 24 in the SIAM Templates book.
//
//   -------------------------------------------------------------
//   return value     indicates
//   ------------     ---------------------
//        0           convergence within max_iter iterations
//        1           no convergence after max_iter iterations
//                    breakdown in:
//        2             rho
//        3             beta
//        4             gamma
//        5             delta
//        6             ep
//        7             xi
//   -------------------------------------------------------------
//   
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax=b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************


#include <math.h>

template < class Matrix, class Vector, class Preconditioner1,
           class Preconditioner2, class Real >
int 
QMR(const Matrix &A, Vector &x, const Vector &b, const Preconditioner1 &M1, 
    const Preconditioner2 &M2, int &max_iter, Real &tol)
{
  Real resid;

  Vector rho(1), rho_1(1), xi(1), gamma(1), gamma_1(1), theta(1), theta_1(1);
  Vector eta(1), delta(1), ep(1), beta(1);

  Vector r, v_tld, y, w_tld, z;
  Vector v, w, y_tld, z_tld;
  Vector p, q, p_tld, d, s;

  Real normb = norm(b);

  r = b - A * x;

  if (normb == 0.0)
    normb = 1;

  if ((resid = norm(r) / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  v_tld = r;
  y = M1.solve(v_tld);
  rho(0) = norm(y);

  w_tld = r;
  z = M2.trans_solve(w_tld);
  xi(0) = norm(z);

  gamma(0) = 1.0;
  eta(0) = -1.0;
  theta(0) = 0.0;

  for (int i = 1; i <= max_iter; i++) {

    if (rho(0) == 0.0)
      return 2;                        // return on breakdown

    if (xi(0) == 0.0)
      return 7;                        // return on breakdown

    v = (1. / rho(0)) * v_tld;
    y = (1. / rho(0)) * y;

    w = (1. / xi(0)) * w_tld;
    z = (1. / xi(0)) * z;

    delta(0) = dot(z, y);
    if (delta(0) == 0.0)
      return 5;                        // return on breakdown

    y_tld = M2.solve(y);               // apply preconditioners
    z_tld = M1.trans_solve(z);

    if (i > 1) {
      p = y_tld - (xi(0) * delta(0) / ep(0)) * p;
      q = z_tld - (rho(0) * delta(0) / ep(0)) * q;
    } else {
      p = y_tld;
      q = z_tld;
    }

    p_tld = A * p;
    ep(0) = dot(q, p_tld);
    if (ep(0) == 0.0)
      return 6;                        // return on breakdown

    beta(0) = ep(0) / delta(0);
    if (beta(0) == 0.0)
      return 3;                        // return on breakdown

    v_tld = p_tld - beta(0) * v;
    y = M1.solve(v_tld);

    rho_1(0) = rho(0);
    rho(0) = norm(y);
    w_tld = A.trans_mult(q) - beta(0) * w;
    z = M2.trans_solve(w_tld);

    xi(0) = norm(z);

    gamma_1(0) = gamma(0);
    theta_1(0) = theta(0);

    theta(0) = rho(0) / (gamma_1(0) * beta(0));
    gamma(0) = 1.0 / sqrt(1.0 + theta(0) * theta(0));

    if (gamma(0) == 0.0)
      return 4;                        // return on breakdown

    eta(0) = -eta(0) * rho_1(0) * gamma(0) * gamma(0) / 
      (beta(0) * gamma_1(0) * gamma_1(0));

    if (i > 1) {
      d = eta(0) * p + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * d;
      s = eta(0) * p_tld + (theta_1(0) * theta_1(0) * gamma(0) * gamma(0)) * s;
    } else {
      d = eta(0) * p;
      s = eta(0) * p_tld;
    }

    x += d;                            // update approximation vector
    r -= s;                            // compute residual

    if ((resid = norm(r) / normb) <= tol) {
      tol = resid;
      max_iter = i;
      return 0;
    }
  }

  tol = resid;
  return 1;                            // no convergence
}
