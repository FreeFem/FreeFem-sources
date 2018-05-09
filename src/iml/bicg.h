/*
 * This file is part of FreeFem++.
 *
 * FreeFem++ is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * FreeFem++ is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
 */

// *****************************************************************
// Iterative template routine -- BiCG
//
// BiCG solves the unsymmetric linear system Ax = b
// using the Preconditioned BiConjugate Gradient method
//
// BiCG follows the algorithm described on p. 22 of the
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
// x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
// tolerance was reached
// tol  --  the residual after the final iteration
//
// *****************************************************************

template<class Matrix, class Vector, class Preconditioner, class Real>
int
BiCG (const Matrix &A, Vector &x, const Vector &b,
      const Preconditioner &M, int &max_iter, Real &tol) {
	Real resid;
	Vector rho_1(1), rho_2(1), alpha(1), beta(1);
	Vector z, ztilde, p, ptilde, q, qtilde;
	Real normb = norm(b);
	Vector r = b - A * x;
	Vector rtilde = r;

	if (normb == 0.0)
		normb = 1;

	if ((resid = norm(r) / normb) <= tol) {
		tol = resid;
		max_iter = 0;
		return 0;
	}

	for (int i = 1; i <= max_iter; i++) {
		z = M.solve(r);
		ztilde = M.trans_solve(rtilde);
		rho_1(0) = dot(z, rtilde);
		if (rho_1(0) == 0) {
			tol = norm(r) / normb;
			max_iter = i;
			return 2;
		}

		if (i == 1) {
			p = z;
			ptilde = ztilde;
		} else {
			beta(0) = rho_1(0) / rho_2(0);
			p = z + beta(0) * p;
			ptilde = ztilde + beta(0) * ptilde;
		}

		q = A * p;
		qtilde = A.trans_mult(ptilde);
		alpha(0) = rho_1(0) / dot(ptilde, q);
		x += alpha(0) * p;
		r -= alpha(0) * q;
		rtilde -= alpha(0) * qtilde;

		rho_2(0) = rho_1(0);
		if ((resid = norm(r) / normb) < tol) {
			tol = resid;
			max_iter = i;
			return 0;
		}
	}

	tol = resid;
	return 1;
}
