//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the 
// Generalized Minimum Residual method
//
// GMRES follows the algorithm described on p. 20 of the 
// SIAM Templates book.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//  
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//  
//*****************************************************************
// modif F. hecht to by compatible with RNM.hpp 
// no dummy vector result 
//  M.solve(xx)  => M*(xx) 
//  dot(u,v) => (u,v)
//  norm(u) => sqrt( (u,u) ) 


template < class Matrix, class Vector >
void 
Update(Vector &x, int k, Matrix &h, Vector &s, Vector v[])
{
  Vector y(s);

  // Backsolve:  
  for (int i = k; i >= 0; i--) {
    y(i) /= h(i,i);
    for (int j = i - 1; j >= 0; j--)
      y(j) -= h(j,i) * y(i);
  }

  for (int j = 0; j <= k; j++)
    x += v[j] * y(j);
}


template < class Real >
Real 
abs(Real x)
{
  return (x > 0 ? x : -x);
}


template < class Operator, class Vector, class Preconditioner,
           class Matrix, class Real >
int 
GMRES(const Operator &A, Vector &x, const Vector &b,
      const Preconditioner &M, Matrix &H, int &m, int &max_iter,
      Real &tol)
{
  Real resid;
  int i, j = 1, k;
  Vector s(m+1), cs(m+1), sn(m+1), w,r,Ax;
  r=M*b;
  Real normb = sqrt((r,r));
  
  Ax=A * x;
  Ax=b-Ax;
  r = M*(Ax);
  Real beta = sqrt((r,r));
  
  if ( abs(normb) < 1.e-30)
    normb = 1;
  
  if ((resid = beta / normb) <= tol) {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  Vector *v = new Vector[m+1];

  while (j <= max_iter) {
    v[0] = r / beta;    
    s = 0.0;
    s(0) = beta;
    
    for (i = 0; i < m && j <= max_iter; i++, j++) {
      w = M*(Ax=A * v[i]);
      for (k = 0; k <= i; k++) {
        H(k, i) = (w, v[k]);
        w -= H(k, i) * v[k];
      }
      H(i+1, i) = sqrt((w,w));
      v[i+1] = w  / H(i+1, i) ; 

      for (k = 0; k < i; k++)
        ApplyPlaneRotation(H(k,i), H(k+1,i), cs(k), sn(k));
      
      GeneratePlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(H(i,i), H(i+1,i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i+1), cs(i), sn(i));
      if(verbosity>5 || (verbosity>2 && j%100==0) )
       cout << "GMRES: " << j << " " << abs(s(i+1)) << " " <<  normb << " " 
           <<  abs(s(i+1)) / normb << " < " << tol << endl;
    
      if ((resid = abs(s(i+1)) / normb) < tol) {
       cout << "GMRES converge: " << j << " " << abs(s(i+1)) << " " <<  normb << " " 
           <<  abs(s(i+1)) / normb << " < " << tol << endl;
      
        Update(x, i, H, s, v);
        tol = resid;
        max_iter = j;
        delete [] v;
        return 0;
      }
    }
    Update(x, m - 1, H, s, v);
    Ax = A*x;
    Ax = b-Ax;
    
    r = M*(Ax);
    beta = sqrt((r,r));
    if(verbosity>4)
      cout << "GMRES: restart" << j << " " << beta << " " <<  normb << " " 
           <<  beta / normb << " < " << tol << endl;
    if ((resid = beta / normb) < tol) {
      tol = resid;
      max_iter = j;
      delete [] v;
      return 0;
    }
  }
  
  tol = resid;
  delete [] v;
  return 1;
}


#include <math.h> 


template<class Real> 
void GeneratePlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  if (dy == 0.0) {
    cs = 1.0;
    sn = 0.0;
  } else if (abs(dy) > abs(dx)) {
    Real temp = dx / dy;
    sn = 1.0 / sqrt( 1.0 + temp*temp );
    cs = temp * sn;
  } else {
    Real temp = dy / dx;
    cs = 1.0 / sqrt( 1.0 + temp*temp );
    sn = temp * cs;
  }
}


template<class Real> 
void ApplyPlaneRotation(Real &dx, Real &dy, Real &cs, Real &sn)
{
  Real temp  =  cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}

