#pragma once

/** holder class for semi-analytical integral */

#include "AnalyticalIntegral.hpp"
#include "Vec.hpp"
#include "Mat3.hpp"

namespace AnalyticalIntegral
{

template <equation E, typename ResultType=double, typename T=double, typename KappaType=double>
void lineInt(ResultType& G, ResultType& dGdn, T z, T x, T v1, T v2, KappaType kappa=0.)
{
  auto theta1 = atan2(v1,x);
  auto theta2 = atan2(v2,x);
  auto dtheta = theta2 - theta1;
  auto thetam = (theta2 + theta1)/2;


  T absZ = fabs(z), signZ;
  if (absZ<1e-10) signZ = 0;
  else            signZ = z/absZ;

  // Loop over gauss points
  ResultType expKr, expKz = exp(-kappa*absZ);
  T thetak, Rtheta, R;

  // define gauss points -- hardcode for now
  const int n_gauss = 5;
  double xk[5] = { -9.06179846e-01, -5.38469310e-01, 1.78162900e-17, 9.06179846e-01,5.38469310e-01 };
  double wk[5] = { 0.23692689, 0.47862867, 0.56888889, 0.23692689, 0.47862867 };
  //double xk[7] = { -9.49107912e-01, -7.41531186e-01, -4.05845151e-01, -5.34227877e-17, 9.49107912e-01, 7.41531186e-01, 4.05845151e-01 };
  //double wk[7] = { 0.12948497, 0.27970539, 0.38183005, 0.41795918, 0.12948497, 0.27970539, 0.38183005 };
  for (int i=0; i<n_gauss; i++)
  {
    thetak = dtheta/2*xk[i] + thetam;
    Rtheta = x/cos(thetak);
    R      = sqrt(Rtheta*Rtheta + z*z);
    expKr  = exp(-kappa*R);
    if (E == HELMHOLTZ) {
      if (std::abs(kappa)>1e-10)
      {
        G += ResultType(0.);
        dGdn += ResultType(0.);
      }
      // devolve to Laplace
      else
      {
        G += wk[i]*(R-absZ) * dtheta/2;
        dGdn += wk[i]*(z/R - signZ) * dtheta/2;
      }
    } else if (E == YUKAWA) {
      if (kappa>1e-10)
      {
        G += -wk[i]*(expKr - expKz)/kappa * dtheta/2;
        dGdn +=  wk[i]*(z/R*expKr - expKz*signZ) * dtheta/2;
      }
      // devolve to Laplace
      else
      {
        G += wk[i]*(R-absZ) * dtheta/2;
        dGdn += wk[i]*(z/R - signZ) * dtheta/2;
      }
    } else {
      // E == LAPLACE
      G += wk[i]*(R-absZ) * dtheta/2;
      dGdn += wk[i]*(z/R - signZ) * dtheta/2;
    }
  }
}

template <typename T>
Vec<3,T> cross(const Vec<3,T>& u, const Vec<3,T>& v) {
  return Vec<3,T>(u[1]*v[2] - u[2]*v[1],
                  u[2]*v[0] - u[0]*v[2],
                  u[0]*v[1] - u[1]*v[0]);
}

template <equation E, typename ResultType=double, typename T=double>
void intSide(ResultType& G, ResultType& dGdn, Vec<3,T>& v1, Vec<3,T>& v2, T p, T Kappa)
{
  typedef Vec<3,T> vec3;
  typedef Mat3<T>  mat3;

  vec3 v21 = v2 - v1;

  auto L21 = norm(v21);
  vec3 v21u = v21/L21;

  vec3 unit(0,0,1);
  vec3 orthog = cross(unit, v21u);

  auto alpha = dot(v21,v1)/(L21*L21);

  vec3 rOrthog = -alpha*v21 + v1;

  mat3 rotateToVertLine;

  for(int i=0; i<3; i++)
  {
    //rotateToVertLine(0,i) = orthog[i];
    //rotateToVertLine(1,i) = v21u[i];
    //rotateToVertLine(2,i) = unit[i];
    rotateToVertLine(i,0) = orthog[i];
    rotateToVertLine(i,1) = v21u[i];
    rotateToVertLine(i,2) = unit[i];
  }

  vec3 v1new = rotateToVertLine.multiply(v1);

  if (v1new[0]<0)
  {
    v21u = -v21u;
    orthog = -orthog;
    rotateToVertLine = -rotateToVertLine;
    rotateToVertLine(2,2) = 1.;
    v1new = rotateToVertLine.multiply(v1);
  }

  vec3 v2new = rotateToVertLine.multiply(v2);
  vec3 rOrthognew = rotateToVertLine.multiply(rOrthog);
  auto x = v1new[0];

  if ((v1new[1]>0 && v2new[1]<0) || (v1new[1]<0 && v2new[1]>0))
  {
    T G1 = 0., dG1dn = 0.;
    T G2 = 0., dG2dn = 0.;

    lineInt<E,ResultType,T>(G1,dG1dn, p, x, 0, v1new[1], Kappa);
    lineInt<E,ResultType,T>(G2,dG2dn, p, x, v2new[1], 0, Kappa);

    G += G1 + G2;
    dGdn += dG1dn + dG2dn;
  }
  else
  {
    T G1 = 0., dG1dn = 0.;
    lineInt<E,ResultType,T>(G1,dG1dn, p, x, v1new[1], v2new[1], Kappa);

    G -= G1;
    dGdn -= dG1dn;
  }

}

template <equation E, typename ResultType=double, typename T=double>
void SemiAnalytical(ResultType& G, ResultType& dGdn, Vec<3,T> y0, Vec<3,T> y1, Vec<3,T> y2, Vec<3,T> x, bool same, T Kappa=0)
{
  typedef Vec<3,T> vec3;
  typedef Mat3<T> mat3;

  // Put first panel at origin
  vec3 x_panel = x-y0;
  vec3 y0_panel;
  vec3 y1_panel = y1-y0;
  vec3 y2_panel = y2-y0;
  vec3 X = y1_panel;

  // Find panel coordinate system X: 0->1
  vec3 Z = cross(y1_panel, y2_panel);
  auto Xnorm = norm(X);
  auto Znorm = norm(Z);
  X /= Xnorm;
  Z /= Znorm;

  vec3 Y = cross(Z,X);

  // Rotate the coordinate system to match panel plane
  mat3 rot_matrix;
  for (auto i=0u; i< 3; i++) {
    rot_matrix(0,i) = X[i];
    rot_matrix(1,i) = Y[i];
    rot_matrix(2,i) = Z[i];
  }

  vec3 panel0_plane = rot_matrix.multiply(y0_panel);
  vec3 panel1_plane = rot_matrix.multiply(y1_panel);
  vec3 panel2_plane = rot_matrix.multiply(y2_panel);
  vec3 x_plane      = rot_matrix.multiply(x_panel);

  // Shift origin so it matches collocation point
  vec3 panel0_final = panel0_plane - x_plane;
  vec3 panel1_final = panel1_plane - x_plane;
  vec3 panel2_final = panel2_plane - x_plane;

  // adjust final value
  panel0_final[2] = panel0_plane[2];
  panel1_final[2] = panel1_plane[2];
  panel2_final[2] = panel2_plane[2];

  // Loop over sides
  intSide<E,ResultType,T>(G, dGdn, panel0_final, panel1_final, x_plane[2], Kappa); // Side 0
  intSide<E,ResultType,T>(G, dGdn, panel1_final, panel2_final, x_plane[2], Kappa); // Side 1
  intSide<E,ResultType,T>(G, dGdn, panel2_final, panel0_final, x_plane[2], Kappa); // Side 2

  if (same)
  {
    if (E == YUKAWA)  dGdn = -2*M_PI;
    if (E == LAPLACE) dGdn = 2*M_PI;
  }
  // printf("G: %.4lg, dGdn: %.4lg\n",G,dGdn);
}

}; // end namespace AnalyticalIntegral
