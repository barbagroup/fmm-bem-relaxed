#pragma once

#include <cmath>

#include "AnalyticalIntegral.hpp"
#include "Vec.hpp"
#include "Mat3.hpp"

namespace AnalyticalIntegral
{

typedef enum { G, dGdn } AnswerType;

/* original extra methods from Fata */
double gammai(double p1, double p2, double qet, double cr1, double cr2, double zn, double zd)
{
  /**********************************************************************************************************************************/
  /***************************************************  Determine gamma_i  **********************************************************/
  /**********************************************************************************************************************************/
  double pix2 = 2*M_PI;

  if ( (p1<0.0) && (p2>0.0) ) // 1-----x----------2 <==> if Projection of collocation point falls in side k
  {
    if ( qet<0.0 )
    {
      if ( (cr1>0.0) && (cr2<0.0) && (zn>0.0) && (zd<0.0) )
        return atan2(zn,zd) - pix2;
      else if ( (cr1<0.0) && (cr2>0.0) && (zn>0.0) && (zd<0.0) )
        return atan2(zn,zd) - pix2;
      else if ( (cr1<0.0) && (cr2<0.0) && (zn>0.0) )
        return atan2(zn,zd) - pix2;
      else
        return atan2(zn,zd);
    }
    else if ( qet>0.0 )
    {
      if ( (cr1>0.0) && (cr2<0.0) && (zn<0.0) && (zd<0.0) )
        return atan2(zn,zd) + pix2;
      else if ( (cr1<0.0) && (cr2>0.0) && (zn<0.0) && (zd<0.0) )
        return atan2(zn,zd) + pix2;
      else if ( (cr1<0.0) && (cr2<0.0) && (zn<0.0) )
        return atan2(zn,zd) + pix2;
      else
        return atan2(zn,zd);
    }
    else
      return 0.0;    // gammai = atan2(zn,zd) = 0.0
  }
  else // x 1----------------2 or 1----------------2 x <==> if Projection of collocation point does not fall in side k
    return atan2(zn,zd);
}

double pt_tria_tst(double w, double u, double v, const double *theta, int *theta_id)
{
  /**********************************************************************************************************************************/
  /*************************************************  Point in Triangle Test  *******************************************************/
  /***********************  Given a point with barycentric coordinates (u,v) on a plane containing a triangle  **********************/
  /************  Determine whether the point (u,v) is outside, or inside, or on an edge, or at a vertex of the triangle  ************/
  /**************************  Set the flag theta_id characterizing the position of the query point (u,v)  **************************/
  /**********************************************************************************************************************************/

  /* Triangle with vertices {y1,y2,y3}
   *
   *           y3
   *           *
   *          / \
   *         /   \
   *        /     \
   *       /       \
   *      /         \
   *     /           \
   *    *-------------*
   *  y1               y2
   *
   * Oriented edge L1 = [y1,y2] <== Edge1
   * Oriented edge L2 = [y2,y3] <== Edge2
   * Oriented edge L3 = [y3,y1] <== Edge3
   */

  /* Note: (u,v) == Barycentric coordinates of the query point on the plane containing the triangle
   * Note: w = 1-(u+v)
   * Note: The length of vector theta is 3, i.e. theta[3] */
  double pi = M_PI, pix2=2*M_PI;
  *theta_id = -1;

  if ( u<0.0 || v<0.0 || w<0.0 )
  {
    *theta_id = 7; return 0.0; // Query point is Outside the triangle
  }
  else if ( u==0.0 )
  {
    if ( v==0.0 )
    {
      *theta_id = 4; return theta[0]; // Query point is at Vertex1
    }
    else if ( v==1.0 )
    {
      *theta_id = 6; return theta[2]; // Query point is at Vertex3
    }
    else if ( v<0.0 || v>1.0)
    {
      *theta_id = 7; return 0.0; // Query point is Outside the triangle
    }
    else // v>0.0 && v<1.0
    {
      *theta_id = 3; return  pi; // Query point is on Edge3 of the triangle
    }
  }
  else if ( v==0.0 )
  {
    if ( u==1.0 )
    {
      *theta_id = 5; return theta[1]; // Query point is at Vertex2
    }
    else if ( u<0.0 || u>1.0 )
    {
      *theta_id = 7; return 0.0; // Query point is Outside the triangle
    }
    else // u>0.0 && u<1.0
    {
      *theta_id = 1; return pi; // Query point is on Edge1 of the triangle
    }
  }
  else if ( w==0.0 )
  {
    if ( u<0.0 || v<0.0 )
    {
      *theta_id = 7; return 0.0; // Query point is Outside the triangle
    }
    else // u>0.0 && v>0.0
    {
      *theta_id = 2; return pi; // Query point is on Edge2 of the triangle
    }
  }
  else // u>0.0 && v>0.0 && w>0.0
  {
    *theta_id = 0; return pix2; // Query point is Inside the triangle
  }
}

inline double asym_chi(double p1, double p2, double d)
{
  /**********************************************************************************************************************************/
  /************************************  return ASYMTOTIC VALUE of chi when q^2 + et^2 = d << 1  ************************************/
  /**********************************************************************************************************************************/
  double z1,z2;

  if ( fabs(p1)>0.0 && fabs(p2)>0.0 )
  {
    z1 = (d/p1)/p1; z2 = (d/p2)/p2;
    return log( fabs(p2/p1)*( (1.0-0.25*z1*(1.0-0.5*z1))/(1.0-0.25*z2*(1.0-0.5*z2)) ) );
  }
  else return 0.0;
}

void ortho_comp_basis(const double *t1, const double *t2, double *e1, double *e2, double *e3)
{
  /**********************************************************************************************************************************/
  /************************************  Given two linearly independent vectors t1 and t2 in 3D  ************************************/
  /** Compute the orthonormal triad {e1,e2,e3} such that e1=t1/norm(t1), norm(e2) = 1 and e2 is orthogonal to e1, and e3 = e1 x e2 **/
  /**********************************************************************************************************************************/

  // Note: The length of vectors t1, t2, e1, e2, e3 is 3
  int i;
  double al,snrm_te,nrm_te,nrm_tx;

  // Orthonormal local basis: e1[3],e2[3],e3[3]
  snrm_te = t1[0]*t1[0] + t1[1]*t1[1] + t1[2]*t1[2];
  nrm_te = sqrt(snrm_te);

  al = (t1[0]*t2[0] + t1[1]*t2[1] + t1[2]*t2[2])/snrm_te;

  for (i=0; i<3; i++) // Loop over all components
    e2[i] = t2[i] - al*t1[i];

  nrm_tx = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);

  for (i=0; i<3; i++) // Loop over all components
  {
    e1[i] = t1[i]/nrm_te; e2[i] = e2[i]/nrm_tx;
  }

  // Unit Normal vector to the plane spanned by {e1,e2}: e3 = e1 x e2
  e3[0] = e1[1]*e2[2] - e2[1]*e1[2];
  e3[1] = e1[2]*e2[0] - e2[2]*e1[0];
  e3[2] = e1[0]*e2[1] - e2[0]*e1[1];
}

/**
 * My additional routines for Stokes / Elasticity integrals
 */

// dyadic product
template <typename T>
Mat3<T> dyadic_prod(T *v1, T *v2)
{
  Mat3<T> m;

  m(0,0) = v1[0]*v2[0]; m(0,1) = v1[0]*v2[1]; m(0,2) = v1[0]*v2[2];
  m(1,0) = v1[1]*v2[0]; m(1,1) = v1[1]*v2[1]; m(1,2) = v1[1]*v2[2];
  m(2,0) = v1[2]*v2[0]; m(2,1) = v1[2]*v2[1]; m(2,2) = v1[2]*v2[2];

  return m;
}

template <typename T>
void add_influence(Mat3<T>& result, T factor, Mat3<T>& m)
{
  int i;
  // result <- result + factor*m (result, m matrices)
  for (i=0; i<9; i++) result.vals_[i] += factor*m.vals_[i];
}

template <typename T, typename Vector>
Vector matvec(Mat3<T>& A, Vector x)
{
  Vector m(0.);

  m[0] = A.vals_[0]*x[0] + A.vals_[1]*x[1] + A.vals_[2]*x[2];
  m[1] = A.vals_[3]*x[0] + A.vals_[4]*x[1] + A.vals_[5]*x[2];
  m[2] = A.vals_[6]*x[0] + A.vals_[7]*x[1] + A.vals_[8]*x[2];

  return m;
}

template <typename T>
void print_matrix(const char *name, Mat3<T>& A)
{
  printf("%s\n",name);
  printf("%.4lg %.4lg %.4lg\n",A.vals_[0],A.vals_[1],A.vals_[2]);
  printf("%.4lg %.4lg %.4lg\n",A.vals_[3],A.vals_[4],A.vals_[5]);
  printf("%.4lg %.4lg %.4lg\n",A.vals_[6],A.vals_[7],A.vals_[8]);
}

// Context to be given to specific integration routines
struct IntegrationContext
{
  // panels self-interacting?
  bool self_interaction;
  // Laplace
  double et, ThetGam, omega;

  // Stokeslet
  double e1[3], e2[3], e3[3]; // orthonormal basis
  double chi[3], alpha[3], q[3], rho[3];

  // Additional for Stresslet
  double p[3][3]; // have 0.0 for unused values

  IntegrationContext() {};
};

template <equation eq>
struct Integration
{
  template <typename Context>
  static void integrate(Context context, AnswerType type) {}
};

template <>
struct Integration<LAPLACE>
{
  template <typename Context>
  static double integrate(Context c, AnswerType type) {
    if (type == G) {
      return c.omega - c.et*c.ThetGam; // G
    } else {
      return -c.ThetGam; // dG/dn
    }
  }
};

template <>
struct Integration<STOKES>
{
  template <typename Context>
  // static Vec<3,double> integrate(Context c, AnswerType type) {
  static Mat3<double> integrate(Context c, AnswerType type) {

    typedef Mat3<double> Matrix;
    Matrix e1e1 = dyadic_prod(c.e1,c.e1);
    Matrix e1e2 = dyadic_prod(c.e1,c.e2);
    Matrix e1e3 = dyadic_prod(c.e1,c.e3);

    Matrix e2e1 = dyadic_prod(c.e2,c.e1);
    Matrix e2e2 = dyadic_prod(c.e2,c.e2);
    Matrix e2e3 = dyadic_prod(c.e2,c.e3);

    Matrix e3e1 = dyadic_prod(c.e3,c.e1);
    Matrix e3e2 = dyadic_prod(c.e3,c.e2);
    Matrix e3e3 = dyadic_prod(c.e3,c.e3);

    // Stokeslet
    if (type == G) {
      double rho_bar[3], I1, etI3, I3_xi, I3_xi_xi, I3_zeta, I3_zeta_zeta, I3_zeta_xi, I3_xi_zeta;

      // simple constants
      //I1 = c.omega - c.et*c.ThetGam;
      I1 = c.omega - c.et*c.ThetGam;
      etI3 = c.ThetGam;

      // single
      I3_xi = c.chi[0]*sin(c.alpha[0]) + c.chi[1]*sin(c.alpha[1]) + c.chi[2]*sin(c.alpha[2]);
      I3_zeta = -(c.chi[0]*cos(c.alpha[0]) + c.chi[1]*cos(c.alpha[1]) + c.chi[2]*cos(c.alpha[2]));

      // double
      rho_bar[0] = c.rho[0]-c.rho[1];
      rho_bar[1] = c.rho[1]-c.rho[2];
      rho_bar[2] = c.rho[2]-c.rho[0];

      I3_xi_xi = 0.;
      I3_zeta_zeta = 0.;
      I3_zeta_xi = 0.;
      for (int i=0; i<3; i++) {
        I3_xi_xi     += (c.q[i]*c.chi[i]*cos(c.alpha[i]) + rho_bar[i]*sin(c.alpha[i]))*cos(c.alpha[i]);
        I3_zeta_zeta += (c.q[i]*c.chi[i]*sin(c.alpha[i]) - rho_bar[i]*cos(c.alpha[i]))*sin(c.alpha[i]);
        I3_zeta_xi   += (c.q[i]*c.chi[i]*cos(c.alpha[i]) + rho_bar[i]*sin(c.alpha[i]))*sin(c.alpha[i]);
      }
      I3_xi_xi -= c.et*c.ThetGam;
      I3_zeta_zeta -= c.et*c.ThetGam;
      I3_xi_zeta = I3_zeta_xi;

      // Assemble final product
      Matrix IU;

      // xi terms
      add_influence(IU,I1+I3_xi_xi,e1e1);
      add_influence(IU,I3_xi_zeta,e1e2);
      add_influence(IU,c.et*I3_xi,e1e3);
      // zeta terms
      add_influence(IU,I3_xi_zeta,e2e1);
      add_influence(IU,I1+I3_zeta_zeta,e2e2);
      add_influence(IU,c.et*I3_zeta,e2e3);
      // eta terms
      add_influence(IU,c.et*I3_xi,e3e1);
      add_influence(IU,c.et*I3_zeta,e3e2);
      add_influence(IU,I1+c.et*etI3,e3e3);

      // Vec<3,double> r = matvec(IU, f);
      return IU;
    }
    // Stresslet
    else if (type == dGdn && c.self_interaction) return Mat3<double>(0.);
    else
    {
      double delta[3], d[3], L[3], et = c.et;

      delta[0] = c.p[0][0]/c.rho[0] - c.p[0][1]/c.rho[1];
      delta[1] = c.p[1][1]/c.rho[1] - c.p[1][2]/c.rho[2];
      delta[2] = c.p[2][2]/c.rho[2] - c.p[2][0]/c.rho[0];

      for (int i=0; i<3; i++) {
        d[i] = c.q[i]*c.q[i] + et*et;
        L[i] = 1./c.rho[i] - 1./c.rho[(i+1)%3];
      }

      // now the integrals
      double I5=0., I5_xi=0., I5_xi_xi=0., I5_zeta=0., I5_zeta_zeta=0., I5_xi_zeta=0.;

      // I5
      for (int i=0; i<3; i++) I5 += 1./(3*et*et)*c.q[i]/d[i]*delta[i];
      I5 += 1./et/et/et/3.*c.ThetGam;

      // single terms
      for (int i=0; i<3; i++) {
        I5_xi   +=  1./3*delta[i]/d[i]*sin(c.alpha[i]);
        I5_zeta += -1./3*delta[i]/d[i]*cos(c.alpha[i]);
      }
      // double and mixed terms
      for (int i=0; i<3; i++) {
        I5_xi_xi     += -1./3*(L[i]*cos(c.alpha[i]) + c.q[i]/d[i]*delta[i]*sin(c.alpha[i]))*sin(c.alpha[i]);
        I5_zeta_zeta +=  1./3*(L[i]*sin(c.alpha[i]) - c.q[i]/d[i]*delta[i]*cos(c.alpha[i]))*cos(c.alpha[i]);
        I5_xi_zeta   += -1./3*(L[i]*sin(c.alpha[i]) - c.q[i]/d[i]*delta[i]*cos(c.alpha[i]))*sin(c.alpha[i]);
      }
      I5_xi_xi     += 1./3/et*c.ThetGam;
      I5_zeta_zeta += 1./3/et*c.ThetGam;

      // assemble the integral
      Matrix IT;

      // adjustment for Stokes -- remove factor of 3
      // xi terms
      add_influence(IT,-et*I5_xi_xi,e1e1);
      add_influence(IT,-et*I5_xi_zeta,e1e2);
      add_influence(IT,-et*et*I5_xi,e1e3);
      // zeta terms
      add_influence(IT,-et*I5_xi_zeta,e2e1);
      add_influence(IT,-et*I5_zeta_zeta,e2e2);
      add_influence(IT,-et*et*I5_zeta,e2e3);
      // eta terms
      add_influence(IT,-et*et*I5_xi,e3e1);
      add_influence(IT,-et*et*I5_zeta,e3e2);
      add_influence(IT,-et*et*et*I5,e3e3);

      // printf("et: %g\n",et);
      // return the matrix
      return IT;
    }
  }
};

template <typename T1, typename T2>
auto multiply(T1 a, T2 b) -> decltype(a*b)
{
  return a*b;
}

template <typename T>
Vec<3,T> multiply(Mat3<T> m, Vec<3,T> x)
{
  return matvec(m,x);
}

template <equation e, typename ChargeType, typename VectorType=Vec<3,double>, typename Context=IntegrationContext>
auto FataAnalytical(VectorType& y1, VectorType& y2, VectorType& y3, ChargeType f,
//                     VectorType& x, bool self_interaction, AnswerType type) -> decltype(multiply(Integration<e>::integrate(Context(),type),f))
                    VectorType& x, bool self_interaction, AnswerType type) -> decltype(Integration<e>::integrate(Context(),type))
{
  (void) f; // for now
  Context context;

  /* code to set up the integration domain goes here
   *
   * e1, e2, e3 etc.
   *
   * Stored in Context
   */
  double pi = M_PI;
  // gm_eps, f_eps = Tolerance on the geometry, Tolerance on the argument of a function
  double gm_eps = 1e-10,f_eps = 1e-10;

  // h, g = potential influence, flux influence
  // double h, g;

  int k, i, th_id;
  double e1[3], e2[3], e3[3]; // Orthonormal companion basis associated with element elQ
  // Geometric parameters of element elQ: bQ = base, aQ = height, cQ = relative position of the 3rd node of element elQ
  double aQ, bQ, cQ;
  double r1[3], xi, zt, eth, et, ets, q[3], qs[3], rho[3], chi[3] = {0.,0.,0.}, gamma[3], gam;
  double vij[2][3], aQs, bmc, theta[3], alpha2, alpha3, Rh[3], Rhs[3], ch11, ch12, ch22, ch23, ch33, ch31;
  double kQ3; // kQ3 = Slope of Edge3 of element elQ
  double shc[3], sncf2, sncf3, cscf2, cscf3, qetx2[3];
  double p11, p12, p22, p23, p33, p31, x3, z3, rh1s, rh2s, rh3s, p11s, p12s, p33s, p11rh1, p31rh1, omega;
  double cr11, cr12, zn1, zd1, cr22, cr23, zn2, zd2, cr31, cr33, zn3, zd3, p12rh2, p22rh2, p23rh3, p33rh3, theta0, ThetGam;

  /* Note: Below, the Projection of the collocation point repnd[i][0:2] means the Projection of the collocation point repnd[i][0:2] on
   * Note: the plane of element elQ in the direction e3. e3 is the unit normal to element elQ */

  for (k=0; k<3; k++) // Loop over all components
  {
    vij[0][k] = y2[k] - y1[k]; // y2-y1 <== Side 1 of element elQ
    vij[1][k] = y3[k] - y1[k]; // y3-y1 <== Side 3 of element elQ
  }

  // Deternime the orthonormal companion basis associated with element
  ortho_comp_basis(vij[0],vij[1],e1,e2,e3);

  // Geometric parameters of element elQ: bQ = base, aQ = height, cQ = relative position of the 3rd node of element elQ
  bQ = vij[0][0]*e1[0] + vij[0][1]*e1[1] + vij[0][2]*e1[2]; // bQ>0.0 always
  aQ = vij[1][0]*e2[0] + vij[1][1]*e2[1] + vij[1][2]*e2[2]; // aQ>0.0 always
  cQ = vij[1][0]*e1[0] + vij[1][1]*e1[1] + vij[1][2]*e1[2];

  bmc = bQ-cQ; aQs = aQ*aQ;

  // theta = inclusion angle at vertex y^i
  theta[0] = acos(cQ/sqrt(cQ*cQ + aQs)); theta[1] = acos(bmc/sqrt(bmc*bmc + aQs)); theta[2] = pi - (theta[0]+theta[1]);

  // eqn (21), p.37
  alpha2 = pi-theta[1]; alpha3 = pi+theta[0]; // alpha1 = 0.0

  kQ3 = cQ/aQ; // kQ3 = Slope of Side 3 of element

  cscf2 = cos(alpha2); sncf2 = sin(alpha2); cscf3 = cos(alpha3); sncf3 = sin(alpha3);

  for (k=0; k<3; k++) // Loop over all components
    r1[k] = x[k] - y1[k]; // r1 = x[0:2]-y1

  /*
  printf("e1: %.12g, %.12g, %.12g\n",e1[0],e1[1],e1[2]);
  printf("e2: %.12g, %.12g, %.12g\n",e2[0],e2[1],e2[2]);
  printf("e3: %.12g, %.12g, %.12g\n",e3[0],e3[1],e3[2]);
  printf("r1: %.12g, %.12g, %.12g\n",r1[0],r1[1],r1[2]);
  */

  xi = r1[0]*e1[0] + r1[1]*e1[1] + r1[2]*e1[2];
  zt = r1[0]*e2[0] + r1[1]*e2[1] + r1[2]*e2[2];
  eth = r1[0]*e3[0] + r1[1]*e3[1] + r1[2]*e3[2]; //  + 1e-9;

  // printf("x1: %.12g, zt: %.12g, eth: %.12g\n",xi,zt,eth);
  p11 = -xi; p12 = bQ-xi; q[0] = -zt; et = -eth;

  if (fabs(et) < gm_eps) et = 0.0;

  x3 = cQ+p11; z3 = aQ+q[0];
  p22 = p12*cscf2 + q[0]*sncf2; p23 = x3*cscf2 + z3*sncf2; q[1] = q[0]*cscf2 - p12*sncf2;
  p31 = p11*cscf3 + q[0]*sncf3; p33 = x3*cscf3 + z3*sncf3; q[2] = q[0]*cscf3 - p11*sncf3;

  ets = et*et; p11s = p11*p11; p12s = p12*p12; p33s = p33*p33; qs[0] = q[0]*q[0]; qs[1] = q[1]*q[1]; qs[2] = q[2]*q[2];

  Rhs[0] = p11s + qs[0]; Rhs[1] = p12s + qs[0]; Rhs[2] = p33s + qs[2]; // Rhs(3) = x3*x3 + z3*z3
  Rh[0] = sqrt(Rhs[0]); Rh[1] = sqrt(Rhs[1]); Rh[2] = sqrt(Rhs[2]);
  rh1s = Rhs[0] + ets; rh2s = Rhs[1] + ets; rh3s = Rhs[2] + ets;
  rho[0] = sqrt(rh1s); rho[1] = sqrt(rh2s); rho[2] = sqrt(rh3s);

  shc[2] = -q[0]/aQ; shc[1] = (kQ3*q[0] - p11)/bQ;
  shc[0] = 1.0 - (shc[2]+shc[1]); // shc[0] = 1.0 + (p11 - kQ2*q[0])/bQ

  if ( fabs(shc[0]) < gm_eps ) { shc[0] = 0.0; q[1] = 0.0; }

  if ( fabs(shc[1]) < gm_eps ) { shc[1] = 0.0; q[2] = 0.0; }
  else if ( fabs(shc[1]-1.0) < gm_eps ) shc[1] = 1.0;

  if ( fabs(shc[2]) < gm_eps ) { shc[2] = 0.0; q[0] = 0.0; }
  else if ( fabs(shc[2]-1.0) < gm_eps ) shc[2] = 1.0;

  // Point in triangle test
  theta0 = pt_tria_tst(shc[0],shc[1],shc[2],theta,&th_id);

  if ( th_id == 4 ) // if the Projection of the collocation point x[0:2] is at vertex y1
  {
    Rh[0] = 0.0;
    if (et == 0.0) rho[0] = 0.0; // x[0:2] = y1
  }
  else if ( th_id == 5 ) // if the Projection of the collocation point x[0:2] is at vertex y2
  {
    Rh[1] = 0.0;
    if (et == 0.0) rho[1] = 0.0; // x[0:2] = y2
  }
  else if ( th_id == 6 ) // if the Projection of the collocation point x[0:2] is at vertex y3
  {
    Rh[2] = 0.0;
    if (et == 0.0) rho[2] = 0.0; // x[0:2] = y3
  }

  qetx2[0] = 2.0*q[0]*et; qetx2[1] = 2.0*q[1]*et; qetx2[2] = 2.0*q[2]*et;

  p11rh1 = p11*rho[0]; p12rh2 = p12*rho[1];
  p22rh2 = p22*rho[1]; p23rh3 = p23*rho[2];
  p33rh3 = p33*rho[2]; p31rh1 = p31*rho[0];

  if ( self_interaction ) // element self-interaction
  {
    gam = 0.0;
    omega = q[0]*log((p11+rho[0])/(p12+rho[1]))+q[1]*log((p22+rho[1])/(p23+rho[2]))+q[2]*log((p33+rho[2])/(p31+rho[0]));
  }
  else if ( rho[0] == 0.0 ) // if The collocation point x[0:2] is at vertex y1, i.e., x[0:2] = y1[0:2]
  {
    gam = 0.0;
    omega = q[1]*log( (p22+rho[1])/(p23+rho[2]) ); // chi[1] = log( (p22+rho[1])/(p23+rho[2]) )
  }
  else if ( rho[1] == 0.0 ) // if The collocation point x[0:2] is at vertex y2, i.e., x[0:2] = y2[0:2]
  {
    gam = 0.0;
    omega = q[2]*log( (p33+rho[2])/(p31+rho[0]) ); // chi[2] = log( (p33+rho[2])/(p31+rho[0]) )
  }
  else if ( rho[2] == 0.0 ) // if The collocation point x[0:2] is at vertex y3, i.e., x[0:2] = y3[0:2]
  {
    gam = 0.0;
    omega = q[0]*log( (p11+rho[0])/(p12+rho[1]) ); // chi[0] = log( (p11+rho[0])/(p12+rho[1]) )
  }
  else
  {
    if ( Rh[0] == 0.0 ) // if the Projection of the collocation point x[0:2] is at vertex y1
    {
      cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
      zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;
      gam = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2); // gam = gamma[1]; gamma[0] = -gamma[2] = +pi or -pi
    }
    else if ( Rh[1] == 0.0 ) // if the Projection of the collocation point x[0:2] is at vertex y2
    {
      cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
      zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;
      gam = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3); // gam = gamma[2]; gamma[0] = -gamma[1] = +pi or -pi
    }
    else if ( Rh[2] == 0.0 ) // if the Projection of the collocation point x[0:2] is at vertex y3
    {
      cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
      zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;
      gam = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1); // gam = gamma[0]; gamma[1] = -gamma[2] = +pi or -pi
    }
    else if ( q[0] == 0.0 ) // if the Projection of the collocation point x[0:2] is on the line parallel to edge [y1,y2]
    {
      cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
      zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

      cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
      zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

      gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
      gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
      gam = gamma[1]+gamma[2];
    }
    else if ( q[1] == 0.0 ) // if the Projection of the collocation point x[0:2] is on the line parallel to edge [y2,y3]
    {
      cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
      zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

      cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
      zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

      gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
      gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
      gam = gamma[0]+gamma[2];
    }
    else if ( q[2] == 0.0 ) // if the Projection of the collocation point x[0:2] is on the line parallel to edge [y3,y1]
    {
      cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
      zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

      cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
      zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

      gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
      gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
      gam = gamma[0]+gamma[1];
    }
    else
    {
      cr11 = qs[0]*rh1s - p11s*ets; cr12 = qs[0]*rh2s - p12s*ets;
      zn1 = qetx2[0]*(p12rh2*cr11 - p11rh1*cr12); zd1 = cr11*cr12 + qetx2[0]*qetx2[0]*p11rh1*p12rh2;

      cr22 = qs[1]*rh2s-p22*p22*ets; cr23 = qs[1]*rh3s-p23*p23*ets;
      zn2 = qetx2[1]*(p23rh3*cr22 - p22rh2*cr23); zd2 = cr22*cr23 + qetx2[1]*qetx2[1]*p22rh2*p23rh3;

      cr33 = qs[2]*rh3s - p33s*ets; cr31 = qs[2]*rh1s - p31*p31*ets;
      zn3 = qetx2[2]*(p31rh1*cr33 - p33rh3*cr31); zd3 = cr31*cr33 + qetx2[2]*qetx2[2]*p31rh1*p33rh3;

      gamma[0] = gammai(p11,p12,qetx2[0],cr11,cr12,zn1,zd1);
      gamma[1] = gammai(p22,p23,qetx2[1],cr22,cr23,zn2,zd2);
      gamma[2] = gammai(p33,p31,qetx2[2],cr33,cr31,zn3,zd3);
      gam = gamma[0]+gamma[1]+gamma[2];
    }

    ch11 = p11+rho[0]; ch12 = p12+rho[1];
    ch22 = p22+rho[1]; ch23 = p23+rho[2];
    ch33 = p33+rho[2]; ch31 = p31+rho[0];

    if ( (ch11 >= f_eps) && (ch12 >= f_eps) )
      chi[0] = log(ch11/ch12);
    else
      chi[0] = asym_chi(p11,p12,qs[0] + ets);

    if ( (ch22 >= f_eps) && (ch23 >= f_eps) )
      chi[1] = log(ch22/ch23);
    else
      chi[1] = asym_chi(p22,p23,qs[1] + ets);

    if ( (ch33 >= f_eps) && (ch31 >= f_eps) )
      chi[2] = log(ch33/ch31);
    else
      chi[2] = asym_chi(p33,p31,qs[2] + ets);
    omega = q[0]*chi[0] + q[1]*chi[1] + q[2]*chi[2];
  }
  double alpha[3] = {0., alpha2, alpha3};
  // ThetGam is eqn (44) in elasticity paper, eqn (27) in Potential
  if ( et > 0.0 )
    ThetGam = 0.5*gam + theta0;
  else
    ThetGam = 0.5*gam - theta0;

  // now we can setup the context
  context.self_interaction = self_interaction;
  context.et = et;
  context.ThetGam = ThetGam;
  context.omega = omega;

  for (i=0; i<3; i++) {
    context.e1[i] = e1[i];
    context.e2[i] = e2[i];
    context.e3[i] = e3[i];
    context.chi[i] = chi[i];
    context.alpha[i] = alpha[i];
    context.q[i] = q[i];
    context.rho[i] = rho[i];
  }
  // set p[][] values
  context.p[0][0] = p11;
  context.p[0][1] = p12;
  context.p[1][1] = p22;
  context.p[1][2] = p23;
  context.p[2][2] = p33;
  context.p[2][0] = p31;

  // now call the actual integration (LAPLACE / STOKES etc.)
  auto r = Integration<e>::integrate(context, type);

  // multiply with the charge to get the final result
  // return multiply(r,f);
  return r;
}

/*
int main(int argc, char **argv)
{
  // triangle vertices
  Vec<3,double> y1(0.,0.,0.), y2(1.,0.,0.), y3(1.,1.,0.);
  // charge
  Vec<3,double> f(0.,0.,1.);
  // target point
  Vec<3,double> x(2.,-1.5,6.);

  AnswerType type = G;

  if (argc > 1) {
    if (strcmp(argv[1],"-dgdn")==0) type = dGdn;
  }

  // auto r = AnalyticalIntegration<STOKES>(y1,y2,y3,f,x,false,type);
  auto r = AnalyticalIntegration<LAPLACE>(y1,y2,y3,1.,x,false,type);
  std::cout << "r: " << r << std::endl;
}*/

}; // end namespace Analytical
