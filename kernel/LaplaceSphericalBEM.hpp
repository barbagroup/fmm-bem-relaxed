#pragma once

#include "LaplaceSpherical.hpp"
#include "SemiAnalytical.hpp"
#include "GaussQuadrature.hpp"

GaussQuadrature<double> GQ;

class LaplaceSphericalBEM : public LaplaceSpherical
{
 // private:
  //! # of quadrature points
 public:
  unsigned K;
  // forward declaration
  struct Panel;
  //! The dimension of the Kernel
  static constexpr unsigned dimension = LaplaceSpherical::dimension;
  //! Point type
  typedef LaplaceSpherical::point_type point_type;
  //! source type
  typedef Panel source_type;
  //! target type
  typedef Panel target_type;
  //! Charge type
  typedef LaplaceSpherical::charge_type charge_type;
  //! The return type of a kernel evaluation
  typedef double kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef double result_type;

  struct Panel
  {
    typedef enum { POTENTIAL, NORMAL_DERIV } BoundaryType;
    struct BoundaryCondition {
      BoundaryType type;
      double value;
      BoundaryCondition(BoundaryType t, double v) : type(t), value(v) {};
      BoundaryCondition() : type(POTENTIAL), value(0.) {};
    };

    point_type center;
    point_type normal;
    std::vector<point_type> vertices;
    //! Stored Quadrature Points
    std::vector<point_type> quad_points;
    double Area;
    // Boundary condition
    BoundaryCondition BC;
    Panel() : center(0), normal(0), Area(0), BC(POTENTIAL,0.) {};
    Panel(const Panel& p) {
      vertices = p.vertices;
      center = p.center;
      normal = p.normal;
      Area = p.Area;
      BC = p.BC;
      quad_points = p.quad_points;
    }
    Panel(point_type p0, point_type p1, point_type p2) : BC(POTENTIAL,0.) {
      vertices.resize(3);
      vertices[0] = p0;
      vertices[1] = p1;
      vertices[2] = p2;

      // get the center
      center = (p0+p1+p2)/3;
      // area & normal

      auto L0 = p2-p0;
      auto L1 = p1-p0;

      auto c = point_type(L0[1]*L1[2]-L0[2]*L1[1],
                          -(L0[0]*L1[2]-L0[2]*L1[0]),
                          L0[0]*L1[1]-L0[1]*L1[0]);
      Area = 0.5*norm(c);
      normal = c/2/Area;

      // generate the quadrature points
      // assume K = 3 for now
      unsigned K = 3;
      quad_points.resize(K);
      // loop over K points performing matvecs
      auto& points = GQ.points(3); // GQ_points_n3;
      for (unsigned i=0; i<K; i++) {
        double x = vertices[0][0]*points[i][0]+vertices[1][0]*points[i][1]+vertices[2][0]*points[i][2];
        double y = vertices[0][1]*points[i][0]+vertices[1][1]*points[i][1]+vertices[2][1]*points[i][2];
        double z = vertices[0][2]*points[i][0]+vertices[1][2]*points[i][1]+vertices[2][2]*points[i][2];
        
        quad_points[i] = point_type(x,y,z);
      }
    }
    // cast to point_type
    operator point_type() const { return center; };

    // copy operator
    Panel& operator=(const Panel& p) {
      center = p.center;
      normal = p.normal;
      vertices.resize(3);
      for (unsigned i=0; i<3; i++) vertices[i] = p.vertices[i];
      Area = p.Area;
      quad_points = p.quad_points;
      return *this;
    }
    // perform gaussian integration from this panel to target
    double eval_dGdn(const point_type target) const {
      double dist = norm(target-center);
      // check for SA / GQ
      if (sqrt(2*Area)/dist >= 0.5) {
        namespace AI = AnalyticalIntegral;
        // semi-analytical integral
        double G = 0., dGdn = 0.;
        AI::SemiAnalytical<AI::LAPLACE>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10);
        return -dGdn;
      } else {
        auto& gauss_weight = GQ.weights(3); // GQ_weight_n3;
        double res=0.;
        for (unsigned i=0; i<quad_points.size(); i++) {
          // auto dx = target-quad_points[i];
          auto dx = quad_points[i]-target;
          // dG/dn = dx.n / |dx|^3
          auto r2 = normSq(dx);
          auto r3 = r2*sqrt(r2);
          // finally accumulate
          res += gauss_weight[i]*Area*(dx[0]*normal[0]+dx[1]*normal[1]+dx[2]*normal[2])/r3;
          // res -= gauss_weight[i]*Area/r2;
        }
        return res;
      }
    };
    double eval_G(const point_type target) const {
      auto dist = norm(target-center);
      // check whether I need SA or GQ
      if (sqrt(2*Area)/dist >= 0.5) {
        namespace AI = AnalyticalIntegral;
        // semi-analytical integral
        double G = 0., dGdn = 0.;
        AI::SemiAnalytical<AI::LAPLACE>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10);
        return G;
      } else {
        // loop over all my quadrature points, accumulating to get \int G
        auto& gauss_weight = GQ.weights(3); // GQ_weight_n3;
        double r = 0.;
        for (unsigned i=0; i<quad_points.size(); i++)
          r += gauss_weight[i]*Area/norm(target-quad_points[i]);
        return r;
      }
    };
    /** flip the boundary condition flag for calculating RHS */
    void switch_BC(void) {
      if (this->BC.type == POTENTIAL) { this->BC.type = NORMAL_DERIV; }
      else                            { this->BC.type = POTENTIAL; }
    };

  };

  //! Multipole expansion type
  typedef LaplaceSpherical::multipole_type multipole_type;
  //! Local expansion type
  typedef LaplaceSpherical::local_type local_type;

  //! Panel type (for BEM kernel(s)
  typedef Panel panel_type;

  //! default constructor - use delegating constructor
  LaplaceSphericalBEM() : LaplaceSphericalBEM(5,4) {};
  //! Constructor
  LaplaceSphericalBEM(int p, unsigned k=4) : LaplaceSpherical(p), K(k) {};

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M, double box_size) const {
    LaplaceSpherical::init_multipole(M,box_size);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L, double box_size) const {
    LaplaceSpherical::init_local(L,box_size);
  }

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   * @result The Laplace potential and force 4-vector on t from s:
   * Potential: 1/|s-t|  Force: (s-t)/|s-t|^3
   */
  kernel_value_type operator()(const source_type& t,
                               const target_type& s) const {
    point_type dist = static_cast<point_type>(s) - static_cast<point_type>(t);         //   Vector from target to source
    real R2 = normSq(dist);          //   R^2
    real invR2 = 1.0 / R2;           //   1 / R^2
    if (R2 < 1e-8 && t.BC.type == Panel::NORMAL_DERIV) {
      // if same panel
      return kernel_value_type(2*M_PI);
    }
    real invR = std::sqrt(invR2);    //   potential
    dist *= invR2 * invR;            //   force
    // return kernel_value_type(invR, dist[0], dist[1], dist[2]);


    if (t.BC.type == Panel::POTENTIAL) {
      // if I know the potential for source panel, need to multiply by dG/dn for RHS, want G for solve
      return kernel_value_type(s.eval_G(static_cast<point_type>(t)));
    } else if (t.BC.type == Panel::NORMAL_DERIV) {
      // I know d(phi)/dn for this panel, need to multiply by G for RHS, want dG/dn for solve
      return kernel_value_type(s.eval_dGdn(static_cast<point_type>(t)));
    } else {
      // should never get here
      return 0.;
    }
  }

  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    auto& gauss_weight = GQ.weights(3); // GQ_weight_n3;
    for (auto i=0u; i<source.quad_points.size(); i++) {
      auto qp = source.quad_points[i];
      point_type dist = static_cast<point_type>(qp) - center;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) { 
        for( int m=0; m<=n; ++m ) { 
          const int nm  = n * n + n + m;
          const int nms = n * (n + 1) / 2 + m;
          if (source.BC.type == Panel::POTENTIAL) {
            // influence of G needed
            M[nms] += charge * gauss_weight[i] * source.Area * Ynm[nm];
          } else { 
            // otherwise influence of dGdn needed
            complex brh = (double)n/rho*Ynm[nm];
            complex bal = YnmTheta[nm];
            complex bbe = -complex(0,1.)*(double)m*Ynm[nm];

            complex bxd = sin(alpha)*cos(beta)*brh + cos(alpha)*cos(beta)/rho*bal - sin(beta)/rho/sin(alpha)*bbe;
            complex byd = sin(alpha)*sin(beta)*brh + cos(alpha)*sin(beta)/rho*bal + cos(beta)/rho/sin(alpha)*bbe;
            complex bzd = cos(alpha)*brh - sin(alpha)/rho*bal;

            auto& normal = source.normal;
            complex mult_term = charge * gauss_weight[i] * source.Area;
            M[nms] += mult_term * normal[0] * bxd;
            M[nms] += mult_term * normal[1] * byd;
            M[nms] += mult_term * normal[2] * bzd;
          }
        }
      }
      M.RMAX = std::max(M.RMAX, norm(dist));
      M.RCRIT = std::min(M.RCRIT, M.RMAX);
    }
  }

  /** Kernel M2M operator
   * M_t += Op(M_s) where M_t is the target and M_s is the source
   *
   * @param[in] source The multipole source at the child level
   * @param[in,out] target The multipole target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Msource includes the influence of all points within its box
   */
  void M2M(const multipole_type& Msource,
           multipole_type& Mtarget,
           const point_type& translation) const {
    LaplaceSpherical::M2M(Msource,Mtarget,translation);
  }

  /** Kernel M2L operation
   * L += Op(M)
   *
   * @param[in] Msource The multpole expansion source
   * @param[in,out] Ltarget The local expansion target
   * @param[in] translation The vector from source to target
   * @pre translation obeys the multipole-acceptance criteria
   * @pre Msource includes the influence of all points within its box
   */
  void M2L(const multipole_type& Msource,
                 local_type& Ltarget,
           const point_type& translation) const {
    LaplaceSpherical::M2L(Msource,Ltarget,translation);
  }

  /** Kernel M2P operation
   * r_i += Op(M)
   *
   * @param[in] M The multpole expansion
   * @param[in] center The center of the box with the multipole expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre M includes the influence of all points within its box
   */
  template <typename PointIter, typename ResultIter>
  void M2P(const multipole_type& M, const point_type& center,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    
    for( ; t_begin != t_end ; ++t_begin, ++r_begin ) { 
      point_type dist = static_cast<point_type>(*t_begin) - center;
      point_type spherical(0);
      point_type cartesian(0);
      double r_temp(0);
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalLocal(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) { 
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        r_temp += std::real(M[nms] * Ynm[nm]);
        spherical[0] -= std::real(M[nms] * Ynm[nm]) / r * (n+1);
        spherical[1] += std::real(M[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) { 
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          r_temp += 2 * std::real(M[nms] * Ynm[nm]);
          spherical[0] -= 2 * std::real(M[nms] *Ynm[nm]) / r * (n+1);
          spherical[1] += 2 * std::real(M[nms] *YnmTheta[nm]);
          spherical[2] += 2 * std::real(M[nms] *Ynm[nm] * CI) * m;
        }
      }
      if   ((*t_begin).BC.type == Panel::POTENTIAL) *r_begin += r_temp;
      else                                          *r_begin -= r_temp;
    }
  }

  /** Kernel L2L operator
   * L_t += Op(L_s) where L_t is the target and L_s is the source
   *
   * @param[in] source The local source at the parent level
   * @param[in,out] target The local target to accumulate into
   * @param[in] translation The vector from source to target
   * @pre Lsource includes the influence of all points outside its box
   */
  void L2L(const local_type& source,
           local_type& target,
           const point_type& translation) const {
    LaplaceSpherical::L2L(source,target,translation);
  }

  /** Kernel L2P operation
   * r_i += Op(L)
   *
   * @param[in] L The local expansion
   * @param[in] center The center of the box with the local expansion
   * @param[in] t_begin,t_end Iterator pair to the target points
   * @param[in] r_begin Iterator to the result accumulator
   * @pre L includes the influence of all points outside its box
   */
  template <typename PointIter, typename ResultIter>
  void L2P(const local_type& L, const point_type& center,
           PointIter t_begin, PointIter t_end,
           ResultIter r_begin) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    double r_temp = 0.;

    for (auto t = t_begin; t != t_end; ++t, ++r_begin) {
      r_temp = 0.;
      point_type dist = static_cast<point_type>(*t) - center;
      point_type spherical(0);
      point_type cartesian(0);
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipole(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) { 
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        r_temp += std::real(L[nms] * Ynm[nm]);
        spherical[0] += std::real(L[nms] * Ynm[nm]) / r * n;
        spherical[1] += std::real(L[nms] * YnmTheta[nm]);
        for( int m=1; m<=n; ++m ) { 
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          r_temp += 2 * std::real(L[nms] * Ynm[nm]);
          spherical[0] += 2 * std::real(L[nms] * Ynm[nm]) / r * n;
          spherical[1] += 2 * std::real(L[nms] * YnmTheta[nm]);
          spherical[2] += 2 * std::real(L[nms] * Ynm[nm] * CI) * m;
        }   
      }
      if   ((*t).BC.type == Panel::POTENTIAL) *r_begin += r_temp;
      else                                    *r_begin -= r_temp;
    }
  }
};


