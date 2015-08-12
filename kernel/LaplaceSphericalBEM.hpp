#pragma once

#include "LaplaceSpherical.hpp"
#include "SemiAnalytical.hpp"
#include "GaussQuadrature.hpp"
#include "BEMConfig.hpp"

// #define USE_ANALYTICAL

#if defined(USE_ANALYTICAL)
#include "FataAnalytical.hpp"
#endif

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

  //! structure for Boundary elements
  struct Panel
  {
    typedef enum { POTENTIAL, NORMAL_DERIV } BoundaryType;

    //! center of the panel
    point_type center;
    //! panel normal
    point_type normal;
    //! vertices of the panel
    std::vector<point_type> vertices;
    //! Stored Quadrature Points
    std::vector<point_type> quad_points;
    double Area;
    //! Boundary condition
    BoundaryType BC;
    Panel() : BC(POTENTIAL) {};
    //! copy constructor
    Panel(const Panel& p) {
      vertices = p.vertices;
      center = p.center;
      normal = p.normal;
      Area = p.Area;
      BC = p.BC;
      quad_points = p.quad_points;
    }
    //! main constructor
    Panel(point_type p0, point_type p1, point_type p2) : BC(POTENTIAL) {
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
      auto Config = BEMConfig::Instance();
      unsigned k = Config->getK();
      auto& points = Config->GaussPoints();
      quad_points.resize(k);
      // loop over K points performing matvecs
      for (unsigned i=0; i<k; i++) {
        double x = vertices[0][0]*points[i][0]+vertices[1][0]*points[i][1]+vertices[2][0]*points[i][2];
        double y = vertices[0][1]*points[i][0]+vertices[1][1]*points[i][1]+vertices[2][1]*points[i][2];
        double z = vertices[0][2]*points[i][0]+vertices[1][2]*points[i][1]+vertices[2][2]*points[i][2];

        quad_points[i] = point_type(x,y,z);
      }
    }
    /** cast to point_type */
    operator point_type() const { return center; };

    /** copy operator */
    Panel& operator=(const Panel& p) {
      center = p.center;
      normal = p.normal;
      vertices.resize(3);
      for (unsigned i=0; i<3; i++) vertices[i] = p.vertices[i];
      Area = p.Area;
      quad_points = p.quad_points;
      BC = p.BC;
      return *this;
    }

    /** flip the boundary condition flag for calculating RHS */
    void switch_BC(void) {
      if (this->BC == POTENTIAL) { this->BC = NORMAL_DERIV; }
      else                       { this->BC = POTENTIAL; }
    };
  };

  //! Multipole expansion type
  typedef std::vector<LaplaceSpherical::multipole_type> multipole_type;
  //! Local expansion type
  typedef std::vector<LaplaceSpherical::local_type> local_type;

  //! Panel type (for BEM kernel(s)
  typedef Panel panel_type;

  //! default constructor - use delegating constructor
  LaplaceSphericalBEM() : LaplaceSphericalBEM(5,3) {};
  //! Constructor
  LaplaceSphericalBEM(int p, unsigned k=3) : LaplaceSpherical(p), K(k) {
    BEMConfig::Init();
    auto Config = BEMConfig::Instance();
    Config->setK(K);
  };

  void set_p(int p)
  {
    LaplaceSpherical::set_p(p);
  }

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned level) const {
    (void) level;
    M.resize(2);
    LaplaceSpherical::init_multipole(M[0], extents, level);
    LaplaceSpherical::init_multipole(M[1], extents, level);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type& extents, unsigned level) const {
    (void) level;
    L.resize(2);
    LaplaceSpherical::init_local(L[0], extents, level);
    LaplaceSpherical::init_local(L[1], extents, level);
  }
  /** perform Gaussian integration over panel to evaluate \int G */
  double eval_G(const source_type& source, const point_type& target) const {
    auto dist = norm(target-source.center);

    // check whether I need Semi-Analytical or Gauss Quadrature
    if (sqrt(2*source.Area)/dist >= 0.5) {
      // perform semi-analytical integral
      namespace AI = AnalyticalIntegral;
      auto& vertices = source.vertices;

      // self-interaction case
      if (true) { // dist < 1e-10) {
#if defined(USE_ANALYTICAL)
        return AI::FataAnalytical<AI::LAPLACE>(vertices[0],vertices[2],vertices[1],1.,target,dist<1e-10,AI::G);
#else
        double G = 0., dGdn = 0.;
        AI::SemiAnalytical<AI::LAPLACE>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10);
        return G;
#endif
      } else {
        // K_fine case
        auto& gauss_weights = BEMConfig::Instance()->GaussWeights(17);
        auto& gauss_points  = BEMConfig::Instance()->GaussPoints(17);

        // loop over gauss points
        double r = 0.;
        for (unsigned i=0; i < gauss_points.size(); i++) {
          auto& gp = gauss_points[i];
          Vec<3,double> point(vertices[0][0]*gp[0] + vertices[1][0]*gp[1] + vertices[2][0]*gp[2],
                              vertices[0][1]*gp[0] + vertices[1][1]*gp[1] + vertices[2][1]*gp[2],
                              vertices[0][2]*gp[0] + vertices[1][2]*gp[1] + vertices[2][2]*gp[2]);

          // accumulate influence
          r += gauss_weights[i]*source.Area/norm(target - point);
        }
        return r;
      }
    }
    else
    {
      // loop over all my quadrature points, accumulating to get \int G
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      double r = 0.;
      for (unsigned i=0; i<K; i++)
        r += gauss_weight[i]*source.Area/norm(target-source.quad_points[i]);
      return r;
    }
  };

  /** Perform gaussian integration from this panel to target */
  double eval_dGdn(const source_type& source, const point_type& target) const {
    double dist = norm(target-source.center);
    // check for self-interaction
    if (dist < 1e-8) {
      return 2*M_PI;
    }
    // check for SA / GQ
    if (sqrt(2*source.Area)/dist >= 0.5) {
#if 0
      namespace AI = AnalyticalIntegral;
      // semi-analytical integral
      auto& vertices = source.vertices;
#if defined(USE_ANALYTICAL)
      return -AI::FataAnalytical<AI::LAPLACE>(vertices[0],vertices[2],vertices[1],1.,target,dist<1e-10,AI::dGdn);
#else
      double G = 0., dGdn = 0.;
      AI::SemiAnalytical<AI::LAPLACE>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10);
      return -dGdn;
#endif
#endif
      // K_fine case
      auto& vertices = source.vertices;
      auto& normal = source.normal;
      auto& gauss_weights = BEMConfig::Instance()->GaussWeights(17);
      auto& gauss_points  = BEMConfig::Instance()->GaussPoints(17);

      // loop over gauss points
      double r = 0.;
      for (unsigned i=0; i < gauss_points.size(); i++) {
        auto& gp = gauss_points[i];
        Vec<3,double> point(vertices[0][0]*gp[0] + vertices[1][0]*gp[1] + vertices[2][0]*gp[2],
                            vertices[0][1]*gp[0] + vertices[1][1]*gp[1] + vertices[2][1]*gp[2],
                            vertices[0][2]*gp[0] + vertices[1][2]*gp[1] + vertices[2][2]*gp[2]);

        point_type dx = point - target;
        double r2 = normSq(dx);
        double r3 = r2*std::sqrt(r2);

        // accumulate influence
        r += gauss_weights[i]*source.Area*(dx[0]*normal[0]+dx[1]*normal[1]+dx[2]*normal[2])/r3;
      }
      return r;
    } else {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights(); // GQ.weights(K);
      double res=0.;
      for (unsigned i=0; i<K; i++) {
        auto dx = source.quad_points[i]-target;
        // dG/dn = dx.n / |dx|^3
        auto r2 = normSq(dx);
        auto r3 = r2*sqrt(r2);
        // finally accumulate
        auto& normal = source.normal;
        res += gauss_weight[i]*source.Area*(dx[0]*normal[0]+dx[1]*normal[1]+dx[2]*normal[2])/r3;
      }
      return res;
    }
  };

  /** Kernel evaluation
   * K(t,s)
   *
   * @param[in] t,s The target and source points to evaluate the kernel
   * @result The Laplace potential and force 4-vector on t from s:
   * Potential: 1/|s-t|  Force: (s-t)/|s-t|^3
   */
  kernel_value_type operator()(const source_type& t,
                               const target_type& s) const {
    point_type dist = static_cast<point_type>(s) - static_cast<point_type>(t);
    // distances
    real R2 = normSq(dist);
    real invR2 = 1.0 / R2;
    real invR = std::sqrt(invR2);
    dist *= invR2 * invR;

    if (t.BC == Panel::POTENTIAL)
    {
      // if I know the potential for source panel, need to multiply by dG/dn for RHS, want G for solve
      return kernel_value_type(eval_G(s, static_cast<point_type>(t)));
    }
    else if (t.BC== Panel::NORMAL_DERIV)
    {
      // I know d(phi)/dn for this panel, need to multiply by G for RHS, want dG/dn for solve
      return kernel_value_type(eval_dGdn(s, static_cast<point_type>(t)));
    }
    else
    {
      // should never get here
      return 0.;
    }
  }

  /** Kernel P2M operation
   * M += Op(s) * c where M is the multipole and s is the source
   *
   * @param[in] source The point source
   * @param[in] charge The source's corresponding charge
   * @param[in] center The center of the box containing the multipole expansion
   * @param[in,out] M The multipole expansion to accumulate into
   */
  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const {
    complex Ynm[4*P*P], YnmTheta[4*P*P];

    auto& gauss_weight = BEMConfig::Instance()->GaussWeights(); // GQ.weights(3);

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
          if (source.BC == Panel::POTENTIAL)
          {
            // influence of G needed
            M[0][nms] += charge * gauss_weight[i] * source.Area * Ynm[nm];
          }
          else
          {
            // otherwise influence of dGdn needed
            complex brh = (double)n/rho*Ynm[nm]; // d(rho)
            complex bal = YnmTheta[nm];          // d(alpha)
            complex bbe = -complex(0,1.)*(double)m*Ynm[nm]; // d(beta)

            complex bxd = sin(alpha)*cos(beta)*brh + cos(alpha)*cos(beta)/rho*bal - sin(beta)/rho/sin(alpha)*bbe; // dx
            complex byd = sin(alpha)*sin(beta)*brh + cos(alpha)*sin(beta)/rho*bal + cos(beta)/rho/sin(alpha)*bbe; // dy
            complex bzd = cos(alpha)*brh - sin(alpha)/rho*bal; // dz

            auto& normal = source.normal;
            complex mult_term = charge * gauss_weight[i] * source.Area;
            M[1][nms] += mult_term * normal[0] * bxd;
            M[1][nms] += mult_term * normal[1] * byd;
            M[1][nms] += mult_term * normal[2] * bzd;
          }
        }
      }
      M[0].RMAX  = std::max(M[0].RMAX, norm(dist));
      M[0].RCRIT = std::min(M[0].RCRIT, M[0].RMAX);
      M[1].RMAX  = std::max(M[1].RMAX, norm(dist));
      M[1].RCRIT = std::min(M[1].RCRIT, M[1].RMAX);
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
    LaplaceSpherical::M2M(Msource[0],Mtarget[0],translation);
    LaplaceSpherical::M2M(Msource[1],Mtarget[1],translation);
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
    LaplaceSpherical::M2L(Msource[0],Ltarget[0],translation);
    LaplaceSpherical::M2L(Msource[1],Ltarget[1],translation);
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
      point_type cartesian;
      double r0_temp(0), r1_temp(0);
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalLocal(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        r0_temp += std::real(M[0][nms] * Ynm[nm]);
        r1_temp += std::real(M[1][nms] * Ynm[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          r0_temp += 2 * std::real(M[0][nms] * Ynm[nm]);
          r1_temp += 2 * std::real(M[1][nms] * Ynm[nm]);
        }
      }
      if   ((*t_begin).BC== Panel::POTENTIAL) *r_begin += r0_temp;
      else                                    *r_begin -= r1_temp;
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
    LaplaceSpherical::L2L(source[0],target[0],translation);
    LaplaceSpherical::L2L(source[1],target[1],translation);
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

    for (auto t = t_begin; t != t_end; ++t, ++r_begin) {
      double r0_temp = 0., r1_temp = 0.;
      point_type dist = static_cast<point_type>(*t) - center;
      point_type cartesian;
      real r, theta, phi;
      cart2sph(r,theta,phi,dist);
      evalMultipole(r,theta,phi,Ynm,YnmTheta);
      for( int n=0; n!=P; ++n ) {
        int nm  = n * n + n;
        int nms = n * (n + 1) / 2;
        r0_temp += std::real(L[0][nms] * Ynm[nm]);
        r1_temp += std::real(L[1][nms] * Ynm[nm]);
        for( int m=1; m<=n; ++m ) {
          nm  = n * n + n + m;
          nms = n * (n + 1) / 2 + m;
          r0_temp += 2 * std::real(L[0][nms] * Ynm[nm]);
          r1_temp += 2 * std::real(L[1][nms] * Ynm[nm]);
        }
      }
      if   ((*t).BC== Panel::POTENTIAL) *r_begin += r0_temp;
      else                              *r_begin -= r1_temp;
    }
  }
};

