#pragma once

#include "YukawaCartesian.hpp"
#include "SemiAnalytical.hpp"
#include "GaussQuadrature.hpp"
#include "BEMConfig.hpp"

class YukawaCartesianBEM : public YukawaCartesian
{
 public:
  //! # of quadrature points
  unsigned K;
  // forward declaration
  struct Panel;
  //! The dimension of the Kernel
  static constexpr unsigned dimension = YukawaCartesian::dimension;
  //! Point type
  typedef YukawaCartesian::point_type point_type;
  //! source type
  typedef Panel source_type;
  //! target type
  typedef Panel target_type;
  //! Charge type
  typedef YukawaCartesian::charge_type charge_type;
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
      point_type L0 = p2-p0;
      point_type L1 = p1-p0;

      point_type c = point_type(L0[1]*L1[2]-L0[2]*L1[1],
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
  typedef std::vector<YukawaCartesian::multipole_type> multipole_type;
  //! Local expansion type
  typedef std::vector<YukawaCartesian::local_type> local_type;

  //! Panel type (for BEM kernel(s)
  typedef Panel panel_type;

  //! default constructor - use delegating constructor
  YukawaCartesianBEM() : YukawaCartesianBEM(5,0.125,3) {};
  //! Constructor
  YukawaCartesianBEM(int p, double kappa, unsigned k=3) : YukawaCartesian(p,kappa), K(k) {
    BEMConfig::Init();
    auto Config = BEMConfig::Instance();
    Config->setK(K);
  };

  /** Initialize a multipole expansion with the size of a box at this level */
  void init_multipole(multipole_type& M,
                      const point_type& extents, unsigned level) const {
    M.resize(2);
    YukawaCartesian::init_multipole(M[0], extents, level);
    YukawaCartesian::init_multipole(M[1], extents, level);
  }
  /** Initialize a local expansion with the size of a box at this level */
  void init_local(local_type& L,
                  const point_type& extents, unsigned level) const {
    L.resize(2);
    YukawaCartesian::init_local(L[0], extents, level);
    YukawaCartesian::init_local(L[1], extents, level);
  }
  /** perform Gaussian integration over panel to evaluate \int G */
  double eval_G(const source_type source, const point_type target) const {
    auto dist = norm(target-source.center);

    // check whether I need Semi-Analytical or Gauss Quadrature
    if (sqrt(2*source.Area)/dist >= 0.5) {
      // perform semi-analytical integral
      namespace AI = AnalyticalIntegral;
      double G = 0., dGdn = 0.;
      auto& vertices = source.vertices;
      AI::SemiAnalytical<AI::YUKAWA>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10, Kappa);
      return G;
    }
    else
    {
      // loop over all my quadrature points, accumulating to get \int G
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      double r = 0.;
      for (unsigned i=0; i<K; i++) {
        auto dist = norm(target-source.quad_points[i]);
        auto inv_dist = 1./dist;
        if (dist < 1e-8) inv_dist = 0.;
        r += gauss_weight[i]*source.Area*exp(-Kappa*dist)*inv_dist;
      }
      return r;
    }
  };

  /** Perform gaussian integration from this panel to target */
  double eval_dGdn(const source_type source, const point_type target) const {
    double dist = norm(target-source.center);
    // check for self-interaction
    if (dist < 1e-8) return 2*M_PI;
    // check for SA / GQ
    if (sqrt(2*source.Area)/dist >= 0.5) {
      namespace AI = AnalyticalIntegral;
      // semi-analytical integral
      auto& vertices = source.vertices;
      double G = 0., dGdn = 0.;
      AI::SemiAnalytical<AI::YUKAWA>(G,dGdn,vertices[0],vertices[1],vertices[2],target,dist < 1e-10, Kappa);
      return -dGdn;
    } else {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights(); // GQ.weights(K);
      double res=0.;
      for (unsigned i=0; i<K; i++) {
        point_type dx = target - source.quad_points[i];
        // dG/dn = dx.n / |dx|^3
        auto r  = norm(dx);
        auto inv_r = 1./r;
        auto inv_r2 = inv_r * inv_r;

        if (r < 1e-8) { inv_r = 0.; inv_r2 = 0.; }
        // finally accumulate
        auto& normal = source.normal;
        auto pot = exp(-Kappa*r)* inv_r;
        dx *= pot * (Kappa * r + 1) * inv_r2;
        res += gauss_weight[i]*source.Area*(-dx[0]*normal[0]-dx[1]*normal[1]-dx[2]*normal[2]);
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
           const point_type& center, multipole_type& M, unsigned p) const {
    auto& gauss_weight = BEMConfig::Instance()->GaussWeights(); // GQ.weights(3);
    auto& I = YukawaCartesian::I;
    auto& J = YukawaCartesian::J;
    auto& K = YukawaCartesian::K;
    auto& fact = YukawaCartesian::fact;

    for (auto j=0u; j<source.quad_points.size(); j++) {
      // quad point specific constants
      auto& qp = source.quad_points[j];
      point_type dX = center - static_cast<point_type>(qp);
      auto mult_term = charge * gauss_weight[j] * source.Area;

      unsigned mterms = (p+1)*(p+2)*(p+3)/6;
      for (unsigned i=0; i < mterms; i++) {
        // Multipole term constants
        auto C = pow(dX[0],I[i]) * pow(dX[1],J[i]) * pow(dX[2],K[i]) / (fact[I[i]]*fact[J[i]]*fact[K[i]]);

        if (source.BC == Panel::POTENTIAL) {
          // influence of G needed
          M[0][i] += mult_term * C;
        }
        else // influence of dG/dn needed
        {
          auto& normal = source.normal;
          M[1][i] -= mult_term * C * I[i] / dX[0] * normal[0];
          M[1][i] -= mult_term * C * J[i] / dX[1] * normal[1];
          M[1][i] -= mult_term * C * K[i] / dX[2] * normal[2];
        }
      }
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
           const point_type& translation, unsigned p) const {
    YukawaCartesian::M2M(Msource[0],Mtarget[0],translation,p);
    YukawaCartesian::M2M(Msource[1],Mtarget[1],translation,p);
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
  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result, unsigned p) const {

    std::vector<real> a_aux(MTERMS,0), ax_aux(MTERMS,0), ay_aux(MTERMS,0), az_aux(MTERMS,0);

    point_type dX = static_cast<point_type>(target) - center;

    // potential, d{x,y,z} coefficients
    YukawaCartesian::getCoeff(a_aux,ax_aux,ay_aux,az_aux,dX);
    auto& I    = YukawaCartesian::I;
    auto& J    = YukawaCartesian::J;
    auto& K    = YukawaCartesian::K;
    auto& fact = YukawaCartesian::fact;

    double result_pot = 0.;
    double result_dn  = 0.;
    // loop over tarSize
    //
    unsigned mterms = (p+1)*(p+2)*(p+3)/6;
    for (unsigned j=0; j<mterms; j++) {
      double fact_term = fact[I[j]]*fact[J[j]]*fact[K[j]];
      result_pot +=  a_aux[j]*M[0][j]*fact_term;
      result_dn  +=  a_aux[j]*M[1][j]*fact_term;
    }

    if   (target.BC== Panel::POTENTIAL) result += result_pot;
    else                                result -= result_dn;
  }

  void M2L(const multipole_type& Msource,
                 local_type& Ltarget,
           const point_type& translation, unsigned p) const
  {
    YukawaCartesian::M2L(Msource[0],Ltarget[0],translation,p);
    YukawaCartesian::M2L(Msource[1],Ltarget[1],translation,p);
  }

  void L2L(const local_type& Lsource,
                 local_type& Ltarget,
           const point_type& translation, unsigned p) const
  {
    YukawaCartesian::L2L(Lsource[0],Ltarget[0],translation,p);
    YukawaCartesian::L2L(Lsource[1],Ltarget[1],translation,p);
  }

  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result, unsigned p) const
  {
    auto dx = static_cast<point_type>(target) - center;
    auto& I    = YukawaCartesian::I;
    auto& J    = YukawaCartesian::J;
    auto& K    = YukawaCartesian::K;
    auto& fact = YukawaCartesian::fact;

    unsigned mterms = (p+1)*(p+2)*(p+3)/6;
    for (unsigned i=0; i<mterms; i++) {
      double mult_term = pow(dx[0],I[i])*pow(dx[1],J[i])*pow(dx[2],K[i]) / (fact[I[i]]*fact[J[i]]*fact[K[i]]);
      double result_pot = L[0][i]*mult_term;
      double result_dn  = L[1][i]*mult_term;

      if (target.BC == Panel::POTENTIAL) result += result_pot;
      else                               result -= result_dn;
    }
  }
};

