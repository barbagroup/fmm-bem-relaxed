#pragma once

#include "StokesSpherical.hpp"
#include "FataAnalytical.hpp"
#include "GaussQuadrature.hpp"
#include "BEMConfig.hpp"

class StokesSphericalBEM : public StokesSpherical
{
 public:
  unsigned K;
  double Mu;
  // forward declaration
  struct Panel;
  static constexpr unsigned dimension = StokesSpherical::dimension;
  //! Point type
  typedef StokesSpherical::point_type point_type;
  //! source type
  typedef Panel source_type;
  //! target type
  typedef Panel target_type;
  //! Charge type
  typedef StokesSpherical::charge_type charge_type;
  //! kernel evaluation type
  typedef StokesSpherical::kernel_value_type kernel_value_type;
  //! result type
  typedef StokesSpherical::result_type result_type;

  //! structure for Boundary elements
  struct Panel
  {
    typedef enum { VELOCITY, TRACTION } BoundaryType;

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
    Panel() : center(0), normal(0), Area(0), BC(VELOCITY) {};
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
    Panel(point_type p0, point_type p1, point_type p2) : BC(VELOCITY) {
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
      if (this->BC == VELOCITY) { this->BC = TRACTION; }
      else                      { this->BC = VELOCITY; }
    };
  };

  //! multipole expansion type
  typedef std::vector<StokesSpherical::multipole_type> multipole_type;
  //! local expansion type
  typedef std::vector<StokesSpherical::local_type> local_type;

  //
  typedef Panel panel_type;

  //! default constructor
  StokesSphericalBEM() : StokesSphericalBEM(5,3,0.1) {};
  //! don't set viscosity
  StokesSphericalBEM(int p, unsigned k) : StokesSphericalBEM(p,k,0.1) {};
  //! main constructor
  StokesSphericalBEM(int p, unsigned k, double mu) : StokesSpherical(p), K(k), Mu(mu) {
    BEMConfig::Init();
    auto Config = BEMConfig::Instance();
    Config->setK(K);
  }

  void init_multipole(multipole_type& M, const point_type& extents, unsigned level) const {
    M.resize(2);
    StokesSpherical::init_multipole(M[0],extents, level);
    StokesSpherical::init_multipole(M[1],extents, level);
  }

  void init_local(local_type& L, const point_type& extents, unsigned level) const {
    L.resize(2);
    StokesSpherical::init_local(L[0],extents,level);
    StokesSpherical::init_local(L[1],extents,level);
  }

  result_type eval_traction_integral(const source_type& source, const target_type& target, const charge_type& g) const
  {
    auto dist = static_cast<point_type>(target) - source.center;
    auto d = norm(dist);

    // check for self-interaction
    bool self = d < 1e-8;

    // check for analytical integral / gauss quadrature
    if (sqrt(2*source.Area)/d>= 0.5)
    {
      namespace AI = AnalyticalIntegral;
      auto& vertices = source.vertices;

      return AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[2],vertices[1],g,target.center,self,AI::dGdn);
    }
    else
    {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      result_type res(0.);

      for (unsigned i=0; i<K; i++) {
        // eval contribution from this gauss point
        auto g2 = g*gauss_weight[i]*source.Area;
        point_type dist_quad = source.quad_points[i] - target.center;
        auto r2 = normSq(dist_quad);
        real R1 = r2;
        // real R2 = R1;
        real invR = 1. / R1;
        if (r2 < 1e-8) invR = 0;

        auto& normal = source.normal;
        auto dxdotn = dist_quad[0]*normal[0] + dist_quad[1]*normal[1] + dist_quad[2]*normal[2];

        //auto invR2 = 1./r2;
        //if (invR2 < 1e-8) invR2 = 0;
        //auto invR = std::sqrt(invR2);
        auto H = std::sqrt(invR) * invR; // 1 / R^3
        H *= dxdotn * invR;  // (dx . n) / r^5

        auto dx0 = dist_quad[0], dx1 = dist_quad[1], dx2 = dist_quad[2];

        res[0] += H * (dx0*dx0*g2[0] + dx0*dx1*g2[1] + dx0*dx2*g2[2]);
        res[1] += H * (dx0*dx1*g2[0] + dx1*dx1*g2[1] + dx1*dx2*g2[2]);
        res[2] += H * (dx0*dx2*g2[0] + dx1*dx2*g2[1] + dx2*dx2*g2[2]);
      }
      return res;
    }
  }

  result_type eval_velocity_integral(const source_type& source, const target_type& target, const charge_type& f) const
  {
    auto dist = static_cast<point_type>(target) - source.center;
    auto d = norm(dist);

    bool self = d < 1e-8;

    // check for analytical / gauss quadrature
    if (sqrt(2*source.Area)/d>= 0.5)
    {
      namespace AI = AnalyticalIntegral;
      auto& vertices = source.vertices;

      return AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[2],vertices[1],f,target.center,self,AI::G);
    }
    else
    {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      auto area = source.Area;
      result_type res(0.);

      for (unsigned i=0; i<K; i++) {
        // eval contribution from this gauss point
        auto dist_quad = source.quad_points[i] - target.center;
        auto r2 = normSq(dist_quad);
        real R1 = r2;
        real R2 = R1;
        real invR = 1. / R1;
        if (r2 < 1e-8) invR = 0;

        auto H = std::sqrt(invR) * invR; // 1 / R^3

        auto f2 = f*gauss_weight[i]*area;
        auto fdx = dist_quad[0]*f[0] + dist_quad[1]*f[1] + dist_quad[2]*f[2];

        res[0] +=  H * (f2[0] * R2 + fdx * dist_quad[0]);
        res[1] +=  H * (f2[1] * R2 + fdx * dist_quad[1]);
        res[2] +=  H * (f2[2] * R2 + fdx * dist_quad[2]);
      }
      return res;
    }
  }

  // Operators
  template <typename SourceIter, typename ChargeIter,
            typename TargetIter, typename ResultIter>
  void P2P(SourceIter s_first, SourceIter s_last, ChargeIter c_first,
           TargetIter t_first, TargetIter t_last, ResultIter r_first) const
  {
    auto ti = t_first;
    auto ri = r_first;

    for ( ; ti != t_last; ++ti, ++ri) {
      auto si = s_first;
      auto ci = c_first;

      result_type res(0.);
      for ( ; si != s_last; ++si, ++ci) {

        // get relevant G / dGdn for this panel
        if (ti->BC == Panel::VELOCITY) {
          res += eval_velocity_integral(*si, *ti, *ci);
        }
        else if (ti->BC == Panel::TRACTION) {
          res -= eval_traction_integral(*si, *ti, *ci);
        }
        else // should never get here
        {
          res += result_type(0.);
        }
      }
    }
  }

  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const
  {
    // setup stuff
    auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
    auto c0 = charge[0], c1 = charge[1], c2 = charge[2];

    complex Ynm[4*P*P], YnmTheta[4*P*P];

    // loop over quadrature points
    for (auto i=0u; i<source.quad_points.size(); i++) {
      auto qp = source.quad_points[i];
      point_type dist = static_cast<point_type>(qp) - center;
      real rho, alpha, beta;
      cart2sph(rho,alpha,beta,dist);
      evalMultipole(rho,alpha,-beta,Ynm,YnmTheta);

      // create expansions
      for (int n=0; n!=P; ++n) {
        for (int m=0; m<=n; ++m) {
          const int nm  = n * (n + 1) + m;
          const int nms = n * (n + 1) / 2 + m;

          // velocity specified
          if (source.BC == Panel::VELOCITY)
          {
            // ease of writing
            auto area = source.Area;
            auto gw = gauss_weight[i];

            auto fdotx = gw*area*c0*qp[0] + gw*area*c1*qp[1] + gw*area*c2*qp[2];
            M[0][0][nms] += gw * area * c0 * Ynm[nm];
            M[0][1][nms] += gw * area * c1 * Ynm[nm];
            M[0][2][nms] += gw * area * c2 * Ynm[nm];
            M[0][3][nms] += fdotx * Ynm[nm];
          }
          // traction specified
          else
          {
            auto& normal = source.normal;
            auto n0 = normal[0], n1 = normal[1], n2 = normal[2];
            auto& area = source.Area;
            auto gw = gauss_weight[i];

            complex brh = (double)n/rho*Ynm[nm]; // d(rho)
            complex bal = YnmTheta[nm];          // d(alpha)
            complex bbe = -complex(0,1.)*(double)m*Ynm[nm]; // d(beta)

            complex bxd = sin(alpha)*cos(beta)*brh + cos(alpha)*cos(beta)/rho*bal - sin(beta)/rho/sin(alpha)*bbe; // dx
            complex byd = sin(alpha)*sin(beta)*brh + cos(alpha)*sin(beta)/rho*bal + cos(beta)/rho/sin(alpha)*bbe; // dy
            complex bzd = cos(alpha)*brh - sin(alpha)/rho*bal; // dz

            // which order should these be in?
            auto rdotn = bxd*n0 + byd*n1 + bzd*n2;
            auto rdotg = bxd*gw*area*c0 + byd*gw*area*c1 + bzd*gw*area*c2;
            M[1][0][nms] += (rdotn * gw*area*c0 + rdotg * n0);
            M[1][1][nms] += (rdotn * gw*area*c1 + rdotg * n1);
            M[1][2][nms] += (rdotn * gw*area*c2 + rdotg * n2);

            auto xdotg = qp[0]*gw*area*c0 + qp[1]*gw*area*c1 + qp[2]*gw*area*c2;
            auto ndotx = n0*qp[0] + n1*qp[1] + n2*qp[2];
            M[1][3][nms] += rdotn * xdotg + rdotg * ndotx;
          }
        }
      }
    }
  }

  void M2M(const multipole_type& Msource,
                 multipole_type& Mtarget,
           const point_type& translation) const {
    StokesSpherical::M2M(Msource[0],Mtarget[0],translation);
    StokesSpherical::M2M(Msource[1],Mtarget[1],translation);
  }

  void M2P(const multipole_type& M, const point_type& center,
           const target_type& target, result_type& result) const
  {
    // temporary result
    result_type r(0.);

    if (target.BC == Panel::VELOCITY)
    {
      StokesSpherical::M2P(M[0],center,target,r);
      result += r;
    }
    else
    {
      StokesSpherical::M2P(M[1],center,target,r);
      result -= 1./6*r;
    }
  }

  void M2L(const multipole_type& Msource,
                 local_type& Ltarget,
           const point_type& translation) const {
    StokesSpherical::M2L(Msource[0],Ltarget[0],translation);
    StokesSpherical::M2L(Msource[1],Ltarget[1],translation);
  }

  void L2L(const local_type& Lsource,
                 local_type& Ltarget,
           const point_type& translation) const {
    StokesSpherical::L2L(Lsource[0],Ltarget[0],translation);
    StokesSpherical::L2L(Lsource[1],Ltarget[1],translation);
  }

  void L2P(const local_type& L, const point_type& center,
           const target_type& target, result_type& result) const
  {
    // temporary result
    result_type r(0.);

    if (target.BC == Panel::VELOCITY)
    {
      StokesSpherical::L2P(L[0],center,target,r);
      result += r;
    }
    else
    {
      StokesSpherical::L2P(L[1],center,target,r);
      result -= 1./6*r;
    }
  }
};
