#pragma once

#include "StokesSpherical.hpp"
#include "FataAnalytical.hpp"
#include "GaussQuadrature.hpp"
#include "BEMConfig.hpp"
#include "Mat3.hpp"

class StokesSphericalBEM : public StokesSpherical
{
 public:
  mutable complex stokeslet_str[4], stresslet_str[4];
  unsigned K;
  unsigned K_fine;
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
  typedef Mat3<real> kernel_value_type;
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

      // check vertices not duplicated
      assert(p0 != p1 && p0 != p2 && p1 != p2);
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

      point_type test = p2;
      assert(Area >= 0.);
      // check for sphere case: norm(center + normal) > R
      // assert(norm(normal+center) >= 1);

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
  StokesSphericalBEM() : StokesSphericalBEM(5,3,1e-3) {};
  //! don't set viscosity
  StokesSphericalBEM(int p, unsigned k) : StokesSphericalBEM(p,k,1e-3) {};
  //! main constructor
  StokesSphericalBEM(int p, unsigned k, double mu) : StokesSpherical(p), K(k), K_fine(25), Mu(mu) {
    BEMConfig::Init();
    auto Config = BEMConfig::Instance();
    Config->setK(K);
    for (unsigned i=0; i<4; i++) {
      stresslet_str[i] = stokeslet_str[i] = 0.;
    }
  }

  void set_Kfine(unsigned k) {
    K_fine = k;
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

  kernel_value_type eval_traction_integral(const source_type& source, const target_type& target) const
  {
    auto dist = static_cast<point_type>(target) - source.center;
    auto d = norm(dist);
    result_type g(0.);

    // check for self-interaction
    bool self = fabs(d) < 1e-8;
    if (fabs(d) < 1e-8) {
      kernel_value_type r(0.);
      r(0,0) = 2*M_PI;
      r(1,1) = 2*M_PI;
      r(2,2) = 2*M_PI;
      return r;
    }

    // check for analytical integral / gauss quadrature
    if (sqrt(2*source.Area)/d>= 0.5)
    {
      namespace AI = AnalyticalIntegral;
      auto& vertices = source.vertices;

      if (self) { // never gets to here -- singular = 2/pi I
        auto M(AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[2],vertices[1],g,target.center,self,AI::dGdn));
        // auto M(AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[1],vertices[2],g,target.center,self,AI::dGdn));
        return M.multiply(-3);
      } else {
        // K_fine case
        auto& gauss_weights = BEMConfig::Instance()->GaussWeights(K_fine);
        auto& gauss_points  = BEMConfig::Instance()->GaussPoints(K_fine);
        real area = source.Area;
        kernel_value_type res(0.), temp(0.);

        // loop over gauss points
        for (unsigned i=0; i<gauss_points.size(); i++) {

          // generate this gauss point
          auto& gp = gauss_points[i];

          point_type point(vertices[0][0]*gp[0] + vertices[1][0]*gp[1] + vertices[2][0]*gp[2],
                           vertices[0][1]*gp[0] + vertices[1][1]*gp[1] + vertices[2][1]*gp[2],
                           vertices[0][2]*gp[0] + vertices[1][2]*gp[1] + vertices[2][2]*gp[2]);

          point_type dist_quad = target.center - point;
          //point_type dist_quad = point - target.center;
          auto r2 = normSq(dist_quad);
          real invR2 = 1. / r2;
          if (r2 < 1e-8) invR2 = 0;
          real invR5 = invR2*invR2*std::sqrt(invR2);

          double dx = dist_quad[0], dy = dist_quad[1], dz = dist_quad[2];

          auto& normal = source.normal;
          real dxdotn = dist_quad[0]*normal[0] + dist_quad[1]*normal[1] + dist_quad[2]*normal[2];

          temp(0,0) = dx*dx; temp(0,1) = dx*dy; temp(0,2) = dx*dz;
          temp(1,0) = dx*dy; temp(1,1) = dy*dy; temp(1,2) = dy*dz;
          temp(2,0) = dx*dz; temp(2,1) = dy*dz; temp(2,2) = dz*dz;

          AnalyticalIntegral::add_influence(res, gauss_weights[i]*area*dxdotn*invR5, temp);

        }
        return res.multiply(-3);
      }
    }
    else
    {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      real area = source.Area;
      kernel_value_type res(0.), temp(0.);

      for (unsigned i=0; i<K; i++) {
        // eval contribution from this gauss point
        // point_type dist_quad = source.quad_points[i] - target.center;
        point_type dist_quad = target.center - source.quad_points[i];
        auto r2 = normSq(dist_quad);
        real invR2 = 1. / r2;
        if (r2 < 1e-8) invR2 = 0;
        real invR5 = invR2*invR2*std::sqrt(invR2);

        double dx = dist_quad[0], dy = dist_quad[1], dz = dist_quad[2];

        auto& normal = source.normal;
        real dxdotn = dist_quad[0]*normal[0] + dist_quad[1]*normal[1] + dist_quad[2]*normal[2];

        temp(0,0) = dx*dx; temp(0,1) = dx*dy; temp(0,2) = dx*dz;
        temp(1,0) = dx*dy; temp(1,1) = dy*dy; temp(1,2) = dy*dz;
        temp(2,0) = dx*dz; temp(2,1) = dy*dz; temp(2,2) = dz*dz;

        AnalyticalIntegral::add_influence(res, gauss_weight[i]*area*dxdotn*invR5, temp);
        /*
        for (int i=0; i<9; i++) {
          assert(!isnan(res.vals_[i]));
        }
        */
      }
      return res.multiply(-3);
    }
  }

  kernel_value_type eval_velocity_integral(const source_type& source, const target_type& target) const
  {
    auto dist = static_cast<point_type>(target) - source.center;
    auto d = norm(dist);
    result_type f(0.);

    bool self = d < 1e-8;

    // check for analytical / gauss quadrature
    if (sqrt(2*source.Area)/d>= 0.5)
    {
      namespace AI = AnalyticalIntegral;
      auto& vertices = source.vertices;

#if 0
      if (false) { // self) {
#else
      if (self) {
#endif
        Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[1],vertices[2],f,target.center,true,AI::G));
        //
        //Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[0],vertices[2],vertices[1],f,target.center,true,AI::G));
        //Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[1],vertices[2],vertices[0],f,target.center,true,AI::G));
        //Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[1],vertices[0],vertices[2],f,target.center,true,AI::G));
        //Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[2],vertices[1],vertices[0],f,target.center,true,AI::G));
        //Mat3<real> M(AI::FataAnalytical<AI::STOKES>(vertices[2],vertices[0],vertices[1],f,target.center,true,AI::G));

        /*
        for (int i=0; i<9; i++) {
          assert(!isnan(1./2/Mu*M.vals_[i]));
        }
        */
        //return M.multiply(1./2/Mu);
        return M.multiply(1./2./Mu);
      } else {
        // K_fine case
        auto& gauss_weights = BEMConfig::Instance()->GaussWeights(K_fine);
        auto& gauss_points  = BEMConfig::Instance()->GaussPoints(K_fine);
        real area = source.Area;
        kernel_value_type res(0.), temp(0.);

        // loop over gauss points
        for (unsigned i=0; i<gauss_points.size(); i++) {

          // generate this gauss point
          auto& gp = gauss_points[i];

          point_type point(vertices[0][0]*gp[0] + vertices[1][0]*gp[1] + vertices[2][0]*gp[2],
                           vertices[0][1]*gp[0] + vertices[1][1]*gp[1] + vertices[2][1]*gp[2],
                           vertices[0][2]*gp[0] + vertices[1][2]*gp[1] + vertices[2][2]*gp[2]);

          auto dist_quad = target.center - point;
          // auto dist_quad = point - target.center;
          auto r2 = normSq(dist_quad);
          real invR2 = 1. / r2;
          if (r2 < 1e-8) invR2 = 0;
          real invR3 = invR2*std::sqrt(invR2);
          double dx = dist_quad[0], dy = dist_quad[1], dz = dist_quad[2];

#if 1
          temp(0,0) = r2 + dx*dx; temp(0,1) = dx*dy; temp(0,2) = dx*dz;
          temp(1,0) = dx*dy; temp(1,1) = r2 + dy*dy; temp(1,2) = dy*dz;
          temp(2,0) = dx*dz; temp(2,1) = dy*dz; temp(2,2) = r2 + dz*dz;

          AnalyticalIntegral::add_influence(res, gauss_weights[i]*area*invR3, temp);
#else
          // regularized version
          auto r = std::sqrt(r2);
          double eps_int = 1e-1;//2e-3; // regularization parameter
          auto lower = 1./std::pow(r2 + eps_int*eps_int,1.5);

          temp(0,0) = r2 + 2*eps_int*eps_int + dx*dx; temp(0,1) = dx*dy; temp(0,2) = dx*dz;
          temp(1,0) = dx*dy; temp(1,1) = r2 + 2*eps_int*eps_int + dy*dy; temp(1,2) = dy*dz;
          temp(2,0) = dx*dz; temp(2,1) = dy*dz; temp(2,2) = r2 + 2*eps_int*eps_int + dz*dz;

          AnalyticalIntegral::add_influence(res, gauss_weights[i]*area*lower, temp);
          for (int i=0; i<9; i++) {
            assert(!isnan(res.val_[i]));
          }
#endif
        }
        return res.multiply(1./2/Mu);
      }
    }
    else
    {
      auto& gauss_weight = BEMConfig::Instance()->GaussWeights();
      auto area = source.Area;
      kernel_value_type res(0.), temp(0.);


      for (unsigned i=0; i<K; i++) {
        // eval contribution from this gauss point
        // auto dist_quad = source.quad_points[i] - target.center;
        auto dist_quad = target.center - source.quad_points[i];
        auto r2 = normSq(dist_quad);
        real invR2 = 1. / r2;
        if (r2 < 1e-8) invR2 = 0;
        real invR3 = invR2*std::sqrt(invR2);
        double dx = dist_quad[0], dy = dist_quad[1], dz = dist_quad[2];

        temp(0,0) = r2 + dx*dx; temp(0,1) = dx*dy; temp(0,2) = dx*dz;
        temp(1,0) = dx*dy; temp(1,1) = r2 + dy*dy; temp(1,2) = dy*dz;
        temp(2,0) = dx*dz; temp(2,1) = dy*dz; temp(2,2) = r2 + dz*dz;

        AnalyticalIntegral::add_influence(res, gauss_weight[i]*area*invR3, temp);
        /*
        for (int i=0; i<9; i++) {
          assert(!isnan(1./2/Mu*res.vals_[i]));
        }
        */

      }
      return res.multiply(1./2/Mu);
    }
  }

  kernel_value_type operator()(const target_type& t, const source_type& s) const
  {
    if (t.BC == Panel::VELOCITY) {
      return eval_velocity_integral(s,t);
    }
    else if (t.BC == Panel::TRACTION) {
      return eval_traction_integral(s,t);
    }
    else {
      printf("[E]: No BC specified on panel\n");
      return kernel_value_type(0);
    }
  }

  void P2M(const source_type& source, const charge_type& charge,
           const point_type& center, multipole_type& M) const
  {
    // setup stuff
    auto& gauss_weight = BEMConfig::Instance()->GaussWeights();

    complex Ynm[4*P*P], YnmTheta[4*P*P];

    // loop over quadrature points
    for (auto i=0u; i<source.quad_points.size(); i++) {
      auto qp = source.quad_points[i];
      auto gw = gauss_weight[i];
      point_type dist = static_cast<point_type>(qp) - center;
      // point_type dist = center - static_cast<point_type>(qp);
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
            real area = source.Area;

            auto f2 = area*gw*charge;
            real fdotx = f2[0]*qp[0] + f2[1]*qp[1] + f2[2]*qp[2];
            stokeslet_str[0] += f2[0];
            stokeslet_str[1] += f2[1];
            stokeslet_str[2] += f2[2];
            stokeslet_str[3] += fdotx;
            M[0][0][nms] += f2[0] * Ynm[nm];
            M[0][1][nms] += f2[1] * Ynm[nm];
            M[0][2][nms] += f2[2] * Ynm[nm];
            M[0][3][nms] += fdotx * Ynm[nm];
          }
          // traction specified
          else
          {
            auto& normal = source.normal;
            real n0 = normal[0], n1 = normal[1], n2 = normal[2];
            real area = source.Area;

            complex brh = (double)n/rho*Ynm[nm]; // d(rho)
            complex bal = YnmTheta[nm];          // d(alpha)
            complex bbe = -complex(0,1.)*(double)m*Ynm[nm]; // d(beta)

            complex bxd = sin(alpha)*cos(beta)*brh + cos(alpha)*cos(beta)/rho*bal - sin(beta)/rho/sin(alpha)*bbe; // dx
            complex byd = sin(alpha)*sin(beta)*brh + cos(alpha)*sin(beta)/rho*bal + cos(beta)/rho/sin(alpha)*bbe; // dy
            complex bzd = cos(alpha)*brh - sin(alpha)/rho*bal; // dz

            auto g2 = area*gw*charge;
            // which order should these be in?
            complex rdotn = bxd*n0 + byd*n1 + bzd*n2;
            complex rdotg = bxd*g2[0] + byd*g2[1] + bzd*g2[2];


            M[1][0][nms] += (rdotn * g2[0] + rdotg * n0);
            M[1][1][nms] += (rdotn * g2[1] + rdotg * n1);
            M[1][2][nms] += (rdotn * g2[2] + rdotg * n2);

            stresslet_str[0] += std::real(rdotn * g2[0] + rdotg * n0);
            stresslet_str[1] += std::real(rdotn * g2[1] + rdotg * n1);
            stresslet_str[2] += std::real(rdotn * g2[2] + rdotg * n2);

            real xdotg = qp[0]*g2[0] + qp[1]*g2[1] + qp[2]*g2[2];
            real ndotx = n0*qp[0] + n1*qp[1] + n2*qp[2];
            M[1][3][nms] += rdotn * xdotg + rdotg * ndotx;
            stresslet_str[3] += xdotg + ndotx;
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
      result += 1./2/Mu*r;
    }
    else
    {
      StokesSpherical::M2P(M[1],center,target,r);
      // need to scale result by -3 (formulation) and 1./ 6 (FMM) = -0.5
      // would then have result -= c*r, c = -0.5 => result += 0.5*r
      result += 0.5*r;
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
      result += 1./2/Mu*r;
    }
    else
    {
      StokesSpherical::L2P(L[1],center,target,r);
      result += 0.5*r;
    }
  }
};
