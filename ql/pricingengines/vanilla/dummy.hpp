#include <ql/types.hpp>
#include <ql/math/matrixutilities/tqreigendecomposition.hpp>
#include <iostream>
#include <boost/function.hpp>

namespace QuantLib {
    class Fj_Helper
        : public std::unary_function<Real, Real>
    {
    public:
        Fj_Helper(Real kappa, Real theta, Real sigma,
            Real v0, Real s0, Real rho,
            Time term,
            Real strike,
            Real ratio,
            Size j);

        Real operator()(Real phi) const;

    private:
        const Size j_;
        //     const VanillaOption::arguments& arg_;
        const Real kappa_, theta_, sigma_, v0_;

        // helper variables
        const Time term_;
        const Real x_, sx_, dd_;
        const Real sigma2_, rsigma_;
        const Real t0_;

        // log branch counter
        mutable int  b_;     // log branch counter
        mutable Real g_km1_; // imag part of last log value

    };

    Fj_Helper::Fj_Helper(Real kappa, Real theta,
        Real sigma, Real v0, Real s0, Real rho,
        Time term,
        Real strike,
        Real ratio,
        Size j)
        :
        j_(j),
        kappa_(kappa),
        theta_(theta),
        sigma_(sigma),
        v0_(v0),
        term_(term),
        x_(std::log(s0)),
        sx_(std::log(strike)),
        dd_(x_-std::log(ratio)),
        sigma2_(sigma_*sigma_),
        rsigma_(rho*sigma_),
        t0_(kappa - ((j== 1)? rho*sigma : 0)),
        b_(0),
        g_km1_(0)
    {
    }

    Real Fj_Helper::operator()(Real phi) const
    {
        const Real rpsig(rsigma_*phi);

        const std::complex<Real> t1 = t0_+std::complex<Real>(0, -rpsig);
        const std::complex<Real> d =
            std::sqrt(t1*t1 - sigma2_*phi
                      *std::complex<Real>(-phi, (j_== 1)? 1 : -1));
        const std::complex<Real> ex = std::exp(-d*term_);
        const std::complex<Real> addOnTerm
            = Real(0.0);

        if (phi != 0.0) {
            if (sigma_ > 1e-5) {
                const std::complex<Real> p = (t1-d)/(t1+d);
                const std::complex<Real> g
                                        = std::log((1.0 - p*ex)/(1.0 - p));

                return
                    std::exp(v0_*(t1-d)*(1.0-ex)/(sigma2_*(1.0-ex*p))
                             + (kappa_*theta_)/sigma2_*((t1-d)*term_-2.0*g)
                             + std::complex<Real>(0.0, phi*(dd_-sx_))
                             + addOnTerm
                             ).imag()/phi;
            }
            else {
                const std::complex<Real> td = phi/(2.0*t1)
                               *std::complex<Real>(-phi, (j_== 1)? 1 : -1);
                const std::complex<Real> p = td*sigma2_/(t1+d);
                const std::complex<Real> g = p*(1.0-ex);

                return
                    std::exp(v0_*td*(1.0-ex)/(1.0-p*ex)
                             + (kappa_*theta_)*(td*term_-2.0*g/sigma2_)
                             + std::complex<Real>(0.0, phi*(dd_-sx_))
                             + addOnTerm
                             ).imag()/phi;
            }
        }
        else {
            // use l'Hospital's rule to get lim_{phi->0}
            if (j_ == 1) {
                const Real kmr = rsigma_-kappa_;
                if (std::fabs(kmr) > 1e-7) {
                    return dd_-sx_
                        + (std::exp(kmr*term_)*kappa_*theta_
                           -kappa_*theta_*(kmr*term_+1.0) ) / (2*kmr*kmr)
                        - v0_*(1.0-std::exp(kmr*term_)) / (2.0*kmr);
                }
                else
                    // \kappa = \rho * \sigma
                    return dd_-sx_ + 0.25*kappa_*theta_*term_*term_
                                   + 0.5*v0_*term_;
            }
            else {
                return dd_-sx_
                    - (std::exp(-kappa_*term_)*kappa_*theta_
                       +kappa_*theta_*(kappa_*term_-1.0))/(2*kappa_*kappa_)
                    - v0_*(1.0-std::exp(-kappa_*term_))/(2*kappa_);
            }
        }
    }

Real integrate( const boost::function1<Real, Real>& f) {
   using std::cout; using std::endl;
   cout << "integrate" << endl;

   //TODO(JAK): this should be an input
   int intOrder = 128; int n = intOrder;

   Array w_(intOrder);
   Array x_(intOrder);

   // set-up matrix to compute the roots and the weights
   Array e(n-1);

   Size i;
   for (i=1; i < n; ++i) {
       x_[i]  = 2*i+1;
       e[i-1] = i;
   }
   x_[0] = 1;

   TqrEigenDecomposition tqr(
                          x_, e,
                          TqrEigenDecomposition::OnlyFirstRowEigenVector,
                          TqrEigenDecomposition::Overrelaxation);

   x_ = tqr.eigenvalues();
   const Matrix& ev = tqr.eigenvectors();
   
   
   Real mu_0 = 1.0;
   for (i=0; i<n; ++i) {
       w_[i] = mu_0*ev[0][i]*ev[0][i] / std::exp(-x_[i]);
   }

   Real sum = 0.0;
   for (Integer i = intOrder-1; i >= 0; --i) {
       sum += w_[i] * f(x_[i]);
   }
   return sum;
}

Real HestonPrice(Real  riskFreeDiscount,
                 Real  dividendDiscount,
                 Real  spotPrice,
                 Real  strikePrice,
                 Real  term,
                 Real  kappa, Real theta, Real sigma, Real v0, Real rho
                 ) {
    using std::cout; using std::endl;
    cout << "HestonPrice" << endl;
    const Real ratio = riskFreeDiscount/dividendDiscount;
    const Real c_inf = std::min(0.2, std::max(0.0001,
        std::sqrt(1.0-rho*rho)/sigma))*(v0 + kappa*theta*term);

    //const Real p1 = integration.calculate(c_inf,
    //    Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
    //              term, strikePrice, ratio, 1))/M_PI;
    const Real p1 = integrate(
        Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
                  term, strikePrice, ratio, 1))/M_PI;
    //int evaluations = 0;
    //evaluations += integration.numberOfEvaluations();

    //const Real p2 = integration.calculate(c_inf,
    //    Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
    //              term, strikePrice, ratio, 2))/M_PI;
    const Real p2 = integrate(
        Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
                  term, strikePrice, ratio, 2))/M_PI;
    //evaluations += integration.numberOfEvaluations();

    Real value = spotPrice*dividendDiscount*(p1+0.5)
                       - strikePrice*riskFreeDiscount*(p2+0.5);
    return value;
}
}
