#include <iostream>
#include <complex>
#include <vector>
#include <cassert>
#include <algorithm>

#include <boost/function.hpp>

namespace Isolated {
    typedef double Real;
    typedef Real Time;
    typedef std::size_t Size;
    typedef int Integer;

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
   //TODO(JAK): this should be an input
   int intOrder = 128; int n = intOrder;

   std::vector<Real> w_(intOrder, 0.0);
   std::vector<Real> x_(intOrder, 0.0);

   // set-up matrix to compute the roots and the weights
   std::vector<Real> e(n-1);

   Size i;
   for (i=1; i < n; ++i) {
       x_[i]  = 2*i+1;
       e[i-1] = i;
   }
   x_[0] = 1;

   std::vector<std::vector<Real> > ev(1, std::vector<Real>(x_.size(), 0));
   Size rows = 1;
   Size columns = x_.size();
   {//----------------------------

     Size iter_(0);
     std::vector<Real> d_(x_);
     const std::vector<Real>& diag(x_);
     const std::vector<Real>& sub(e);
     std::vector<std::vector<Real> >& ev_(ev);

     Size n = diag.size();

     assert(n == sub.size()+1); //, "Wrong dimensions");

     std::vector<Real> e(n, 0.0);
     std::copy(sub.begin(),sub.end(),e.begin()+1);
     for (Size i=0; i < rows; ++i) {
         ev_[i][i] = 1.0;
     }

     for (Size k=n-1; k >=1; --k) {
         while (!(std::fabs(d_[k-1])+std::fabs(d_[k]) == std::fabs(d_[k-1])+std::fabs(d_[k])+std::fabs(e[k]))) {
             Size l = k;
             while (--l > 0 && !(std::fabs(d_[l-1])+std::fabs(d_[l]) == std::fabs(d_[l-1])+std::fabs(d_[l])+std::fabs(e[l])));
             iter_++;

             Real q = d_[l];
             if (true /*strategy != NoShift*/) {
                 // calculated eigenvalue of 2x2 sub matrix of
                 // [ d_[k-1] e_[k] ]
                 // [  e_[k]  d_[k] ]
                 // which is closer to d_[k+1].
                 // FLOATING_POINT_EXCEPTION
                 const Real t1 = std::sqrt(
                                       0.25*(d_[k]*d_[k] + d_[k-1]*d_[k-1])
                                       - 0.5*d_[k-1]*d_[k] + e[k]*e[k]);
                 const Real t2 = 0.5*(d_[k]+d_[k-1]);

                 const Real lambda =
                     (std::fabs(t2+t1 - d_[k]) < std::fabs(t2-t1 - d_[k]))?
                     t2+t1 : t2-t1;

                 q-=((k==n-1)? 1.25 : 1.0)*lambda;
             }

             // the QR transformation
             Real sine = 1.0;
             Real cosine = 1.0;
             Real u = 0.0;

             bool recoverUnderflow = false;
             for (Size i=l+1; i <= k && !recoverUnderflow; ++i) {
                 const Real h = cosine*e[i];
                 const Real p = sine*e[i];

                 e[i-1] = std::sqrt(p*p+q*q);
                 if (e[i-1] != 0.0) {
                     sine = p/e[i-1];
                     cosine = q/e[i-1];

                     const Real g = d_[i-1]-u;
                     const Real t = (d_[i]-g)*sine+2*cosine*h;

                     u = sine*t;
                     d_[i-1] = g + u;
                     q = cosine*t - h;

                     for (Size j=0; j < rows; ++j) {
                         const Real tmp = ev_[j][i-1];
                         ev_[j][i-1] = sine*ev_[j][i] + cosine*tmp;
                         ev_[j][i] = cosine*ev_[j][i] - sine*tmp;
                     }
                 } else {
                     // recover from underflow
                     d_[i-1] -= u;
                     e[l] = 0.0;
                     recoverUnderflow = true;
                 }
             }

             if (!recoverUnderflow) {
                 d_[k] -= u;
                 e[k] = q;
                 e[l] = 0.0;
             }
         }
    }
 
    // sort (eigenvalues, eigenvectors),
    // code taken from symmetricSchureDecomposition.cpp
    std::vector<std::pair<Real, std::vector<Real> > > temp(n);
    std::vector<Real> eigenVector(rows);
    for (Size i=0; i<n; i++) {
        // warning hardcoded:
        //if (rows > 0)
        //    std::copy(ev_[i].begin(),
        //              ev_[i].end(), eigenVector.begin());
        eigenVector[0] = ev_[0][i];              
        temp[i] = std::make_pair(d_[i], eigenVector);
    }
    std::sort(temp.begin(), temp.end(),
              std::greater<std::pair<Real, std::vector<Real> > >());
    // first element is positive
    for (Size i=0; i<n; i++) {
        d_[i] = temp[i].first;
        Real sign = 1.0;
        if (rows > 0 && temp[i].second[0]<0.0)
            sign = -1.0;
        for (Size j=0; j<rows; ++j) {
            ev_[j][i] = sign * temp[i].second[j];
        }
    }

    x_ = d_;

   }//---------------------------
   
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
                 Real  Tdays,
                 Real  kappa, Real theta, Real sigma, Real v0, Real rho
                 ) {
    Real term = Tdays/365.0;
    const Real ratio = riskFreeDiscount/dividendDiscount;
    const Real p1 = integrate(
        Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
                  term, strikePrice, ratio, 1))/M_PI;
    const Real p2 = integrate(
        Fj_Helper(kappa, theta, sigma, v0, spotPrice, rho,
                  term, strikePrice, ratio, 2))/M_PI;

    Real value = spotPrice*dividendDiscount*(p1+0.5)
                       - strikePrice*riskFreeDiscount*(p2+0.5);
    return value;
}
}
