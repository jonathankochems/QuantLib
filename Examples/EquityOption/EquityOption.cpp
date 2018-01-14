/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*!
 Copyright (C) 2005, 2006, 2007, 2009 StatPro Italia srl

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

#include <ql/pricingengines/vanilla/dummy.hpp>
#include <iomanip>

//#include <ql/qldefines.hpp>
//#ifdef BOOST_MSVC
//#  include <ql/auto_link.hpp>
//#endif
//#include <ql/instruments/vanillaoption.hpp>
//#include <ql/pricingengines/vanilla/binomialengine.hpp>
//#include <ql/pricingengines/vanilla/analyticeuropeanengine.hpp>
//#include <ql/pricingengines/vanilla/analytichestonengine.hpp>
//#include <ql/pricingengines/vanilla/baroneadesiwhaleyengine.hpp>
//#include <ql/pricingengines/vanilla/bjerksundstenslandengine.hpp>
//#include <ql/pricingengines/vanilla/batesengine.hpp>
//#include <ql/pricingengines/vanilla/integralengine.hpp>
//#include <ql/pricingengines/vanilla/fdeuropeanengine.hpp>
//#include <ql/pricingengines/vanilla/fdbermudanengine.hpp>
//#include <ql/pricingengines/vanilla/fdamericanengine.hpp>
//#include <ql/pricingengines/vanilla/mceuropeanengine.hpp>
//#include <ql/pricingengines/vanilla/mcamericanengine.hpp>
//#include <ql/time/calendars/target.hpp>
//#include <ql/utilities/dataformatters.hpp>
//
//#include <boost/timer.hpp>
//#include <iostream>
//#include <iomanip>
//
//using namespace QuantLib;
//
//#if defined(QL_ENABLE_SESSIONS)
//namespace QuantLib {
//
//    Integer sessionId() { return 0; }
//
//}
//#endif


int main(int, char* []) {
        double riskFreeDiscount = 1.0 ;
        double dividendDiscount = 1.0;
        double spotPrice        = 100.0;
        double strikePrice      = 100.0;
        double Tdays            = 30.0;
        double kappa            = 1.0;
        double theta            = 0.04;
        double sigma            = 2.0;
        double v0               = 0.04;
        double rho              = -0.7; 
        std::cout << "price: " << std::setprecision(12)
           << Isolated::HestonPrice(riskFreeDiscount,
                                    dividendDiscount,
                                    spotPrice,
                                    strikePrice,
                                    Tdays,
                                    kappa, 
                                    theta, 
                                    sigma, 
                                    v0, 
                                    rho) 
           << std::endl;

}
