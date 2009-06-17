/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2009 Klaus Spanderen

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

/*! \file fdhestonhullwhitevanillaengine.hpp
    \brief Finite-Differences Heston Hull-White vanilla option engine
*/

#ifndef quantlib_fd_heston_hull_white_vanilla_engine_hpp
#define quantlib_fd_heston_hull_white_vanilla_engine_hpp

#include <ql/instruments/dividendvanillaoption.hpp>
#include <ql/models/equity/hestonmodel.hpp>
#include <ql/processes/hullwhiteprocess.hpp>
#include <ql/pricingengines/genericmodelengine.hpp>
#include <ql/experimental/finitedifferences/fdmhestonhullwhitesolver.hpp>

namespace QuantLib {

    //! Finite-Differences Heston Hull-White Vanilla Option engine

    /*! \ingroup vanillaengines

        \test the correctness of the returned value is tested by
              reproducing results available in web/literature
              and comparison with Black/Heston pricing.
    */
    class FdHestonHullWhiteVanillaEngine
        : public GenericModelEngine<HestonModel,
                                    DividendVanillaOption::arguments,
                                    DividendVanillaOption::results> {
      public:
        // Constructor
        FdHestonHullWhiteVanillaEngine(
            const boost::shared_ptr<HestonModel>& model,
            const boost::shared_ptr<HullWhiteProcess>& hwProcess,
            Real corrEquityShortRate,
            Size tGrid = 50, Size xGrid = 100, 
            Size vGrid = 40, Size rGrid = 20,
            bool controlVariate = true,
            FdmHestonHullWhiteSolver::FdmSchemeType type 
                                = FdmHestonHullWhiteSolver::HundsdorferScheme,
            Real theta = 0.3, Real mu = 0.5);

        void calculate() const;

      private:
        const boost::shared_ptr<HullWhiteProcess> hwProcess_;
        const Real corrEquityShortRate_;
        const Size tGrid_, xGrid_, vGrid_, rGrid_;
        const bool controlVariate_;
        const FdmHestonHullWhiteSolver::FdmSchemeType type_;
        const Real theta_, mu_;
    };
}
#endif