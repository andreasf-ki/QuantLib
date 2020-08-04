/*
 Copyright (C) 2017 Quaternion Risk Management Ltd
 All rights reserved.

 This file is part of ORE, a free-software/open-source library
 for transparent pricing and risk analysis - http://opensourcerisk.org

 ORE is free software: you can redistribute it and/or modify it
 under the terms of the Modified BSD License.  You should have received a
 copy of the license along with this program.
 The license is also available online at <http://opensourcerisk.org>

 This program is distributed on the basis that it will form a useful
 contribution to risk analytics and model standardisation, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 FITNESS FOR A PARTICULAR PURPOSE. See the license for more details.
*/

/*! \file subperiodsswap.hpp
    \brief Single currency sub periods swap helper
    \ingroup termstructures
*/

#ifndef quantext_sub_period_swap_helper_hpp_
#define quantext_sub_period_swap_helper_hpp_

#include <ql/termstructures/yield/ratehelpers.hpp>

#include <ql/ext/instruments/subperiodsswap.hpp>

//using namespace QuantLib;

namespace QuantLib {

//! Rate helper for bootstrapping using Sub Periods Swaps
/*! \ingroup termstructures
*/
class SubPeriodsSwapHelper : public RelativeDateRateHelper {
public:
    SubPeriodsSwapHelper(Handle<Quote> spread, const Period& swapTenor, const Period& fixedTenor,
                         const Calendar& fixedCalendar, const DayCounter& fixedDayCount,
                         BusinessDayConvention fixedConvention, const Period& floatPayTenor,
                         const ext::shared_ptr<IborIndex>& iborIndex, const DayCounter& floatDayCount,
                         const Handle<YieldTermStructure>& discountingCurve = Handle<YieldTermStructure>(),
                         SubPeriodsCoupon_ext::Type type = SubPeriodsCoupon_ext::Compounding);

    //! \name RateHelper interface
    //@{
    Real impliedQuote() const;
    void setTermStructure(YieldTermStructure*);
    //@}
    //! \name SubPeriodsSwapHelper inspectors
    //@{
    ext::shared_ptr<SubPeriodsSwap> swap() const { return swap_; }
    //@}
    //! \name Visitability
    //@{
    void accept(AcyclicVisitor&);
    //@}

protected:
    void initializeDates();

private:
    ext::shared_ptr<SubPeriodsSwap> swap_;
    ext::shared_ptr<IborIndex> iborIndex_;
    Period swapTenor_;
    Period fixedTenor_;
    Calendar fixedCalendar_;
    DayCounter fixedDayCount_;
    BusinessDayConvention fixedConvention_;
    Period floatPayTenor_;
    DayCounter floatDayCount_;
    SubPeriodsCoupon_ext::Type type_;

    RelinkableHandle<YieldTermStructure> termStructureHandle_;
    Handle<YieldTermStructure> discountHandle_;
    RelinkableHandle<YieldTermStructure> discountRelinkableHandle_;
};
}

#endif
