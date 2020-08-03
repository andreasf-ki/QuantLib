/*
 Copyright (C) 2016 Quaternion Risk Management Ltd
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

/*! \file qle/termstructures/oibasisswaphelper.hpp
    \brief Overnight Indexed Basis Swap rate helpers
    \ingroup termstructures
*/

#ifndef quantlib_oisbasisswaphelper_hpp
#define quantlib_oisbasisswaphelper_hpp

#include <ql/termstructures/yield/ratehelpers.hpp>
#include <ql/instruments/overnightindexedswap.hpp>
#include <ql/ext/instruments/oibasisswap.hpp>

namespace QuantLib {

//! Rate helper for bootstrapping over Overnight Indexed Basis Swap Spreads
/*! \ingroup termstructures
*/
class OIBSHelper : public RelativeDateRateHelper {
public:
    OIBSHelper(Natural settlementDays,
               const Period& tenor, // swap maturity
               const Handle<Quote>& oisSpread, const ext::shared_ptr<OvernightIndex>& overnightIndex,
               const ext::shared_ptr<IborIndex>& iborIndex);
    //! \name RateHelper interface
    //@{
    Real impliedQuote() const;
    void setTermStructure(YieldTermStructure*);
    //@}
    //! \name inspectors
    //@{
    ext::shared_ptr<OvernightIndexedBasisSwap> swap() const { return swap_; }
    //@}
    //! \name Visitability
    //@{
    void accept(AcyclicVisitor&);
    //@}
protected:
    void initializeDates();

    Natural settlementDays_;
    Period tenor_;
    ext::shared_ptr<OvernightIndex> overnightIndex_;
    ext::shared_ptr<IborIndex> iborIndex_;

    ext::shared_ptr<OvernightIndexedBasisSwap> swap_;
    RelinkableHandle<YieldTermStructure> termStructureHandle_;
};
}

#endif
