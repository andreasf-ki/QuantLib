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

#include <ql/indexes/indexmanager.hpp>
#include <ql/ext/cashflows/fxlinkedcashflow.hpp>

namespace QuantLib {

FXLinked::FXLinked(const Date& fxFixingDate, Real foreignAmount, ext::shared_ptr<FxIndex> fxIndex, bool invertIndex)
    : fxFixingDate_(fxFixingDate), foreignAmount_(foreignAmount), fxIndex_(fxIndex), invertIndex_(invertIndex) {}

Real FXLinked::fxRate() const {
    Real fixing = fxIndex_->fixing(fxFixingDate_);
    return invertIndex_ ? 1.0 / fixing : fixing;
}

FXLinkedCashFlow::FXLinkedCashFlow(const Date& cashFlowDate, const Date& fxFixingDate, Real foreignAmount,
                                   ext::shared_ptr<FxIndex> fxIndex, bool invertIndex)
    : FXLinked(fxFixingDate, foreignAmount, fxIndex, invertIndex), cashFlowDate_(cashFlowDate) {
    registerWith(FXLinked::fxIndex());
}

} // namespace QuantExt
