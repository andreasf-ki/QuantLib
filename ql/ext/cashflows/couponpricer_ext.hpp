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

/*! \file ql/ext/cashflows/couponpricer.hpp
    \brief Utility functions for setting coupon pricers on legs
        \ingroup cashflows
*/

#ifndef quantext_coupon_pricer_hpp
#define quantext_coupon_pricer_hpp

#include <boost/shared_ptr.hpp>

#include <ql/cashflows/couponpricer.hpp>

#include <ql/ext/cashflows/averageonindexedcouponpricer.hpp>
#include <ql/ext/cashflows/subperiodscouponpricer.hpp>



namespace QuantLib {
/*!	\addtogroup cashflows
    @{
*/
//! Set Coupon Pricer
void setCouponPricer_ext(const Leg& leg, const ext::shared_ptr<FloatingRateCouponPricer>&);
//! Set Coupon Pricers
void setCouponPricers_ext(const Leg& leg, const std::vector<ext::shared_ptr<FloatingRateCouponPricer> >&);

// @}
}

#endif
