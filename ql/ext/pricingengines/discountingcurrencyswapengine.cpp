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

#include <ql/ext/pricingengines/discountingcurrencyswapengine.hpp>

#include <ql/cashflows/cashflows.hpp>
#include <ql/cashflows/floatingratecoupon.hpp>
#include <ql/utilities/dataformatters.hpp>
#include <ql/exchangerate.hpp>

#include <ql/errors.hpp>
#include <ql/cashflows/cashflows.hpp>

namespace QuantLib {

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
    const std::vector<Handle<YieldTermStructure> >& discountCurves, const std::vector<Handle<Quote> >& fxQuotes,
    const std::vector<Currency>& currencies, const Currency& npvCurrency,
    boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
    : discountCurves_(discountCurves), fxQuotes_(fxQuotes), currencies_(currencies), npvCurrency_(npvCurrency),
      includeSettlementDateFlows_(includeSettlementDateFlows), settlementDate_(settlementDate), npvDate_(npvDate) {

    QL_REQUIRE(discountCurves_.size() == currencies_.size(), "Number of "
                                                             "currencies does not match number of discount curves.");
    QL_REQUIRE(fxQuotes_.size() == currencies_.size(), "Number of "
                                                       "currencies does not match number of FX quotes.");

    for (Size i = 0; i < discountCurves_.size(); i++) {
        registerWith(discountCurves_[i]);
        registerWith(fxQuotes_[i]);
    }

	spot_calendar_ = NullCalendar();
}

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
	const Handle<YieldTermStructure>& discountCurve1, const Handle<YieldTermStructure>& discountCurve2, 
	const Handle<Quote>& fxQuote1, const Handle<Quote>& fxQuote2,
	const Currency& currency1, const Currency& currency2, const Currency& npvCurrency,
	boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
	: npvCurrency_(npvCurrency),
	includeSettlementDateFlows_(includeSettlementDateFlows), settlementDate_(settlementDate), npvDate_(npvDate) {

	QL_REQUIRE(discountCurves_.size() == currencies_.size(), "Number of "
		"currencies does not match number of discount curves.");
	QL_REQUIRE(fxQuotes_.size() == currencies_.size(), "Number of "
		"currencies does not match number of FX quotes.");

	discountCurves_.resize(2);
	discountCurves_[0] = discountCurve1;
	discountCurves_[1] = discountCurve2;
	fxQuotes_.resize(2);
	fxQuotes_[0] = fxQuote1;
	fxQuotes_[1] = fxQuote2;
	currencies_.resize(2);
	currencies_[0] = currency1;
	currencies_[1] = currency2;

	spot_calendar_ = NullCalendar();

	for (Size i = 0; i < discountCurves_.size(); i++) {
		registerWith(discountCurves_[i]);
		registerWith(fxQuotes_[i]);
	}
}

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
	const Handle<YieldTermStructure>& discountCurve1, const Handle<YieldTermStructure>& discountCurve2, 
	const Currency & currency1, const Currency & currency2, const Handle<Quote>& fxQuote2,
	boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
	: npvCurrency_(currency1),
	includeSettlementDateFlows_(includeSettlementDateFlows), settlementDate_(settlementDate), npvDate_(npvDate) {

	QL_REQUIRE(discountCurves_.size() == currencies_.size(), "Number of "
		"currencies does not match number of discount curves.");
	QL_REQUIRE(fxQuotes_.size() == currencies_.size(), "Number of "
		"currencies does not match number of FX quotes.");

	discountCurves_.resize(2);
	discountCurves_[0] = discountCurve1;
	discountCurves_[1] = discountCurve2;
	fxQuotes_.resize(2);
	ext::shared_ptr<SimpleQuote> fxQuote1 = ext::shared_ptr<SimpleQuote>( new SimpleQuote(1.0) );
	fxQuotes_[0] = Handle<Quote>( fxQuote1 );
	fxQuotes_[1] = fxQuote2;
	currencies_.resize(2);
	currencies_[0] = currency1;
	currencies_[1] = currency2;

	spot_calendar_ = NullCalendar();

	for (Size i = 0; i < discountCurves_.size(); i++) {
		registerWith(discountCurves_[i]);
		registerWith(fxQuotes_[i]);
	}
}

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
	const Handle<YieldTermStructure>& discountCurve1, const Handle<YieldTermStructure>& discountCurve2, 
	const FxIndex & fxIndex1, const FxIndex & fxIndex2, 
	const Currency & currency1, const Currency & currency2, const Currency & npvCurrency, 
	boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
	: DiscountingCurrencySwapEngine(
		discountCurve1, discountCurve2, fxIndex1.fxQuote(), fxIndex2.fxQuote(), 
		currency1, currency2, npvCurrency,
		includeSettlementDateFlows, settlementDate, npvDate)
{
	QL_REQUIRE(fxIndex1.fixingDays() == fxIndex1.fixingDays(), "Number of "
		"fixing days does not match for FxIndices in PricingEngine.");
	quote_spot_lag_ = fxIndex1.fixingDays();
	spot_calendar_  = JointCalendar(fxIndex1.fixingCalendar(), fxIndex2.fixingCalendar());
}

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
	const Handle<YieldTermStructure>& discountCurveDom, const Handle<YieldTermStructure>& discountCurveFor, 
	const ext::shared_ptr<FxIndex>& fxIndex,
	boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
	: DiscountingCurrencySwapEngine(
		discountCurveDom, discountCurveFor, 
		Handle<Quote>(ext::shared_ptr<Quote>(new SimpleQuote(1.0))), fxIndex->fxQuote(),
		fxIndex->targetCurrency(), fxIndex->sourceCurrency(), fxIndex->targetCurrency(),
		includeSettlementDateFlows, settlementDate, npvDate)
{
	quote_spot_lag_ = fxIndex->fixingDays();
	spot_calendar_  = fxIndex->fixingCalendar();
}

DiscountingCurrencySwapEngine::DiscountingCurrencySwapEngine(
	const Handle<YieldTermStructure>& discountCurve1, const Handle<YieldTermStructure>& discountCurve2, 
	const Handle<Quote>& fxQuote1, const Handle<Quote>& fxQuote2, 
	const Currency & currency1, const Currency & currency2, const Currency & npvCurrency, 
	const Calendar & spot_calendar, const int quote_spot_lag, 
	boost::optional<bool> includeSettlementDateFlows, Date settlementDate, Date npvDate)
	: DiscountingCurrencySwapEngine(
		discountCurve1, discountCurve2, fxQuote1, fxQuote2,
		currency1, currency2, npvCurrency,
		includeSettlementDateFlows, settlementDate, npvDate)
{
	quote_spot_lag_ = quote_spot_lag;
	spot_calendar_  = spot_calendar;
}

Handle<YieldTermStructure> DiscountingCurrencySwapEngine::fetchTS(Currency ccy) const {
    std::vector<Currency>::const_iterator i = std::find(currencies_.begin(), currencies_.end(), ccy);
    if (i == currencies_.end())
        return Handle<YieldTermStructure>();
    else
        return discountCurves_[i - currencies_.begin()];
}

Handle<Quote> DiscountingCurrencySwapEngine::fetchFX(Currency ccy) const {
    std::vector<Currency>::const_iterator i = std::find(currencies_.begin(), currencies_.end(), ccy);
    if (i == currencies_.end())
        return Handle<Quote>();
    else
        return fxQuotes_[i - currencies_.begin()];
}

void DiscountingCurrencySwapEngine::calculate() const {
	
    for (Size i = 0; i < arguments_.currency.size(); i++) {
        Currency ccy = arguments_.currency[i];
        Handle<YieldTermStructure> yts = fetchTS(ccy);
        QL_REQUIRE(!yts.empty(), "Discounting term structure is "
                                 "empty for "
                                     << ccy.name());
        Handle<Quote> fxQuote = fetchFX(ccy);
        QL_REQUIRE(!fxQuote.empty(), "FX quote is empty "
                                     "for "
                                         << ccy.name());
    }

    Handle<YieldTermStructure> npvCcyYts = fetchTS(npvCurrency_);

    // Instrument settlement date
    Date referenceDate = npvCcyYts->referenceDate();
    Date settlementDate = settlementDate_;
    if (settlementDate_ == Date()) {
        settlementDate = referenceDate;
    } else {
        QL_REQUIRE(settlementDate >= referenceDate, "Settlement date (" << settlementDate
                                                                        << ") cannot be before discount curve "
                                                                           "reference date (" << referenceDate << ")");
    }

    // Prepare the results containers
    Size numLegs = arguments_.legs.size();

    // - Instrument::results
    if (npvDate_ == Date()) {
        results_.valuationDate = referenceDate;
    } else {
        QL_REQUIRE(npvDate_ >= referenceDate, "NPV date (" << npvDate_ << ") cannot be before "
                                                                          "discount curve reference date ("
                                                           << referenceDate << ")");
        results_.valuationDate = npvDate_;
    }
    results_.value = 0.0;
    results_.errorEstimate = Null<Real>();

    // - CurrencySwap::results
    results_.legNPV.resize(numLegs);
    results_.legBPS.resize(numLegs);
    results_.inCcyLegNPV.resize(numLegs);
    results_.inCcyLegBPS.resize(numLegs);
    results_.startDiscounts.resize(numLegs);
    results_.endDiscounts.resize(numLegs);

    bool includeRefDateFlows =
        includeSettlementDateFlows_ ? *includeSettlementDateFlows_ : Settings::instance().includeReferenceDateEvents();

    results_.npvDateDiscount = npvCcyYts->discount(results_.valuationDate);

    for (Size i = 0; i < numLegs; ++i) {
        try {
            Currency ccy = arguments_.currency[i];
            Handle<YieldTermStructure> yts = fetchTS(ccy);

            QuantLib::CashFlows::npvbps(arguments_.legs[i], **yts, includeRefDateFlows, settlementDate,
                                        results_.valuationDate, results_.inCcyLegNPV[i], results_.inCcyLegBPS[i]);

            results_.inCcyLegNPV[i] *= arguments_.payer[i];
            if (results_.inCcyLegBPS[i] != Null<Real>()) {
                results_.inCcyLegBPS[i] *= arguments_.payer[i];
            }

            // Converts into base currency and adds.
            Handle<Quote> fx = fetchFX(ccy); 

			// Adjust the Fx Quote from being spot to instantaneous Fx rate if quote_spot_lag_>0.
			Date spot_date = spot_calendar_.advance(results_.valuationDate, quote_spot_lag_, Days);
			Real fx_conv_factor = quote_spot_lag_>0 ? npvCcyYts->discount(spot_date) / yts->discount(spot_date) : 1;

            results_.legNPV[i] = results_.inCcyLegNPV[i] * fx->value() * fx_conv_factor;
            if (results_.inCcyLegBPS[i] != Null<Real>()) {
                results_.legBPS[i] = results_.inCcyLegBPS[i] * fx->value() * fx_conv_factor;
            } else {
                results_.legBPS[i] = Null<Real>();
            }

            results_.value += results_.legNPV[i];

            if (!arguments_.legs[i].empty()) {
                Date d1 = CashFlows::startDate(arguments_.legs[i]);
                if (d1 >= referenceDate)
                    results_.startDiscounts[i] = yts->discount(d1);
                else
                    results_.startDiscounts[i] = Null<DiscountFactor>();

                Date d2 = CashFlows::maturityDate(arguments_.legs[i]);
                if (d2 >= referenceDate)
                    results_.endDiscounts[i] = yts->discount(d2);
                else
                    results_.endDiscounts[i] = Null<DiscountFactor>();
            } else {
                results_.startDiscounts[i] = Null<DiscountFactor>();
                results_.endDiscounts[i] = Null<DiscountFactor>();
            }

        } catch (std::exception& e) {
            QL_FAIL("leg " << i << ": " << e.what());
        }
    }
}
}
