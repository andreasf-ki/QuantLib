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

#include <ql/ext/instruments/currencyswap.hpp>

#include <ql/cashflows/cashflows.hpp>
#include <ql/cashflows/coupon.hpp>
#include <ql/cashflows/iborcoupon.hpp>
#include <ql/cashflows/cashflowvectors.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>



namespace QuantLib {

CurrencySwap::CurrencySwap(Size nLegs) {
    legs_.resize(nLegs);
    payer_.resize(nLegs);
    currency_.resize(nLegs);
    legNPV_.resize(nLegs);
    inCcyLegNPV_.resize(nLegs);
    legBPS_.resize(nLegs);
    inCcyLegBPS_.resize(nLegs);
    startDiscounts_.resize(nLegs);
    endDiscounts_.resize(nLegs);
}

CurrencySwap::CurrencySwap(const std::vector<Leg>& legs, const std::vector<bool>& payer,
                           const std::vector<Currency>& currency)
    : legs_(legs), payer_(legs.size(), 1.0), currency_(currency), legNPV_(legs.size(), 0.0),
      inCcyLegNPV_(legs.size(), 0.0), legBPS_(legs.size(), 0.0), inCcyLegBPS_(legs.size(), 0.0),
      startDiscounts_(legs.size(), 0.0), endDiscounts_(legs.size(), 0.0), npvDateDiscount_(0.0) {
    QL_REQUIRE(payer.size() == legs_.size(), "size mismatch between payer (" << payer.size() << ") and legs ("
                                                                             << legs_.size() << ")");
    QL_REQUIRE(currency.size() == legs_.size(), "size mismatch between currency (" << currency.size() << ") and legs ("
                                                                                   << legs_.size() << ")");
    for (Size j = 0; j < legs_.size(); ++j) {
        if (payer[j])
            payer_[j] = -1.0;
        for (Leg::iterator i = legs_[j].begin(); i != legs_[j].end(); ++i)
            registerWith(*i);
    }
}

bool CurrencySwap::isExpired() const {
    for (Size j = 0; j < legs_.size(); ++j) {
        Leg::const_iterator i;
        for (i = legs_[j].begin(); i != legs_[j].end(); ++i)
            if (!(*i)->hasOccurred())
                return false;
    }
    return true;
}

void CurrencySwap::setupExpired() const {
    Instrument::setupExpired();
    std::fill(legBPS_.begin(), legBPS_.end(), 0.0);
    std::fill(legNPV_.begin(), legNPV_.end(), 0.0);
    std::fill(inCcyLegBPS_.begin(), inCcyLegBPS_.end(), 0.0);
    std::fill(inCcyLegNPV_.begin(), inCcyLegNPV_.end(), 0.0);
    std::fill(startDiscounts_.begin(), startDiscounts_.end(), 0.0);
    std::fill(endDiscounts_.begin(), endDiscounts_.end(), 0.0);
    npvDateDiscount_ = 0.0;
}

void CurrencySwap::setupArguments(PricingEngine::arguments* args) const {
    CurrencySwap::arguments* arguments = dynamic_cast<CurrencySwap::arguments*>(args);
    QL_REQUIRE(arguments != 0, "wrong argument type");

    arguments->legs = legs_;
    arguments->payer = payer_;
    arguments->currency = currency_;
}

void CurrencySwap::fetchResults(const PricingEngine::results* r) const {
    Instrument::fetchResults(r);

    const CurrencySwap::results* results = dynamic_cast<const CurrencySwap::results*>(r);
    QL_REQUIRE(results != 0, "wrong result type");

    if (!results->legNPV.empty()) {
        QL_REQUIRE(results->legNPV.size() == legNPV_.size(), "wrong number of leg NPV returned");
        legNPV_ = results->legNPV;
    } else {
        std::fill(legNPV_.begin(), legNPV_.end(), Null<Real>());
    }

    if (!results->legBPS.empty()) {
        QL_REQUIRE(results->legBPS.size() == legBPS_.size(), "wrong number of leg BPS returned");
        legBPS_ = results->legBPS;
    } else {
        std::fill(legBPS_.begin(), legBPS_.end(), Null<Real>());
    }

    if (!results->inCcyLegNPV.empty()) {
        QL_REQUIRE(results->inCcyLegNPV.size() == inCcyLegNPV_.size(), "wrong number of leg NPV returned");
        inCcyLegNPV_ = results->inCcyLegNPV;
    } else {
        std::fill(inCcyLegNPV_.begin(), inCcyLegNPV_.end(), Null<Real>());
    }

    if (!results->inCcyLegBPS.empty()) {
        QL_REQUIRE(results->inCcyLegBPS.size() == inCcyLegBPS_.size(), "wrong number of leg BPS returned");
        inCcyLegBPS_ = results->inCcyLegBPS;
    } else {
        std::fill(inCcyLegBPS_.begin(), inCcyLegBPS_.end(), Null<Real>());
    }

    if (!results->startDiscounts.empty()) {
        QL_REQUIRE(results->startDiscounts.size() == startDiscounts_.size(),
                   "wrong number of leg start discounts returned");
        startDiscounts_ = results->startDiscounts;
    } else {
        std::fill(startDiscounts_.begin(), startDiscounts_.end(), Null<DiscountFactor>());
    }

    if (!results->endDiscounts.empty()) {
        QL_REQUIRE(results->endDiscounts.size() == endDiscounts_.size(), "wrong number of leg end discounts returned");
        endDiscounts_ = results->endDiscounts;
    } else {
        std::fill(endDiscounts_.begin(), endDiscounts_.end(), Null<DiscountFactor>());
    }

    if (results->npvDateDiscount != Null<DiscountFactor>()) {
        npvDateDiscount_ = results->npvDateDiscount;
    } else {
        npvDateDiscount_ = Null<DiscountFactor>();
    }
}

Date CurrencySwap::startDate() const {
    QL_REQUIRE(!legs_.empty(), "no legs given");
    Date d = CashFlows::startDate(legs_[0]);
    for (Size j = 1; j < legs_.size(); ++j)
        d = std::min(d, CashFlows::startDate(legs_[j]));
    return d;
}

Date CurrencySwap::maturityDate() const {
    QL_REQUIRE(!legs_.empty(), "no legs given");
    Date d = CashFlows::maturityDate(legs_[0]);
    for (Size j = 1; j < legs_.size(); ++j)
        d = std::max(d, CashFlows::maturityDate(legs_[j]));
    return d;
}

void CurrencySwap::arguments::validate() const {
    QL_REQUIRE(legs.size() == payer.size(), "number of legs and multipliers differ");
    QL_REQUIRE(currency.size() == legs.size(), "number of legs and currencies differ");
}

void CurrencySwap::results::reset() {
    Instrument::results::reset();
    legNPV.clear();
    legBPS.clear();
    inCcyLegNPV.clear();
    inCcyLegBPS.clear();
    startDiscounts.clear();
    endDiscounts.clear();
    npvDateDiscount = Null<DiscountFactor>();
}

//=========================================================================
// Constructors for specialised currency swaps
//=========================================================================

VanillaCrossCurrencySwap::VanillaCrossCurrencySwap(bool payFixed, const Currency& fixedCcy, Real fixedNominal,
                                                   const Schedule& fixedSchedule, Rate fixedRate,
                                                   const DayCounter& fixedDayCount, const Currency& floatCcy,
                                                   Real floatNominal, const Schedule& floatSchedule,
                                                   const ext::shared_ptr<IborIndex>& iborIndex, Rate floatSpread,
                                                   boost::optional<BusinessDayConvention> paymentConvention)
    : CurrencySwap(4) {

    BusinessDayConvention convention;
    if (paymentConvention)
        convention = *paymentConvention;
    else
        convention = floatSchedule.businessDayConvention();

    // fixed leg
    currency_[0] = fixedCcy;
    payer_[0] = (payFixed ? -1 : +1);
    legs_[0] = FixedRateLeg(fixedSchedule)
                   .withNotionals(fixedNominal)
                   .withCouponRates(fixedRate, fixedDayCount)
                   .withPaymentAdjustment(convention);

    // add initial and final notional exchange
    currency_[1] = fixedCcy;
    payer_[1] = payer_[0];
    legs_[1].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-fixedNominal, fixedSchedule.calendar().adjust(fixedSchedule.dates().front(), convention))));
    legs_[1].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(fixedNominal, fixedSchedule.calendar().adjust(fixedSchedule.dates().back(), convention))));

    // floating leg
    currency_[2] = floatCcy;
    payer_[2] = (payFixed ? +1 : -1);
    legs_[2] = IborLeg(floatSchedule, iborIndex)
                   .withNotionals(floatNominal)
                   .withPaymentDayCounter(iborIndex->dayCounter())
                   .withPaymentAdjustment(convention)
                   .withSpreads(floatSpread);
    for (Leg::const_iterator i = legs_[2].begin(); i < legs_[2].end(); ++i)
        registerWith(*i);

    // add initial and final notional exchange
    currency_[3] = floatCcy;
    payer_[3] = payer_[2];
    legs_[3].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-floatNominal, floatSchedule.calendar().adjust(floatSchedule.dates().front(), convention))));
    legs_[3].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(floatNominal, floatSchedule.calendar().adjust(floatSchedule.dates().back(), convention))));
}

//-------------------------------------------------------------------------
CrossCurrencySwap::CrossCurrencySwap(bool payFixed, const Currency& fixedCcy, std::vector<Real> fixedNominals,
                                     const Schedule& fixedSchedule, std::vector<Rate> fixedRates,
                                     const DayCounter& fixedDayCount, const Currency& floatCcy,
                                     std::vector<Real> floatNominals, const Schedule& floatSchedule,
                                     const ext::shared_ptr<IborIndex>& iborIndex, std::vector<Rate> floatSpreads,
                                     boost::optional<BusinessDayConvention> paymentConvention)
    : CurrencySwap(4) {

    BusinessDayConvention convention;
    if (paymentConvention)
        convention = *paymentConvention;
    else
        convention = floatSchedule.businessDayConvention();

    // fixed leg
    currency_[0] = fixedCcy;
    payer_[0] = (payFixed ? -1 : +1);
    legs_[0] = FixedRateLeg(fixedSchedule)
                   .withNotionals(fixedNominals)
                   .withCouponRates(fixedRates, fixedDayCount)
                   .withPaymentAdjustment(convention);

    // add initial, interim and final notional flows
    currency_[1] = fixedCcy;
    payer_[1] = payer_[0];
    legs_[1].push_back(
        ext::shared_ptr<CashFlow>(new SimpleCashFlow(-fixedNominals[0], fixedSchedule.dates().front())));
    QL_REQUIRE(fixedNominals.size() < fixedSchedule.size(), "too many fixed nominals provided");
    for (Size i = 1; i < fixedNominals.size(); i++) {
        Real flow = fixedNominals[i - 1] - fixedNominals[i];
        legs_[1].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, fixedSchedule.calendar().adjust(fixedSchedule[i], convention))));
    }
    if (fixedNominals.back() > 0)
        legs_[1].push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(
            fixedNominals.back(), fixedSchedule.calendar().adjust(fixedSchedule.dates().back(), convention))));

    // floating leg
    currency_[2] = floatCcy;
    payer_[2] = (payFixed ? +1 : -1);
    legs_[2] = IborLeg(floatSchedule, iborIndex)
                   .withNotionals(floatNominals)
                   .withPaymentDayCounter(iborIndex->dayCounter())
                   .withPaymentAdjustment(convention)
                   .withSpreads(floatSpreads);
    for (Leg::const_iterator i = legs_[2].begin(); i < legs_[2].end(); ++i)
        registerWith(*i);

    // add initial, interim and final notional flows
    currency_[3] = floatCcy;
    payer_[3] = payer_[2];
    legs_[3].push_back(
        ext::shared_ptr<CashFlow>(new SimpleCashFlow(-floatNominals[0], floatSchedule.dates().front())));
    QL_REQUIRE(floatNominals.size() < floatSchedule.size(), "too many float nominals provided");
    for (Size i = 1; i < floatNominals.size(); i++) {
        Real flow = floatNominals[i - 1] - floatNominals[i];
        legs_[3].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, floatSchedule.calendar().adjust(floatSchedule[i], convention))));
    }
    if (floatNominals.back() > 0)
        legs_[3].push_back(ext::shared_ptr<CashFlow>(new SimpleCashFlow(
            floatNominals.back(), floatSchedule.calendar().adjust(floatSchedule.dates().back(), convention))));
}

//-------------------------------------------------------------------------
CrossCurrencySwap::CrossCurrencySwap(bool pay1, const Currency& ccy1, std::vector<Real> nominals1, const Schedule& schedule1,
                                     std::vector<Rate> rates1, const DayCounter& dayCount1, const Currency& ccy2,
                                     std::vector<Real> nominals2, const Schedule& schedule2, std::vector<Rate> rates2,
                                     const DayCounter& dayCount2,
                                     boost::optional<BusinessDayConvention> paymentConvention)
    : CurrencySwap(4) {

    BusinessDayConvention convention;
    if (paymentConvention)
        convention = *paymentConvention;
    else
        convention = schedule1.businessDayConvention();

    // fixed leg 1
    currency_[0] = ccy1;
    payer_[0] = (pay1 ? -1 : +1);
    legs_[0] = FixedRateLeg(schedule1)
                   .withNotionals(nominals1)
                   .withCouponRates(rates1, dayCount1)
                   .withPaymentAdjustment(convention);

    // add initial, interim and final notional flows
    currency_[1] = ccy1;
    payer_[1] = payer_[0];
    legs_[1].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-nominals1[0], schedule1.calendar().adjust(schedule1.dates().front(), convention))));
    QL_REQUIRE(nominals1.size() < schedule1.size(), "too many fixed nominals provided, leg 1");
    for (Size i = 1; i < nominals1.size(); i++) {
        Real flow = nominals1[i - 1] - nominals1[i];
        legs_[1].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, schedule1.calendar().adjust(schedule1[i], convention))));
    }
    if (nominals1.back() > 0)
        legs_[1].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(nominals1.back(), schedule1.calendar().adjust(schedule1.dates().back(), convention))));

    // fixed leg 2
    currency_[2] = ccy2;
    payer_[2] = (pay1 ? +1 : -1);
    legs_[2] = FixedRateLeg(schedule2)
                   .withNotionals(nominals2)
                   .withCouponRates(rates2, dayCount2)
                   .withPaymentAdjustment(convention);

    // add initial, interim and final notional flows
    currency_[3] = ccy2;
    payer_[3] = payer_[2];
    legs_[3].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-nominals2[0], schedule2.calendar().adjust(schedule2.dates().front(), convention))));
    QL_REQUIRE(nominals2.size() < schedule2.size(), "too many fixed nominals provided, leg 2");
    for (Size i = 1; i < nominals2.size(); i++) {
        Real flow = nominals2[i - 1] - nominals2[i];
        legs_[3].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, schedule2.calendar().adjust(schedule2[i], convention))));
    }
    if (nominals2.back() > 0)
        legs_[3].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(nominals2.back(), schedule2.calendar().adjust(schedule2.dates().back(), convention))));
}

//-------------------------------------------------------------------------
CrossCurrencySwap::CrossCurrencySwap(bool pay1, const Currency& ccy1, std::vector<Real> nominals1, const Schedule& schedule1,
                                     const ext::shared_ptr<IborIndex>& iborIndex1, std::vector<Rate> spreads1,
                                     const Currency& ccy2, std::vector<Real> nominals2, const Schedule& schedule2,
                                     const ext::shared_ptr<IborIndex>& iborIndex2, std::vector<Rate> spreads2,
                                     boost::optional<BusinessDayConvention> paymentConvention)
    : CurrencySwap(4) {

    BusinessDayConvention convention;
    if (paymentConvention)
        convention = *paymentConvention;
    else
        convention = schedule1.businessDayConvention();

    // floating leg 1
    currency_[0] = ccy1;
    payer_[0] = (pay1 ? -1 : +1);
    legs_[0] = IborLeg(schedule1, iborIndex1)
                   .withNotionals(nominals1)
                   .withPaymentDayCounter(iborIndex1->dayCounter())
                   .withPaymentAdjustment(convention)
                   .withSpreads(spreads1);
    for (Leg::const_iterator i = legs_[0].begin(); i < legs_[0].end(); ++i)
        registerWith(*i);

    // add initial, interim and final notional flows
    currency_[1] = ccy1;
    payer_[1] = payer_[0];
    legs_[1].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-nominals1[0], schedule1.calendar().adjust(schedule1.dates().front(), convention))));
    QL_REQUIRE(nominals1.size() < schedule1.size(), "too many float nominals provided");
    for (Size i = 1; i < nominals1.size(); i++) {
        Real flow = nominals1[i - 1] - nominals1[i];
        legs_[1].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, schedule1.calendar().adjust(schedule1[i], convention))));
    }
    if (nominals1.back() > 0)
        legs_[1].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(nominals1.back(), schedule1.calendar().adjust(schedule1.dates().back(), convention))));

    // floating leg 2
    currency_[2] = ccy2;
    payer_[2] = (pay1 ? +1 : -1);
    legs_[2] = IborLeg(schedule2, iborIndex2)
                   .withNotionals(nominals2)
                   .withPaymentDayCounter(iborIndex2->dayCounter())
                   .withPaymentAdjustment(convention)
                   .withSpreads(spreads2);
    for (Leg::const_iterator i = legs_[2].begin(); i < legs_[2].end(); ++i)
        registerWith(*i);

    // add initial, interim and final notional flows
    currency_[3] = ccy2;
    payer_[3] = payer_[2];
    legs_[3].push_back(ext::shared_ptr<CashFlow>(
        new SimpleCashFlow(-nominals2[0], schedule2.calendar().adjust(schedule2.dates().front(), convention))));
    QL_REQUIRE(nominals2.size() < schedule2.size(), "too many float nominals provided");
    for (Size i = 1; i < nominals2.size(); i++) {
        Real flow = nominals2[i - 1] - nominals2[i];
        legs_[3].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(flow, schedule2.calendar().adjust(schedule2[i], convention))));
    }
    if (nominals2.back() > 0)
        legs_[3].push_back(ext::shared_ptr<CashFlow>(
            new SimpleCashFlow(nominals2.back(), schedule2.calendar().adjust(schedule2.dates().back(), convention))));
}

//-------------------------------------------------------------------------
ResetableCrossCurrencySwap::ResetableCrossCurrencySwap(bool payDom, const Currency& ccyDom, Real nominalDomInitial, const Schedule& scheduleDom,
	const ext::shared_ptr<IborIndex>& iborIndexDom, std::vector<Rate> spreadsDom,
	const Currency& ccyFor, std::vector<Real> nominalsFor, const Schedule& scheduleFor,
	const ext::shared_ptr<IborIndex>& iborIndexFor, std::vector<Rate> spreadsFor,
	const ext::shared_ptr<FxIndex>& fxIndex, bool forecastFxToday, bool fixedNominalDomInitial,
	boost::optional<BusinessDayConvention> paymentConvention)
	: CurrencySwap(4), scheduleFor_(scheduleFor), scheduleDom_(scheduleDom), 
	  fxIndex_(fxIndex), iborIndexDom_(iborIndexDom), iborIndexFor_(iborIndexFor),
	  spreadsFor_(spreadsFor), spreadsDom_(spreadsDom),
	  nominalsFor_(nominalsFor), nominalsDom_(scheduleDom.size()-1, nominalDomInitial),
	  forecastFxToday_(forecastFxToday), fixedNominalDomInitial_(fixedNominalDomInitial) {

	if (nominalsFor_.size() < (scheduleFor_.size() - 1) && nominalsFor_.size() == 1) {
		nominalsFor_.resize(scheduleFor_.size() - 1);
		std::fill(nominalsFor_.begin(), nominalsFor_.end(), nominalsFor[0]);
	}

	if (paymentConvention)
		convention_ = *paymentConvention;
	else
		convention_ = scheduleDom.businessDayConvention();

	// floating leg 1
	currency_[0] = ccyDom;
	payer_[0] = (payDom ? -1 : +1);

	// add initial, interim and final notional flows
	currency_[1] = ccyDom;
	payer_[1] = payer_[0];
		
	// floating leg 2
	currency_[2] = ccyFor;
	payer_[2] = (payDom ? +1 : -1);
	
	// add initial, interim and final notional flows
	currency_[3] = ccyFor;
	payer_[3] = payer_[2];

	init();
	
	}
	
ResetableCrossCurrencySwap::ResetableCrossCurrencySwap(
	bool payDom, const Currency& ccyDom, Real nominalDomInitial, const Schedule& scheduleDom,
	const ext::shared_ptr<IborIndex>& iborIndexDom, Spread spreadDom,
	const Currency& ccyFor, Real nominalFor, const Schedule& scheduleFor,
	const ext::shared_ptr<IborIndex>& iborIndexFor, Spread spreadFor,
	const ext::shared_ptr<FxIndex>& fxIndex, bool forecastFxToday, bool fixedNominalDomInitial,
	boost::optional<BusinessDayConvention> paymentConvention)
	: CurrencySwap(4), scheduleFor_(scheduleFor), scheduleDom_(scheduleDom), 
	  fxIndex_(fxIndex), iborIndexDom_(iborIndexDom), iborIndexFor_(iborIndexFor),
	  spreadsFor_(scheduleFor.size()-1, spreadFor), 
	  spreadsDom_(scheduleDom.size()-1, spreadDom),
	  nominalsFor_(scheduleFor.size()-1, nominalFor), 
	  nominalsDom_(scheduleDom.size()-1, nominalDomInitial),
	  forecastFxToday_(forecastFxToday), fixedNominalDomInitial_(fixedNominalDomInitial) {

	if (paymentConvention)
		convention_ = *paymentConvention;
	else
		convention_ = scheduleDom.businessDayConvention();

	// floating leg 1
	currency_[0] = ccyDom;
	payer_[0] = (payDom ? -1 : +1);

	// add initial, interim and final notional flows
	currency_[1] = ccyDom;
	payer_[1] = payer_[0];

	// floating leg 2
	currency_[2] = ccyFor;
	payer_[2] = (payDom ? +1 : -1);

	// add initial, interim and final notional flows
	currency_[3] = ccyFor;
	payer_[3] = payer_[2];

	init();
	
	}

ResetableCrossCurrencySwap::ResetableCrossCurrencySwap(
	bool payDom, const Currency & ccyDom, const Schedule & scheduleDom, const ext::shared_ptr<IborIndex>& iborIndexDom, 
	const Currency & ccyFor, Real nominalFor, const Schedule & scheduleFor, const ext::shared_ptr<IborIndex>& iborIndexFor, 
	Spread spreadFor, const ext::shared_ptr<FxIndex>& fxIndex)
	: ResetableCrossCurrencySwap(payDom, ccyDom, -1.0, scheduleDom, iborIndexDom, 0.0,
		ccyFor, nominalFor, scheduleFor, iborIndexFor,spreadFor, fxIndex) {}
	
void ResetableCrossCurrencySwap::init() {

	registerWith(fxIndex_);
	registerWith(iborIndexDom_);
    registerWith(iborIndexFor_);

	QL_REQUIRE(scheduleFor_.size() == scheduleDom_.size(), "Schedules are not aligned. Same payment frequency required.");
	
	updateForLegFlows();
	updateDomLegFlows();
	
	}

/* foreign leg #2 #3 */
void ResetableCrossCurrencySwap::updateForLegFlows() const {
	
	legs_[3].resize(scheduleFor_.size());

	// add initial, interim and final notional flows on leg #3
	legs_[3].front() = ext::shared_ptr<CashFlow>(
		new SimpleCashFlow(-nominalsFor_[0], scheduleFor_.calendar().adjust(scheduleFor_.dates().front(), convention_)));

	for (Size i = 1; i < nominalsFor_.size(); i++) {
		Real flow = (nominalsFor_[i - 1] - nominalsFor_[i]);
		legs_[3].at(i) = ext::shared_ptr<CashFlow>(
			new SimpleCashFlow(flow, scheduleFor_.calendar().adjust(scheduleFor_[i], convention_)));
	}

	if (nominalsFor_.back() > 0)
		legs_[3].back() = ext::shared_ptr<CashFlow>(
			new SimpleCashFlow(nominalsFor_.back(), scheduleFor_.calendar().adjust(scheduleFor_.dates().back(), convention_)));

	// floating leg #2
	for (Leg::const_iterator i = legs_[2].begin(); i < legs_[2].end(); ++i)
		const_cast<ResetableCrossCurrencySwap*>(this)->unregisterWith(*i);

	legs_[2] = IborLeg(scheduleFor_, iborIndexFor_)
		.withNotionals(nominalsFor_)
		.withPaymentDayCounter(iborIndexFor_->dayCounter())
		.withPaymentAdjustment(convention_)
		.withSpreads(spreadsFor_);
	for (Leg::const_iterator i = legs_[2].begin(); i < legs_[2].end(); ++i)
		const_cast<ResetableCrossCurrencySwap*>(this)->registerWith(*i);

} 

/* domestic leg #0 #1 */
void ResetableCrossCurrencySwap::updateDomLegFlows() const  {

	legs_[1].resize(scheduleDom_.size());

	try {
		// update nominals
		if (nominalsDom_[0] < 0.0 || !fixedNominalDomInitial_)
			nominalsDom_[0] = nominalsFor_[0] * fxIndex_->fixing(fxIndex_->fixingDate(scheduleDom_[0]), forecastFxToday_);

		for (Size i = 1; i < nominalsDom_.size(); ++i)
			nominalsDom_[i] = nominalsFor_[i] * fxIndex_->fixing(fxIndex_->fixingDate(scheduleDom_[i]), forecastFxToday_);
	} catch (...) {
		try { QL_FAIL("Domestic leg nominals could not be set - set to 1."); }
		catch (...) { std::fill(nominalsDom_.begin(), nominalsDom_.end(), 1.0); }
	}

	legs_[1].front() = ext::shared_ptr<CashFlow>(
		new SimpleCashFlow(-nominalsDom_[0], scheduleDom_.calendar().adjust(scheduleDom_.dates().front(), convention_)));

	for (Size i = 1; i < nominalsDom_.size(); i++) {
		Real flow = (nominalsDom_[i - 1] - nominalsDom_[i]);
		legs_[1].at(i) = ext::shared_ptr<CashFlow>(
			new SimpleCashFlow(flow, scheduleDom_.calendar().adjust(scheduleDom_[i], convention_)));
	}

	if (nominalsDom_.back() > 0)
		legs_[1].back() = ext::shared_ptr<CashFlow>(
			new SimpleCashFlow(nominalsDom_.back(), scheduleDom_.calendar().adjust(scheduleDom_.dates().back(), convention_)));

	// floating leg - update
	for (Leg::const_iterator i = legs_[0].begin(); i < legs_[0].end(); ++i)
		const_cast<ResetableCrossCurrencySwap*>(this)->unregisterWith(*i) ;

	legs_[0] = IborLeg(scheduleDom_, iborIndexDom_)
		.withNotionals(nominalsDom_)
		.withPaymentDayCounter(iborIndexDom_->dayCounter())
		.withPaymentAdjustment(convention_)
		.withSpreads(spreadsDom_);

	for (Leg::const_iterator i = legs_[0].begin(); i < legs_[0].end(); ++i)
		const_cast<ResetableCrossCurrencySwap*>(this)->registerWith(*i) ;

}

void ResetableCrossCurrencySwap::setupArguments(PricingEngine::arguments* args) const {

	CurrencySwap::setupArguments(args);

	ResetableCrossCurrencySwap::arguments* arguments = dynamic_cast<ResetableCrossCurrencySwap::arguments*>(args);

	/* Returns here if e.g. args is CrossCcySwap::arguments which
	is the case if PricingEngine is a CrossCcySwap::engine. */
	if (!arguments)
		return;

	arguments->spreadsFor = spreadsFor_;
	arguments->spreadsDom = spreadsDom_;
}

void ResetableCrossCurrencySwap::fetchResults(const PricingEngine::results* r) const {

	CurrencySwap::fetchResults(r);

	const ResetableCrossCurrencySwap::results* results = dynamic_cast<const ResetableCrossCurrencySwap::results*>(r);
	if (results) {
		/* If PricingEngine::results are of type
		CrossCcyBasisSwap::results */
		fairForSpread_ = results->fairForSpread;
		fairDomSpread_ = results->fairDomSpread;
	}
	else {
		/* If not, e.g. if the engine is a CrossCcySwap::engine */
		fairForSpread_ = std::vector<Rate>(spreadsFor_.size(), Null<Rate>());
		fairDomSpread_ = std::vector<Rate>(spreadsDom_.size(), Null<Rate>());
	}

	/* Calculate the fair pay and receive spreads if they are null */
	static Rate basisPoint = 1.0e-4;
	if (fairForSpread_.front() == Null<Rate>()) {
		if (legBPS_[2] != Null<Real>()) // leg[2] is foreign IborLeg
			for (Size i = 0; i < spreadsFor_.size(); i++) 
				fairForSpread_[i] = spreadsFor_[i] - NPV_ / (legBPS_[2] / basisPoint);
	}
	if (fairDomSpread_.front() == Null<Rate>()) {
		if (legBPS_[0] != Null<Real>()) // leg[0] is domestic IborLeg
			for (Size i = 0; i < spreadsDom_.size(); i++)
				fairDomSpread_[i] = spreadsDom_[i] - NPV_ / (legBPS_[0] / basisPoint);
	}

}

void ResetableCrossCurrencySwap::results::reset() {
	CurrencySwap::results::reset();
	fairForSpread.clear();
	fairDomSpread.clear();
}

}
