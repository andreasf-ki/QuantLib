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

/*! \file ql/ext/instruments/currencyswap.hpp
    \brief Interest rate swap with extended interface

        \ingroup instruments
*/

#ifndef quantext_currencyswap_hpp
#define quantext_currencyswap_hpp

#include <ql/instruments/swap.hpp>
#include <ql/time/daycounter.hpp>
#include <ql/time/schedule.hpp>
#include <ql/cashflows/cashflows.hpp>
#include <ql/cashflows/simplecashflow.hpp>
#include <ql/indexes/iborindex.hpp>
#include <ql/currency.hpp>
#include <ql/currencies/europe.hpp>
#include <ql/money.hpp>
#include <ql/ext/indexes/fxindex.hpp>
// #include <ql/ext/cashflows/fxlinkedcashflow.hpp>

#include <iostream>

namespace QuantLib {

//! %Currency Interest Rate Swap
/*!
  This instrument generalizes the QuantLib Swap instrument in that
  it allows multiple legs with different currencies (one per leg)

  \ingroup instruments
*/
class CurrencySwap : public Instrument {
public:
    class arguments;
    class results;
    class engine;
    //! \name Constructors
    //@{
    /*! Multi leg constructor. */
    CurrencySwap(const std::vector<Leg>& legs, const std::vector<bool>& payer, const std::vector<Currency>& currency);
    //@}
    //! \name Instrument interface
    //@{
    bool isExpired() const;
    void setupArguments(PricingEngine::arguments*) const;
    void fetchResults(const PricingEngine::results*) const;
    //@}
    //! \name Additional interface
    //@{
    Date startDate() const;
    Date maturityDate() const;
    Real legBPS(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg# " << j << " doesn't exist!");
        calculate();
        return legBPS_[j];
    }
    Real legNPV(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        calculate();
        return legNPV_[j];
    }
    Real inCcyLegBPS(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg# " << j << " doesn't exist!");
        calculate();
        return inCcyLegBPS_[j];
    }
    Real inCcyLegNPV(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        calculate();
        return inCcyLegNPV_[j];
    }
    DiscountFactor startDiscounts(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        calculate();
        return startDiscounts_[j];
    }
    DiscountFactor endDiscounts(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        calculate();
        return endDiscounts_[j];
    }
    DiscountFactor npvDateDiscount() const {
        calculate();
        return npvDateDiscount_;
    }
    const Leg& leg(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        return legs_[j];
    }
    const Currency& legCurrency(Size j) const {
        QL_REQUIRE(j < legs_.size(), "leg #" << j << " doesn't exist!");
        return currency_[j];
    }
    std::vector<Leg> legs() { return legs_; }
    std::vector<Currency> currencies() { return currency_; }
    ext::shared_ptr<PricingEngine> engine() { return engine_; }
    //@}
protected:
    //! \name Constructors
    //@{
    /*! This constructor can be used by derived classes that will
        build their legs themselves.
    */
    CurrencySwap(Size legs);
    //! \name Instrument interface
    //@{
    void setupExpired() const;
    //@}
    // data members
    mutable std::vector<Leg> legs_;
    std::vector<Real> payer_;
    std::vector<Currency> currency_;
    mutable std::vector<Real> legNPV_, inCcyLegNPV_;
    mutable std::vector<Real> legBPS_, inCcyLegBPS_;
    mutable std::vector<DiscountFactor> startDiscounts_, endDiscounts_;
    mutable DiscountFactor npvDateDiscount_;
};

class CurrencySwap::arguments : public virtual PricingEngine::arguments {
public:
    std::vector<Leg> legs;
    std::vector<Real> payer;
    std::vector<Currency> currency;
    void validate() const;
};

class CurrencySwap::results : public Instrument::results {
public:
    std::vector<Real> legNPV, inCcyLegNPV;
    std::vector<Real> legBPS, inCcyLegBPS;
    std::vector<DiscountFactor> startDiscounts, endDiscounts;
    DiscountFactor npvDateDiscount;
    void reset();
};

class CurrencySwap::engine : public GenericEngine<CurrencySwap::arguments, CurrencySwap::results> {};

//! Vanilla cross currency interest rate swap
/*! Specialised CurrencySwap: Two currencies, fixed vs. floating,
  constant notionals, rate and spread.

      \ingroup instruments{
*/
class VanillaCrossCurrencySwap : public CurrencySwap {
public:
    VanillaCrossCurrencySwap(bool payFixed, const Currency& fixedCcy, Real fixedNominal, const Schedule& fixedSchedule,
                             Rate fixedRate, const DayCounter& fixedDayCount, const Currency& floatCcy, Real floatNominal,
                             const Schedule& floatSchedule, const ext::shared_ptr<IborIndex>& iborIndex,
                             Rate floatSpread, boost::optional<BusinessDayConvention> paymentConvention = boost::none);
};

//! Cross currency swap
/*! Specialised CurrencySwap: Two currencies, variable notionals,
    rates and spreads; flavours fix/float, fix/fix, float/float

        \ingroup instruments
*/
class CrossCurrencySwap : public CurrencySwap {
public:
	// fixed/floating
	CrossCurrencySwap(bool payFixed, const Currency& fixedCcy, std::vector<Real> fixedNominals, const Schedule& fixedSchedule,
		std::vector<Rate> fixedRates, const DayCounter& fixedDayCount, const Currency& floatCcy,
		std::vector<Real> floatNominals, const Schedule& floatSchedule,
		const ext::shared_ptr<IborIndex>& iborIndex, std::vector<Rate> floatSpreads,
		boost::optional<BusinessDayConvention> paymentConvention = boost::none);

	// fixed/fixed
	CrossCurrencySwap(bool pay1, const Currency& ccy1, std::vector<Real> nominals1, const Schedule& schedule1,
		std::vector<Rate> fixedRates1, const DayCounter& fixedDayCount1, const Currency& ccy2,
		std::vector<Real> nominals2, const Schedule& schedule2, std::vector<Rate> fixedRates2,
		const DayCounter& fixedDayCount2,
		boost::optional<BusinessDayConvention> paymentConvention = boost::none);

	// floating/floating
	CrossCurrencySwap(bool pay1, const Currency& ccy1, std::vector<Real> nominals1, const Schedule& schedule1,
		const ext::shared_ptr<IborIndex>& iborIndex1, std::vector<Rate> spreads1, const Currency& ccy2,
		std::vector<Real> nominals2, const Schedule& schedule2,
		const ext::shared_ptr<IborIndex>& iborIndex2, std::vector<Rate> spreads2,
		boost::optional<BusinessDayConvention> paymentConvention = boost::none);
};

class ResetableCrossCurrencySwap : public CurrencySwap {
public:
	class arguments;
	class results;

	// floating/floating w Mark-to-Market notional resets on the domestic (Dom) leg
	ResetableCrossCurrencySwap(bool payDom, const Currency& ccyDom, Real nominalDomInitial, const Schedule& scheduleDom,
		const ext::shared_ptr<IborIndex>& iborIndexDom, std::vector<Rate> spreadsDom, const Currency& ccyFor,
		std::vector<Real> nominalsFor, const Schedule& scheduleFor,
		const ext::shared_ptr<IborIndex>& iborIndexFor, std::vector<Rate> spreadsFor,
		const ext::shared_ptr<FxIndex>& fxIndex, bool forecastFxToday = false, bool fixedNominalDomInitial = false,
		boost::optional<BusinessDayConvention> paymentConvention = boost::none);
	ResetableCrossCurrencySwap(bool payDom, const Currency& ccyDom, Real nominalDomInitial, const Schedule& scheduleDom,
		const ext::shared_ptr<IborIndex>& iborIndexDom, Spread spreadDom, const Currency& ccyFor,
		Real nominalFor, const Schedule& scheduleFor,
		const ext::shared_ptr<IborIndex>& iborIndexFor, Spread spreadFor,
		const ext::shared_ptr<FxIndex>& fxIndex, bool forecastFxToday = false, bool fixedNominalDomInitial = false,
		boost::optional<BusinessDayConvention> paymentConvention = boost::none);
	ResetableCrossCurrencySwap(bool payDom, const Currency& ccyDom, const Schedule& scheduleDom,
		const ext::shared_ptr<IborIndex>& iborIndexDom, const Currency& ccyFor,
		Real nominalFor, const Schedule& scheduleFor,
		const ext::shared_ptr<IborIndex>& iborIndexFor, Spread spreadFor,
		const ext::shared_ptr<FxIndex>& fxIndex);

	void setupArguments(PricingEngine::arguments* args) const;
	void fetchResults(const PricingEngine::results*) const;

	void performCalculations() const { /*updateForLegFlows();*/ updateDomLegFlows(); CurrencySwap::performCalculations(); }

	void setFxForecast(bool forecast) { forecastFxToday_ = forecast; update(); }
	bool forecastFxToday() { return forecastFxToday_; }

	const Schedule& domesticSchedule() const { return scheduleDom_; }
	const Schedule& foreignSchedule() const { return scheduleFor_; }
	
	const ext::shared_ptr<FxIndex>& fxIndex() { return fxIndex_; }
	const ext::shared_ptr<IborIndex>& foreignIndex() const { return iborIndexFor_; }
	const ext::shared_ptr<IborIndex>& domesticIndex() const { return iborIndexDom_; }
	
	std::vector<Real> domesticNominals() { updateDomLegFlows(); return nominalsDom_; }
	std::vector<Real> foreignNominals() { updateForLegFlows(); return nominalsFor_; }
	
	bool payDom() { return (payer_[0] == +1 ? false : true); }
	
	void setForeignSpread(std::vector<Rate> spreads) { spreadsFor_ = spreads; updateForLegFlows(); };
	
	Real inCcyDomLegNPV() {
		calculate();
		return inCcyLegNPV_[0] + inCcyLegNPV_[1];
	}
	Real inCcyDomLegBPS() {
		calculate();
		return inCcyLegBPS_[0] + inCcyLegBPS_[1];
	}
	Real inCcyForLegNPV() {
		calculate();
		return inCcyLegNPV_[2] + inCcyLegNPV_[3];
	}
	Real inCcyForLegBPS() {
		calculate();
		return inCcyLegBPS_[2] + inCcyLegBPS_[3];
	}
	Real domLegNPV() {
		calculate();
		return legNPV_[0] + legNPV_[1];
	}
	Real forLegNPV() {
		calculate();
		return legNPV_[2] + legNPV_[3];
	}
	std::vector<Rate> fairForSpreads() {
		calculate();
		QL_REQUIRE(!fairForSpread_.empty(), "Fair foreign leg spread is not available");
		return fairForSpread_;
	}
	std::vector<Rate> fairDomSpreads() {
		calculate();
		QL_REQUIRE(!fairDomSpread_.empty(), "Fair domestic leg spread is not available");
		return fairDomSpread_;
	}
	Spread fairForSpread() {
		return fairForSpreads()[0];
	}
	Spread fairDomSpread() {
		return fairDomSpreads()[0];
	}

private:
	void init();
	void updateForLegFlows() const; // bool hasUpdateForeignFlow_ = false;
	void updateDomLegFlows() const; // bool hasUpdateDomesticFlow_ = false;

	const Schedule scheduleFor_;
	const Schedule scheduleDom_;
	
	ext::shared_ptr<IborIndex> iborIndexFor_;
	ext::shared_ptr<IborIndex> iborIndexDom_;
	ext::shared_ptr<FxIndex> fxIndex_;
	
	bool forecastFxToday_ = true;
	bool fixedNominalDomInitial_ = true;
	
	mutable std::vector<Real> nominalsFor_;
	mutable std::vector<Real> nominalsDom_;
	
	std::vector<Rate> spreadsFor_;
	std::vector<Rate> spreadsDom_;

	BusinessDayConvention convention_;

	mutable std::vector<Rate> fairForSpread_;
	mutable std::vector<Rate> fairDomSpread_;
};

class ResetableCrossCurrencySwap::arguments : public CurrencySwap::arguments {
public:
	std::vector<Rate> spreadsFor;
	std::vector<Rate> spreadsDom;
	//void validate() const;
};

class ResetableCrossCurrencySwap::results : public CurrencySwap::results {
public:
	std::vector<Rate> fairForSpread;
	std::vector<Rate> fairDomSpread;
	void reset();
};

};

#endif
