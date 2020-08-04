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

/*! \file ql/ext/pricingengines/discountingcurrencyswapengine.hpp
    \brief discounting currency swap engine

        \ingroup engines
*/

#ifndef quantext_discounting_currencyswap_engine_hpp
#define quantext_discounting_currencyswap_engine_hpp

#include <ql/currency.hpp>
#include <ql/termstructures/yieldtermstructure.hpp>
#include <ql/handle.hpp>
#include <ql/time/calendars/all.hpp>
#include <ql/quotes/simplequote.hpp>

#include <ql/ext/instruments/currencyswap.hpp>

namespace QuantLib {

//! %Discounting %CurrencySwap Engine

/*! This class generalizes QuantLib's DiscountingSwapEngine. It takes
    leg currencies into account and converts into the provided "npv
    currency", which must be one of the leg currencies. The evaluation
    date is the reference date of either of the discounting curves (which
    must be equal).

            \ingroup engines
*/
class DiscountingCurrencySwapEngine : public CurrencySwap::engine {
public:
    /*! The FX spots must be given as units of npvCurrency per respective
      currency. The spots must be given w.r.t. a settlement date equal
      to the npv date. */
    DiscountingCurrencySwapEngine(const std::vector<Handle<YieldTermStructure> >& discountCurves,
                                  const std::vector<Handle<Quote> >& fxQuotes, const std::vector<Currency>& currencies,
                                  const Currency& npvCurrency,
                                  boost::optional<bool> includeSettlementDateFlows = boost::none,
                                  Date settlementDate = Date(), Date npvDate = Date());
	DiscountingCurrencySwapEngine(const Handle<YieldTermStructure>& discountCurve1,
		                          const Handle<YieldTermStructure>& discountCurve2,
		                          const Handle<Quote>& fxQuote1,
		                          const Handle<Quote>& fxQuote2,
		                          const Currency& currency1,
		                          const Currency& currency2,
		                          const Currency& npvCurrency,
		                          boost::optional<bool> includeSettlementDateFlows = boost::none,
		                          Date settlementDate = Date(), Date npvDate = Date());
	DiscountingCurrencySwapEngine(const Handle<YieldTermStructure>& discountCurve1,
								  const Handle<YieldTermStructure>& discountCurve2,
								  const Currency& currency1,
							      const Currency& currency2,
		                          const Handle<Quote>& fxQuote2,
								  boost::optional<bool> includeSettlementDateFlows = boost::none,
								  Date settlementDate = Date(), Date npvDate = Date());
	DiscountingCurrencySwapEngine(const Handle<YieldTermStructure>& discountCurve1,
								  const Handle<YieldTermStructure>& discountCurve2,
								  const FxIndex& fxIndex1,
		                          const FxIndex& fxIndex2,
								  const Currency& currency1,
								  const Currency& currency2,
								  const Currency& npvCurrency,
								  boost::optional<bool> includeSettlementDateFlows = boost::none,
								  Date settlementDate = Date(), Date npvDate = Date());
	DiscountingCurrencySwapEngine(const Handle<YieldTermStructure>& discountCurveDom,
								  const Handle<YieldTermStructure>& discountCurveFor,
								  const ext::shared_ptr<FxIndex>& fxIndex,
								  boost::optional<bool> includeSettlementDateFlows = boost::none,
								  Date settlementDate = Date(), Date npvDate = Date());
	DiscountingCurrencySwapEngine(const Handle<YieldTermStructure>& discountCurve1,
		                          const Handle<YieldTermStructure>& discountCurve2,
								  const Handle<Quote>& fxQuote1,
								  const Handle<Quote>& fxQuote2,
								  const Currency& currency1,
								  const Currency& currency2,
								  const Currency& npvCurrency,
								  const Calendar& spot_calendar, const int quote_spot_lag,
								  boost::optional<bool> includeSettlementDateFlows = boost::none,
								  Date settlementDate = Date(), Date npvDate = Date());


    void calculate() const;
    std::vector<Handle<YieldTermStructure> > discountCurves() { return discountCurves_; }
    std::vector<Currency> currencies() { return currencies_; }
    Currency npvCurrency() { return npvCurrency_; }
	Calendar spotCalendar() { return spot_calendar_; }

private:
    Handle<YieldTermStructure> fetchTS(Currency ccy) const;
    Handle<Quote> fetchFX(Currency ccy) const;

    std::vector<Handle<YieldTermStructure> > discountCurves_;
    std::vector<Handle<Quote> > fxQuotes_;
    std::vector<Currency> currencies_;
	Calendar spot_calendar_;
    Currency npvCurrency_;
    boost::optional<bool> includeSettlementDateFlows_;
    Date settlementDate_;
    Date npvDate_;
	int quote_spot_lag_ = 0;
};
}

#endif
