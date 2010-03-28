/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2006, 2007 Ferdinando Ametrano
 Copyright (C) 2005 Aurelien Chanudet
 Copyright (C) 2006, 2007 Cristina Duminuco
 Copyright (C) 2006 Giorgio Facchinetti

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

#ifndef qla_leg_hpp
#define qla_leg_hpp

#include <oh/libraryobject.hpp>

#include <ql/types.hpp>
#include <ql/cashflow.hpp>
#include <ql/compounding.hpp>
#include <ql/cashflows/cashflows.hpp>

namespace QuantLib {
    class FloatingRateCouponPricer;
}

namespace QuantLibAddin {
    class FloatingRateCouponPricer;
}

namespace QuantLibAddin {

    class Leg : public ObjectHandler::Object {
      public:
        QuantLib::Date startDate() const;
        QuantLib::Date maturityDate() const;
        bool isExpired(bool includeSettlementDateFlows,
                       const QuantLib::Date& refDate) const;

        QuantLib::Date previousCashFlowDate(bool includeSettlementDateFlows,
                                            const QuantLib::Date& refDate) const;
        QuantLib::Date nextCashFlowDate(bool includeSettlementDateFlows,
                                        const QuantLib::Date& refDate) const;
        QuantLib::Real previousCashFlowAmount(bool includeSettlementDateFlows,
                                              const QuantLib::Date& refDate) const;
        QuantLib::Real nextCashFlowAmount(bool includeSettlementDateFlows,
                                          const QuantLib::Date& refDate) const;

        QuantLib::Rate previousCouponRate(bool includeSettlementDateFlows,
                                          const QuantLib::Date& settlementDate) const;
        QuantLib::Rate nextCouponRate(bool includeSettlementDateFlows,
                                      const QuantLib::Date& settlementDate) const;
        QuantLib::Real accruedAmount(bool includeSettlementDateFlow,
                                     const QuantLib::Date& settlementDates) const;

        QuantLib::Real npv(const QuantLib::YieldTermStructure&,
                           bool includeSettlementDateFlows,
                           QuantLib::Date settlementDate,
                           const QuantLib::Date& npvDate) const;
        QuantLib::Real bps(const QuantLib::YieldTermStructure&,
                           bool includeSettlementDateFlows,
                           QuantLib::Date settlementDate,
                           const QuantLib::Date& npvDate) const;
        QuantLib::Rate atmRate(const QuantLib::YieldTermStructure&,
                               bool includeSettlementDateFlows,
                               QuantLib::Date settlementDate,
                               const QuantLib::Date& npvDate,
                               QuantLib::Real npv) const;

        QuantLib::Real npv(QuantLib::Rate y,
                           const QuantLib::DayCounter& dayCounter,
                           QuantLib::Compounding compounding,
                           QuantLib::Frequency frequency,
                           bool includeSettlementDateFlows,
                           QuantLib::Date settlementDate,
                           const QuantLib::Date& npvDate) const;
        QuantLib::Real bps(QuantLib::Rate y,
                           const QuantLib::DayCounter& dayCounter,
                           QuantLib::Compounding compounding,
                           QuantLib::Frequency frequency,
                           bool includeSettlementDateFlows,
                           QuantLib::Date settlementDate,
                           const QuantLib::Date& npvDate) const;
        QuantLib::Rate yield(QuantLib::Real npv,
                             const QuantLib::DayCounter& dayCounter,
                             QuantLib::Compounding compounding,
                             QuantLib::Frequency frequency,
                             bool includeSettlementDateFlows,
                             QuantLib::Date settlementDate,
                             const QuantLib::Date& npvDate,
                             QuantLib::Real accuracy,
                             QuantLib::Size maxIterations,
                             QuantLib::Rate guess) const;
        QuantLib::Time duration(QuantLib::Rate y,
                                const QuantLib::DayCounter& dayCounter,
                                QuantLib::Compounding compounding,
                                QuantLib::Frequency frequency,
                                QuantLib::Duration::Type type,
                                bool includeSettlementDateFlows,
                                QuantLib::Date settlementDate,
                                const QuantLib::Date& npvDate) const;
        QuantLib::Real convexity(QuantLib::Rate y,
                                 const QuantLib::DayCounter& dayCounter,
                                 QuantLib::Compounding compounding,
                                 QuantLib::Frequency frequency,
                                 bool includeSettlementDateFlows,
                                 QuantLib::Date settlementDate,
                                 const QuantLib::Date& npvDate) const;
        QuantLib::Real basisPointValue(QuantLib::Rate y,
                                       const QuantLib::DayCounter& dayCounter,
                                       QuantLib::Compounding compounding,
                                       QuantLib::Frequency frequency,
                                       bool includeSettlementDateFlows,
                                       QuantLib::Date settlementDate,
                                       const QuantLib::Date& npvDate) const;
        QuantLib::Real yieldValueBasisPoint(QuantLib::Rate y,
                                            const QuantLib::DayCounter& dayCounter,
                                            QuantLib::Compounding compounding,
                                            QuantLib::Frequency frequency,
                                            bool includeSettlementDateFlows,
                                            QuantLib::Date settlementDate,
                                            const QuantLib::Date& npvDate) const;

        QuantLib::Real npv(const boost::shared_ptr<QuantLib::YieldTermStructure>&,
                           QuantLib::Spread zSpread,
                           const QuantLib::DayCounter& dayCounter,
                           QuantLib::Compounding compounding,
                           QuantLib::Frequency frequency,
                           bool includeSettlementDateFlows,
                           QuantLib::Date settlementDate,
                           const QuantLib::Date& npvDate) const;
        QuantLib::Spread zSpread(QuantLib::Real npv,
                                 const boost::shared_ptr<QuantLib::YieldTermStructure>&,
                                 const QuantLib::DayCounter& dayCounter,
                                 QuantLib::Compounding compounding,
                                 QuantLib::Frequency frequency,
                                 bool includeSettlementDateFlows,
                                 QuantLib::Date settlementDate,
                                 const QuantLib::Date& npvDate,
                                 QuantLib::Real accuracy,
                                 QuantLib::Size maxIterations,
                                 QuantLib::Rate guess) const;

        void setCouponPricers(const std::vector<boost::shared_ptr<QuantLibAddin::FloatingRateCouponPricer> >& p);

        std::vector<std::vector<ObjectHandler::property_t> > analysis() const;
        const QuantLib::Leg& getQuantLibLeg();
      protected:
        OH_OBJ_CTOR(Leg, ObjectHandler::Object);
        // copy or shared_ptr?
        QuantLib::Leg leg_;
    };

    class MultiPhaseLeg : public Leg {
      public:
        MultiPhaseLeg(const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
                      const std::vector<boost::shared_ptr<Leg> >& legs,
                      bool toBeSorted,
                      bool permanent);
    };

    class SimpleCashFlowVector : public Leg {
      public:
        SimpleCashFlowVector(const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
                             const std::vector<QuantLib::Real>& amounts,
                             const std::vector<QuantLib::Date>& dates,
                             bool permanent);
    };

    class InterestRate : public ObjectHandler::LibraryObject<QuantLib::InterestRate> {
      public:
      InterestRate(const boost::shared_ptr<ObjectHandler::ValueObject>& properties,
                   QuantLib::Rate r,
                   const QuantLib::DayCounter& dc,
                   QuantLib::Compounding comp,
                   QuantLib::Frequency freq,
                   bool permanent);
    };
}

#endif