// pda_powersiterator.h
// Probability DeltaDistribution Arithmetic
// Header for basic classes PowersIterator and PDV
// (c) Markus Olbrich

#ifndef PDA_POWERS_ITERATOR_H
#define PDA_POWERS_ITERATOR_H

#include "pda_pda.h"
#include <functional>

namespace Pda {
    namespace Util {

/** 
 * @class PowersIterator
 * Iterator for calculating the product of several PDVs.
 * \f[\prod_{i=0}^{\mbox{nNumberOfFactors}}X_i\f]
 * \f[X_i \in PDV \f]
 * \f[X = \sum\limits_{\vec{i}\in D} 
 \left(x_{i_1,\dots,i_n} \prod\limits_{j=1}^n \Delta_j^{i_j}\right)\f]
 */
        class PowersIterator {
        private:
            PowersIterator &operator=(const PowersIterator &);

            const PDA &m_pda;

            /** Number of factors which are to iterate. */
            size_t m_nNumberOfFactors;

            /** Delta powers.
             * m_aaPowers[nFactor][nDelta] */
            std::vector<std::vector<size_t>> m_aaPowers;

            /** Array of m_nNumberOfDeltas power sums in dependence of delta index.
             * aFactorsPowersSums[nDelta] is the sum of the powers of delta nDelta in all factors.
             * \f[ \mbox{aFactorsPowersSums}[nDelta]=
             \sum_{nFactor=0}^{\mbox{nNumberOfFactors}}
             \mbox{aaPowers}[nFactor][nDelta] \f] */
            std::vector<size_t> m_aFactorsPowersSums;

            /** Array of m_nNumberOfFactors power sums in dependence of factor index.
             * aDeltasPowersSums[nFactor] is the sum of the delta powers in factor nFactor.
             * \f[ \mbox{aDeltasPowersSums}[nFactor]=
             \sum_{nDelta=0}^{\mbox{nNumberOfDeltas}}
             \mbox{aaPowers}[nFactor][nDelta] \f] */
            std::vector<size_t> m_aDeltasPowersSums;

            /** Sum of all delta powers in the resulting product.
             * \f[ \mbox{nPowersPowersSum} = \sum_{nFactor=0}^{\mbox{nNumberOfFactors}} \mbox{aDeltasPowersSums}[nFactor] = \sum_{nDelta=0}^{\mbox{nNumberOfDeltas}} \mbox{aFactorsPowersSums}[nDelta] \f] */
            size_t m_nTotalPowersSum;

            /** Offset inside m_aCoeff for each factor */
            std::vector<size_t> m_aPositions;

            /** Maximum total deltas powers sum for resulting product of all factors.
             * \f[ \mbox{nTotalPowersSum} \leq \mbox{nMaxTotalPowersSum}\f]
             */
            const size_t m_nMaxTotalPowersSum;

            /** Maximum delta power for each delta in the resulting product of all factors.
             * \f[ \forall  \mbox{nDelta}\, :\, \mbox{aFactorsPowersSums}[\mbox{nDelta}] \leq \mbox{nMaxTotalDeltaPower}\f]
             */
            const size_t m_nMaxTotalDeltaPower;

            /** Maximum delta power for each delta in each factor.
             * \f[ \forall  \mbox{nDelta}, \mbox{nFactor}\, :\, \mbox{aaPowers}[\mbox{nFactor}][\mbox{nDelta}] \leq \mbox{nMaxDeltaPower}\f]
             */
            const size_t m_nMaxDeltaPower;
            void iterateDelta(size_t nDelta, size_t nFactor, const std::function<void(void)> &callback);
            void iterateFactor(size_t nFactor, const std::function<void(void)> &callback);

        public:
            PowersIterator(const PDA &pda, size_t nNumberOfFactors,
                           size_t nMaxTotalPowersSum = 0, size_t nMaxTotalDeltaPower = 0, size_t nMaxDeltaPower = 0);

            bool next();
            void iterate(const std::function<void()> &callback);
            const PDA &getPDA() { return m_pda; };
            size_t getPosition() const;
            std::vector<size_t> &getPositions();
            std::vector<size_t> &getFactorDeltasPowers(size_t nFactor);
            std::vector<size_t> &getFactorsPowersSum();
            std::vector<size_t> &getDeltasPowersSums();

            size_t getTotelPowersSum() const { return m_nTotalPowersSum; };
            size_t getNumberOfFactors() const { return m_nNumberOfFactors; };
            size_t getMaxTotalPowersSum() const { return m_nMaxTotalPowersSum; };
            size_t getMaxTotalDeltaPower() const { return m_nMaxTotalDeltaPower; };
            size_t getMaxDeltaPower() const { return m_nMaxDeltaPower; };
            void dump() const;
        };
    }
}
#endif
