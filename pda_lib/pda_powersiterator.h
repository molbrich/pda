/**
 * pda_powersiterator.h
 * Probability Distribution Arithmetic
 * Header for basic classes PowersIterator and PDV
 * (c) Markus Olbrich
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

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
            std::vector<size_t> m_aTemp;

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

        public:
            PowersIterator(const PDA &pda, size_t nNumberOfFactors, size_t nMaxTotalPowersSum = 0);

            inline bool next();
            const PDA &getPDA() const { return m_pda; }
            inline size_t getPosition();
            inline const std::vector<size_t>& getPositions() const;
            inline std::vector<size_t> getFactorDeltasPowers(size_t nFactor) const;
            inline const std::vector<size_t>& getFactorsPowersSum() const;
            inline const std::vector<size_t>& getDeltasPowersSums() const;

            size_t getTotelPowersSum() const { return m_nTotalPowersSum; }
            size_t getNumberOfFactors() const { return m_nNumberOfFactors; }
            size_t getMaxTotalPowersSum() const { return m_nMaxTotalPowersSum; }

            friend std::ostream& operator<<(std::ostream&, const PowersIterator&);
            void dump(std::ostream&) const;
        };

        /**
         * Returns the coefficient position of the result
         * @returns Position
         */
        size_t PowersIterator::getPosition() {
            assert(m_nTotalPowersSum <= m_pda.getOrder());
            m_aTemp.clear();
            for(const auto& aPowerIndices : m_aaPowers)
                for(const size_t& deltaIdx : aPowerIndices)
                    m_aTemp.emplace(std::upper_bound(m_aTemp.begin(), m_aTemp.end(), deltaIdx, std::greater<>()),
                                    deltaIdx);
            size_t i = 0;
            size_t deltaWeight = m_pda.m_nOrder;
            for (const size_t& deltaIdx : m_aTemp)
            {
                i += m_pda.getBinCoeff(deltaIdx + deltaWeight, deltaWeight);
                --deltaWeight;
            }
            return i;
        }

        /**
         * Returns an array of nNumberOfFactors coefficient positions.
         * @returns Position array
         */
        const std::vector<size_t>& PowersIterator::getPositions() const {
            return m_aPositions;
        }

        /**
         * Returns an array of nNumberOfDeltas power sums.
         * @see PowersIterator::m_aFactorsPowersSums
         * @returns Power sum array
         */
        const std::vector<size_t>& PowersIterator::getFactorsPowersSum() const {
            return m_aFactorsPowersSums;
        }

        /**
         * Returns an array of nNumberOfDeltas delta powers for a given factor.
         * Costly method which should only be used for debug purposes.
         * @param nFactor Factor number
         * @returns Delta power array
         */
        std::vector<size_t> PowersIterator::getFactorDeltasPowers(size_t nFactor) const {
            std::vector<size_t> aPowers(m_pda.getNumberOfDeltas(), 0);
            for (const size_t& nDelta : m_aaPowers[nFactor])
                ++aPowers[nDelta];
            return aPowers;
        }

        /**
         * Returns an array of nNumberOfFactors power sums.
         * @returns Power sum array
         */
        const std::vector<size_t>& PowersIterator::getDeltasPowersSums() const {
            return m_aDeltasPowersSums;
        }

        /**
         * Go to next combination.
         */
        bool PowersIterator::next() {
            for (size_t factorIdx = 0; factorIdx < m_nNumberOfFactors; ++factorIdx) {
                std::vector<size_t>& factorPowerIndices = m_aaPowers[factorIdx];
                for (auto powIdxIt = factorPowerIndices.begin(); powIdxIt < factorPowerIndices.end(); ++powIdxIt) {
                    const size_t nextDeltaIndex = *powIdxIt + 1;
                    if (nextDeltaIndex < m_pda.getNumberOfDeltas()) {
                        size_t deltaPowWgt = m_pda.getOrder();
                        for (auto powIdxIt2 = factorPowerIndices.begin(); powIdxIt2 <= powIdxIt; ++powIdxIt2) {
                            m_aPositions[factorIdx] -= m_pda.getBinCoeff(*powIdxIt2 + deltaPowWgt, deltaPowWgt);
                            m_aPositions[factorIdx] += m_pda.getBinCoeff(nextDeltaIndex + deltaPowWgt, deltaPowWgt);

                            --m_aFactorsPowersSums[*powIdxIt2];
                            ++m_aFactorsPowersSums[nextDeltaIndex];

                            *powIdxIt2 = nextDeltaIndex;

                            --deltaPowWgt;
                        }

                        return true;
                    }
                }
                if (factorPowerIndices.size() < m_pda.getOrder() && m_nTotalPowersSum < m_nMaxTotalPowersSum) {
                    m_aFactorsPowersSums.back() -= factorPowerIndices.size();

                    std::fill(factorPowerIndices.begin(), factorPowerIndices.end(), 0);
                    factorPowerIndices.push_back(0);

                    m_aFactorsPowersSums[0] += factorPowerIndices.size();
                    m_aPositions[factorIdx] = factorPowerIndices.size();

                    ++m_aDeltasPowersSums[factorIdx];
                    ++m_nTotalPowersSum;

                    return true;
                } else {
                    m_aFactorsPowersSums.back() -= factorPowerIndices.size();
                    m_nTotalPowersSum -= factorPowerIndices.size();

                    factorPowerIndices.clear();

                    m_aDeltasPowersSums[factorIdx] = 0;
                    m_aPositions[factorIdx] = 0;
                }
            }
            return false;
        }
    }
}
#endif
