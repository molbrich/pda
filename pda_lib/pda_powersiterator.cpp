// pda_powersiterator.cpp
// Probability DeltaDistribution Arithmetic
// Implementation of classes PowersIterator
// (c) Markus Olbrich

#include "pda_powersiterator.h"
#include "pda_pda.h"
#include <cassert>
#include <iostream>

namespace Pda {
    namespace Util {

//----------------------------------------------------------------
// PowersIterator
//----------------------------------------------------------------
/**
 * Constructor
 * Intitializes the iterator
 * @param pda
 * @param nNumberOfFactors
 * @param nMaxDeltaPower Maximum delta power for each delta in each factor
 * @param nMaxTotalDeltasPowerSum Maximum sum of delta powers in the product
 * @param nMaxDeltaPowersSum Maximum Delta
 * Postcondition: The iterator points to nominal values.
 */
        PowersIterator::PowersIterator(const PDA &pda, size_t nNumberOfFactors,
                                       const size_t nMaxTotalPowersSum, const size_t nMaxTotalDeltaPower,
                                       const size_t nMaxDeltaPower) :
                m_pda(pda),
                m_nNumberOfFactors(nNumberOfFactors),
                m_aaPowers(m_nNumberOfFactors, std::vector<size_t>(pda.getNumberOfDeltas(), 0)),
                m_aFactorsPowersSums(pda.getNumberOfDeltas(), 0),
                m_aDeltasPowersSums(m_nNumberOfFactors, 0),
                m_nTotalPowersSum(0),
                m_aPositions(m_nNumberOfFactors, 0),
                m_nMaxTotalPowersSum(
                        nMaxTotalPowersSum == 0 ? nNumberOfFactors * pda.getNumberOfDeltas() * pda.getOrder()
                                                : nMaxTotalPowersSum),
                m_nMaxTotalDeltaPower(nMaxTotalDeltaPower == 0 ?
                                      std::min(nNumberOfFactors * pda.getOrder(), m_nMaxTotalPowersSum) :
                                      std::min(nMaxTotalDeltaPower, m_nMaxTotalPowersSum)),
                m_nMaxDeltaPower(nMaxDeltaPower == 0 ?
                                 std::min(pda.getOrder(), m_nMaxTotalDeltaPower) :
                                 std::min(nMaxDeltaPower, m_nMaxTotalDeltaPower)) {
            assert(m_nMaxDeltaPower <= pda.getOrder());
            assert(m_nMaxDeltaPower <= m_nMaxTotalDeltaPower);
            assert(m_nTotalPowersSum <= m_nMaxTotalPowersSum);
            assert(nNumberOfFactors > 0);
            assert(nNumberOfFactors <= 4);
        }

/** 
 * Go to next combination.
 */
        bool PowersIterator::next() {
            // NOTE: m_nMaxTotalDeltaPower is not considered
            for (size_t nFactor = 0; nFactor < m_nNumberOfFactors; ++nFactor) {
                size_t nDelta;
                do {
                    for (nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                        if (m_aaPowers[nFactor][nDelta] < m_nMaxDeltaPower
                            && m_aDeltasPowersSums[nFactor] < m_pda.getOrder()
                            && m_nTotalPowersSum < m_nMaxTotalPowersSum) {
                            ++(m_aaPowers[nFactor][nDelta]);
                            ++m_aFactorsPowersSums[nDelta];
                            ++m_aDeltasPowersSums[nFactor];
                            ++m_nTotalPowersSum;
                            break;
                        } else {
                            m_aFactorsPowersSums[nDelta] -= m_aaPowers[nFactor][nDelta];
                            m_aDeltasPowersSums[nFactor] -= m_aaPowers[nFactor][nDelta];
                            m_nTotalPowersSum -= m_aaPowers[nFactor][nDelta];
                            m_aaPowers[nFactor][nDelta] = 0;
                        }
                    }
                    m_aPositions[nFactor] = m_pda.calcCoeffPos(m_aaPowers[nFactor]);
                } while (m_nTotalPowersSum > m_nMaxTotalPowersSum && nDelta != m_pda.getNumberOfDeltas());

                if (nDelta != m_pda.getNumberOfDeltas()) {
                    // nFactor is not at end.
                    assert (m_nTotalPowersSum <= m_nMaxTotalPowersSum);
                    return true;
                } else {
                    // nFactor is at end.
                    assert (m_aDeltasPowersSums[nFactor] == 0);
                    m_aPositions[nFactor] = 0;
                }
            }
            return false;
        }

        void PowersIterator::iterateFactor(size_t nFactor, const std::function<void(void)> &callback) {
            iterateDelta(0, nFactor, callback);
            m_aPositions[nFactor] = 0;
        }

        void PowersIterator::iterateDelta(size_t nDelta, size_t nFactor, const std::function<void(void)> &callback) {
            bool incremented;
            do {
                if (nDelta < m_pda.getNumberOfDeltas() - 1)
                    iterateDelta(nDelta + 1, nFactor, callback);
                else {
                    m_aPositions[nFactor] = m_pda.calcCoeffPos(m_aaPowers[nFactor]);
                    if (m_nNumberOfFactors > 0 && nFactor < m_nNumberOfFactors - 1)
                        iterateFactor(nFactor + 1, callback);
                    else
                        callback();
                }
                incremented = false;
                if (m_aaPowers[nFactor][nDelta] < m_nMaxDeltaPower
                    && m_aDeltasPowersSums[nFactor] < m_pda.getOrder()
                    && m_nTotalPowersSum < m_nMaxTotalPowersSum) {
                    ++(m_aaPowers[nFactor][nDelta]);
                    ++m_aFactorsPowersSums[nDelta];
                    ++m_aDeltasPowersSums[nFactor];
                    ++m_nTotalPowersSum;
                    incremented = true;
                }
            } while (incremented);
            // Set this Delta power to 0:
            m_aFactorsPowersSums[nDelta] -= m_aaPowers[nFactor][nDelta];
            m_aDeltasPowersSums[nFactor] -= m_aaPowers[nFactor][nDelta];
            m_nTotalPowersSum -= m_aaPowers[nFactor][nDelta];
            m_aaPowers[nFactor][nDelta] = 0;
        }

/**
 * Iterates all combinations
 * @param callback
 */
        void PowersIterator::iterate(const std::function<void(void)> &callback) {
            size_t nNumberOfDeltas = m_pda.getNumberOfDeltas();
            std::function<void(size_t)> iterateFactor = [&](size_t nFactor) -> void {
                std::function<void(size_t)> iterateDelta = [&](size_t nDelta) -> void {
                    bool incremented;
                    do {
                        if (nDelta < nNumberOfDeltas - 1)
                            iterateDelta(nDelta + 1);
                        else {
                            m_aPositions[nFactor] = m_pda.calcCoeffPos(m_aaPowers[nFactor]);
                            if (m_nNumberOfFactors > 0 && nFactor < m_nNumberOfFactors - 1)
                                iterateFactor(nFactor + 1);
                            else
                                callback();
                        }
                        incremented = false;
                        if (m_aaPowers[nFactor][nDelta] < m_nMaxDeltaPower
                            && m_aDeltasPowersSums[nFactor] < m_pda.getOrder()
                            && m_nTotalPowersSum < m_nMaxTotalPowersSum) {
                            ++(m_aaPowers[nFactor][nDelta]);
                            ++m_aFactorsPowersSums[nDelta];
                            ++m_aDeltasPowersSums[nFactor];
                            ++m_nTotalPowersSum;
                            incremented = true;
                        }
                    } while (incremented);
                    // Set this Delta power to 0:
                    m_aFactorsPowersSums[nDelta] -= m_aaPowers[nFactor][nDelta];
                    m_aDeltasPowersSums[nFactor] -= m_aaPowers[nFactor][nDelta];
                    m_nTotalPowersSum -= m_aaPowers[nFactor][nDelta];
                    m_aaPowers[nFactor][nDelta] = 0;
                };
                iterateDelta(0);
                m_aPositions[nFactor] = 0;
            };
            if (m_pda.getNumberOfDeltas() == 0)
                callback();
            else
                this->iterateFactor(0, callback);
        }


/** 
 * Returns the coefficient position of the result
 * @returns Position
 */
        size_t PowersIterator::getPosition() const {
            assert(m_nTotalPowersSum <= m_pda.getOrder());
            return m_pda.calcCoeffPos(m_aFactorsPowersSums);
        }

/** 
 * Returns an array of nNumberOfFactors coefficient positions.
 * @returns Position array
 */
        std::vector<size_t> &PowersIterator::getPositions() {
            return m_aPositions;
        }

/**
 * Returns an array of nNumberOfDeltas power sums.
 * @see PowersIterator::m_aFactorsPowersSums
 * @returns Power sum array
 */
        std::vector<size_t> &PowersIterator::getFactorsPowersSum() {
            return m_aFactorsPowersSums;
        }

/** 
 * Returns an array of nNumberOfDeltas delta powers for a given factor.
 * @param nFactor Factor number
 * @returns Delta power array
 */
        std::vector<size_t> &PowersIterator::getFactorDeltasPowers(size_t nFactor) {
            return m_aaPowers[nFactor];
        }

/** 
 * Returns an array of nNumberOfFactors power sums.
 * @returns Power sum array
 */
        std::vector<size_t> &PowersIterator::getDeltasPowersSums() {
            return m_aDeltasPowersSums;
        }

/**
 * Dumps current combination
 */
        void PowersIterator::dump() const {
            for (size_t nFactor = 0; nFactor < m_nNumberOfFactors; ++nFactor) {
                std::cout << "(Factor " << nFactor << ":";
                std::cout << "Powers:";
                for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                    std::cout << m_aaPowers[nFactor][nDelta];
                }
                std::cout << " ";
                std::cout << "Pos: " << m_aPositions[nFactor];
                std::cout << ")";
            }

            // PowerSums:
            std::cout << "PowerSums:(";
            for (size_t nFactor = 0; nFactor < m_nNumberOfFactors; ++nFactor) {
                std::cout << m_aDeltasPowersSums[nFactor];
            }
            std::cout << ")";

            // TotalPowers:
            std::cout << "TotalPowers:(";
            for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                std::cout << m_aFactorsPowersSums[nDelta];
            }
            std::cout << ")";

            std::cout << std::endl;
        }

    }
}