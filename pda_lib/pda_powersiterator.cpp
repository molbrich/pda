/**
 * pda_powersiterator.cpp
 * Probability Distribution Arithmetic
 * Implementation of classes PowersIterator
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
                                       const size_t nMaxTotalPowersSum) :
                m_pda(pda),
                m_nNumberOfFactors(pda.getNumberOfDeltas() == 0 ? 0 : nNumberOfFactors),
                m_aaPowers(nNumberOfFactors, std::vector<size_t>(pda.getOrder(), 0)),
                m_aFactorsPowersSums(pda.getNumberOfDeltas(), 0),
                m_aDeltasPowersSums(nNumberOfFactors, 0),
                m_nTotalPowersSum(0),
                m_aPositions(nNumberOfFactors, 0),
                m_nMaxTotalPowersSum(nMaxTotalPowersSum == 0 ? nNumberOfFactors * pda.getNumberOfDeltas() * pda.getOrder()
                                                             : nMaxTotalPowersSum)
        {
            assert(m_nTotalPowersSum <= m_nMaxTotalPowersSum);
            assert(nNumberOfFactors > 0);
            assert(nNumberOfFactors <= 4);
            for (auto& vec : m_aaPowers)
                vec.clear();
        }

        /**
         * Operator [ostream] << [PowersIterator]
         */
        std::ostream& operator<<(std::ostream& os, const PowersIterator& pi) {
            pi.dump(os);
            return os;
        }

        /**
         * Dumps current combination
         */
        void PowersIterator::dump(std::ostream& os) const {
            for (size_t nFactor = 0; nFactor < m_nNumberOfFactors; ++nFactor) {
                os << "(Factor " << nFactor << ":";
                os << "Powers:";

                std::vector<size_t> aPowers = getFactorDeltasPowers(nFactor);

                for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                    os << aPowers[nDelta];
                }
                os << " ";
                os << "Pos: " << m_aPositions[nFactor];
                os << ")";
            }

            // PowerSums:
            os << "PowerSums:(";
            for (size_t nFactor = 0; nFactor < m_nNumberOfFactors; ++nFactor) {
                os << m_aDeltasPowersSums[nFactor];
            }
            os << ")";

            // TotalPowers:
            os << "TotalPowers:(";
            for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                os << m_aFactorsPowersSums[nDelta];
            }
            os << ")";

            os << std::endl;
        }
    }
}
