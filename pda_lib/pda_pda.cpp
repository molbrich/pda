/**
 * pda_pda.cpp
 * Probability Distribution Arithmetic
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

#include "pda_pda.h"
#include <cassert>
#include <cstring>
#include <cmath>
#include <random>
#include <stdexcept>

namespace Pda {

    PDA::PDA(size_t order, size_t numberOfDeltas) :
            m_nNumberOfDeltas(numberOfDeltas),
            m_nOrder(order),
            m_nNumberOfPDVInstances(0),
            m_nSizeInner(order + 1),
            m_nSizeOuter(m_nSizeInner * m_nSizeInner),
            m_aCoeffIndexingData(numberOfDeltas * m_nSizeOuter),
            m_deltaDistributions(numberOfDeltas),
            m_deltaNames(numberOfDeltas) {
        if (order > 4)
            throw std::runtime_error("Tried to initialize PDA with order "
                                     + std::to_string(order)
                                     + ". Order must be between 0 and 4.");

        // Fill helper array m_aBinCoeffs:
        m_nMaxBinCoeffN = m_nNumberOfDeltas + m_nOrder;
        m_nMaxBinCoeffK = m_nOrder;
        m_aBinCoeffs.reserve(m_nMaxBinCoeffN + 1);
        for (size_t n = 0; n < m_nMaxBinCoeffN + 1; ++n) {
            m_aBinCoeffs.emplace_back(std::min(n + 1, m_nMaxBinCoeffK + 1), 1);
            for (size_t k = 1; k < std::min(n, m_nMaxBinCoeffK + 1); ++k)
                m_aBinCoeffs[n][k] = m_aBinCoeffs[n - 1][k] + m_aBinCoeffs[n - 1][k - 1];
        }

        // m_nNumberOfCoeffs
        m_nNumberOfCoeffs = getBinCoeff(m_nNumberOfDeltas + m_nOrder, m_nOrder);

        // Fill helper array m_aCoeffIndexingData:
        for (size_t i = 0; i < this->getNumberOfDeltas(); ++i) { // nDelta
            for (size_t j = 1; j <= this->getOrder(); ++j) {     // l
                for (size_t k = 1; k <= j; ++k) {                // aPowers[nDelta]
                    size_t l = j;
                    size_t nPos = 0;
                    for(size_t m = 0; m < k; ++m) {
                        nPos += this->getBinCoeff(i + l, l);
                        --l;
                    }
                    m_aCoeffIndexingData[i * this->m_nSizeOuter + j * this->m_nSizeInner + k] = nPos;
                }
            }
        }

        // Init Delta moments:
        //m_aDeltaMoments.resize(m_nNumberOfDeltas * (m_nMaxMoment + 1));

        // Initialize all Delta symbols as normal distributed:
        for (size_t nDelta = 0; nDelta < m_nNumberOfDeltas; ++nDelta)
            setDeltaAsNormal(nDelta, 1.0);

        // Initialize delta names:
        for (size_t nDelta = 0; nDelta < m_nNumberOfDeltas; ++nDelta)
            setDeltaName(nDelta, "D" + std::to_string(nDelta));
    }

    PDA::~PDA() {
        assert(m_nNumberOfPDVInstances == 0);
    }

    /**
     * Sets the moments of a Delta symbol according to a
     * standard normal distribution with given deviation.
     * @param nDelta Index of the Delta symbol.
     * The index of the first symbol is 0.
     * @param variance Given deviation of distribution
     */
    void PDA::setDeltaAsNormal(size_t nDelta,
                               pdaValueType variance) {
        setDeltaDistribution(nDelta, Util::NormalDistribution(sqrt(variance)));
    }

    /**
     * Determines a single raw moment of a Delta symbol.
     * @param nDelta Index of the delta symbol.
     * The index of the first symbol is 0.
     * @param nMoment Order of the central moment
     * @returns Value of the moment
     */
    pdaValueType PDA::getDeltaMoment(size_t nDelta,
                                     size_t nMoment) {
        return m_deltaDistributions[nDelta]->getMoment(nMoment);
    }

    /**
     * Sets the name of a delta symbol.
     * @param nDelta
     * @param name
     */
    void PDA::setDeltaName(const size_t nDelta, const std::string &name) {
        m_deltaNames[nDelta] = name;
    }

    /**
     * Gets the name of a delta symbol.
     * @param nDelta
     * @return Name of the delta symbol
     */
    std::string PDA::getDeltaName(const size_t nDelta) {
        return m_deltaNames[nDelta];
    }

    namespace Util {

    /**
     * Calculates raw moments of the normal distribution
     * @return Vector of moments
     */
        std::vector<pdaValueType> NormalDistribution::calculateMoments() const {
            size_t nMaxMoment = 16;
            std::vector<pdaValueType> moments(nMaxMoment + 1);
            auto doubleFaculty = [](signed long n) -> unsigned long {
                if (n < 0)
                    return 1;
                unsigned long f = 1;
                for (unsigned long i = n; i > 1; i -= 2)
                    f *= i;
                return f;
            };
            for (size_t k = 0; k <= nMaxMoment; ++k) {
                pdaValueType moment = 0;
                for (size_t i = 0; i <= k / 2; ++i) {
                    if (i == 0)
                        moment += pow(m, static_cast<pdaValueType >(k));
                    else
                        moment += static_cast<pdaValueType>(Util::getBinCoeff(k, 2 * i)) *
                                  static_cast<pdaValueType>(doubleFaculty(2 * i - 1)) *
                                  pow(s, 2 * static_cast<pdaValueType >(i)) *
                                  pow(m, static_cast<pdaValueType >(k) - 2 * static_cast<pdaValueType >(i));
                }
                moments[k] = moment;
            }
#ifndef NDEBUG
            std::vector<pdaValueType> check = {
                    1,
                    m,
                    pow(m,2)+pow(s,2),
                    pow(m,3)+3*m*pow(s,2),
                    pow(m,4)+6*pow(m,2)+3*pow(s, 4),
                    pow(m,5)+10*pow(m,3)*pow(s,2)+15*m*pow(s,4),
                    pow(m,6)+15*pow(m,4)*pow(s,2)+45*pow(m,2)*pow(s,4)+15*pow(s,6),
                    pow(m,7)+21*pow(m,5)*pow(s,2)+105*pow(m,3)*pow(s,4)+105*m*pow(s,6),
                    pow(m,8)+28*pow(m,6)*pow(s,2)+210*pow(m,4)*pow(s,4)+420*pow(m,2)*pow(s,6)+105*pow(s, 8)};
            for (size_t i=0; i<check.size(); ++i)
                assert(moments[i] == check[i]);
#if false
            std::cout << "Moments of normal distribution N(" << this->m << ", " << this->s << ")" << std::endl;
            for (size_t i = 0; i<moments.size(); ++i)
                std::cout << i << ": " << moments[i] << std::endl;
#endif
#endif
            return moments;
        }

        /**
         * Calculate raw moments of the log normal distribution
         * @return Vector of moments
         */
        std::vector<pdaValueType> LogNormalDistribution::calculateMoments() const {
            size_t nMaxMoment = 16;
            std::vector<pdaValueType> moments(nMaxMoment + 1);
            std::vector<pdaValueType> shiftedMoments(nMaxMoment + 1);
            // non shifted distribution:
            for (size_t n = 0; n <= nMaxMoment; ++n)
                moments[n] = exp(static_cast<pdaValueType>(n) * m + static_cast<pdaValueType >(n * n) * s * s / 2.0);
            auto intPow = [] (pdaValueType x, size_t n) -> pdaValueType {
                pdaValueType r = 1;
                for (size_t i=0; i<n; ++i)
                    r *= x;
                return r;
            };

            // shifted distribution:
            for (size_t n = 0; n <= nMaxMoment; ++n) {
                shiftedMoments[n] = 0;
                for (size_t k = 0; k <= n; ++k) {
                    shiftedMoments[n] += static_cast<pdaValueType>(Util::getBinCoeff(n, k)) * intPow(-getOffset(), n-k)  * moments[k];
                }
            }
#if false
            std::cout << "Moments of log normal distribution LN(" << this->m << ", " << this->s << ")" << std::endl;
            for (size_t i = 0; i<moments.size(); ++i)
                std::cout << i << ": " << shiftedMoments[i] << " (non-shifted: " << moments[i] << ")" << std::endl;
#endif
            return shiftedMoments;
        }

    } // end of namespace Util
} // end of namespace Pda
