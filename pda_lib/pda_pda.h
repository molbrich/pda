// pda_pda.h
// Probability DeltaDistribution Arithmetic
// (c) Markus Olbrich

#ifndef PDA_PDA_H
#define PDA_PDA_H

#include <iosfwd>
#include <vector>
#include <iostream>
#include <memory>
#include <cmath>
#include <random>
#include <cassert>
#include <bits/unique_ptr.h>

namespace Pda {

/** Type of real numbers inside PDV */
    using pdaValueType = double;
/** Forward declarations */
    class PDA;
    class PDV;
    namespace Util {
        struct DeltaDistribution;
        class PowersIterator;
        size_t getBinCoeff(size_t n, size_t k);
        std::vector<pdaValueType> monteCarlo(
                const std::vector<PDV*>& aArguments,
                PDV (*function)(PDA& pda),
                size_t nSampleNumber);
    }

/**
 * @class PDA
 * @brief Probability DeltaDistribution Arithmetic
 * @author Markus Olbrich
 *
 * This class implements a distribution arithmetic.
 *
 * There is a fixed number of Delta symbols which represent
 * the source of any variation.
 * Delta symbols are random variables.
 * Their distribution is given by the central moments
 * of their PDF (probability density function).
 * The mean (or expectancy value) of each Delta symbol is zero.
 *
 * Each PDV is a sum of weighted Delta symbols (BRVs).
 * PDVs can be used like normal float or double values
 * including arithmetic operations.
 * At any time, the nominal value, the mean value (expectancy value)
 * and central moments of the PDF of a PDV can be calculated.
 */
    class PDA
    {
    private:
        /** Number of Delta symbols */
        const size_t m_nNumberOfDeltas;

        /** Order of Taylor series expansion */
        size_t m_nOrder;

        /** Number of PDV instances */
        size_t m_nNumberOfPDVInstances;

        /**
         * Helper array of binomeal coefficients
         * n over k is m_aBinCoeffs[n][k];
         */
        size_t m_nMaxBinCoeffN;
        size_t m_nMaxBinCoeffK;
        std::vector<std::vector<size_t>> m_aBinCoeffs;

        /** Number of coefficients of each PDV value */
        size_t m_nNumberOfCoeffs;

        /** Delta distributions */
        std::vector<std::unique_ptr<Util::DeltaDistribution>> m_deltaDistributions;

        /** Delta names */
        std::vector<std::string> m_deltaNames;

    public:
        explicit PDA(size_t order=0, size_t numberOfDeltas=0);
        ~PDA();

        size_t getNumberOfDeltas() const { return m_nNumberOfDeltas; }
        size_t getOrder() const { return m_nOrder; }
        size_t getNumberOfCoeffs() const { return m_nNumberOfCoeffs; }
        inline size_t calcCoeffPos(const std::vector<size_t>& aPowers) const;
        size_t getBinCoeff(size_t n, size_t k) const { assert(n <= m_nMaxBinCoeffN && k <= m_nMaxBinCoeffK); return m_aBinCoeffs[n][k]; }
        template<class Distribution>
        void setDeltaDistribution(size_t nDelta, const Distribution& distribution) {
            this->m_deltaDistributions[nDelta] = std::make_unique<Distribution>(distribution);
        }
        Util::DeltaDistribution* getDeltaDistribution(size_t deltaNumber) { return m_deltaDistributions[deltaNumber].get(); }
        void setDeltaAsNormal(size_t nDelta, pdaValueType variance);
        pdaValueType getDeltaMoment(size_t nDelta, size_t nMoment);
        void setDeltaName(size_t nDelta, const std::string& name);
        std::string getDeltaName(size_t nDelta);

        friend PDV;
        friend Util::PowersIterator;
        friend std::vector<pdaValueType> Pda::Util::monteCarlo(
                const std::vector<PDV*>& aArguments,
                PDV (*function)(PDA& pda),
                size_t nSampleNumber);
    };

/**
 * Calculates the position inside the m_aCoeff array
 * depending on the Delta symbol powers
 * @param aPowers Array of Delta symbol powers
 * @return Corresponding index inside m_aCoeff
 */
    size_t PDA::calcCoeffPos(const std::vector<size_t>& aPowers) const {
        size_t nPos = 0;
        size_t l = this->getOrder();
        for (size_t nDelta = this->getNumberOfDeltas(); nDelta > 0; --nDelta) {
            size_t actualNDelta = nDelta - 1;
            for (size_t m = 0; m < aPowers[actualNDelta]; ++m){
                nPos += this->getBinCoeff(actualNDelta + l, l);
                --l;
            }
        }
        return nPos;
    }

    namespace Util {
        static std::random_device PDA_randomDevice;
        static std::mt19937 PDA_randomGen(Util::PDA_randomDevice());

/**
 * Base class for BRV distributions
 */
        class DeltaDistribution {
        protected:
            std::vector<pdaValueType> moments;
        public:
            /** Draw a random number according to this distribution */
            virtual pdaValueType drawSample() = 0;
            /** Raw moments */
            pdaValueType getMoment(size_t nOrder) const { assert(nOrder <= moments.size()); return moments[nOrder]; };
            /** Virtual destructor for base class */
            virtual ~DeltaDistribution() = default;
        };


/**
 * Normal distribution
 * \f[\mbox{PDF}=\frac{1}{\sigma \sqrt{2\pi}}e^{-\frac{1}{2}(\frac{x-\mu}{\sigma})^2}\f]
 * \f[\mu = \mbox{m}, \sigma = \mbox{s}\f]
 */
        class NormalDistribution: public DeltaDistribution {
        private:
            pdaValueType m = 0; // mu
            pdaValueType s = 1; // sigma
            std::normal_distribution<> dist;
        public:
            pdaValueType drawSample() override { return dist(Util::PDA_randomGen); };
            std::vector<pdaValueType> calculateMoments() const;
            explicit NormalDistribution(pdaValueType s): s(s), dist(m, s) {
                moments = calculateMoments();  };
        };

/**
 * Log normal distribution
 * Let \f$Z\f$ be standard normal distributed.
 * Then the distribution of \f$X=\exp^{\mu + \sigma Z}\f$ is log normal distributed.
 * The parameters of the distribution are \f$\mbox{m} = \mu\f$ and \f$\mbox{s} = \sigma\f$.
 * Since the expected value of a BRV distribution is expected to be zero,
 * the log normal distribution is moved by the offset \f$\exp^{\mu + \frac{\sigma^2}{2}}\f$ to the left.
 */
        class LogNormalDistribution: public DeltaDistribution {
        private:
            pdaValueType m = 0; // mu
            pdaValueType s = 1; // sigma
            pdaValueType offset;
            std::lognormal_distribution<> dist;
        public:
            pdaValueType drawSample() override { return dist(Util::PDA_randomGen)-offset; };
            std::vector<pdaValueType> calculateMoments() const;
            LogNormalDistribution(pdaValueType m, pdaValueType s): m(m), s(s), offset(::exp(m+s*s/2)), dist(m, s) { moments = calculateMoments(); };
            pdaValueType getOffset() const { return offset; };
        };

    } // end of namespace Util
} // end of namespace Pda

#endif // PDA_PDA_H