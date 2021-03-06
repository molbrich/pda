// pda_pdv.h
// Probability DeltaDistribution Arithmetic Variable
// (c) Markus Olbrich

#ifndef PDA_PDV_H
#define PDA_PDV_H

#include "pda_pda.h"
#include <iosfwd>
#include <vector>
#include <iostream>
#include <cassert>

namespace Pda {

/**
 * Type of moment calculation
 * Full: The following optional parameter is the maximal moment order that will be used
 * MonteCarlo: The following optional parameter is the number of samples
 * MonteCarloTime: The following optional parameter ist the time in clock ticks (seconds*CLOCKS_PER_SEC)
 */
    enum class MomentMethod {
        Auto,
        Full,
        PowerLimit,
        MonteCarloSamples,
        MonteCarloTime
    };

/**
 * @class PDV
 * @brief Probability DeltaDistribution Variable
 * @author Markus Olbrich
 */
    class PDV {

    private:
        /** PDA, this variable is part of */
        PDA& m_pda;

        // private object members:
        /** Array of coefficients before powers of Delta symbols */
        std::vector<pdaValueType> m_aCoeff;

        PDV unaryOperation(const std::vector<pdaValueType>& aDerivatives) const;
        pdaValueType getRawMomentUnscaled(size_t nOrder, MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;

    public:
        // Constructors:
        PDV() = delete; // Default constructor not allowed
        PDV(PDA& pda);
        PDV(PDA& pda, pdaValueType nomValue);
        PDV(const PDV & P);

        // Destructor
        ~PDV();

        // Compare two values for similarity
        static bool similar(pdaValueType, pdaValueType, pdaValueType = 1e-10);
        bool similar(const PDV&, pdaValueType = 1e-14) const;
        friend bool similar(const PDV&, const PDV&);

        PDA& getPDA() const { return m_pda; }

        // Operators:
        bool operator<(const PDV&) const;
        bool operator<=(const PDV&) const;
        bool operator>(const PDV&) const;
        bool operator>=(const PDV&) const;
        bool operator==(const PDV&) const;
        bool operator!=(const PDV&) const;

        PDV& operator=(const PDV&);
        PDV& operator=(pdaValueType d);

        friend PDV operator+(const PDV&, const PDV&);
        friend PDV operator+(PDV&&, const PDV&);
        friend PDV operator+(const PDV&, PDV&&);
        friend PDV operator+(PDV&&, PDV&&);

        friend PDV operator+(const PDV&, pdaValueType);
        friend PDV operator+(PDV&&, pdaValueType);
        friend PDV operator+(pdaValueType, const PDV&);
        friend PDV operator+(pdaValueType, PDV&&);

        PDV& operator+=(const PDV& );
        PDV& operator+=(pdaValueType d);

        PDV operator-() const &;
        PDV operator-() &&;

        friend PDV operator-(const PDV&, const PDV&);
        friend PDV operator-(PDV&&, const PDV&);
        friend PDV operator-(const PDV&, PDV&&);
        friend PDV operator-(PDV&&, PDV&&);

        friend PDV operator-(const PDV&, pdaValueType);
        friend PDV operator-(PDV&&, pdaValueType);
        friend PDV operator-(pdaValueType, const PDV&);
        friend PDV operator-(pdaValueType, PDV&&);

        PDV& operator-=(const PDV&);
        PDV& operator-=(pdaValueType);

        friend PDV operator*(const PDV&, const PDV&);

        friend PDV operator*(const PDV&, pdaValueType);
        friend PDV operator*(PDV&&, pdaValueType);
        friend PDV operator*(pdaValueType, const PDV&);
        friend PDV operator*(pdaValueType, PDV&&);
        PDV& operator*=(const PDV&);
        PDV& operator*=(pdaValueType);

        friend PDV operator/(const PDV&, const PDV&);
        PDV operator/(pdaValueType d) const;
        friend PDV operator /(pdaValueType, const PDV&);
        PDV& operator/=(const PDV&);
        PDV& operator/=(pdaValueType);

        PDV operator^(const PDV&) const;
        PDV operator^(pdaValueType) const;
        friend PDV operator^(pdaValueType, const PDV&);

        // Math functions:
        friend PDV log(const PDV&);
        friend PDV exp(const PDV&);
        friend PDV sqrt(const PDV&);
        friend PDV isqrt(const PDV&);
        friend PDV logexp(const PDV&);
        friend PDV inv(const PDV&);
        friend PDV sin(const PDV&);
        friend PDV asin(const PDV&);
        friend PDV cos(const PDV&);
        friend PDV acos(const PDV&);
        friend PDV tan(const PDV&);
        friend PDV atan(const PDV&);
        friend PDV tanh(const PDV&);
        friend PDV pow(const PDV&,
                       const PDV&);
        friend PDV pow(const PDV&, pdaValueType);

        // Set/Get values:
        void setNom(pdaValueType value);
        pdaValueType getNom() const;
        void setDeltaCoeff(size_t nDeltaNumber,
                           pdaValueType factor);
        pdaValueType getDeltaCoeff(size_t nDeltaNumber) const;
        void setCoeff(const std::vector<size_t>& aPowers,
                      pdaValueType value);
        pdaValueType getCoeff(const std::vector<size_t>& aPowers) const;

        // Calculate moments:
        pdaValueType getMean(MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        pdaValueType getStandardDeviation(MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        pdaValueType getVariance(MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        pdaValueType getSkewness(MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        pdaValueType getExcessKurtosis(MomentMethod method= MomentMethod::Auto, size_t nMax= 0) const;
        friend pdaValueType E(const PDV& value, MomentMethod method, size_t nMax);
        friend pdaValueType Var(const PDV& value, MomentMethod method, size_t nMax);
        pdaValueType getRawMoment(size_t nOrder, MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        std::vector<pdaValueType> getRawMoments(size_t nMaxOrder, MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        pdaValueType getCentralMoment(size_t nOrder, MomentMethod method=MomentMethod::Auto, size_t nMax=0) const;
        std::vector<pdaValueType> getCentralMoments(size_t nMaxOrder, MomentMethod method= MomentMethod::Auto, size_t nMax= 0) const;

        // Sensitivities:
        pdaValueType getSensitivity(size_t nDeltaNumber);

        // Monte Carlo estimation:
        std::vector<pdaValueType> drawSamples(size_t sampleCount = 100);
        std::vector<std::pair<pdaValueType, pdaValueType>> estimateCDF(size_t pointCount = 50, size_t averagingCount = 500);
        std::vector<std::pair<pdaValueType, pdaValueType>> estimatePDF(size_t pointCount = 50, size_t averagingCount = 2000);

        // Calculate covariance:
        friend pdaValueType Cov(const PDV& value1, const PDV& value2, MomentMethod method, size_t nMax);
        friend pdaValueType Cor(const PDV& value1, const PDV& value2, MomentMethod method, size_t nMax);

        // Output:
        friend std::ostream& operator<<(std::ostream& os,
                                        const PDV& value);
        void dump(std::ostream &os) const;

        // Cleanup:
        void clean();

        // Debug:
        void check();

        friend void swap(PDV& a, PDV& b) { a.m_aCoeff.swap(b.m_aCoeff); };
    };

    pdaValueType E(const PDV& value, MomentMethod method=MomentMethod::Auto, size_t nMax=0);
    pdaValueType Var(const PDV& value, MomentMethod method=MomentMethod::Auto, size_t nMax=0);
    pdaValueType Cov(const PDV& value1, const PDV& value2, MomentMethod method=MomentMethod::Auto, size_t nMax=0);
    pdaValueType Cor(const PDV& value1, const PDV& value2, MomentMethod method=MomentMethod::Auto, size_t nMax=0);

} // end of namespace Pda

#endif // PDA_PDV_H
