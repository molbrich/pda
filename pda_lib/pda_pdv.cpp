// pda_pdv.cpp
// Probability DeltaDistribution Arithmetic
// Implementation of classes PowersIterator and PDV
// (c) Markus Olbrich

#include "pda_pda.h"
#include "pda_pdv.h"
#include "pda_powersiterator.h"
#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <cmath>

namespace Pda {

//#define USE_CALLBACK_ITERATOR

//----------------------------------------------------------------
// private methods:
//----------------------------------------------------------------

/** 
 * Performs a unary numerical operation on (*this). 
 * \f[\sum\limits_{i=0}^{l} \frac{1}{i!}\frac{\partial^{i}}{\partial X^i}T(x_0)
 (X-x_0)^i\f]
 * Sets the coefficients according to given derivatives
 * using a Taylor series.
 * @param aDerivatives First four derivatives of the 
 *        operation.
 * @returns The result of the operation
 * @note Deprecated
 */
    PDV PDV::unaryOperation(const std::vector<pdaValueType>& aDerivatives) const {
        PDV result(m_pda, aDerivatives[0]);
        pdaValueType dProduct;
        // Loop over Taylor series terms:
        pdaValueType dDivisor = 1;
        for (size_t nTerm = 1; nTerm <= m_pda.getOrder(); ++nTerm){
            if (aDerivatives.size() <= m_pda.getOrder()) continue;
            dDivisor *= static_cast<pdaValueType>(nTerm);
            Util::PowersIterator pi(m_pda, nTerm, m_pda.getOrder(), m_pda.getOrder());
#ifdef USE_CALLBACK_ITERATOR
            pi.iterate([&]()->void {
                dProduct = 1;
                std::vector<size_t>& aPositions = pi.getPositions();
                for (size_t nFactor = 0; nFactor < nTerm; ++nFactor)
                    if (aPositions[nFactor] == 0) {
                        dProduct = 0;
                        break;
                    }
                    else dProduct *= m_aCoeff[aPositions[nFactor]];
                result.m_aCoeff[pi.getPosition()] += aDerivatives[nTerm] / dDivisor * dProduct;
            });
#else
            do {
                dProduct = 1;
                std::vector<size_t>& aPositions = pi.getPositions();
                for (size_t nFactor = 0; nFactor < nTerm; ++nFactor)
                    if (aPositions[nFactor] == 0) {
                        dProduct = 0;
                        break;
                    }
                    else dProduct *= m_aCoeff[aPositions[nFactor]];
                result.m_aCoeff[pi.getPosition()] += aDerivatives[nTerm] / dDivisor * dProduct;
            } while (pi.next());
#endif
        }
        return result;

#if false
        if (m_pda.getOrder() == 3){
        PDV dX = *this;
        dX.m_aCoeff[0] = 0;
        PDV dX2 = dX * dX;
        PDV dX3 = dX * dX2;
        return aDerivatives[0] + aDerivatives[1] * dX + (aDerivatives[2] / 2) * dX2 + (aDerivatives[3] / 6) * dX3;
    }

    if (m_pda.getOrder() == 4){
        PDV dX = *this;
        dX.m_aCoeff[0] = 0;
        PDV dX2 = dX * dX;
        PDV dX3 = dX * dX2;
        PDV dX4 = dX2 * dX2;
        return aDerivatives[0] + aDerivatives[1] * dX + (aDerivatives[2] / 2) * dX2 + (aDerivatives[3] / 6) * dX3 + (aDerivatives[4] / 24) * dX4;
    }

    PDV result(m_pda, aDerivatives[0]);
    return result;
#endif
    }

//----------------------------------------------------------------
// public methods:
//----------------------------------------------------------------

// Constructors:
/**
 * Simple Constructor.
 */
    PDV::PDV(PDA& pda) :
            m_pda(pda),
            m_aCoeff(pda.getNumberOfCoeffs(), 0.0) {
        ++pda.m_nNumberOfPDVInstances;
    }

/** 
 * Constructor that sets the nominal value.
 * @param nomValue Nominal value
 */
    PDV::PDV(PDA& pda, pdaValueType nomValue) :
            m_pda(pda),
            m_aCoeff(pda.getNumberOfCoeffs(), 0.0) {
        m_aCoeff[0] = nomValue;
        ++pda.m_nNumberOfPDVInstances;
    }

/** 
 * Copy constructor.
 */
    PDV::PDV(const PDV& x) : m_pda(x.m_pda),  m_aCoeff(x.m_aCoeff) {
        ++m_pda.m_nNumberOfPDVInstances;
    }

/** 
 * Destructor.
 */
    PDV::~PDV() {
        --m_pda.m_nNumberOfPDVInstances;
    }

/**
 * Similarity test.
 * @param value, value2
 * @returns True if values are nearly equal
 */
    bool PDV::similar(pdaValueType value1, pdaValueType value2, pdaValueType epsilon) {
        pdaValueType diff = fabs(value1 - value2);
        if (diff < 1e-20)
            return true;
        pdaValueType sum = fabs(value1 + value2);
        if (sum < 1e-20)
            return true;
        if(diff / sum < epsilon)
            return true;
        return false;
    }

/**
 * Similarity test
 * @param epsilon Relative threashold value
 * @returns True if all coefficients are equal
 */
    bool PDV::similar(const PDV& x, const pdaValueType epsilon) const {
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            if (!PDV::similar(m_aCoeff[i], x.m_aCoeff[i])) {
                if (fabs((m_aCoeff[i] - x.m_aCoeff[i]) / (m_aCoeff[i] + x.m_aCoeff[i])) < epsilon)
                    continue;
                if (!PDV::similar(m_aCoeff[0], 0) && i != 0)
                    if (fabs(m_aCoeff[i] - x.m_aCoeff[i]) / m_aCoeff[0] < epsilon)
                        continue;
                return false;
            }
        return true;
    }
    bool similar(const PDV& x, const PDV& y) {
        return x.similar(y);
    }

/** 
 * Comparison
 * @param x Right value
 * @returns True if *this.getNom()<x.getNom()
 */
    bool PDV::operator<(const PDV& x) const {
        return m_aCoeff[0] < x.m_aCoeff[0];
    }

/** 
 * Comparison
 * @param x Right value
 * @returns True if *this.getNom()<=x.getNom()
 */
    bool PDV::operator<=(const PDV& x) const {
        return m_aCoeff[0] <= x.m_aCoeff[0];
    }

/** 
 * Comparison
 * @param x Right value
 * @returns True if *this.getNom()>x.getNom()
 */
    bool PDV::operator>(const PDV& x) const {
        return m_aCoeff[0] > x.m_aCoeff[0];
    }

/** 
 * Comparison
 * @param x Right value
 * @returns True if *this.getNom()>=x.getNom()
 */
    bool PDV::operator>=(const PDV& x) const {
        return m_aCoeff[0] >= x.m_aCoeff[0];
    }

/** 
 * Equality test
 * @param x Right value
 * @returns True if all coefficients are equal
 */
    bool PDV::operator==(const PDV &P) const {
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            if (!similar(m_aCoeff[i], P.m_aCoeff[i]))
                return false;
        return true;
    }

/** 
 * Not equal test
 * @param x Right value
 * @returns True if any coefficient is not equal
 */
    bool PDV::operator!=(const PDV& x) const {
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            if (!similar(m_aCoeff[i], x.m_aCoeff[i])) {
                return true;
            }
        return false;
    }

/** 
 * Assignment operator. [PDV] = [PDV]
 * @param x Right value
 * @returns *this
 */
    PDV& PDV::operator=(const PDV& x) {
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            m_aCoeff[i] = x.m_aCoeff[i];
        return *this;
    }

/** 
 * Assignment operator. [PDV] = [pdaValueType]
 * @param d Right value
 * @returns *this
 */
    PDV& PDV::operator=(const pdaValueType d) {
        m_aCoeff[0] = d;
        for (size_t i = 1; i < m_pda.getNumberOfCoeffs(); ++i)
            m_aCoeff[i] = 0;
        return *this;
    }

/** 
 * Operator [PDV] + [PDV]
 * @param x&
 * @param y&
 * @returns x + y
 */
    PDV operator+(const PDV& x, const PDV& y) {
        PDV temp(x.m_pda);
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] = x.m_aCoeff[i] + y.m_aCoeff[i];
        return temp;
    }
    PDV operator+(PDV&& x, const PDV& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            x.m_aCoeff[i] += y.m_aCoeff[i];
        return x;
    }
    PDV operator+(const PDV& x, PDV&& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            y.m_aCoeff[i] += x.m_aCoeff[i];
        return y;
    }
    PDV operator+(PDV&& x, PDV&& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            x.m_aCoeff[i] += y.m_aCoeff[i];
        return x;
    }

/** 
 * Operator [PDV] + [pdaValueType]
 * @param x Left value
 * @param d Right value
 * @returns x + d
 */
    PDV operator+(const PDV& x, const pdaValueType d) {
        PDV temp(x);
        temp.m_aCoeff[0] += d;
        return temp;
    }
    PDV operator+(PDV&& x, const pdaValueType d) {
        x.m_aCoeff[0] += d;
        return x;
    }

/** 
 * Operator [pdaValueType] + [PDV]
 * @param d Left value
 * @param x Right value
 * @returns d + P
 */
    PDV operator+(const pdaValueType d, const PDV& x) {
        PDV temp(x);
        temp.m_aCoeff[0] += d;
        return temp;
    }
    PDV operator+(const pdaValueType d, PDV&& x) {
        x.m_aCoeff[0] += d;
        return x;
    }

/** 
 * Calculating assignment operator. [PDV] += [PDV]
 * @param P Right value
 * @returns *this
 */
    PDV& PDV::operator+=(const PDV& P) {
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            this->m_aCoeff[i] += P.m_aCoeff[i];
        return *this;
    }

/** 
 * Calculating assignment operator. [PDV] += [pdaValueType]
 * @param d Right value
 * @returns *this
 */
    PDV& PDV::operator+=(const pdaValueType d) {
        m_aCoeff[0] += d;
        return *this;
    }

/** 
 * Operator - [PDV&]
 * @returns -(*this)
 */
    PDV PDV::operator-() const & {
        PDV temp{m_pda};
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] = -m_aCoeff[i];
        return temp;
    }
/**
 * Operator - [PDV&&]
 * @returns -(*this)
 */
    PDV PDV::operator-() && {
        assert(m_aCoeff.size() == m_pda.getNumberOfCoeffs());
        for (auto& coeff: m_aCoeff)
            coeff = -coeff;
        return *this;
    }

/** 
 * Operator [PDV] - [PDV]
 * @param x Lef value
 * @param y Right value
 * @returns x - y
 */
    PDV operator-(const PDV& x, const PDV& y) {
        PDV temp(x.m_pda);
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] = x.m_aCoeff[i] - y.m_aCoeff[i];
        return temp;
    }
    PDV operator-(PDV&& x, const PDV& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            x.m_aCoeff[i] -= y.m_aCoeff[i];
        return x;
    }
    PDV operator-(const PDV& x, PDV&& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            y.m_aCoeff[i] = x.m_aCoeff[i] - y.m_aCoeff[i];
        return y;
    }
    PDV operator-(PDV&& x, PDV&& y) {
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            x.m_aCoeff[i] -= y.m_aCoeff[i];
        return x;
    }

/** 
 * Operator [PDV] - [pdaValueType]
 */
    PDV operator-(const PDV& x, const pdaValueType d) {
        PDV temp(x);
        temp.m_aCoeff[0] -= d;
        return temp;
    }
    PDV operator-(PDV&& x, const pdaValueType d) {
        x.m_aCoeff[0] -= d;
        return x;
    }

/** 
 * Operator [pdaValueType] - [PDV]
 */
    PDV operator-(const pdaValueType d, const PDV& x) {
        PDV temp(-x);
        temp.m_aCoeff[0] += d;
        return temp;
    }
/**
 * Operator [pdaValueType] - [PDV]
 * */
    PDV operator-(const pdaValueType d, PDV&& x) {
        x=-x;
        x.m_aCoeff[0] += d;
        return x;
    }

/** 
 * Calculating assignment operator. [PDV] -= [PDV]
 */
    PDV & PDV::operator-=(const PDV & P){
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            m_aCoeff[i] -= P.m_aCoeff[i];
        return *this;
    }

/** 
 * Calculating assignment operator. [PDV] -= [pdaValueType]
 */
    PDV & PDV::operator-=(const pdaValueType d){
        m_aCoeff[0] -= d;
        return *this;
    }

/** 
 * Operator [PDV] * [PDV]
 */
    PDV operator*(const PDV& x, const PDV& y) {
        assert(&x.m_pda == &y.m_pda);
        PDV result(x.m_pda, 0);
#ifdef USE_CALLBACK_ITERATOR
        PowersIterator pi(x.m_pda, 2, x.m_pda.getOrder()) ;
        pi.iterate([&]()->void{
            result.m_aCoeff[pi.getPosition()] += x.m_aCoeff[pi.getPositions()[0]] * y.m_aCoeff[pi.getPositions()[1]];
        });
#else
        Util::PowersIterator pi(x.m_pda, 2, x.m_pda.getOrder()) ;
        do {
            result.m_aCoeff[pi.getPosition()] += x.m_aCoeff[pi.getPositions()[0]] * y.m_aCoeff[pi.getPositions()[1]];
        } while (pi.next());
#endif
        return result;
    }

/** 
 * Operator [PDV] * [pdaValueType]
 */
    PDV operator*(const PDV& x, pdaValueType d) {
        PDV temp{x.m_pda};
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] = x.m_aCoeff[i] * d;
        return temp;
    }
/**
 * Operator [PDV] * [pdaValueType]
 */
    PDV operator*(PDV&& x, pdaValueType d) {
        assert(x.m_aCoeff.size() == x.m_pda.getNumberOfCoeffs());
        for (auto& coeff: x.m_aCoeff)
            coeff *= d;
        return x;
    }

/** 
 * Operator [pdaValueType] * [PDV]
 */
    PDV operator *(const pdaValueType d, const PDV& x) {
        PDV temp{x.getPDA()};
        for (size_t i = 0; i < x.m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] = d * x.m_aCoeff[i];
        return temp;
    }
    PDV operator *(const pdaValueType d, PDV&& x) {
        assert(x.m_aCoeff.size() == x.m_pda.getNumberOfCoeffs());
        for (auto& coeff: x.m_aCoeff)
            coeff *= d;
        return x;
    }

/** 
 * Calculating assignment operator. [PDV] *= [PDV]
 */
    PDV& PDV::operator*=(const PDV& x) {
        *this = *this * x;
        return *this;
    }

/** 
 * Calculating assignment operator. [PDV] *= [pdaValueType]
 */
    PDV& PDV::operator*=(const pdaValueType d) {
        assert (m_aCoeff.size() == m_pda.getNumberOfCoeffs());
        for (auto& coeff: m_aCoeff)
            coeff *= d;
        return *this;
    }

/** 
 * Operator [PDV] / [PDV]
 */
    PDV operator/(const PDV& x, const PDV& y) {
        return x * inv(y);
    }

/** 
 * Operator [PDV] / [pdaValueType]
 */
    PDV PDV::operator/(pdaValueType d) const {
        PDV temp(*this);
        for (size_t i = 0; i < m_pda.getNumberOfCoeffs(); ++i)
            temp.m_aCoeff[i] /= d;
        return temp;
    }

/** 
 * Operator [pdaValueType] / [PDV]
 */
    PDV operator/(const pdaValueType d, const PDV& x){
        return d * inv(x);
    }

/** 
 * Calculating assignment operator. [PDV] /= [PDV]
 */
    PDV & PDV::operator/=(const PDV & P) {
        *this = *this / P;
        return *this;
    }

/** 
 * Calculating assignment operator. [PDV] /= [pdaValueType]
 * @param d Right value
 * @returns *this
 */
    PDV & PDV::operator/=(const pdaValueType d) {
        for (auto& coeff: m_aCoeff)
            coeff /= d;
        return *this;
    }

/** 
 * Operator [PDV] ^ [PDV]
 * @param P Right value
 * @returns (*this) ^ P
 */
    PDV PDV::operator^(const PDV & P) const {
        return exp(P * log(*this));
#if false
        pdaValueType aDerivatives[5][5];
     pdaValueType x = getNom();
     pdaValueType y = P.getNom();
     pdaValueType xpy = pow(x,y);
     pdaValueType lnx = log(x);
     aDerivatives[0][0] = xpy;
     aDerivatives[1][0] = xpy*x/y;
     aDerivatives[0][1] = xpy*lnx;
     aDerivatives[1][1] = (y*lnx+1)*xpy/x;
     if (n_mOrder > 1) {
     aDerivatives[2][0] = (y-1)*y*xpy/(x*x);
     aDerivatives[2][1] = ((y-1)*y*lnx+2*y-1)*xpy/(x*x);
     aDerivatives[0][2] = xpy*lnx*lnx;
     aDerivatives[1][2] = (y*lnx+2)*xpy*lnx/x;
     aDerivatives[2][2] = (((y*lnx-lnx+4)*y-2)*lnx+2)*xpy/(x*x);
     }
     if (n_mOrder > 2) {
     aDerivatives[3][0] = ((y-3)*y-2)*y*xpy/(x*x*x);
     aDerivatives[3][1] = xpy/(x*x*x);
     aDerivatives[3][2] = xpy/(x*x*x);
     aDerivatives[0][3] = xpy/(x*x*x);
     aDerivatives[1][3] = 1;
     aDerivatives[2][3] = xpy/(x*x*x);
     aDerivatives[3][3] = xpy/(x*x*x);
     }
     return binaryOperation(P, aDerivatives);
#endif
    }

/** 
 * Operator [PDV] ^ [pdaValueType]
 * @param d Right value
 * @returns *this ^ d
 */
    PDV PDV::operator ^(pdaValueType d) const {
        return exp(d * log(*this));
    }

/** 
 * Operator [pdaValueType] ^ [PDV]
 * @param d Left value
 * @param P Right value
 * @returns d ^ P
 */
    PDV operator ^(const pdaValueType d,
                   const PDV & P){
        return exp(P * log(d));
    }

/**
 * Natural logarithm.
 * @param x Argument
 * @returns The natural logarithm of P.
 */
    PDV log(const PDV& x) {
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        assert (x.getNom() > 0.0);
#if true
        pdaValueType xNom = x.getNom();
        aDerivatives[0] = std::log(xNom);
        if (x.m_pda.getOrder() >= 1) {
            aDerivatives[1] = 1.0 / xNom;
            if (x.m_pda.getOrder() >= 2) {
                pdaValueType xNom2 = xNom * xNom;
                aDerivatives[2] = -1.0 / xNom2;
                if (x.m_pda.getOrder() >= 3) {
                    pdaValueType xNom3 = xNom2 * xNom;
                    aDerivatives[3] = 2.0 / xNom3;
                    if (x.m_pda.getOrder() == 4) {
                        pdaValueType xNom4 = xNom3 * xNom;
                        aDerivatives[4] = -6.0 / xNom4;
                    }
                }
            }
        }
#else
        pdaValueType xNom = x.getNom();
        pdaValueType xNom2 = xNom * xNom;
        pdaValueType xNom3 = xNom2 * xNom;
        pdaValueType xNom4 = xNom3 * xNom;
        aDerivatives[0] = std::log(xNom);
        aDerivatives[1] = 1.0 / xNom;
        aDerivatives[2] = -1.0 / xNom2;
        aDerivatives[3] = 2.0 / xNom3;
        aDerivatives[4] = -6.0 / xNom4;
#endif
        return x.unaryOperation(aDerivatives);
    }

/**
 * Natural exponential.
 * @param x Argument
 * @returns e to the power of x.
 */
    PDV exp(const PDV& x) {
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        aDerivatives[4] = aDerivatives[3] = aDerivatives[2] = aDerivatives[1] = aDerivatives[0] = std::exp(x.getNom());
        return x.unaryOperation(aDerivatives);
    }

/**
 * Square root
 * @param x Argument
 * @returns Square root of x.
 */
    PDV sqrt(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType xNom = x.getNom();
        assert (xNom >= 0);
        pdaValueType s = std::sqrt(xNom);
        aDerivatives[0] = s;
#if true
        if (x.m_pda.getOrder() >= 1) {
            aDerivatives[1] = 1.0 / (2.0 * s);
            if (x.m_pda.getOrder() >= 2) {
                aDerivatives[2] = -1.0 / (4.0 * xNom * s);
                if (x.m_pda.getOrder() >= 3) {
                    pdaValueType xNom2 = xNom * xNom;
                    aDerivatives[3] = 3.0 / (8.0 * xNom2 * s);
                    if (x.m_pda.getOrder() == 4) {
                        pdaValueType xNom3 = xNom2 * xNom;
                        aDerivatives[4] = -15.0 / (16.0 * xNom3 * s);
                    }
                }
            }
        }
#else
        aDerivatives[1] = 1.0 / (2.0 * s);
        aDerivatives[2] = -1.0 / (4.0 * xNom * s);
        pdaValueType xNom2 = xNom * xNom;
        aDerivatives[3] = 3.0 / (8.0 * xNom2 * s);
        pdaValueType xNom3 = xNom2 * xNom;
        aDerivatives[4] = -15.0 / (16.0 * xNom3 * s);
#endif
        return x.unaryOperation(aDerivatives);
    }

/**
 * Inverse of square root
 * @param x Argument
 * @returns Square root of x.
 */
    PDV isqrt(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType xNom = x.getNom();
        pdaValueType xNom2 = xNom * xNom;
        pdaValueType xNom3 = xNom2 * xNom;
        pdaValueType xNom4 = xNom3 * xNom;
        assert (xNom >= 0);
        pdaValueType s = std::sqrt(xNom);
        aDerivatives[0] = 1 / s;
        aDerivatives[1] = -1 / (2 * xNom * s);
        aDerivatives[2] = 3 / (4 * xNom2 * s);
        aDerivatives[3] = -15 / (8 * xNom3 * s);
        aDerivatives[4] = 105 / (16 * xNom4 * s);
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes log(1.0 + exp(x))
 * @param x Argument
 * @returns Result.
 */
    PDV logexp(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType ex = ::exp(x.getNom());
        pdaValueType fx = ex / (1.0 + ex);
        pdaValueType fx2 = fx * fx;
        pdaValueType fx3 = fx2 * fx;
        pdaValueType fx4 = fx3 * fx;
        aDerivatives[0] = std::log(1.0 + ex);
        aDerivatives[1] = fx;
        aDerivatives[2] = fx - fx2;
        aDerivatives[3] = fx - 3 * fx2 + 2 * fx3;
        aDerivatives[4] = fx - 7 * fx2 + 12 * fx3 - 6 * fx4;
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes arcus tangens
 * @param x Argument
 * @returns Result.
 */
    PDV atan(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType xNom = x.getNom();
        pdaValueType d = 1 + xNom * xNom;
        pdaValueType d2 = d * d;
        pdaValueType d3 = d2 * d;
        pdaValueType d4 = d3 * d;
        aDerivatives[0] = std::atan(xNom);
        aDerivatives[1] = 1 / d;
        aDerivatives[2] = -2 * xNom / d2;
        aDerivatives[3] = 8 * xNom * xNom / d3 - 2 / d2;
        aDerivatives[4] = -48 * xNom * xNom * xNom / d4 + 24 * xNom / d3;
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes tangens hyperbolicus
 * @param x Argument
 * @returns Result.
 */
    PDV tanh(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType t = std::tanh(x.getNom());
        pdaValueType t2 = t * t;
        pdaValueType t3 = t2 * t;
        pdaValueType et = 1 - t2;
        pdaValueType et2 = et * et;
        aDerivatives[0] = t;
        aDerivatives[1] = 1 - t2;
        aDerivatives[2] = -2 * t * et;
        aDerivatives[3] = -2 * et2 + 4 * t2 * et;
        aDerivatives[4] = 16 * et2 * t - 8 * t3 * et;
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes the inverse
 * @param x Argument
 * @returns 1/x
 */
    PDV inv(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType xNom = x.getNom();
        pdaValueType xNom2 = xNom * xNom;
        pdaValueType xNom3 = xNom2 * xNom;
        pdaValueType xNom4 = xNom3 * xNom;
        pdaValueType xNom5 = xNom4 * xNom;
        aDerivatives[0] = 1 / xNom;
        aDerivatives[1] = -1 / xNom2;
        aDerivatives[2] = 2 / xNom3;
        aDerivatives[3] = -6 / xNom4;
        aDerivatives[4] = 24 / xNom5;
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes sine
 * @param x Argument
 * @returns sin(x).
 */
    PDV sin(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType s = std::sin(x.getNom());
        pdaValueType c = std::cos(x.getNom());
        aDerivatives[0] = s;
        aDerivatives[1] = c;
        aDerivatives[2] = -s;
        aDerivatives[3] = -c;
        aDerivatives[4] = s;
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes cosine
 * @param x Argument
 * @returns cos(x).
 */
    PDV cos(const PDV& x){
        static std::vector<pdaValueType> aDerivatives(5); // The first four derivatives
        pdaValueType s = std::sin(x.getNom());
        pdaValueType c = std::cos(x.getNom());
        aDerivatives[0] = c;
        if (x.m_pda.getOrder() <= 1) {
            aDerivatives[1] = -s;
            if (x.m_pda.getOrder() <= 2) {
                aDerivatives[2] = -c;
                if (x.m_pda.getOrder() <= 3) {
                    aDerivatives[3] = s;
                    if (x.m_pda.getOrder() <= 4) {
                        aDerivatives[4] = c;
                    }
                }
            }
        }
        return x.unaryOperation(aDerivatives);
    }

/**
 * Computes power
 * @param x Base
 * @param y Exponent
 * @returns x ^ y.
 */
    PDV pow(const PDV& x, const PDV& y){
        return exp(y * log(x));
    }

/**
 * Computes power
 * @param x Base
 * @param y Exponent
 * @returns x ^ y.
 */
    PDV pow(const PDV& x, const pdaValueType y){
        return exp(y * log(x));
    }

/** 
 * Sets the nominal value. 
 * Leaves all other coefficients unchanged.
 * @param value New nominal value
 */
    void PDV::setNom(pdaValueType value){
        m_aCoeff[0] = value;
    }

/** 
 * Sets the factor (coeffecient) of a single Delta symbol. 
 * Leaves all other coefficients unchanged
 * @param nDeltaNumber Index of the Delta Symbol.
 * The index of the first symbol is 0.
 * @param factor New coefficient
 */
    void PDV::setDeltaCoeff(size_t nDeltaNumber,
                            pdaValueType factor){
        assert(nDeltaNumber < m_pda.getNumberOfDeltas());
        if (m_pda.getOrder() == 0) return;
        assert(nDeltaNumber < m_pda.getNumberOfDeltas());
        std::vector<size_t> aPowers(m_pda.getNumberOfDeltas());
        for (size_t i = 0; i < m_pda.getNumberOfDeltas(); ++i)
            aPowers[i] = 0;
        aPowers[nDeltaNumber] = 1;
        setCoeff(aPowers, factor);
    }

/** 
 * Determines the factor (coeffecient) of a single Delta symbol. 
 * @param nDeltaNumber Index of the Delta Symbol.
 * The index of the first symbol is 0.
 * @returns factor
 */
    pdaValueType PDV::getDeltaCoeff(size_t nDeltaNumber){
        assert(nDeltaNumber < m_pda.getNumberOfDeltas());
        std::vector<size_t> aPowers(m_pda.getNumberOfDeltas());
        for (size_t i = 0; i < m_pda.getNumberOfDeltas(); ++i)
            aPowers[i] = 0;
        aPowers[nDeltaNumber] = 1;
        pdaValueType r;
        r = getCoeff(aPowers);
        return r;
    }

/** 
 * Sets the coeffecient of an arbitrary Delta symbol power combination. 
 * Leaves all other coefficients unchanged
 * @param aPowers
 *        Array of m_nNumberOfDeltas powers
 * @param value New coefficient value
 */
    void PDV::setCoeff(const std::vector<size_t>& aPowers,
                       pdaValueType value) {
        m_aCoeff[m_pda.calcCoeffPos(aPowers)] = value;
    }

/** 
 * Gets the coeffecient of an arbitrary Delta symbol power combination. 
 * Leaves all other coefficients unchanged
 * @param aPowers
 *        Array of m_nNumberOfDeltas powers
 * @returns Coefficient
 */
    pdaValueType PDV::getCoeff(const std::vector<size_t>& aPowers) {
        return m_aCoeff[m_pda.calcCoeffPos(aPowers)];
    }

/** 
 * Determins the nominal value.
 * @returns Nominal value
 */
    pdaValueType PDV::getNom(void) const {
        return m_aCoeff[0];
    }

/** 
 * Calculates the expected (mean) value.
 * @returns  Mean value
 * @see E()
 */
    pdaValueType PDV::getMean(const MomentMethod method, const size_t nMax) const {
        return getRawMoment(1, method, nMax);
    }

/**
 * Calculates the variance.
 * @returns Variance = getCentralMoment(2)
 */
    pdaValueType PDV::getVariance(const MomentMethod method, const size_t nMax) const {
        return getCentralMoment(2, method, nMax);
    }

/**
 * Determins the moment coefficient of skewness.
 * @returns Skewness = getCentralMoment(3)/getCentralMoment(2)^(3/2)
 */
    pdaValueType PDV::getSkewness(const MomentMethod method, const size_t nMax) const {
        auto centralMoments = this->getCentralMoments(3, method, nMax);
        if (centralMoments[2] == 0)
            return 0;
        return centralMoments[3] / ::pow(centralMoments[2], 1.5);
    }

/**
 * Determins the kurtosis.
 * @returns Kurtosis = getCentralMoment(4)/getCentralMoment(2)^2-3
 */
    pdaValueType PDV::getExcessKurtosis(MomentMethod method, size_t nMax) const {
        auto centralMoments = this->getCentralMoments(4, method, nMax);
        if (centralMoments[2] == 0)
            return 0;
        return centralMoments[4] / ::pow(centralMoments[2], 2) - 3;
    }

/**
 * Determins the mean (expactancy) value.
 * Same as value.getMean().
 * @returns  Mean value
 * @see getMean()
 */
    pdaValueType E(const PDV& value, const MomentMethod method, const size_t nMax){
        return value.getMean(method, nMax);
    }

/**
 * Determins the variance.
 * Same as value.getVariance().
 * @returns Variance = getCentralMoment(2)
 * @see getVariance()
 */
    pdaValueType Var(const PDV& value, const MomentMethod method, const size_t nMax){
        return value.getVariance(method, nMax);
    }

/** 
 * Calculates the nOrder'th raw moment.
 * @param nOrder Order of the raw moment to calculate
 * @param method
 * @param nMax FullMoments method: Maximum sum of the orders of the used moments, MonteCarlo method: Sample number
 * @returns Value of the raw moment
 */
    pdaValueType PDV::getRawMoment(const size_t nOrder, const MomentMethod method, const size_t nMax) const {
        pdaValueType dMoment = 0;
        if(method == MomentMethod::Auto) {
            size_t nMaxTotalPowersSum = (nMax == 0 ? std::max(m_pda.getOrder() * nOrder, nOrder) : nMax);
            return getRawMoment(nOrder, MomentMethod::FullMoments, nMaxTotalPowersSum);
        }
        if(method == MomentMethod::FullMoments) {
            if (nOrder == 0)
                return 1;
            // Make sure no moment is calculated wrongly:
            assert(nOrder * m_pda.getOrder() <= nMax);
            Util::PowersIterator pi(m_pda, nOrder, nMax);
            do {
                pdaValueType dAddend = 1;
                std::vector<size_t>& aPositions = pi.getPositions();
                for (size_t nFactor = 0; nFactor < nOrder; ++nFactor) {
                    dAddend *= m_aCoeff[aPositions[nFactor]];
                }
                if (similar(dAddend, 0)) continue;
                std::vector<size_t>& aFactorsPowersSum = pi.getFactorsPowersSum();
                for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                    dAddend *= m_pda.getDeltaMoment(nDelta, aFactorsPowersSum[nDelta]);
                }
                dMoment += dAddend;
            } while (pi.next());
        }
        else {
            pdaValueType dMomentSum = 0;
            size_t nMaxDeltaPower = nOrder * m_pda.getOrder(); // Full power combinations
            size_t count = 0;

            auto processSample = [&] {
                std::vector<std::vector<pdaValueType>> powers(m_pda.getNumberOfDeltas(),
                                                              std::vector<pdaValueType>(m_pda.getOrder() + 1));
                for (size_t delta = 0; delta < m_pda.getNumberOfDeltas(); ++delta) {
                    pdaValueType deltaValue = m_pda.m_deltaDistributions[delta]->drawSample();
                    for (size_t order = 0; order <= m_pda.getOrder(); ++order)
                        powers[delta][order] = ::pow(deltaValue, static_cast<double>(order));
                }
                pdaValueType dValue = 0;
                pdaValueType dAddend;
                Util::PowersIterator pi(m_pda, 1, nMaxDeltaPower);
                do {
                    dAddend = m_aCoeff[pi.getPositions()[0]];
                    if (similar(dAddend, 0)) continue;
                    std::vector<size_t> &aFactorsPowersSum = pi.getFactorsPowersSum();
                    for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                        if (aFactorsPowersSum[nDelta] > nMaxDeltaPower) {
                            dAddend = 0;
                            break;
                        }
                        dAddend *= powers[nDelta][aFactorsPowersSum[nDelta]];
                    }
                    dValue += dAddend;
                } while (pi.next());
                switch (nOrder) {
                    case 0:
                        dMomentSum += 1;
                        break;
                    case 1:
                        dMomentSum += dValue;
                        break;
                    case 2:
                        dMomentSum += dValue * dValue;
                        break;
                    case 3:
                        dMomentSum += dValue * dValue * dValue;
                        break;
                    case 4:
                        dMomentSum += dValue * dValue * dValue * dValue;
                        break;
                    default:
                        throw std::runtime_error("Wrong moment order.");
                }

            };
            if (method == MomentMethod::MonteCarloSamples) {
                count = nMax;
                if (count == 0)
                    count = 10000 * nOrder * nOrder * nOrder * m_pda.m_nOrder + 1;
                for (size_t i = 0; i < count; ++i) {
                    processSample();
                }
            }
            else {
                assert(method == MomentMethod::MonteCarloTime);
                clock_t nMaxClocks = nMax; //(nMax == 0 ? 1000 : nMax);
                clock_t nStarttime;
                clock_t nClocks = 0;
                do {
                    nStarttime = clock();
                    processSample();
                    nClocks += clock() - nStarttime;
                    ++count;
                } while (nClocks < nMaxClocks);
            }
            dMoment = dMomentSum / static_cast<pdaValueType>(count);
        }
        return dMoment;
    }

/** 
 * Calculates the first raw moments.
 * @param nMaxOrder Order of the highest raw moment to determine.
 * @param method
 * @param nMaxT Meaning depends on MomentMethod
 */
    std::vector<pdaValueType> PDV::getRawMoments(const size_t nMaxOrder, const MomentMethod method, const size_t nMax) const {
        assert(nMaxOrder <= 4);
        std::vector<pdaValueType> aMoments(nMaxOrder+1);
        if (method == MomentMethod::FullMoments || method == MomentMethod::Auto)
            for (size_t nOrder = 0; nOrder <= nMaxOrder; ++nOrder) {
                size_t nMaxTotalPowersSum = (nMax == 0 ? std::max(m_pda.getOrder() * nOrder, nOrder) : nMax);
                aMoments[nOrder] = getRawMoment(nOrder, method, nMaxTotalPowersSum);
            }
        else {
            assert(method == MomentMethod::MonteCarloSamples || method == MomentMethod::MonteCarloTime);
            size_t count=nMax;
            size_t nMaxTotalPowersSum = nMaxOrder * m_pda.getOrder(); // Full power combinations
            std::vector<std::vector<pdaValueType>> powers(m_pda.getNumberOfDeltas(),
                                                          std::vector<pdaValueType>(m_pda.getOrder() + 1));

            auto processSample = [&] {
                for (size_t delta = 0; delta < m_pda.getNumberOfDeltas(); ++delta) {
                    pdaValueType deltaValue = m_pda.m_deltaDistributions[delta]->drawSample();
                    for (size_t order = 0; order <= m_pda.getOrder(); ++order)
                        powers[delta][order] = ::pow(deltaValue, static_cast<double>(order));
                }
                pdaValueType dValue = 0;
                pdaValueType dAddend;
                Util::PowersIterator pi(m_pda, 1, nMaxTotalPowersSum);
                do {
                    dAddend = m_aCoeff[pi.getPositions()[0]];
                    if (similar(dAddend, 0)) continue;
                    std::vector<size_t>& aTotalPowers = pi.getFactorsPowersSum();
                    for (size_t nDelta = 0; nDelta < m_pda.getNumberOfDeltas(); ++nDelta) {
                        if (aTotalPowers[nDelta] > m_pda.getOrder()) {
                            dAddend = 0;
                            break;
                        }
                        dAddend *= powers[nDelta][aTotalPowers[nDelta]];
                    }
                    dValue += dAddend;
                } while (pi.next());
                aMoments[0] += 1;
                assert(aMoments.size() == nMaxOrder+1);
                if (nMaxOrder >= 1)
                    aMoments[1] += dValue;
                if (nMaxOrder >= 2)
                    aMoments[2] += dValue*dValue;
                if (nMaxOrder >= 3)
                    aMoments[3] += dValue*dValue*dValue;
                if (nMaxOrder >= 4)
                    aMoments[4] += dValue*dValue*dValue*dValue;
                //dMomentSum += pow(dValue, nOrder);
            };

            if (method == MomentMethod::MonteCarloSamples) {
                if (count == 0)
                    count = 10000 * nMaxOrder * nMaxOrder * nMaxOrder * m_pda.m_nOrder + 1;
                for (size_t i = 0; i < count; ++i) {
                    processSample();
                }
            } else {
                assert(method == MomentMethod::MonteCarloTime);
                clock_t nMaxClocks = (nMax == 0 ? 1000 : nMax);
                clock_t nStarttime;
                clock_t nClocks = 0;
                count = 0;
                do {
                    nStarttime = clock();
                    processSample();
                    nClocks += clock() - nStarttime;
                    ++count;
                } while (nClocks < nMaxClocks);
            }
            for (pdaValueType& moment: aMoments)
                moment /= static_cast<pdaValueType>(count);
        }
        return aMoments;
    }

/** 
 * Calculates the nOrder'th central moment.
 * @param nOrder Order of the central moment to calculate. Order <= 4.
 * @returns Value of the central moment
 */
    pdaValueType PDV::getCentralMoment(const size_t nOrder, const MomentMethod method, const size_t nMax) const {
        assert(nOrder <= 4);

        std::vector<pdaValueType> aRawMoments;
        if (nOrder > 1)
            aRawMoments = getRawMoments(nOrder, method, nMax);

        switch (nOrder) {
            case 0:
                return 1;
                break;
            case 1:
                return 0;
                break;
            case 2:
                return aRawMoments[2] - ::pow(aRawMoments[1], 2);
                break;
            case 3:
                return aRawMoments[3] - 3 * aRawMoments[1] * aRawMoments[2] + 2 * ::pow(aRawMoments[1], 3);
                break;
            case 4:
                return aRawMoments[4] - 4 * aRawMoments[1] * aRawMoments[3] + 6 * ::pow(aRawMoments[1], 2) * aRawMoments[2] - 3 * ::pow(aRawMoments[1], 4);
                break;
            default:
                break;
        }
        return 0;
    }

/** 
 * Calculates the first central moments.
 * Exceptionally, the first moment is set to the first raw moment (expected value), not zero!
 * @param nMaxOrder Order of the highest raw moment to determine.
 */
    std::vector<pdaValueType> PDV::getCentralMoments(const size_t nMaxOrder, const MomentMethod method, const size_t nMax) const {
        assert(nMaxOrder <= 4);
        std::vector<pdaValueType> aMoments(nMaxOrder+1);

        std::vector<pdaValueType> aRawMoments = getRawMoments(nMaxOrder, method, nMax);

        aMoments[0] = 1;
        if (nMaxOrder >= 1)
            aMoments[1] = aRawMoments[1]; // Not 0 !;
        if (nMaxOrder >= 2)
            aMoments[2] = aRawMoments[2] - ::pow(aRawMoments[1], 2);
        if (nMaxOrder >= 3)
            aMoments[3] = aRawMoments[3] - 3 * aRawMoments[1] * aRawMoments[2] + 2 * ::pow(aRawMoments[1], 3);
        if (nMaxOrder >= 4)
            aMoments[4] = aRawMoments[4] - 4 * aRawMoments[1] * aRawMoments[3] + 6 * ::pow(aRawMoments[1], 2) * aRawMoments[2] - 3 * ::pow(aRawMoments[1], 4);
        return aMoments;
    }

/** 
 * Determins the sensitivity or partial derivative 
 * with respect to \f[\Delta_{\mbox{nDeltaNumber}}\f].
 * @param nDeltaNumber Index of the Delta Symbol.
 * The index of the first symbol is 0.
 * @returns Sensitivity value
 */
    pdaValueType PDV::getSensitivity(const size_t nDeltaNumber){
        std::vector<size_t> aPowers(m_pda.getNumberOfDeltas());
        for (size_t i = 0; i < m_pda.getNumberOfDeltas(); ++i)
            if (i == nDeltaNumber) aPowers[i] = 1;
            else aPowers[i] = 0;
        pdaValueType sensitivity = m_aCoeff[m_pda.calcCoeffPos(aPowers)];
        return sensitivity;
    }

/**
 * drawSamples
 * @param sampleCount
 * @return
 */
    std::vector<pdaValueType> PDV::drawSamples(size_t sampleCount) {
        assert (m_pda.getOrder() > 0);
        std::vector<pdaValueType> r(sampleCount);
        std::vector<std::vector<pdaValueType>> deltaPowers(m_pda.getNumberOfDeltas(),
                                                           std::vector<pdaValueType>(m_pda.getOrder()+1));
        for (size_t i=0; i<sampleCount; ++i) {
            // Draw Deltas:
            for (size_t d=0; d<m_pda.getNumberOfDeltas(); ++d) {
                deltaPowers[d][0] = 1.0;
                deltaPowers[d][1] = m_pda.getDeltaDistribution(d)->drawSample();
                for (size_t p=2; p<m_pda.getOrder()+1; ++p)
                    deltaPowers[d][p] = deltaPowers[d][p-1] * deltaPowers[d][1];
            }

            // Evaluate PDV:
            pdaValueType v = 0.0;
            Util::PowersIterator pi(m_pda, 1);
            do {
                pdaValueType term = getCoeff(pi.getFactorsPowersSum());
                for (size_t d = 0; d < m_pda.getNumberOfDeltas(); ++d)
                    term *= deltaPowers[d][pi.getFactorsPowersSum()[d]];
                v += term;
            } while (pi.next());
            r[i] = v;
        }
        std::sort(r.begin(), r.end());
        return r;
    }

/**
 * Returns a vector of CDF function points
 * @param pointCount
 * @param averagingCount
 * @return Vector of x/CDF(x) pairs
 */
    std::vector<std::pair<pdaValueType, pdaValueType>> PDV::estimateCDF(size_t pointCount, size_t averagingCount) {
        std::vector<std::pair<pdaValueType, pdaValueType>> r(pointCount, {0,0});
        std::vector<pdaValueType> s = this->drawSamples(pointCount * averagingCount);
        for (size_t i=0; i<pointCount; ++i) {
            pdaValueType average = 0;
            for (size_t j=0; j<averagingCount; ++j)
                average += s[i*averagingCount+j];
            average /= static_cast<pdaValueType>(averagingCount);
            pdaValueType percentage = (static_cast<pdaValueType>(i*averagingCount) + static_cast<pdaValueType>(averagingCount)/2)
                                      / static_cast<pdaValueType>(pointCount*averagingCount);
            r[i] = std::make_pair(average, percentage);
        }
        return r;
    }

/**
 * Returns a vector of PDF function points
 * @param x
 * @param pointCount
 * @param averagingCount
 * @return Vector of x/PDF(x) pairs
 */
    std::vector<std::pair<pdaValueType, pdaValueType>> PDV::estimatePDF(size_t pointCount, size_t averagingCount) {
        std::vector<std::pair<pdaValueType, pdaValueType>> r(pointCount, {0,0});
        auto cdf = this->estimateCDF(pointCount+1, averagingCount);
        for (size_t i=0; i<pointCount; ++i) {
            pdaValueType value = (cdf[i+1].first + cdf[i].first) / 2;
            pdaValueType p = (cdf[i+1].second-cdf[i].second) / (cdf[i+1].first-cdf[i].first);
            r[i] = std::make_pair(value, p);
        }
        return r;
    }

 /**
 * Determins the covariance.
 * \f[Cov(X,Y)=\langle (X-\langle X\rangle)(Y-\langle Y\rangle) \rangle\f]
 * @returns Covariance of value1 and value2
 */
    pdaValueType Cov(const PDV& value1,
                     const PDV& value2,
                     const MomentMethod method, const size_t nMax){
        return E(value1 * value2, method, nMax) - E(value1) * E(value2, method, nMax);
    }

/** 
 * Determins the statistical correlation.
 * \f[Cor(X,Y)=\frac{Cov(X,Y)}{sqrt(Var(X))sqrt(Var(Y))}\f]
 * @returns Correlation of value1 and value2
 */
    pdaValueType Cor(const PDV& value1,
                     const PDV& value2,
                     const MomentMethod method, const size_t nMax){
        return Cov(value1, value2, method, nMax) / ::sqrt(Var(value1, method, nMax) * Var(value2, method, nMax));
    }

/** 
 * Operator [ostream] << [PDV]
 */
    std::ostream& operator<<(std::ostream& os,
                             const PDV& value){
        value.dump(os);
        return os;
    }

/** 
 * Dumps the PDV as a sum of all Delta symbol products
 * @param os Stream to wich the dump is done.
 */
    void PDV::dump(std::ostream &os) const {
        Util::PowersIterator pi(m_pda, 1, m_pda.getOrder());
        do{
            if (pi.getPosition() != 0) {
                if (similar(m_aCoeff[pi.getPosition()], 0)) continue;
                os << " + ";
            }
            os << m_aCoeff[pi.getPosition()];
            for (size_t i = 0; i < m_pda.getNumberOfDeltas(); ++i){
                if (pi.getFactorsPowersSum()[i] > 0) {
#if false
                    os << " Delta_" << i;
#else
                    os << " " << m_pda.getDeltaName(i);
#endif
                }
                if (pi.getFactorsPowersSum()[i] > 1)
                    os << "^" << pi.getFactorsPowersSum()[i];
            }
        } while (pi.next());
        //os << endl;
    }

/** 
 * Sets small coefficient to zero
 */
    void PDV::clean() {
        static const pdaValueType epsilon = 1e-14;
        Util::PowersIterator pi(m_pda, 1);
        do{
            if (fabs(m_aCoeff[pi.getPosition()]) < epsilon)
                m_aCoeff[pi.getPosition()] = 0;
        } while (pi.next());
    }

/** 
 * Check
 */
    void PDV::check() {
        Util::PowersIterator pi(m_pda, 1);
        do {
            assert(m_pda.calcCoeffPos(pi.getFactorsPowersSum()) == pi.getPosition());
            std::cout << "Pos:" << pi.getPosition() << std::endl;
        } while (pi.next());

        Util::PowersIterator pi2(m_pda, 2);
        do{
            assert(m_pda.calcCoeffPos(pi2.getFactorsPowersSum()) == pi2.getPosition());
            std::cout << "Pos:" << pi2.getPosition() << std::endl;
        } while (pi2.next());
    }

}