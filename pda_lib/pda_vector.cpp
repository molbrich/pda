/**
 * pda_vector.cpp
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

#include "pda_vector.h"
#include "pda_matrix.h"
#include <cassert>
#include <iostream>
#include <cmath>

namespace Pda {
    namespace Util {

//----------------------------------------------------------------
// Vector:
//----------------------------------------------------------------

/**
 * Constructor that makes a vector of a given size. 
 */
        Vector::Vector(PDA &pda, size_t size) :
                m_pda(pda),
                m_n(size),
                m_e(m_n, PDV(pda)) {}

/**
 * Gets size
 */
        size_t Vector::getSize() const {
            return m_n;
        }

/** 
 * Operator []
 */
        PDV const &Vector::operator[](const size_t i) const {
            return m_e[i];
        }

        PDV &Vector::operator[](const size_t i) {
            return m_e[i];
        }

/**
 * Assignment
 */
        Vector &Vector::operator=(const Vector &a) {
            m_e = a.m_e;
            return *this;
        }

/**
 * Add a
 */
        Vector &Vector::operator+=(const Vector &a) {
            for (size_t i = 0; i < a.m_n; ++i)
                (*this)[i] += a[i];
            return *this;
        }

/**
 * Return -*this
 */
        Vector Vector::operator -() const {
            Vector temp(this->m_pda, this->m_n);
            for (size_t i = 0; i < m_n; ++i)
                temp.m_e[i] = -m_e[i];
            return temp;
        }

        /**
 * *this = -*this
 */
        void Vector::neg() {
            for (size_t i = 0; i < m_n; ++i)
                m_e[i] = -m_e[i];
        }

/**
 * Subtract b
 */
        Vector &Vector::operator-=(const Vector &a) {
            for (size_t i = 0; i < a.m_n; ++i)
                (*this)[i] -= a[i];
            return *this;
        }

/**
 * Scalar Multiplication (with scalar result)
 */
        PDV Vector::operator*(const Vector &a) const {
            PDV res(m_pda);
            for (size_t i = 0; i < m_n; ++i)
                res += (*this)[i] * a[i];
            return res;
        }

/**
 * Multiplication with scalar
 */
        Vector &Vector::operator*=(const PDV &a) {
            for (size_t i = 0; i < m_n; ++i)
                (*this)[i] *= a;
            return (*this);
        }

/**
 * Norm
 */
        PDV Vector::norm() {
            PDV normValue(m_pda, 0);
            pdaValueType c;
            for (size_t i = 0; i < m_n; ++i) {
                c = m_e[i].getMean();
                normValue += c * c;
            }
            return sqrt(normValue);
        }

/** 
 * Operator [ostream] << [Vector]
 */
        std::ostream &operator<<(std::ostream &os, const Vector &value) {
            value.dump(os);
            return os;
        }

/** 
 * Dumps the Vector as a sum of all Delta symbol products
 * @param os Stream to wich the dump is done.
 */
        void Vector::dump(std::ostream &os) const {
            os << m_e[0];
            for (size_t i = 1; i < m_n; ++i)
                os << ", " << m_e[i];
        }

    } // end of Util
} // end of Pda