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

#include "pda_matrix.h"
#include <cassert>
#include <iostream>
#include <cmath>

namespace Pda {
    namespace Util {

//----------------------------------------------------------------
// Matrix:
//----------------------------------------------------------------

        /**
         * Constructor that makes a symmetric matrix of a given size.
         */
        Matrix::Matrix(PDA &pda, size_t size) :
                m_pda{pda},
                m_n{size},
                m_e{size, std::vector<PDV>{size, PDV(pda)}} {}

        /**
         * Copy constructor
         */
        Matrix::Matrix(const Util::Matrix &M) :
                m_pda{M.m_pda},
                m_n(M.m_n),
                m_e(m_n, std::vector<PDV>{m_n, PDV(m_pda)}) {
            for (size_t i = 0; i < m_n; ++i)
                m_e[i] = M.m_e[i];
        }

        /**
         * Assignment
         */
        Matrix &Matrix::operator=(const Util::Matrix &M) {
            for (size_t i = 0; i < m_n; ++i)
                m_e[i] = M[i];
            return *this;
        }

        /**
         * Plus Assignment
         */
        Matrix &Matrix::operator+=(const Util::Matrix &M) {
            for (size_t i = 0; i < m_n; ++i)
                for (size_t j = 0; j < m_n; ++j)
                    m_e[i][j] += M[i][j];
            return *this;
        }

        /**
         * Matrix Vector multiplication
         */
        Vector Matrix::operator*(const Vector &b) const {
            Vector r(m_pda, this->getSize());
            for (size_t i = 0; i < m_n; ++i) {
                r[i] = PDV{m_pda, 0.0};
                for (size_t j = 0; j < m_n; ++j)
                    r[i] += (*this)[i][j] * b[j];
            }
            return r;
        }

        /**
         * Operator [ostream] << [Matrix]
         */
        std::ostream &operator<<(std::ostream &os, const Matrix &value) {
            value.dump(os);
            return os;
        }

        /**
         * Dumps the Matrix as a sum of all Delta symbol products
         * @param os Stream to wich the dump is done.
         */
        void Matrix::dump(std::ostream &os) const {
            for (size_t i = 0; i < m_n; ++i) {
                //os <<  m_e[i];
                for (size_t j = 1; j < m_n; ++j)
                    os << "\t" << m_e[i][j];
                os << std::endl;
            }
        }

    }
}