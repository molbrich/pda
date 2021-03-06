/**
 * pda_matrix.h
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

#ifndef PDA_MATRIX_H
#define PDA_MATRIX_H

#include "pda_pda.h"
#include "pda_pdv.h"
#include "pda_vector.h"

namespace Pda {
    namespace Util {

/** @class Matrix
 *  @brief Square Matrix of PDV variables
 *  @author Markus Olbrich
 */
        class Matrix {
        private:
            PDA& m_pda;
            /** Matrix size */
            size_t m_n;
            /** Vector of row vectors */
            std::vector<std::vector<PDV>> m_e;
        public:
            Matrix(PDA& pda, size_t size);
            Matrix(const Matrix & M);
            size_t getSize() const { return m_n; };
            PDA& getPDA() const { return m_pda; }

            /** Return i-th row vector */
            std::vector<PDV>& operator[](size_t i) { return m_e[i]; };
            const std::vector<PDV>& operator[](size_t i) const { return m_e[i]; };

            Matrix& operator = (const Matrix& M);
            Matrix& operator += (const Matrix& M);
            Vector operator * (const Vector& b) const ;

            // output:
            friend std::ostream& operator<<(std::ostream& os, const Matrix& value);
            void dump(std::ostream &os) const;
        };

    }
}

#endif // PDA_MATRIX_H
