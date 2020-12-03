// pda_matrix.h
// Probability DeltaDistribution Arithmetic
// (c) Markus Olbrich

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
