// pda_vector.h
// Probability DeltaDistribution Arithmetic
// (c) Markus Olbrich

#ifndef PDA_VECTOR_H
#define PDA_VECTOR_H

#include "pda_pda.h"

namespace Pda {
    namespace Util {

        class Matrix;

/** @class Vector
 *  @brief Vector of PDV variables
 *  @author Markus Olbrich
 */
        class Vector {
        private:
            PDA &m_pda;
            /** Vector size */
            size_t m_n;
            /** Vector of components */
            std::vector<PDV> m_e;
        public:
            Vector(PDA &pda, size_t size);
            Vector(const Vector &v) = default;
            size_t getSize();
            PDA &getPDA() const { return m_pda; }
            PDV const &operator[](size_t i) const;
            PDV &operator[](size_t i);
            Vector &operator=(const Vector &a);
            Vector &operator+=(const Vector &a);
            Vector operator+(const Vector &a) const {
                Vector r(m_pda, this->m_n);
                r += a;
                return r;
            }
            Vector operator-() const;
            void neg();
            Vector &operator-=(const Vector &a);
            Vector operator-(const Vector &a) const {
                Vector r(*this);
                r -= a;
                return r;
            }
            PDV operator*(const Vector &a) const;
            Vector &operator*=(const PDV &a);
            Vector operator*(const PDV &a) const {
                Vector r(m_pda, this->m_n);
                r *= a;
                return r;
            }
            friend Vector operator*(const PDV &a, const Vector &b) {
                Vector r(b);
                r *= a;
                return r;
            }
            PDV norm();

            // output:
            friend std::ostream &operator<<(std::ostream &os, const Vector &value);
            void dump(std::ostream &os) const;
        };

    }
}

#endif // PDA_VECTOR_H
