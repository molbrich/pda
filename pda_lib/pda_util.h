// pda_util.h
// Probability DeltaDistribution Arithmetic
// (c) Markus Olbrich

#ifndef PDA_UTIL_H
#define PDA_UTIL_H

#include "pda_pda.h"
#include "pda_pdv.h"
#include "pda_powersiterator.h"
#include "pda_matrix.h"
//#include <array>

namespace Pda {
    namespace Util {
        size_t getBinCoeff(size_t n, size_t k);
        pdaValueType getNormQuantile(pdaValueType a);
        pdaValueType getCFQuantile(pdaValueType a, std::vector<pdaValueType> aMoments);

        // Monte Carlo:
        std::vector<pdaValueType> monteCarlo(
                PDV x,
                pdaValueType (*function)(pdaValueType argument),
                size_t nSampleNumber);
        std::vector<pdaValueType> monteCarlo(
                std::function<PDV(std::vector<PDV>)> function,
                size_t nSampleNumber,
                time_t maxTime);

        // Linear Equation Solver:
        void solve_LES_CG(const Matrix& A, Vector& x, const Vector& b);
        void solve_LES_GRAD(const Matrix& A, Vector& x, const Vector& b);
        bool solve_LES_LU(const Matrix& A, Vector& x, const Vector& b);
        bool LU_decomposition(Matrix& A, std::vector<size_t>& pivot);
        bool LU_solve(const Matrix& LU, Vector b, const std::vector<size_t>& pivot, Vector& x);
        Matrix inv_GE(Matrix& A_in); // @warning Under construction

        // Nonlinear Equation Solver
        std::pair<size_t, bool> newton_raphson(Vector & x, void (*load)(const Vector& x,
                                                                        Vector & f,
                                                                        Matrix& J),
                                               size_t nMaxIteration = 60);
    }
}

#endif // PDA_UTIL_H
