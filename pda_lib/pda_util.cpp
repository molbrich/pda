// pda_util.cpp
// Probability DeltaDistribution Arithmetic
// (c) Markus Olbrich

#include "pda_pda.h"
#include "pda_powersiterator.h"
#include "pda_util.h"
#include "pda_vector.h"
#include <cassert>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <random>

//----------------------------------------------------------------
// Pda::Util
//----------------------------------------------------------------
namespace Pda {

    size_t Util::getBinCoeff(size_t n, size_t k) {
        assert(n >= k);
        size_t res = 1;
        if (k > n-k)
            k = n-k;
        for (size_t i=0; i < k; ++i) {
            assert(res*(n-1)/(n-1) == res);
            res *= n-i;
            res /= i+1;
        }
#ifndef NDEBUG
        double b = 1;
        for (size_t i = 1; i <= k; ++i)
            b *= static_cast<double>(n + 1 - i) / static_cast<double>(i);
        b = round(b);
        assert(static_cast<double>(static_cast<size_t>(b)) == b);
        assert(static_cast<size_t>(b) == res);
#endif
        return res;
    }

/**
 * Monte Carlo simulation
 * @param x Argument
 * @param aMoments Place for moments to return
 * @param function Function that is to simulate
 * @param nSampleNumber Number of samples
 * @returns aMoments[2..4] The first central moments. 
 *          aMoments[0] is the 1st raw moment (mean).
 * @note Assumes standard normal distributed Delta symbols.
 */
    std::vector<pdaValueType> Util::monteCarlo(PDV x,
                                               pdaValueType (*function)(pdaValueType argument),
                                               size_t nSampleNumber) {
        assert(nSampleNumber != 0);
        PDA& pda = x.getPDA();
        size_t nNumberOfDeltas = x.getPDA().getNumberOfDeltas();
        std::vector<pdaValueType> aDeltaValue(nNumberOfDeltas);
        std::vector<pdaValueType> aSums(5);
        std::vector<pdaValueType> aMoments(5);

        // Initialize aSums:
        for (size_t j = 0; j <= 4; ++j)
            aSums[j] = 0;
        // Loop over samples:
        size_t nRejected = 0;
        for (size_t nSample = 0; nSample < nSampleNumber; ++nSample) {
            // Generate random deltas:
            for (size_t i = 0; i < nNumberOfDeltas; ++i)
                aDeltaValue[i] = pda.getDeltaDistribution(i)->drawSample() * sqrt(x.getPDA().getDeltaMoment(i, 2));

            // Generate argument according to m_aCoeff:
            pdaValueType argument = 0;
            PowersIterator pi(x.getPDA(), 1);
            do {
                pdaValueType factor = 1;
                for (size_t nDelta = 0; nDelta < nNumberOfDeltas; ++nDelta)
                    for (size_t nPower = 0; nPower < pi.getFactorsPowersSum()[nDelta]; ++nPower)
                        factor *= aDeltaValue[nDelta];
                argument += x.getCoeff(pi.getFactorsPowersSum()) * factor;
            } while (pi.next());

            // Calculate function:
            pdaValueType result;
            //std::cout << "argument=" << argument << std::endl;

            if (argument <= 0) {
                ++nRejected;
                continue;
            }

            result = (*function)(argument);
            // Consider result:
            aSums[0] += 1;
            aSums[1] += result;
            aSums[2] += result * result;
            aSums[3] += result * result * result;
            aSums[4] += result * result * result * result;
        } // Loop over samples

        if (nRejected != 0)
            std::cout << nRejected << "samples rejected." << std::endl;
        auto nSamples = static_cast<pdaValueType>(nSampleNumber - nRejected);
        pdaValueType aRawMoments[5];
        // Raw moments:
        aRawMoments[1] = aSums[1] / nSamples;
        aRawMoments[2] = aSums[2] / nSamples;
        aRawMoments[3] = aSums[3] / nSamples;
        aRawMoments[4] = aSums[4] / nSamples;
        // Powers of 1st moment:
        pdaValueType m1 = aRawMoments[1];
        pdaValueType m1_2 = m1 * m1;
        pdaValueType m1_3 = m1_2 * m1;
        pdaValueType m1_4 = m1_3 * m1;
        // Resulting moments:
        aMoments[0] = aRawMoments[1];
        aMoments[1] = aRawMoments[1];
        // Central moments:
        aMoments[2] = aRawMoments[2] - m1_2;
        aMoments[3] = aRawMoments[3] - 3 * aRawMoments[1] * aRawMoments[2] + 2 * m1_3;
        aMoments[4] = aRawMoments[4] - 4 * aRawMoments[1] * aRawMoments[3]
                      + 6 * aRawMoments[1] * aRawMoments[1] * aRawMoments[2] - 3 * m1_4;
        return aMoments;
    }

/**
 * Monte Carlo Simulation
 * Performs a Monte Carlo Simulation of a given function which is
 * evaluated nominally.
 * Determins the first four moments of the resulting distribution.
 * The function can have an arbitrary number of arguments.
 * @param aArguments Vector of pointers to PDV arguments of the function
 * @param function Function that is to simulate (evaluated nominally)
 * @param nSampleNumber Number of samples
 * @returns aMoments[1..4] The first central moments.
 *          aMoments[0] is the 1st raw moment (mean).
 * @note Assumes standard normal distributed Delta symbols.
 */
    std::vector<pdaValueType> Util::monteCarlo(
            const std::vector<PDV *> &aArguments,
            PDV (*function)(PDA &pda),
            size_t nSampleNumber) {
        PDA &pda = aArguments[0]->getPDA();
        size_t nNumberOfDeltas = pda.getNumberOfDeltas();
        std::vector<pdaValueType> aDeltaValue(nNumberOfDeltas);
        std::vector<pdaValueType> aSums(5);
        std::vector<size_t> aPowers(nNumberOfDeltas);
        std::vector<pdaValueType> aOriginalNomValues(aArguments.size());
        std::vector<pdaValueType> aMoments(5);

        // Save original nominal values:
        for (size_t i = 0; i < aArguments.size(); ++i)
            aOriginalNomValues[i] = (aArguments[i])->getNom();

        // Initialize aSums:
        for (size_t i = 0; i <= 4; ++i)
            aSums[i] = 0;

        // Loop over samples:
        for (size_t nSample = 0; nSample < nSampleNumber; ++nSample) {
            // Generate random deltas:
            for (size_t i = 0; i < nNumberOfDeltas; ++i)
                aDeltaValue[i] = pda.getDeltaDistribution(i)->drawSample(); //normalRand();

            // Generate arguments according to m_aCoeffs:
#if false
            for (size_t i=0; i<aArguments.size(); ++i)
            {
                pdaValueType argument = aOriginalNomValues[i];
                for (size_t d=0; d<nNumberOfDeltas; ++d)
                    argument += pda.getBinCoeff(d,1) * aDeltaValue[i];
                aArguments[i]->setNom(argument);
            }
#else
            for (size_t i = 0; i < aArguments.size(); ++i) {
                pdaValueType argument = aOriginalNomValues[i];
                PowersIterator pi(pda, 1);
                do {
                    if (pi.getTotelPowersSum()) {
                        pdaValueType factor = 1;
                        auto &factorsPowersSum = pi.getFactorsPowersSum();
                        for (size_t nDelta = 0; nDelta < nNumberOfDeltas; ++nDelta) {
                            for (size_t nPower = 0; nPower < factorsPowersSum[nDelta]; ++nPower)
                                factor *= aDeltaValue[nDelta];
                        }
                        argument += aArguments[i]->getCoeff(pi.getFactorsPowersSum()) * factor;
                    }
                } while (pi.next());
                aArguments[i]->setNom(argument);
                //std::cout << "Argument #" << i << ": " << argument << endl;
            }
#endif

            // Calculate function nominally:
            size_t oldOrder = pda.getOrder();
            pda.m_nOrder = 0; // Hack for performing nominal evaluation
            pdaValueType result = (*function)(pda).getNom();
            pda.m_nOrder = oldOrder;
            //cout << "Result: " << result << endl;

            // Consider result:
            aSums[0] += 1;
            aSums[1] += result;
            aSums[2] += result * result;
            aSums[3] += result * result * result;
            aSums[4] += result * result * result * result;
        } // Loop over samples

        // Restore original nominal values:
        for (size_t i = 0; i < aArguments.size(); ++i)
            aArguments[i]->setNom(aOriginalNomValues[i]);

        auto nSampleNumberDivisor = static_cast<pdaValueType>(nSampleNumber);
        pdaValueType aRawMoments[5];
        // Raw moments:
        aRawMoments[1] = aSums[1] / nSampleNumberDivisor;
        aRawMoments[2] = aSums[2] / nSampleNumberDivisor;
        aRawMoments[3] = aSums[3] / nSampleNumberDivisor;
        aRawMoments[4] = aSums[4] / nSampleNumberDivisor;
        // Powers of 1st moment:
        pdaValueType m1 = aRawMoments[1];
        pdaValueType m1_2 = m1 * m1;
        pdaValueType m1_3 = m1_2 * m1;
        pdaValueType m1_4 = m1_3 * m1;
        // Resulting moments:
        aMoments[0] = aRawMoments[1];
        aMoments[1] = aRawMoments[1];
        // Central moments:
        aMoments[2] = aRawMoments[2] - m1_2;
        aMoments[3] = aRawMoments[3] - 3 * aRawMoments[1] * aRawMoments[2] + 2 * m1_3;
        aMoments[4] = aRawMoments[4] - 4 * aRawMoments[1] * aRawMoments[3]
                      + 6 * aRawMoments[1] * aRawMoments[1] * aRawMoments[2] - 3 * m1_4;
        return aMoments;
    }

/**
 * Solves Ax=b using the Conjugated Gradient Method.
 * Converges for symmetric, positiv definite A.
 * Precondition: x is the start value for iterative solver.
 */
    void Util::solve_LES_CG(const Matrix &A,
                            Vector &x,
                            const Vector &b) {
        PDA &pda = A.getPDA();
        Vector g(pda, A.getSize());
        Vector d(pda, A.getSize());
        Vector Ad(pda, A.getSize());
        PDV alpha(pda), beta(pda), temp(pda), g2(pda);

        // g = b-Ax
        g = b - A * x;
        d = g;
        d.neg();

        g2 = g * g;

        for (size_t i = 0; i < 5 * A.getSize(); ++i) {
            Ad = A * d;
            temp = d * Ad;
            if (temp == PDV(pda, 0))
                break;
            alpha = g2 / temp;
            x -= alpha * d;
            g += alpha * Ad;
            beta = 1 / g2;
            g2 = g * g;
            if (g2 == PDV(pda, 0))
                break;
            beta = g2 * beta;
            d = beta * d - g;
        }
    }

/**
 * Solves Ax=b using a simple gradient method step.
 * Converges for symmetric, positive definite A.
 * Precondition: x is the start value for iterative solver.
 */
    void Util::solve_LES_GRAD(const Matrix &A,
                              Vector &x,
                              const Vector &b) {
        PDA &pda = A.getPDA();
        Vector r(pda, A.getSize());
        Vector p(pda, A.getSize());
        Vector h(pda, A.getSize());
        Vector tempVec(pda, A.getSize());
        PDV t(pda);

        r = A * x -b;
        //cout << "Residuum: " << r << endl;
        for (size_t i = 0; i < 5 * A.getSize(); ++i) {
            p = -r;
            h = A * p;
            t = -(p * r) / (p * h);
            x += t * p;
            r += t * h;
            //cout << "Residuum: " << r << endl;
        }
    }

/**
 * Solves Ax=b using a Crout LU decomposition.
 * This solver is recommended for general purposes.
 * It is used for the Newton-Raphson nonlinear equation solver.
 * Precondition: x is the start value for iterative solver.
 * @see // http://mymathlib.webtrellis.net/matrices/c_source/linearsystems/crout_pivot.c
 */
    bool Util::solve_LES_LU(const Matrix &A,
                            Vector &x,
                            const Vector &b) {
        PDA &pda = A.getPDA();
        size_t n = A.getSize();
        Matrix LU(pda, n);
        std::vector<size_t> pivot(n);

        LU = A;
        if (!LU_decomposition(LU, pivot)) {
            //std::cerr << "Matrix singular!" << std::endl;
            return false;
        }
        else
            return LU_solve(LU, b, pivot, x);
    }

/**
 * Does an LU decomposition using Crouts algorithm.
 * @returns if solved
 */
    bool Util::LU_decomposition(Matrix &A, std::vector<size_t> &pivot) {
        size_t n = A.getSize();
        assert(pivot.size() == n);
        pdaValueType dMax;

        for (size_t k = 0; k < n; k++) {

            // Find the pivot row:
            pivot[k] = k;

            dMax = fabs(A[k][k].getNom());
            for (size_t j = k + 1; j < n; j++) {
                if (dMax < fabs(A[j][k].getNom())) {
                    dMax = fabs(A[j][k].getNom());
                    pivot[k] = j;
                }
            }

            // If it differs from the current row, interchange the two rows.
            if (pivot[k] != k) {
                A[k].swap(A[pivot[k]]);
            }

            // Break if matrix is singular
            if (PDV::similar(A[k][k].getNom(), 0.0))
                return false;

            // Find the upper triangular matrix elements for row k.
            for (size_t j = k + 1; j < n; j++) {
                A[k][j] /= A[k][k];
            }

            // Find the lower triangular matrix elements for column k
            for (size_t i = k + 1; i < n; i++) {
                for (size_t j = k + 1; j < n; j++)
                    A[i][j] -= A[i][k] * A[k][j];
            }
        }
        return true;
    }

/**
 * Solves x on the basis of an LU decomposition (LUx=b)
 * @param LU Decomposition (@see Util::LU_decomposition)
 * @param b 
 * @param pivot Exchange vector (Result of the LU decomposition)
 * @param x Resulting vector
 * @param n Size
 * @returns if solved
 */
    bool Util::LU_solve(const Matrix &LU, Vector b, const std::vector<size_t> &pivot, Vector &x) {
        size_t n = LU.getSize();
        assert(pivot.size() == n);
        PDA &pda = LU.getPDA();
        PDV dum(pda);

        // Solve the linear equation Lx = B for x, where L is a lower
        // triangular matrix.
        for (size_t k = 0; k < n; ++k) {
            if (pivot[k] != k) {
                swap(b[k], b[pivot[k]]);
            }
            x[k] = b[k];
            for (size_t i = 0; i < k; i++) {
                x[k] -= x[i] * LU[k][i];
            }
            x[k] /= LU[k][k];
        }

        // Solve the linear equation Ux = y, where y is the solution
        // obtained above of Lx = B and U is an upper triangular matrix.
        // The diagonal part of the upper triangular part of the matrix is
        // assumed to be 1.0.
        size_t k = n - 1;
        do {
            if (pivot[k] != k) {
                swap(b[k], b[pivot[k]]);
            }
            for (size_t i = k + 1; i < n; i++)
                x[k] -= x[i] * LU[k][i];
            if (PDV::similar(LU[k][k].getNom(), 0.0))
                return false;
        } while ((k--) > 0);
        return true;
    }

/**
 * Inverts A_in using Gauss elimination
 */
    Util::Matrix Util::inv_GE(Util::Matrix &A_in) {
        PDA &pda = A_in.getPDA();
        PDV temp(pda, 0.0);
        size_t n = A_in.getSize();
        Util::Matrix A(pda, n);
        A = A_in;
        Util::Matrix I(pda, n);
        for (size_t i = 0; i < n; i++)
            I[i][i] = 1;
        for (size_t i = 0; i < n; i++) {
            PDV main_element(pda, 0.0);
            size_t max_row = 0;
            for (size_t j = i; j < n; j++) {
                if (main_element.getNom() < fabs(A[j][i].getNom())) {
                    main_element = fabs(A[j][i].getNom());
                    max_row = j;
                }
            }
            if (main_element.getNom() == 0) {
                std::cerr << "A is singular!" << std::endl;
                assert(0);
            } else {
                for (size_t j = 0; j < n; j++) {
                    temp = A[i][j];
                    A[i][j] = A[max_row][j];
                    A[max_row][j] = temp;
                    temp = I[i][j];
                    I[i][j] = I[max_row][j];
                    I[max_row][j] = temp;
                }
                temp = A[i][i];
                for (size_t j = 0; j < n; j++) {
                    A[i][j] = A[i][j] / temp;
                    I[i][j] = I[i][j] / temp;
                }
                for (size_t j = 0; j < n; j++) {
                    if (i != j) {
                        temp = A[j][i];
                        for (size_t k = 0; k < n; k++) {
                            A[j][k] = A[j][k] - temp * A[i][k];
                            I[j][k] = I[j][k] - temp * I[i][k];
                        }
                    }
                }
            }
        }
        return I;
    }

/**
 * Nonlinear equation solver
 * @param x Resulting vector
 * @param load Loading function. Sets Jacobian matrix and f(x).
 * @param nMaxIterations Maximum number of iterations.
 * @returns Number of iteration steps, and if converged.
 */
    std::pair<size_t, bool> Util::newton_raphson(Vector &x,
                                                 void (*load)(const Vector &x,
                                           Vector &f,
                                           Matrix &J),
                                                 size_t nMaxIterations) {
        size_t n = x.getSize();
        PDA &pda = x.getPDA();
        Vector r(pda, n);
        Vector dx(pda, n);
        Matrix J(pda, n);
        Vector xNeu(pda, n);
        PDV zero(pda, 0);
        size_t i;
        for (i = 0; i < nMaxIterations; ++i) {
            (*load)(x, r, J);
            bool success = solve_LES_LU(J, dx, r);
            if (!success)
                return std::make_pair(i, false);
            xNeu = x - dx;
            size_t j;
            for (j = 0; j < n; ++j)
                if (!xNeu[j].similar(x[j], 1e-6)) {
                    break;
                }
            if (j == n) {
                break;
            }
            x = xNeu;
        }
        return std::make_pair(i, i != nMaxIterations);
    }

}