/**
 * pda_util.h
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
        std::vector<pdaValueType> centralFromRawMoments(std::vector<pdaValueType>&);

        std::vector<pdaValueType> monteCarlo(
                PDV x,
                pdaValueType (*function)(pdaValueType argument),
                size_t nSampleNumber,
                clock_t maxClocks = 24*60*60*CLOCKS_PER_SEC);
        std::vector<pdaValueType> monteCarlo(
                const std::vector<PDV *> &aArguments,
                PDV (*function)(PDA &pda),
                size_t nSampleNumber,
                clock_t maxClocks = 24*60*60*CLOCKS_PER_SEC);

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
