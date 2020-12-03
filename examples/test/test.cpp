// test.cpp
// Checks the functions of the distribution arithmetics

//#include "pda_lib.h"
#include "pda_pda.h"
#include "pda_pdv.h"
#include "pda_matrix.h"
#include "pda_util.h"
#include <iostream>
#include <assert.h>
#include <math.h>

using namespace Pda;

pdaValueType myfabs(pdaValueType x) {
    if (x<0)
        return -x;
    else
        return x;
}

/**
 * f = x^2 - p = 0
 */
void function_01(const Util::Vector & x, Util::Vector & f, Util::Matrix & J) {
    PDV p(x[0].getPDA(), 3);
    p.setDeltaCoeff(0, 0.5);
    p.setDeltaCoeff(1, 0.3);
    f[0] = x[0]*x[0]-p;
    J[0][0] = 2*x[0];
}

/**
 * f = x - sqrt(p) = 0
 */
void function_02(const Util::Vector& x, Util::Vector& f, Util::Matrix& J) {
    PDV p(x[0].getPDA(), 3);
    p.setDeltaCoeff(0, 0.5);
    p.setDeltaCoeff(1, 0.3);
    PDV one(x[0].getPDA(), 1);
    f[0] = x[0]-sqrt(p);
    J[0][0] = one;
}

void checkPowersIterator(Util::PowersIterator& pi, long nExpectedCombinationCount=0) {
    long nCombinationCounter = 0;

    auto check = [&]()->void {
        const PDA& pda = pi.getPDA();
        size_t nTotalPowerSum = 0;
        for (size_t nDelta=0; nDelta < pda.getNumberOfDeltas(); ++nDelta) {
            size_t nPowerSumForFactor = 0;
            for (size_t nFactor=0; nFactor < pi.getNumberOfFactors(); ++nFactor)
                nPowerSumForFactor += pi.getFactorDeltasPowers(nFactor)[nDelta];
            assert(nPowerSumForFactor == pi.getFactorsPowersSum()[nDelta]);
            nTotalPowerSum += nPowerSumForFactor;
        }
        assert(pi.getTotelPowersSum() == nTotalPowerSum);
        nTotalPowerSum = 0;
        for (size_t nFactor=0; nFactor < pi.getNumberOfFactors(); ++nFactor) {
            size_t  nPowerSumForDelta = 0;
            for (size_t nDelta=0; nDelta < pda.getNumberOfDeltas(); ++nDelta) {
                assert(pi.getFactorDeltasPowers(nFactor)[nDelta] <= pi.getMaxDeltaPower());
                nPowerSumForDelta += pi.getFactorDeltasPowers(nFactor)[nDelta];
            }
            assert(nPowerSumForDelta == pi.getDeltasPowersSums()[nFactor]);
            assert(pi.getPositions()[nFactor] == pda.calcCoeffPos(pi.getFactorDeltasPowers(nFactor)));
            nTotalPowerSum += nPowerSumForDelta;
        }
        assert(pi.getTotelPowersSum() == nTotalPowerSum);
        ++nCombinationCounter;
    };

    //pi.iterate([&]()->void{pi.dump();});
    //std::cout << std::endl;

    nCombinationCounter = 0;
    do {
        //pi.dump();
        check();
    } while(pi.next());
    if (nExpectedCombinationCount)
        assert(nExpectedCombinationCount == nCombinationCounter);
    //std::cout << nCombinationCounter << std::endl;

    nCombinationCounter = 0;
    pi.iterate(check);
    if (nExpectedCombinationCount)
        assert(nExpectedCombinationCount == nCombinationCounter);
    //std::cout << nCombinationCounter << std::endl;

    for (size_t nFactor=0; nFactor < pi.getNumberOfFactors(); ++nFactor) {
        assert(pi.getPositions()[nFactor] == 0);
        assert(pi.getDeltasPowersSums()[nFactor] == 0);
        for (size_t nDelta=0; nDelta < pi.getPDA().getNumberOfDeltas(); ++nDelta) {
            assert(pi.getFactorsPowersSum()[nDelta] == 0);
            assert(pi.getFactorDeltasPowers(nFactor)[nDelta] == 0);
        }
    }

    for (size_t nFactor=0; nFactor < pi.getNumberOfFactors(); ++nFactor) {
        assert(pi.getPositions()[nFactor] == 0);
        assert(pi.getDeltasPowersSums()[nFactor] == 0);
        for (size_t nDelta=0; nDelta < pi.getPDA().getNumberOfDeltas(); ++nDelta) {
            assert(pi.getFactorsPowersSum()[nDelta] == 0);
            assert(pi.getFactorDeltasPowers(nFactor)[nDelta] == 0);
        }
    }
}

void da_checks() {
    PDA pda(4,1);
    pda.setDeltaAsNormal(0,2);
    PDV x(pda, 2);
    x.setDeltaCoeff(0, 1);
    x = log(x);
    //std::cout << x << std::endl;

    const pdaValueType epsilon = 1E-9;

    /*
    std::cout << "Delta symbol moments:" << std::endl;
    PDV::dumpMoments(cout);
    std::cout << std::endl;
    */

    // Powers Iterator
    {
        {
            PDA pda(0, 0);
            Util::PowersIterator pi(pda, 2);
            checkPowersIterator(pi, 1);
        }
        {
            PDA pda(2, 10);
            Util::PowersIterator pi(pda, 2);
            checkPowersIterator(pi, 4356);
        }
        {
            PDA pda(2, 10);
            Util::PowersIterator pi(pda, 2, 1);
            checkPowersIterator(pi, 2*10+1);
        }
        {
            PDA pda(2, 10);
            Util::PowersIterator pi(pda, 2, 2);
            checkPowersIterator(pi, Util::getBinCoeff(2*10,2)+2*10+2*10+1);
        }
        {
            PDA pda(4, 1);
            Util::PowersIterator pi(pda, 3);
            checkPowersIterator(pi);
        }
        {
            PDA pda(4, 1);
            Util::PowersIterator pi(pda, 4);
            checkPowersIterator(pi);
        }
        {
            PDA pda(2, 2);
            Util::PowersIterator pi(pda, 1, 2, 2);
            checkPowersIterator(pi, 6);
        }
        {
            PDA pda(2, 3);
            Util::PowersIterator pi(pda, 2, 2, 4);
            checkPowersIterator(pi, 28);
        }
        {
            PDA pda(3, 2);
            Util::PowersIterator pi(pda, 3);
            checkPowersIterator(pi);
        }
        {
            PDA pda(4, 1);
            Util::PowersIterator pi(pda, 4);
            checkPowersIterator(pi);
        }
    }

    // Variable definition and init:
    {
        PDA pda(1,2);
        PDV x(pda);
        assert(x.getNom() == 0);

        x.setNom(3);
        assert(x.getNom() == 3);

        x.setDeltaCoeff(1,1.5);
        assert(x.getDeltaCoeff(1) == 1.5);

        PDV y(x);
        assert(x==y);

        PDV z(pda);
        assert(x!=z);
        z = x;
        assert(x==z);
    }

    {
        PDV y{pda,4};
        assert(y.getNom() == 4);
    }
    // Comparison
    {
        PDV x{pda, 3};
        PDV y{pda,2};
        PDV z{pda,2};

        assert(x > y);
        assert(x >= y);
        assert(z <= y);
        assert(z == y);
        assert(x != y);
        assert(z == y);
    }
    // Simple +
    {
        PDV x{pda, 3};
        PDV y{pda, 2};
        assert(y.getNom() == 2);
        x = x+y;
        assert(x.getNom() == 5);
        assert(y.getNom() == 2);
        x = y+4;
        assert(x.getNom() == 6);
        assert(y.getNom() == 2);
        x = 4+y;
        assert(x.getNom() == 6);
        assert(y.getNom() == 2);
        PDV z(x+y);
        assert(x+y+z == y+z+x);
        assert(y+z+x == z+y+x);
        assert(y+z+x == z+(y+x));
        assert((x+y)+(x+z) == x+y+x+z);
    }
    // Simple -
    {
        PDV x{pda, 3};
        PDV y{pda, 2};

        assert(y.getNom() == 2);
        x = x-y;
        assert(x.getNom() == 1);
        assert(y.getNom() == 2);
        x = y-4;
        assert(x.getNom() == -2);
        assert(y.getNom() == 2);
        x = 4-y;
        assert(x.getNom() == 2);
        assert(y.getNom() == 2);
        PDV z(pda);
        z = y-4;
        assert(z.getNom() == -2);
        x = -x;
        x = 1.0 - (x+x);
        assert(x-(y-z) == x-y+z);
        assert((x-y)-z == x-y-z);
        assert((x-y)-(x-z) == x-y-x+z);
    }
    // Calculating assignments:
    {
        PDV x(pda, 3);
        x.setDeltaCoeff(0,2);

        {
            PDV y(x);
            assert(x==y);
            y += x;
            assert(2*x == y);
        }
        {
            PDV y(x);
            y -= x;
            assert(y == PDV(pda, 0));
        }
        {
            PDV y(x);
            y *= x;
            assert(x*x == y);
        }
        {
            PDV y(x);
            //std::cout << x << std::endl;
            //std::cout << y << std::endl;
            y /= x;
            //std::cout << y << std::endl;
            assert(y == PDV(pda, 1));
        }
    }
    {
        pda.setDeltaAsNormal(0,2);
        PDV x(pda);
        x.setNom(3);
        x.setDeltaCoeff(0,1.5);
        assert(x.getMean() == 3);
        //assert(x.getCentralMoment(2) == 1.5);
    }
    // Unary first order operation:
    {
        PDA pda(1, 2);

        pdaValueType x0 = 2;
        pdaValueType a = 1;
        pdaValueType b = 2;
        PDV x(pda, x0);
        x.setDeltaCoeff(0,a);
        x.setDeltaCoeff(1,b);
        //std::cout << x << std::endl;
        //std::cout << log(x) << std::endl;
        //std::cout << exp(log(x)) << std::endl;
        PDV result = log(x);
        std::vector<size_t> r00 = {0, 0};
        std::vector<size_t> r10 = {1, 0};
        std::vector<size_t> r01 = {0, 1};
        assert (result.getCoeff(r00) == log(x0));
        assert (result.getCoeff(r10) == a/x0);
        assert (result.getCoeff(r01) == b/x0);
    }
    // Unary fourth order operation:
    {
        PDA pda(4, 2);
        pdaValueType x0 = 2;
        pdaValueType a = 1;
        pdaValueType b = 2;
        pdaValueType c = 3;
        pdaValueType d = 4;
        PDV x_1(pda, 0);
        x_1.setDeltaCoeff(0, 1);
        PDV x_2(pda, 0);
        x_2.setDeltaCoeff(1, 1);
        PDV x(pda, x0);

        //std::cout<< x_1 << endl;
        //std::cout<< x_2 << endl;
        //std::cout<< x_1*x_2 << endl;

        x = x + x_1 * (b * x_1 + d * x_2);
        x.setDeltaCoeff(0, a);
        x.setDeltaCoeff(1, c);

        PDV result(pda);
        result = log(x);
        //std::cout << std::endl << x << std::endl;
        //std::cout << result << std::endl;
        std::vector<size_t> r00 = {0, 0};
        std::vector<size_t> r10 = {1, 0};
        std::vector<size_t> r20 = {2, 0};
        std::vector<size_t> r01 = {0, 1};
        std::vector<size_t> r02 = {0, 2};
        std::vector<size_t> r11 = {1, 1};
        assert (result.getCoeff(r00) == log(x0));
        assert (result.getCoeff(r10) == a / x0);
        assert (result.getCoeff(r20) == b / x0 - a * a / (2 * x0 * x0));
        assert (result.getCoeff(r01) == c / x0);
        assert (result.getCoeff(r02) == -c * c / (2 * x0 * x0));
        assert (result.getCoeff(r11) == d / x0 - a * c / (x0 * x0));
    }
    {
        PDA pda(4, 20);
        pdaValueType x0 = 2;
        PDV x(pda, x0);
        //std::cout << "x= " << x << std::endl;
        //std::cout << "log(exp(x))= " << log(exp(x)) << std::endl;
        assert(similar(log(exp(x)), x));
        assert(similar(exp(log(x)), x));
        assert(similar(sqrt(x)*sqrt(x), x));
        assert(similar(pow(sqrt(x), 2.0) , x));
    }
    // Math functions:
    {
        for (int nOrder = 0; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);
            //PDV::use(pda);

            //std::cout<< "nOrder:" << nOrder << endl;

            pda.setDeltaAsNormal(0,0.2);
            pda.setDeltaAsNormal(1,0.1);
            pda.setDeltaAsNormal(2,0.1);

            PDV x{pda}, y{pda};
            x.setNom(3);
            x.setDeltaCoeff(0,1);
            x.setDeltaCoeff(1,1);
            x.setDeltaCoeff(2,1.2);

            PDV v{2 * x + 2};
            PDV w{v / 2 - 1};
            assert(w==x);

            //std::cout<<"x   = " << x << endl<< flush;
            y=log(x);
            //std::cout<<"y   = " << y << endl;
            y=exp(x);
            //std::cout<<"x   = " << x << endl;
            assert(similar(exp(log(x)), x));
            y=3*x;
            //std::cout<<"3*x = " << y << endl;
            y=x*3;
            //std::cout<<"x*3 = " << y << endl;

            //std::cout<<"x   = " << x << endl;
            //std::cout<<"y   = " << y << endl;

            x = y + 3;

            //std::cout<<"x   = " << x << endl;
            //std::cout<<"y   = " << y << endl;

            y = y -2.2;

            //std::cout<<"x   = " << x << endl;
            //std::cout<<"y   = " << y << endl;
            PDV z=x*y;
            //std::cout<<"z=x*y = " << z << endl;
            //std::cout<<"x*z = " << x*z << endl;
            PDV a=z*z;
            //std::cout<<"a=z*z = " << a << endl;
            //std::cout<<"x*a = " << x*a << endl;
            //std::cout<<"x*z*z = " << x*z*z << endl;

            //std::cout <<"x   = " << x << std::endl;
            //std::cout<<"E[x] = " << x.getRawMoment(1) << endl;
            //std::cout<<"x*x = " << x*x << endl;
            //std::cout<<"N[x*x] = " << (x*x).getNom() << endl;
            //std::cout<<"E[x*x] = " << (x*x).getRawMoment(1) << endl;
            //std::cout<<"V[x*x] = " << (x*x).getCentralMoment(2) << endl;

            //std::cout<<"log(exp(x)) = " << log(exp(x)) << endl << endl;
            //std::cout<<"E[log(exp(x))] = " << (log(exp(x))).getRawMoment(1) << endl << endl;
            //std::cout<<"(log(exp(x))).getRawMoment(1)-x.getRawMoment(1)=" << (log(exp(x))).getRawMoment(1)-x.getRawMoment(1) << endl;
            assert (myfabs((log(exp(x))).getRawMoment(1)-x.getRawMoment(1)) < epsilon);
            assert (myfabs((log(exp(x))).getCentralMoment(2)-x.getCentralMoment(2))
                    < epsilon);

            //std::cout<<"sqrt(x*x) = " << sqrt(x*x) << endl << endl;
            //std::cout<<"E[sqrt(x*x)] = " << (sqrt(x*x)).getRawMoment(1) << endl << endl;
            assert (myfabs((sqrt(x*x)).getRawMoment(1)-x.getRawMoment(1))
                    < epsilon);
            assert (myfabs((sqrt(x*x)).getCentralMoment(2)-x.getCentralMoment(2))
                    < epsilon);

            //std::cout<<"-x = " << -x << endl << endl;
            assert (-(-x)==x);

            //std::cout<< "1/x = " << 1/x << endl << endl;
            //std::cout<< "inv(x) = " << inv(x) << endl << endl;
            //std::cout<< "1/x-inv(x) = " << 1/x-inv(x) << endl << endl;
            //std::cout<< "inv(inv(x)) = " << inv(inv(x)) << endl << endl;
            //std::cout<< "x-inv(inv(x)) = " << x-inv(inv(x)) << endl << endl;
            assert (1/x == inv(x));
            //std::cout<< "(x-inv(inv(x))).getMean()" << (x-inv(inv(x))).getMean() << endl;
            assert (myfabs((x-inv(inv(x))).getMean()) < epsilon);
            assert ((x-inv(inv(x))).getNom() == 0);

            //std::cout <<"x*inv(x) = " << x*inv(x) << std::endl << std::endl;
            assert ((x*inv(x)).getMean() == 1);
            assert ((x*inv(x)).getNom() == 1);

            //std::cout<<"x^2.0 = " << (x^2.0) << endl << endl;
            //std::cout<<"pow(x,2.0) = " << pow(x,2.0) << endl << endl;
            assert ((x^2.0) == pow(x,2.0));
            assert ((2.0^x) == pow(PDV{pda, 2.0},x));
            assert ((x^y) == pow(x,y));
            //std::cout<<"x^2 = " << (x^2) << endl << endl;
            //std::cout<<"x^2.0 = " << (x^2.0) << endl << endl;
            //std::cout<<"x*x = " << x*x << endl << endl;
            //std::cout<<"x^2-x*x = " << (x^2)-x*x << endl << endl;
            assert (myfabs(((x^2)-x*x).getRawMoment(1)) < epsilon);
            assert (myfabs(((x^2)-x*x).getCentralMoment(2)) < epsilon);

            // trigonometry:
            //assert (myfabs(sin(x).getMean()-sin(x)) < epsilon);
        }
    }
    // FullMoments:
    {
        {
            PDA pda(4, 1);
            PDV x(pda,0);
            x.setDeltaCoeff(0, 1);
            auto rawMoments = x.getRawMoments(4);
            assert(rawMoments[0] == 1);
            auto centralMoments = x.getCentralMoments(4);
            assert(centralMoments[0] == 1);
        }
        for (size_t nOrder = 0; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 2);
            //std::cout << "nOrder:" << nOrder << std::endl;

            PDV x{pda, 1};
            x.setDeltaCoeff(0, 1);
            x.setDeltaCoeff(1, 0.5);
            assert(x.getNom() == 1);
            assert(x.getRawMoment(1) == 1);

            assert(x.getVariance() == x.getCentralMoment(2));
            auto centralMoments = x.getCentralMoments(4);
            if (centralMoments[2]) {
                //std::cout << x.getSkewness() << " " << std::flush;
                //std::cout << x.getCentralMoment(3) / pow(x.getCentralMoment(2), 1.5) << std::endl;
                assert(x.getSkewness() == x.getCentralMoment(3) / pow(x.getCentralMoment(2), 1.5));
                //std::cout << x.getExcessKurtosis() << " ";
                //std::cout << x.getCentralMoment(4) / pow(x.getCentralMoment(2), 2) - 3 << std::endl;
                assert(myfabs(x.getExcessKurtosis() - (x.getCentralMoment(4) / pow(x.getCentralMoment(2), 2) - 3)) <
                       epsilon);
            }

            auto cm = x.getCentralMoments(4);
            auto rm = x.getRawMoments(4);
            for (size_t n = 0; n <= 4; ++n) {
                if (n == 1)
                    assert(cm[n] == x.getRawMoment(n));
                else
                    assert(cm[n] == x.getCentralMoment(n));
                assert(rm[n] == x.getRawMoment(n));
            }
            // Normal distribution
            {
                //std::cout << "Order " << nOrder << std::endl;
                PDA pda(nOrder, 1);
                PDV x(pda, 0);
                x.setDeltaCoeff(0, 1.0);
                pda.setDeltaDistribution(0, Util::NormalDistribution(1));
                auto m = x.getRawMoments(nOrder, MomentMethod::FullMoments);
                for (size_t i = 0; i < nOrder; ++i) {
                    //std::cout << m[i] << "\t" << pda.getDeltaMoment(0, i) << std::endl;
                    assert(m[i] == pda.getDeltaMoment(0, i));
                }
            }
            // Log normal distribution
            {
                //std::cout << "Order " << nOrder << std::endl;
                PDA pda(nOrder, 1);
                PDV x(pda, 0);
                x.setDeltaCoeff(0, 1.0);
                pda.setDeltaDistribution(0, Util::LogNormalDistribution(1,1));
                auto m = x.getRawMoments(nOrder, MomentMethod::FullMoments);
                for (size_t i = 1; i < nOrder; ++i) {
                    //std::cout << m[i] << "\t" << pda.getDeltaMoment(0, i) << std::endl;
                    assert(m[i] == pda.getDeltaMoment(0, i));
                }
            }
        }
    }
    // Covariance:
    {
        for (size_t nOrder = 2; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);

            PDV x{pda, 1};
            x.setDeltaCoeff(0, 1);
            x.setDeltaCoeff(1, 0.5);
            PDV y{pda, 5};
            y.setDeltaCoeff(0, 5);
            y.setDeltaCoeff(2, 0.3);

            //std::cout << "Cov(x,y)=" << Cov(x,y) << std::endl;
            //std::cout << Var(x+y) << "=";
            //std::cout << Var(x)+Var(y)+2*Cov(x,y) << std::endl;
            assert(myfabs(Var(x+y) - Var(x) - Var(y) - 2*Cov(x,y)) < epsilon);
            //std::cout << Cor(x,y)*sqrt(Var(x)*Var(y)) << "=";
            //std::cout << Cov(x,y) << std::endl;
            assert(myfabs(Cor(x,y)*sqrt(Var(x)*Var(y)) - Cov(x,y)) < epsilon);
        }
    }
    // Sensitivity:
    {
        for (int nOrder = 1; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);

            PDV x{pda, 1};
            x.setDeltaCoeff(0, 1);
            x.setDeltaCoeff(1, 0.5);
            PDV y={pda, 5};
            y.setDeltaCoeff(0, 5);
            y.setDeltaCoeff(2, 0.3);

            assert(x.getSensitivity(0) == 1);
            assert(x.getSensitivity(1) == 0.5);

            assert(y.getSensitivity(0) == 5);
            assert(y.getSensitivity(1) == 0);
            assert(y.getSensitivity(2) == 0.3);

            //std::cout << "dXY/dD0=" << (x*y).getSensitivity(0) << std::endl;
            //std::cout << "dXY/dD1=" << (x*y).getSensitivity(1) << std::endl;
            //std::cout << "dXY/dD2=" << (x*y).getSensitivity(2) << std::endl;
            assert((x*y).getSensitivity(0) == 10 );
            assert((x*y).getSensitivity(1) == 2.5);
            assert((x*y).getSensitivity(2) == 0.3);

            //std::cout << exp(x).getSensitivity(0) << std::endl;
        }
    }
    // Vector, Matrics:
    {
        for (int nOrder = 1; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);

            PDV x{pda,1};
            x.setDeltaCoeff(0, 1);
            x.setDeltaCoeff(1, 0.5);
            PDV y{pda, 5};
            y.setDeltaCoeff(0, 5);
            y.setDeltaCoeff(2, 0.3);

            size_t size = 5;

            Util::Vector a{pda, size};
            Util::Vector b{pda, size};
            a[0] = x;
            //std::cout<< a << endl;
            b[0] = y;
            //std::cout<< b << endl;
            Util::Vector c{pda, size};
            c = a+b;
            //std::cout<< c << endl;
            //std::cout<< a.mul(&b) << endl;

            Util::Matrix A{pda, size};
            A[0][0] = 11;
            A[0][1] = 12;
            A[1][0] = 21;
            //std::cout<< A;

            //std::cout<< endl;
        }
    }
    // Linear Equation Solvers:
    {
        for (int nOrder = 1; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);
            //PDV::use(pda);
            //std::cout<< "Order: " << nOrder << endl;

            PDV xd{pda, 1};
            xd.setDeltaCoeff(0, 0.5);
            xd.setDeltaCoeff(1, 1);
            PDV yd{pda, 4};
            yd.setDeltaCoeff(0, 5);
            yd.setDeltaCoeff(2, 0.3);

            size_t size = 2;

            Util::Matrix A{pda, size};
            A[0][0] = xd;
            assert(A[0][0] == xd);
            A[0][1] = 4;
            A[1][0] = 1;
            A[1][1] = 0.5;

            Util::Vector b{pda, size};
            b[0] = 1;
            b[1] = -0.307;

            Util::Vector x{pda, size};
            x[0] = -0.493714;
            x[1] = 0.373429;
            x[0] = -1;
            x[1] = 1;

            // LU decomposition
            Util::solve_LES_LU(A,x,b);
            //std::cout<< "x=" << x << endl;
            Util::Vector r = A * x;
            r = r - b;
            //std::cout << "Residuum LU: " << r << std::endl;

            // Inversion
            Util::Matrix I = Util::inv_GE(A);
            r = I * b;
            r = r - x;
            //std::cout << "Residuum GE: " << r << std::endl;
        }
    }
    // Nonlinear Equation Solver:
    {
        for (int nOrder = 1; nOrder <= 4; ++nOrder) {
            PDA pda(nOrder, 3);
            //PDV::use(pda);
            //std::cout << "Order: " << nOrder << std::endl;

            PDV xd={pda, 1};
            xd.setDeltaCoeff(0, 0.5);
            xd.setDeltaCoeff(1, 1);
            PDV yd={pda, 4};
            yd.setDeltaCoeff(0, 5);
            yd.setDeltaCoeff(2, 0.3);

            size_t size = 1;

            Util::Vector x{pda, size};
            x[0] = 1.732;

            Util::newton_raphson(x, function_01);
            //std::cout << "x=" << x << std::endl;

            Util::Vector r{pda, size};
            Util::Matrix J{pda, size};
            function_01(x, r, J);
            //std::cout << "Residuum: " << r << std::endl;

            PDV p{pda, 3};
            p.setDeltaCoeff(0, 0.5);
            p.setDeltaCoeff(1, 0.3);
            //std::cout << "Direct solution:" << sqrt(p) << std::endl;

            x[0] = 1.732;

            Util::newton_raphson(x, function_02);
            //std::cout << "x=" << x << std::endl;

            //std::cout << std::endl;
        }
        // Monte Carlo:
        {
            for (int nOrder = 1; nOrder <= 4; ++nOrder) {
                PDA pda(nOrder, 3);
                std::function<PDV(std::vector<PDV>)> iterateFactor = [&] (std::vector<PDV> x) -> PDV {
                    return x[0];
                };
            }
        }
    }
}

int main (void) {
    da_checks();
    std::cout << "Test done." << std::endl;
    return(0);
}
