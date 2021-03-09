/* 
 * hello-world.cpp
 * Simple demonstration example for the probability distribution arithmetic
 */

#include <pda_lib.h>
#include <iostream>

using namespace Pda; // This imports PDA and PDV

int main () {
    PDA pda(2,1);  // order: 2, #BRVs: 1
    PDV x(pda, 2); // nominal: 2
    x.setDeltaCoeff(0, 0.5); // set deviation of x to 0.5
    std::cout << log(x*x).getMean() << std::endl; // calculate mean of log(x*x)
    return 0;
}
