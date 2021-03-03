/* rc.cpp
* Simple demonstration example for the probability distribution arithmetic
*/

#include <pda_lib.h>
#include <iostream>

using namespace Pda; // This imports PDA and PDV

int main () {
    PDA pda(2,2);  // Order: 2, : Independent variables: 2
    PDV r(pda, 200); // Mean of r is 200
    PDV c(pda, 0.01); // Mean of c is 0.01
    r.setDeltaCoeff(0, 50);   // Set standard deviation of r to 50 (using Delta_0)
    c.setDeltaCoeff(1, 0.001); // Set standard deviation of c to 0.01 (using Delta_1)
    PDV t(pda);
    t = r * c;
    std::cout << "Mean of t: " << t.getMean() << std::endl;
    std::cout << "Standard deviation of t: " << t.getStandardDeviation() << std::endl;
    return 0;
}
