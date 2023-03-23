#include <iostream>
#include <vector>
#include <algorithm>
#include "FiniteDifferences.h"
#include "FDExplicitBond.h"
#include "FDExplicitTRS.h"

int main()
{
    std::vector<double> couponSchedule{ 25, 50, 75 };
    auto Bond = std::make_shared<FDExplicitBond>(0.05, 2, 0.1, 0.3, 0.01, 0.0, 0.2, 100, 100, couponSchedule, 0.025);
    auto TRS = FDExplicitTRS(0.05, 1, 0.1, 0.3, 0.01, 0.0, 0.2, 100, 80, Bond, 1.0, 0.1, true);
    std::cout << Bond->price() << std::endl;
    std::cout << TRS.price() << std::endl;

    std::vector<double> couponSchedule2{ 30, 70, 85 };
    auto Bond2 = std::make_shared<FDExplicitBond>(0.05, 2, 0.1, 0.3, 0.01, 0.0, 0.2, 100, 100, couponSchedule2, 0.035);
    TRS.substitution(Bond2);
    std::cout << Bond2->price() << std::endl;
    std::cout << TRS.price() << std::endl;

    return 0;
}

