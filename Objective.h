#ifndef OBJECTIVE_H
#define OBJECTIVE_H

// CALL AUXILIARY .H FILES

#include "Data.h"

// CALL C++ LIBRARIES 

#include <algorithm>
#include <iostream>
#include <vector>

// AUXILIARY FUNCTION TO PRINT A SOLUTION 

void PrintSolution(const std::vector<Bin>& sol) {
    for (size_t i = 0; i < sol.size(); ++i) {
        std::cout << "[";
        for (size_t j = 0; j < sol[i].items.size(); ++j) {
            std::cout << sol[i].items[j];
            if (j < sol[i].items.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "] size: " << sol[i].size;
        if (i < sol.size() - 1) {
            std::cout << ", ";
        }
    }
    std::cout << std::endl;
}

// COMPUTE THE OBJECTIVE FUNCTION VALUE

double compute_objective(const std::vector<Bin>& sol, const ProblemInstance& inst) {

    double ct = 0.0;
    double Tmax = 0.0;
    int len = sol.size();

    for (int i = 0; i < len; ++i) {
        ct += inst.s;

        for (int j = 0; j < sol[i].items.size(); ++j) {
            ct += inst.t;
            double T = std::max(0.0, ct - inst.d[sol[i].items[j]]);

            if (T > Tmax) {
                Tmax = T;
            }
        }
    }

    double of = (inst.w1 * len) + (inst.w2 * Tmax);

    return of;
}

#endif