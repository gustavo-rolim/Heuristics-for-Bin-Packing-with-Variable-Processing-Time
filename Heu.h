#ifndef HEU_H
#define HEU_H

// Include C++ libraries

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

// Include auxiliary header files

#include "Data.h"
#include "Objective.h"

double main_constructive(std::vector<Bin>& sol,
    const ProblemInstance& inst) {
    // Define the set of unscheduled items 

    std::vector<int> U(inst.n);
    std::iota(U.begin(), U.end(), 0);

    // Define the best objective value 

    double best_objective = 0.0;

    while (!U.empty()) {
        // Auxiliary variables for each iteration

        bool open_bin = true;
        double min_objective = std::numeric_limits<double>::max();
        int best_pos_schedule = 0;
        int best_pos_bin = 0;
        int len = sol.size();
        int itm = U[0];

        for (int i = 0; i <= len; ++i) {
            // Evaluate the creation of a new 

            Bin new_bin;
            new_bin.items.push_back(itm);
            new_bin.size += inst.l[itm];
            sol.insert(sol.begin() + i, new_bin);

            // Compute the objective function

            double of = compute_objective(sol, inst);

            if (of < min_objective) {
                min_objective = of;
                best_pos_schedule = i;
            }

            // Return the current solution 

            sol.erase(sol.begin() + i);

            // Evaluate the insertion of the candidate item in an existing bin

            if (len != 0 && i != len) {
                // Only checks the insertion of a new item if it fits the existing bin 

                if (sol[i].size + inst.l[itm] <= inst.lmax) {
                    // Checks the insertion of the existing item into every possible position of the existing bin  

                    int len_bin = sol[i].items.size();

                    // Insert the item into every possible position of the existing bin

                    for (int j = 0; j <= len_bin; ++j) {
                        sol[i].items.insert(sol[i].items.begin() + j, itm);

                        // Compute the objective function            

                        of = compute_objective(sol, inst);

                        if (of < min_objective) {
                            min_objective = of;
                            best_pos_schedule = i;
                            best_pos_bin = j;
                            open_bin = false;
                        }

                        // Return to the current solution 

                        sol[i].items.erase(sol[i].items.begin() + j);
                    }
                }
            }
        }

        // Accept the best movement 

        if (open_bin == true) {
            // Create a new bin and insert it in the best possible position

            Bin new_bin;
            new_bin.items.push_back(itm);
            new_bin.size += inst.l[itm];
            sol.insert(sol.begin() + best_pos_schedule, new_bin);

        }
        else {
            // Insert the candidate item in the best possible position of an existing bin

            sol[best_pos_schedule].items.insert(sol[best_pos_schedule].items.begin() + best_pos_bin, itm);
            sol[best_pos_schedule].size += inst.l[itm];
        }

        best_objective = min_objective;
        U.erase(U.begin());
    }

    return best_objective;
}

#endif
