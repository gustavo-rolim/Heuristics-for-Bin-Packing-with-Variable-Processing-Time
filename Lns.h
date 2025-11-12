#ifndef LNS_H
#define LNS_H

// Include C++ libraries

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <vector>

// Include auxiliary header files

#include "Data.h"
#include "Objective.h"
#include "Rng.h"
#include "Vnd.h"

void random_removal(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const int nD,
    std::vector<int>& R,
    std::mt19937& gen) {
    // RANDOM REMOVAL OF SOLUTION ELEMENTS

    int count = nD;

    while (count--) {
        // Select a random bin 

        int rand_pos_bin = random_int(0, sol.size() - 1, gen);

        // Select a random item

        int len_bin = sol[rand_pos_bin].items.size();
        int rand_pos_item = random_int(0, len_bin - 1, gen);
        int itm = sol[rand_pos_bin].items[rand_pos_item];

        if (len_bin == 1) {
            // Remove the existing bin from the solution 

            sol.erase(sol.begin() + rand_pos_bin);

            // Insert the item in the removed items list

            R.push_back(itm);
        }
        else {
            // Remove the item from an existing bin 

            sol[rand_pos_bin].items.erase(sol[rand_pos_bin].items.begin() + rand_pos_item);
            sol[rand_pos_bin].size -= inst.l[itm];

            // Insert the item in the removed items list

            R.push_back(itm);
        }
    }
}

void length_based_removal(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const int nD,
    std::vector<int>& R,
    std::mt19937& gen) {
    // REMOVAL OF THE LARGEST ITEM IN A BIN

    int count = nD;

    while (count--) {
        // Select a random bin

        int rand_pos_bin = random_int(0, sol.size() - 1, gen);
        int len_bin = sol[rand_pos_bin].items.size();

        if (len_bin > 1) {
            // Find the index that corresponds to the item with largest size

            int max_size = inst.l[sol[rand_pos_bin].items[0]];
            int indx_max_size = 0;

            for (int i = 0; i < sol[rand_pos_bin].items.size(); ++i) {
                if (inst.l[sol[rand_pos_bin].items[i]] > max_size) {
                    max_size = inst.l[sol[rand_pos_bin].items[i]];
                    indx_max_size = i;
                }
            }

            // Break ties arbitrarely 

            std::vector<int> idx_vec;

            for (int i = 0; i < sol[rand_pos_bin].items.size(); ++i) {
                if (inst.l[sol[rand_pos_bin].items[i]] == max_size) {
                    idx_vec.push_back(i);
                }
            }

            // Select item 

            int itm = 0;
            int len_indx_vec = idx_vec.size();

            if (len_indx_vec == 1) {
                itm = sol[rand_pos_bin].items[indx_max_size];
            }
            else if (len_indx_vec > 1) {
                int pos_vec = random_int(0, len_indx_vec - 1, gen);
                indx_max_size = idx_vec[pos_vec];
                itm = sol[rand_pos_bin].items[indx_max_size];
            }
            else {
                std::cerr << "Error, empty index vector." << std::endl;
                std::exit(EXIT_FAILURE);
            }

            // Remove the item from current solution

            sol[rand_pos_bin].items.erase(sol[rand_pos_bin].items.begin() + indx_max_size);
            sol[rand_pos_bin].size -= inst.l[itm];

            // Insert the item with the largest size of the bin on the removed items list

            R.push_back(itm);
        }
        else {
            // Remove the existing bin from the solution 

            int itm = sol[rand_pos_bin].items[0];
            sol.erase(sol.begin() + rand_pos_bin);

            // Insert the item in the removed items list

            R.push_back(itm);
        }
    }
}

void tardiness_removal(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const int nD,
    std::vector<int>& R,
    std::mt19937& gen) {
    // REMOVAL OF THE ITEMS THAT CONTRIBUTE THE MOST TO THE SOLUTION'S TARDINESS

    int count = nD;

    while (count--) {
        int ct = 0;
        int bin_pos = -1;
        int item_pos = -1;
        int Tmax = -1;

        // Compute the item representing the maximum tardiness of a solution

        for (int i = 0; i < sol.size(); ++i) {
            ct += inst.s;

            for (int j = 0; j < sol[i].items.size(); ++j) {
                ct += inst.t;
                int itm = sol[i].items[j];
                double T = std::max(0, ct - inst.d[itm]);

                if (T > Tmax) {
                    Tmax = T;
                    bin_pos = i;
                    item_pos = j;
                }
            }
        }

        int len_bin = sol[bin_pos].items.size();
        int itm = sol[bin_pos].items[item_pos];

        if (len_bin == 1) {
            sol.erase(sol.begin() + bin_pos);
        }
        else if (len_bin > 1) {
            sol[bin_pos].items.erase(sol[bin_pos].items.begin() + item_pos);
            sol[bin_pos].size -= inst.l[itm];
        }
        else {
            std::cerr << "Error, empty index vector." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        R.push_back(itm);
    }

    std::shuffle(R.begin(), R.end(), gen);
}

void early_placement_removal(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const int nD,
    std::vector<int>& R,
    std::mt19937& gen) {
    // REMOVAL OF ITEMS WITH LARGE DUE DATES FROM THE BEGINNING OF THE SEQUENCE 

    int count = nD;
    int max_d = *std::max_element(inst.d.begin(), inst.d.end());

    while (count--) {
        double max_score = -1.0;
        int len_sol = sol.size();
        int bin_pos = -1;
        int item_pos = -1;

        for (int i = 0; i < len_sol; ++i) {
            for (int j = 0; j < sol[i].items.size(); ++j) {
                double score = (inst.d[sol[i].items[j]] / max_d) + (1 - i / len_sol);

                if (max_score < score) {
                    max_score = score;
                    bin_pos = i;
                    item_pos = j;
                }
            }
        }

        int len_bin = sol[bin_pos].items.size();
        int itm = sol[bin_pos].items[item_pos];

        if (len_bin == 1) {
            sol.erase(sol.begin() + bin_pos);
        }
        else if (len_bin > 1) {
            sol[bin_pos].items.erase(sol[bin_pos].items.begin() + item_pos);
            sol[bin_pos].size -= inst.l[itm];
        }
        else {
            std::cerr << "Error, empty index vector." << std::endl;
            std::exit(EXIT_FAILURE);
        }

        R.push_back(itm);
    }

    std::shuffle(R.begin(), R.end(), gen);
}

void apply_removal_heuristic(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const int nD,
    std::vector<int>& R,
    std::mt19937& gen) {
    // APPLY A REMOVAL HEURISTIC

    int N = random_int(0, 3, gen);

    if (N == 0) {
        random_removal(sol, inst, nD, R, gen);
    }
    else if (N == 1) {
        length_based_removal(sol, inst, nD, R, gen);
    }
    else if (N == 2) {
        tardiness_removal(sol, inst, nD, R, gen);
    }
    else if (N == 3) {
        early_placement_removal(sol, inst, nD, R, gen);
    }
}

double greedy_insertion(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    std::vector<int>& R) {
    // GREEDY INSERTION

    double best_objective = 0.0f;

    while (!R.empty()) {
        // Auxiliary variables for each iteration

        bool open_bin = true;
        double min_objective = std::numeric_limits<double>::max();
        int best_pos_schedule = 0;
        int best_pos_bin = 0;
        int len = sol.size();
        int itm = R[0];

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

                // Only checks the insertion of a new item if it fits inside the existing bin 

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
        R.erase(R.begin());
    }

    return best_objective;
}

double apply_insertion_heuristic(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    std::vector<int>& R,
    std::mt19937& gen) {
    // SELECTS THE GREEDY STRATEGY 

    int r = random_int(0, 1, gen);

    if (r == 0) {
        // Sorts R according to the EDD rule

        std::sort(R.begin(), R.end(),
            [&inst](int a, int b) {
                return inst.d[a] < inst.d[b];
            });
    }

    double best_objective = greedy_insertion(sol, inst, R);

    return best_objective;
}

double main_lns(std::vector<Bin>& sol,
    const ProblemInstance& inst,
    const double of,
    std::mt19937& gen) {
    // PARAMETER DECLARATION

    const double q = 0.35;
    int nD = std::ceil(q * inst.n);

    // SET CURRENT SOLUTION 

    std::vector<Bin> curr_sol = sol;
    std::vector<Bin> inc_sol = sol;
    std::vector<Bin> best_sol = sol;

    // SET CURRENT OBJECTIVE

    double curr_of = of;
    double inc_of = of;
    double best_of = of;

    // SET INITIAL TEMPERATURE

    double T_init = (0.1 * of) / log(2);
    double T = T_init;

    // LIST OF REMOVED JOBS

    std::vector<int> R;

    // SET ITERATION COUNTER 

    int iter = 0;
    int max_iter = 1000;

    while (true) {
        // UPDATE ITERATION COUNTER 
        
        iter++;

        // APPLY REMOVAL HEURISTIC 

        apply_removal_heuristic(curr_sol, inst, nD, R, gen);

        // APPLY INSERTION HEURISTIC

        curr_of = apply_insertion_heuristic(curr_sol, inst, R, gen);

        // SOLUTION ACCEPTANCE

        double delta = curr_of - inc_of;

        if (delta < 0) {
            inc_sol = curr_sol;

            if (curr_of < best_of) {
                best_of = curr_of;
                best_sol = curr_sol;
            }
        }
        else {
            double r = random_double(0, 1, gen);

            if (r < (exp(-delta / T))) {
                inc_sol = curr_sol;
            }
            else {
                curr_sol = inc_sol;
            }
        }

        if (iter == max_iter) {
            break;
        }

        // UPDATE TEMPERATURE 

        T = T_init * (1 - iter / max_iter);
    }

    sol = best_sol;

    return best_of;
}

#endif 
