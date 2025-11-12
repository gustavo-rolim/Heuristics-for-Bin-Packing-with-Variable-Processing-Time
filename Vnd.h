#ifndef VND_H
#define VND_H

// CALL AUXILIARY .H FILES

#include "Data.h"
#include "Objective.h"

// CALL C++ LIBRARIES

#include <iostream>
#include <limits>
#include <vector>

double item_insertion(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	const double of) {

	double best_of = of;
	int len = sol.size();

	// Auxiliary variables to mark the best item position

	int rmv_item_bin = -1;
	int rmv_item_pos = -1;

	// Auxiliary variable to mark the best insetion position

	int ins_item_bin = -1;
	int ins_item_pos = -1;

	// Auxiliary binary variables

	bool open_new_bin = false;

	for (int i = 0; i < len; ++i) {

		int bin_len = sol[i].items.size();
		bool singleton = false;

		for (int j = 0; j < bin_len; ++j) {

			int itm = sol[i].items[j];

			// Remove the current item from solution 

			if (bin_len > 1) {
				sol[i].items.erase(sol[i].items.begin() + j);
				sol[i].size -= inst.l[itm];
			}
			else {
				sol.erase(sol.begin() + i);
				singleton = true;
			}

			// Execute insertion movement 

			for (int k = 0; k <= sol.size(); ++k) {

				// Create a new bin 

				Bin new_bin;
				new_bin.items.push_back(itm);
				new_bin.size += inst.l[itm];
				sol.insert(sol.begin() + k, new_bin);

				// Compute the objective function

				double new_of = compute_objective(sol, inst);

				if (new_of < best_of) {
					best_of = new_of;
					open_new_bin = true;
					rmv_item_bin = i;
					rmv_item_pos = j;
					ins_item_bin = k;
				}

				// Return the current solution

				sol.erase(sol.begin() + k);

				// Evaluate the insertion of the candidate item in an existing bin

				if (sol.size() != 0 && k != sol.size()) {
					// Only checks the insertion of a new item if it fits the existing bin 

					if (sol[k].size + inst.l[itm] <= inst.lmax) {

						for (int x = 0; x <= sol[k].items.size(); ++x) {
							// Insert the item into every possible position of the existing bin

							sol[k].items.insert(sol[k].items.begin() + x, itm);

							// Compute the objective function

							new_of = compute_objective(sol, inst);

							if (new_of < best_of) {
								best_of = new_of;
								open_new_bin = false;
								rmv_item_bin = i;
								rmv_item_pos = j;
								ins_item_bin = k;
								ins_item_pos = x;
							}

							// Return to the current solution

							sol[k].items.erase(sol[k].items.begin() + x);
						}
					}
				}
			}

			// Reinsert the candidate item into the solution 

			if (singleton == true) {
				Bin new_bin;
				new_bin.items.push_back(itm);
				new_bin.size += inst.l[itm];
				sol.insert(sol.begin() + i, new_bin);
			}
			else {
				sol[i].items.insert(sol[i].items.begin() + j, itm);
				sol[i].size += inst.l[itm];
			}
		}
	}

	// Accept the movement if there is an improvement

	if (best_of < of) {
		// Remove the selected item

		int itm = sol[rmv_item_bin].items[rmv_item_pos];

		if (sol[rmv_item_bin].items.size() > 1) {
			sol[rmv_item_bin].items.erase(sol[rmv_item_bin].items.begin() + rmv_item_pos);
			sol[rmv_item_bin].size -= inst.l[itm];
		}
		else {
			sol.erase(sol.begin() + rmv_item_bin);
		}

		if (open_new_bin == true) {
			// Create a new bin and insert it in the sequence

			Bin insert_bin;
			insert_bin.items.push_back(itm);
			insert_bin.size += inst.l[itm];
			sol.insert(sol.begin() + ins_item_bin, insert_bin);
		}
		else {
			// Insert in the best position of the sequence

			sol[ins_item_bin].items.insert(sol[ins_item_bin].items.begin() + ins_item_pos, itm);
			sol[ins_item_bin].size += inst.l[itm];
		}
	}

	return best_of;
}

double item_swap(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	const double of) {

	double best_of = of;
	int len = sol.size();

	// Auxiliary variables to execute the swap 

	int best_bin1 = -1;
	int best_bin2 = -1;
	int best_pos1 = -1;
	int best_pos2 = -1;

	for (int i = 0; i < len; ++i) {

		int bin_len = sol[i].items.size();

		for (int j = 0; j < bin_len; ++j) {

			int itm1 = sol[i].items[j];

			for (int k = 0; k < len; ++k) {

				// Only swap items from different bins 

				if (i != k) {
					for (int x = 0; x < sol[k].items.size(); ++x) {

						int itm2 = sol[k].items[x];

						if (sol[i].size - inst.l[itm1] + inst.l[itm2] <= inst.lmax && sol[k].size - inst.l[itm2] + inst.l[itm1] <= inst.lmax) {

							// Execute the swap

							sol[i].items[j] = itm2;
							sol[k].items[x] = itm1;

							double new_of = compute_objective(sol, inst);

							if (new_of < best_of) {
								best_of = new_of;
								best_bin1 = i;
								best_pos1 = j;
								best_bin2 = k;
								best_pos2 = x;
							}

							// Return to the original sequence

							sol[i].items[j] = itm1;
							sol[k].items[x] = itm2;
						}
					}
				}
			}
		}
	}

	// Accept the movement if there is an improvement

	if (best_of < of) {
		// Execute the swap 

		int itm1 = sol[best_bin1].items[best_pos1];
		int itm2 = sol[best_bin2].items[best_pos2];
		sol[best_bin1].items[best_pos1] = itm2;
		sol[best_bin2].items[best_pos2] = itm1;

		// Update bins' lengths

		sol[best_bin1].size = sol[best_bin1].size - inst.l[itm1] + inst.l[itm2];
		sol[best_bin2].size = sol[best_bin2].size - inst.l[itm2] + inst.l[itm1];
	}

	return best_of;
}

double bin_insertion(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	const double of) {

	double best_of = of;
	int len = sol.size();

	// Auxiliary variables to mark the best insetion position

	int rmv_bin_pos = -1;
	int ins_bin_pos = -1;

	for (int i = 0; i < len; ++i) {
		// Remove the bin from its current position 

		Bin rmv_bin = sol[i];
		sol.erase(sol.begin() + i);

		for (int j = 0; j < len; ++j) {
			if (i != j) {
				// Insert the removed bin in a new position

				sol.insert(sol.begin() + j, rmv_bin);

				// Compute the objective function

				double new_of = compute_objective(sol, inst);

				if (new_of < best_of) {
					best_of = new_of;
					rmv_bin_pos = i;
					ins_bin_pos = j;
				}

				// Return to the current solution

				sol.erase(sol.begin() + j);
			}
		}

		// Return the current solution 

		sol.insert(sol.begin() + i, rmv_bin);
	}

	// Accept the movement if there is an improvement

	if (best_of < of) {
		// Select and remove the best bin from its current position 

		Bin selected_bin = sol[rmv_bin_pos];
		sol.erase(sol.begin() + rmv_bin_pos);

		// Insert in the new best position

		sol.insert(sol.begin() + ins_bin_pos, selected_bin);
	}

	return best_of;
}

double bin_swap(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	const double of) {

	double best_of = of;
	int len = sol.size();

	// Auxiliary variables to mark the best swap

	int swap_pos1 = -1;
	int swap_pos2 = -1;

	for (int i = 0; i < len; ++i) {
		for (int j = 0; j < len; ++j) {
			if (i != j) {
				// Swap two bins 

				std::swap(sol[i], sol[j]);

				// Compute the objective function

				double new_of = compute_objective(sol, inst);

				if (new_of < best_of) {
					best_of = new_of;
					swap_pos1 = i;
					swap_pos2 = j;
				}

				std::swap(sol[i], sol[j]);
			}
		}
	}

	// Accept the movement if there is an improvement

	if (best_of < of) std::swap(sol[swap_pos1], sol[swap_pos2]);

	return best_of;
}

double bin_merge(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	const double of) {

	double best_of = of;
	int len = sol.size();

	// Auxiliary variables to mark the best bin merge 

	int bin_merge_pos_i = -1;
	int bin_merge_pos_j = -1;

	for (int i = 0; i < len; ++i) {
		for (int j = i + 1; j < len; ++j) {
			if (sol[i].size + sol[j].size <= inst.lmax) {
				// Merge the bins

				int size_bin_i = sol[i].items.size();
				Bin backup = sol[j];

				sol[i].items.insert(sol[i].items.end(), sol[j].items.begin(), sol[j].items.end());

				// Remove the right most bin from its original position

				sol.erase(sol.begin() + j);

				// Compute the objective function

				double new_of = compute_objective(sol, inst);

				if (new_of < best_of) {
					best_of = new_of;
					bin_merge_pos_i = i;
					bin_merge_pos_j = j;
				}

				// Return to the original solution

				sol[i].items.resize(size_bin_i);
				sol.insert(sol.begin() + j, backup);
			}
		}
	}

	// Accept the movement if there is an improvement

	if (best_of < of) {
		// Merge the bins in the best position 

		sol[bin_merge_pos_i].items.insert(sol[bin_merge_pos_i].items.end(),
			sol[bin_merge_pos_j].items.begin(), sol[bin_merge_pos_j].items.end());

		// Update capacity

		sol[bin_merge_pos_i].size += sol[bin_merge_pos_j].size;

		// Remove the right most bin from its original position

		sol.erase(sol.begin() + bin_merge_pos_j);
	}

	return best_of;
}

double main_vnd(std::vector<Bin>& sol, const ProblemInstance& inst, double of) {

	bool improvement;

	do {
		improvement = false;

		// Try item insertion (N1)

		double new_of = item_insertion(sol, inst, of);

		if (new_of < of) {
			of = new_of;
			improvement = true;
			continue;
		}

		// Try item swap (N2)

		new_of = item_swap(sol, inst, of);

		if (new_of < of) {
			of = new_of;
			improvement = true;
			continue;
		}

		// Try bin insertion (N3)

		new_of = bin_insertion(sol, inst, of);

		if (new_of < of) {
			of = new_of;
			improvement = true;
			continue;
		}

		// Try bin swap (N4)

		new_of = bin_swap(sol, inst, of);

		if (new_of < of) {
			of = new_of;
			improvement = true;
			continue;
		}

		// Try bin merge (N5)

		new_of = bin_merge(sol, inst, of);

		if (new_of < of) {
			of = new_of;
			improvement = true;
			continue;
		}

	} while (improvement);

	return of;
}

#endif