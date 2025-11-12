#ifndef SVC_H
#define SVC_H

// Include C++ libraries

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <vector>
#include <utility>

// Include auxiliary header files

#include "Data.h"
#include "Objective.h"
#include "Rng.h"
#include "Vnd.h"

void initialize_pseudo_prices(std::vector<std::pair<int, double>>& m,
	const ProblemInstance& inst) {
	// Initializes the vector of pseudo-prices

	m.reserve(inst.n);

	// Obtain the maximum due date value to calculate pseudo-prices 

	int d_max = *std::max_element(inst.d.begin(), inst.d.end());

	// Generate pseudo-prices

	for (int i = 0; i < inst.n; ++i) {
		double pseudo_price = d_max / static_cast<double>(inst.d[i]) + inst.l[i] / static_cast<double>(inst.lmax);
		m.push_back(std::pair<int, double>(i, pseudo_price));
	}
}

void knapsack_solver(Bin& new_bin,
	std::vector<std::pair<int, double>>& m,
	double**& dp,
	const ProblemInstance& inst) {
	// Solves the 0-1 knapsack problem by maximizing pseudo-prices

	const int capacity = inst.lmax;
	const int n = m.size();

	// Construct dp table 

	for (int i = 1; i <= n; ++i) {
		int item_index = m[i - 1].first;
		double value = m[i - 1].second;
		int weight = inst.l[item_index];

		for (int w = 0; w <= capacity; ++w) {
			if (weight <= w) {
				dp[i][w] = std::max(dp[i - 1][w], dp[i - 1][w - weight] + value);
			}
			else {
				dp[i][w] = dp[i - 1][w];
			}
		}
	}

	// Backtrack to find selected items 

	int remaining_capacity = capacity;

	for (int i = n; i > 0 && remaining_capacity > 0; --i) {
		int item_index = m[i - 1].first;
		int weight = inst.l[item_index];

		if (dp[i][remaining_capacity] != dp[i - 1][remaining_capacity]) {
			new_bin.items.push_back(item_index);
			new_bin.size += weight;
			remaining_capacity -= weight;
		}
	}

	// Sort items

	std::sort(new_bin.items.begin(), new_bin.items.end());

	// Reset dp table 

	for (int i = 0; i <= inst.n; ++i) {
		std::fill(dp[i], dp[i] + capacity + 1, 0.0);
	}
}

double build_solution(std::vector<Bin>& current_sol,
	double**& dp,
	const ProblemInstance& inst,
	std::vector<std::pair<int, double>>& m) {

	while (!m.empty()) {
		// Obtain patterns by maximizing the sum of pseudo-prices using dynamic programming 

		Bin new_bin;
		knapsack_solver(new_bin, m, dp, inst);

		// Add the new bin to the solution

		current_sol.push_back(new_bin);

		// Delete the selected items from the list of cadidate pseudo-prices 

		std::erase_if(m, [&](const std::pair<int, double>& p) {
			return std::find(new_bin.items.begin(),
				new_bin.items.end(),
				p.first) != new_bin.items.end();
			});
	}

	double of = compute_objective(current_sol, inst);

	return of;
}

void critical_item(std::pair<int, double>& tardy_item,
	const ProblemInstance& inst,
	const std::vector<Bin>& sol) {
	// Identify the critical item for pseudo-price update

	double ct = 0.0f; // Completion time
	double Tmax = 0.0f; // Maximum Tardiness
	int len = sol.size(); // Solution vector length 

	for (int i = 0; i < len; ++i) {
		ct += inst.s;

		for (int j = 0; j < sol[i].items.size(); ++j) {
			ct += inst.t;
			double T = std::max(0.0, ct - inst.d[sol[i].items[j]]);

			if (T > Tmax) {
				Tmax = T;
				tardy_item.first = sol[i].items[j];
				tardy_item.second = ct;
			}
		}
	}
}

void update_pseudo_prices(const double eta1,
	const double eta2,
	const double eta3,
	const double epsilon,
	const ProblemInstance& inst,
	std::vector<Bin>& sol,
	std::vector<std::pair<int, double>>& m,
	std::mt19937& gen) {
	// Identify the critical item 

	std::pair<int, double> tardy_item;
	critical_item(tardy_item, inst, sol);

	// Create a distribution for eta calculation

	const double left_limit = 1.0; 
	const double right_limit = 1.0 + eta2;
	double eta = random_double(left_limit, right_limit, gen);

	// Update the pseudo-prices 

	for (int i = 0; i < sol.size(); i++) {
		int slack = inst.lmax - sol[i].size;

		for (int j = 0; j < sol[i].items.size(); ++j) {
			int itm = sol[i].items[j];

			if (itm != tardy_item.first) {
				m[itm].second *= (1.0 + static_cast<double>(slack) / (eta1 * inst.lmax)); // (36)
			}
			else {
				m[itm].second *= (1.0 + eta * std::pow(tardy_item.second / inst.d[itm], eta3)); // (37)
			}

			// Prevent pseudo-prices to reach zero 

			if (m[itm].second < epsilon) {
				m[itm].second = epsilon;
			}
		}
	}
}

double main_svc(std::vector<Bin>& sol,
	const ProblemInstance& inst,
	std::mt19937& gen) {
	// Parameter declaration 

	double epsilon = 1.0;
	double eta1 = 50.0;
	double eta2 = 0.3;
	double eta3 = 1.0;
	int iter_max = 1000;

	// Initialize the vector of pseudo-prices 

	std::vector<std::pair<int, double>> m;

	initialize_pseudo_prices(m, inst);

	std::vector<std::pair<int, double>> m_reference = m;

	// Initialize dynamic programming table 

	double** dp = new double* [inst.n + 1];
	for (int i = 0; i <= inst.n; ++i) {
		dp[i] = new double[inst.lmax + 1]();
	}

	// Declare solution vector 

	std::vector<Bin> current_sol = sol;

	// Store best objective value

	double best_of = std::numeric_limits<double>::max();

	// Iteratevely construct new solutions

	int iter = 0;

	while (iter <= iter_max) {
		// Increment counter 

		iter++;

		// Build a new solution

		double of = build_solution(current_sol, dp, inst, m);

		// Update pseudo-prices 

		update_pseudo_prices(eta1, eta2, eta3, epsilon,
			inst, current_sol, m_reference, gen);

		m = m_reference;

		// Solution acceptance 

		if (of < best_of) {
			best_of = of;
			sol = current_sol;
		}

		current_sol.clear();
	}

	// Free memory 

	for (int i = 0; i <= inst.n; ++i) {
		delete[] dp[i];
	}
	delete[] dp;

	return best_of;
}

#endif