// Include C++ libraries

#include <chrono>
#include <iostream>
#include <random>
#include <vector>

// Include auxiliary header files

#include "Data.h"
#include "Heu.h"
#include "Io.h"
#include "Lns.h"
#include "Svc.h"
#include "Vnd.h"

int main(int argc,
	char* argv[]) {
	// ======================
	// 1. COMMAND LINE VALIDATION
	// ======================
	
	if (argc != 4) {
		std::cerr << "Please, use the following format:\n<executable_file> <instance_file> <parameter_file> <solution_method>" << std::endl;
		return EXIT_FAILURE;
	}

	// ======================
    // 2. INITIALIZATION
    // ======================

	const std::string instance_path = argv[1];
	const std::string param_path = argv[2];
	const std::string method = argv[3];
	ProblemInstance inst;

	try {
		// Load parameters

		read_parameter_file(param_path, inst);
		std::cout << "Parameters loaded from: " << param_path << std::endl;

		// Load problem instance

		read_problem_instance(instance_path, inst);
		std::cout << "Instance loaded: " << instance_path
			<< " (" << inst.n << " items, capacity " << inst.lmax << ")\n";
	}
	catch (const std::exception& e) {
		std::cerr << "Error during initialization: " << e.what() << std::endl;
		return EXIT_FAILURE;
	}

	// ======================
	// 3. START RANDOM GENERATOR
	// ======================

	std::random_device rd;
	std::mt19937 gen(rd()); 

	// ======================
	// 4. SOLUTION METHODS
	// ======================

	std::vector<Bin> sol; 
	double of = 0.0; 
	auto start_time = std::chrono::high_resolution_clock::now(); 

	if (method == "heu") {
		of = main_constructive(sol, inst);
		of = main_vnd(sol, inst, of);
	}
	else if (method == "lns") {
		of = main_constructive(sol, inst);
		of = main_vnd(sol, inst, of);
		of = main_lns(sol, inst, of, gen);
	}
	else if (method == "svc") {
		of = main_svc(sol, inst, gen);
	}
	else if (method == "svc_vnd") {
		of = main_svc(sol, inst, gen);
		of = main_vnd(sol, inst, of);
	}

	auto end_time = std::chrono::high_resolution_clock::now(); 
	double running_time = std::chrono::duration<double>(end_time - start_time).count();

	// ======================
	// 5. WRITE RESULTS
	// ======================

	write_results(inst, sol, running_time, instance_path, of, method);
}