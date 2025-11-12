#ifndef DATA_H
#define DATA_H

#include<vector>

// DEFINE THE DATA STRUCTURE FOR PROBLEM INSTANCE

struct ProblemInstance {
    int n;                    // Number of items in the instance 
    int lmax;                 // Bin capacity
    double s;                 // Setup time 
    double t;                 // Processing time
    double w1;                // Weight for number of open bins 
    double w2;                // Weight for maximum tardiness
    std::vector<int> d;       // Vector of due dates 
    std::vector<int> l;       // Vector of items' lengths 
};

// DEFINE THE DATA STRUCTURE TO REPRESENT AN OPEN BIN 

struct Bin {
    std::vector<int> items;
    int size = 0;
};

#endif