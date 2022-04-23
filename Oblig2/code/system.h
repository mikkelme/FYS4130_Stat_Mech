#ifndef SYSTEM_H
#define SYSTEM_H

#include <iostream>
#include <fstream>
#include <iomanip> 
#include <complex>
#include <random>
#include <math.h>
using namespace std;

class System
{
public:

    // Class variables
    string system_name;
    size_t dim;             // dimension
    size_t N;               // Number of spin sites
    size_t L;               // Lattice length 
    size_t *lattice;        // Lattice 
    size_t q;               // number of spin states
    size_t *M;              // Spin count
    complex<double> *W_M;   // Weights for m_j
    complex<double> *C_r;   // Correlation function
    double T;               // Temperature
    double pconnect;        // Probability for equal sizes to create bonds

    size_t cluster_per_cycle;
    size_t MC_equil;
    size_t MC_measure;
    size_t bins;

    // Measurable values
    complex<double> m;
    double m1, m2, m4;

    // Handle neighbours
    enum dirs {RIGHT, LEFT, UP, DOWN};
    int num_neighbours;

    // Class methods
    System(string system_name_input, size_t dim_input, size_t L_input, size_t q_input, double T_input, size_t cluster_per_cycle_input, size_t MC_equil_cycles_input, size_t MC_measure_cycles_input, size_t bins); // Constructor
    void Intialize(string mode);
    int get_neighbour(size_t s, size_t dir);
    complex<double> get_m();
    void FlipandBuildFrom(size_t s);
    void Equilibriate();
    void Sample();
    void Temp_sample(double start_temp, double end_temp, double dt, size_t equil_each_temp);
    void Write_to_file(size_t create_file);
    void print_lattice();
    };

#endif // SYSTEM_H

