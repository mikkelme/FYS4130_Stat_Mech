#include "system.h"


System::System(string system_name_input, size_t dim_input, size_t L_input, size_t q_input, double T_input, size_t cluster_per_cycle_input, size_t MC_equil_input, size_t MC_measure_input, size_t bins_input)
{
    /* Constructor */

    srand(4130); // seed
    // srand((unsigned) time(NULL)); // RN seed

    //
    system_name = system_name_input;
    dim = dim_input;
    L = L_input;
    N = pow(L_input, dim);
    num_neighbours = pow(2, dim);
    T = T_input;
    pconnect = 1 - exp(-1/T);  

    q = q_input;
    M = (size_t *)malloc(q * sizeof(size_t));
    W_M = (complex<double> *)malloc(q * sizeof(complex<double>));

    cluster_per_cycle = cluster_per_cycle_input;
    MC_equil = MC_equil_input;
    MC_measure = MC_measure_input;
    bins = bins_input;

    for (size_t s = 0; s < q; s++)
    {
        M[s] = 0;
        W_M[s] = complex<double>(cos(2*M_PI * s/q), sin(2*M_PI * s/q));
    }
    

    if (dim == 1){ // 1D Chain
        lattice = (size_t *)malloc(L * sizeof(size_t));
    }

    else if (dim == 2){ // 2D square lattice
        lattice = (size_t *)malloc(L * L * sizeof(size_t)); 
    }

    else {
        cout << dim << endl;
        cout << "Lattice dimension must be either 1 or 2" << endl;
        exit(1);
    }

}


void System::Intialize(string mode){
    /*  Intialize system: either ranbdom or 
        all aligned in some direction  */

    if (mode == "RN") // Random configuratyion 
    {
        for (size_t i = 0; i < N; i++)
        {
            lattice[i] =  rand()%3; // random spin {0, 1, 2}
            M[lattice[i]] += 1;
        }
    }

    if (mode == "0" or mode == "1" or mode == "2"){ // All alligned
        int start_val = stoi(mode);
        for (size_t i = 0; i < N; i++)
        {
            lattice[i] = start_val;
        }
    M[start_val] = N;
    }

}


int System::get_neighbour(size_t s, size_t dir){
    /* Returns lattice index for the desired neighbour to spin s */

    size_t x = s%L;
    size_t y = s/(int) L;

    switch(dir)
    {
        case RIGHT: return (s+1)%L + y*L;
        case LEFT:  return (s-1+L)%L + y*L;
        case UP:    return x + (y+1)%L * L;
        case DOWN:  return x + (y-1+L)%L * L;
    }

    cout << "Illegal direction passed in System::Neighbours" << endl;
    exit(1);
}

void System::FlipandBuildFrom(size_t s){
    /*  Iterative algorithm to develop clusters 
        for the wolf algorithm  */

    size_t old_state = lattice[s];
    size_t new_state = (lattice[s] + 1)%q;

    // Flip spin
    lattice[s] = new_state;

    // Update spin count
    M[old_state] -= 1; M[new_state] += 1;


    size_t j;
    for (size_t dir = 0; dir < num_neighbours; dir++) // Go through neighbours
    {
        j = get_neighbour(s, dir);
        if (lattice[j] == old_state) {
            if (rand()/(RAND_MAX + 1.) < pconnect ){ FlipandBuildFrom(j); }
        }
    }
}

void System::Equilibriate(){
    /* Equilibriate the system before sampling */

    cout << "Equilibrating | " << MC_equil << " MC cycles" << endl;  
    pconnect = 1 - exp(-1 / T);
    for (size_t MC_cycle = 0; MC_cycle < MC_equil; MC_cycle++){
        for (size_t i = 0; i < cluster_per_cycle; i++)
        {
            FlipandBuildFrom(rand()%N); // Flip and build from random site
        }
    }
}

complex<double> System::get_m()
{
    complex<double> tmp_m(0., 0.);
    for (size_t s = 0; s < q; s++) // Calculate order parameter
    {
        tmp_m += W_M[s] * (double) M[s];
    }
    return tmp_m;
}

void System::Sample() {
    /* Sample system using the wolf algorithm */

    pconnect = 1 - exp(-1 / T);
    complex<double> m_cycle, m0conj, m0conj_cycle;
    double m1_cycle, m2_cycle;

    complex<double> *m_r =          (complex<double> *)calloc(N, sizeof(complex<double>));
    complex<double> *m_0m_r =       (complex<double> *)calloc(N, sizeof(complex<double>));
    complex<double> *m_r_cycle =    (complex<double> *)calloc(N, sizeof(complex<double>));
    C_r = (complex<double> *)calloc(N, sizeof(complex<double>));

    m = complex<double> (0., 0.); m0conj = complex<double>(0., 0.);
    m1 = 0; m2 = 0; m4 = 0;


    cout << "\33[2K\rSampling (MC cycles = " << MC_measure <<") | T = " << T << ", bin = " << 0 << "/" << bins << flush;
    for (size_t bin = 0; bin < bins; bin++)
    {
        for (size_t MC_cycle = 0; MC_cycle < MC_measure; MC_cycle++)
        {
            for (size_t i = 0; i < cluster_per_cycle; i++)
            {
                FlipandBuildFrom(rand() % N);
            } // Flip and build from random site
            m_cycle = get_m();
            m_cycle /= N;

            m1_cycle = abs(m_cycle);
            m2_cycle = m1_cycle * m1_cycle;

            m += m_cycle;              //  m
            m1 += m1_cycle;            // |m|
            m2 += m2_cycle;            // |m|^2
            m4 += m2_cycle * m2_cycle; // |m|^4

            m0conj_cycle = conj(W_M[lattice[0]]);
            m0conj += m0conj_cycle;
            for (size_t r = 0; r < N; r++)
            {
                m_r_cycle[r] = W_M[lattice[r]];
                m_r[r] += m_r_cycle[r];
                m_0m_r[r] += m0conj_cycle * m_r_cycle[r];
            }
        } // end of MC_measure-loop

        m /= MC_measure;
        m0conj /= MC_measure;
        m1 /= MC_measure;
        m2 /= MC_measure;
        m4 /= MC_measure;

        for (size_t r = 0; r < N; r++)
        {
            m_0m_r[r] /= MC_measure;
            m_r[r] /= MC_measure;
            C_r[r] = m_0m_r[r] - m0conj * m_r[r];
        }
    cout << "\33[2K\rSampling (MC cycles = " << MC_measure << ") | T = " << T << ", bin = " << bin + 1 << "/" << bins << flush;
    Write_to_file(0); // Append results

    } // end of bin-loop
}



void System::Temp_sample(double start_temp, double end_temp, double dt, size_t equil_each_temp){
    /* Sample system for different temperatures */

    double temp = start_temp;
    double eps = 0.001;
    Write_to_file(1); // Create outfile

    cout << "Temperature sampling | T in [" << start_temp << ", " << end_temp << "], dt = " << dt << endl;
    while (temp <= end_temp + eps)
    {
        if (equil_each_temp){
            cout << endl;
            Equilibriate();}
            
        System::T = temp;
        Sample();
        temp += dt;
    }
    cout << endl;
    
}

void System::Write_to_file(size_t create_file){
    /* Write current results to file */

    ofstream outfile;
    string output_file = system_name + ".txt";

    if (create_file)
    {
        outfile.open(output_file, ios::out);
        outfile << setw(15) << setprecision(8) << "dim " << dim << endl;
        outfile << setw(15) << setprecision(8) << "L " << L << endl;
        outfile << setw(15) << setprecision(8) << "q " << q << endl;
        outfile << setw(15) << setprecision(8) << "cluster_per_cycle " << cluster_per_cycle << endl;
        outfile << setw(15) << setprecision(8) << "MC_equil " << MC_equil << endl;
        outfile << setw(15) << setprecision(8) << "MC_measure " << MC_measure << endl;
        outfile << setw(15) << setprecision(8) << "bins " << bins << endl;

        outfile << setw(15) << setprecision(8) << "T";
        outfile << setw(15) << setprecision(8) << "m_real";
        outfile << setw(15) << setprecision(8) << "m1";
        outfile << setw(15) << setprecision(8) << "m2";
        outfile << setw(15) << setprecision(8) << "m4";
        outfile << setw(15) << setprecision(8) << "C_r.real";
    }

    else
    {
        outfile.open(output_file, ios::out | ios::app);
        outfile << setw(15) << setprecision(8) << T;
        outfile << setw(15) << setprecision(8) << m.real();
        outfile << setw(15) << setprecision(8) << m1;
        outfile << setw(15) << setprecision(8) << m2;
        outfile << setw(15) << setprecision(8) << m4;
        for (size_t r = 0; r < N; r++)
        {
            outfile << setw(15) << setprecision(8) << C_r[r].real();
        }
    }

    outfile << endl;
    outfile.close();
}


void System::print_lattice(){ 
    /* Print lattice: 1D or 2D */

    if (dim == 1){
        for (size_t i = 0; i < N; i++)
        {
            cout << lattice[i] << " ";
        }
        cout << endl;
    }

    if (dim == 2){
        for (size_t i = 0; i < L; i++){
            for (size_t j = 0; j < L; j++){
                cout << lattice[i*L + j] << " ";
            }
        cout << endl;
        }
        
    }
}
