#include "system.h"


int main(){
    string system_name = "2D_chain_XXX"; 
    size_t dim = 2;
    size_t L = 32;
    size_t q = 3;

    double T = 0.9900; // [J]
    double end_temp = 1.0000; 
    double dt = 0.00025;
    
    size_t cluster_per_cycle = 1;
    size_t MC_equil = 1e6;
    size_t MC_measure = 1e6;
    size_t bins = 10;

    System system = System(system_name, dim, L, q, T, cluster_per_cycle, MC_equil, MC_measure, bins);
    system.Intialize("RN");
    system.Equilibriate();
    system.Temp_sample(T, end_temp, dt, 0);

    return 0;
}
 

// Commands to run on ssh server
// CUDA_VISIBLE_DEVICES=x make
// CUDA_VISIBLE_DEVICES=x  nohup make > out1.log 2> error1.log &
// ps -u mikkelme | grep -i python
// nvidia-smi
