#include <random>
#include <cstdlib>
#include <ctime>
#include "./random_utilities.h"

double* generate_obs_coordinates(double* fig_center, double radius){
    double* set = new double [3];
    const char sign[2] = {-1, 1};
    
    std::random_device rd;
    std::default_random_engine eng(rd());
    srand(time(NULL));
    
    for (int i = 0; i < 3; i++){
        int MIN_ort = 100 + int(fig_center[i]);
        std::uniform_real_distribution<double> distr(MIN_ort, radius);
        int RandIndex = rand() % 2;
        set[i] = distr(eng) * sign[RandIndex];
    }
    return set;
}