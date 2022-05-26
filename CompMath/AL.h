#ifndef AL_H
#define AL_H
#include "model.h"

class MVN
{

public:
    double (MVN::*mean)();
    double (MVN::*kernel)();
    MVN(/*double (MVN::*kernel_func)() = nullptr, double (MVN::*mean_func)() = nullptr*/)
    {
        //if (kernel_func == nullptr)
        //    kernel = &MVN::ARD;
    };

    int const_mean(double rotate, double rayY, double rayZ)
    {
        return 0;
    }

    double ARD(double rotate_a, double rayY_a, double rayZ_a, double rotate_b, double rayY_b, double rayZ_b, double sigma = 1.5 )//vector& l = vector(3, 100.0))
    {   /* Automatic Relevance Determination */
        vector l;
        l.push_back(10.0);
        l.push_back(10.0);
        l.push_back(10.0);
        double sigma_suqared = std::pow(sigma, 2);
        double cumulative = 0;
        cumulative += std::pow( (rotate_a - rotate_b) / l[0], 2);
        cumulative += std::pow( (rayY_a - rayY_b) / l[1], 2);
        cumulative += std::pow( (rayZ_a - rayZ_b) / l[2], 2);
        return sigma_suqared * std::exp(-0.5 * cumulative);
    }

    Matrix get_C_block(Matrix& data_points)
    {
        int rows = data_points.size1();
        Matrix C_block(rows, rows);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < rows; j++)
                C_block(i, j) = ARD(data_points(i, 0), data_points(i, 1), data_points(i, 2),
                                    data_points(j, 0), data_points(j, 1), data_points(j, 2));
        return C_block;
    }

    Matrix get_K_block(Matrix& data_points, double rotate, double rayY, double rayZ)
    {
        Matrix K_block(1, data_points.size1());
        for (size_t i = 0; i < data_points.size1(); i++)
            K_block(0, i) = ARD(data_points(i, 0), data_points(i, 1), data_points(i, 2),
                                    rotate, rayY, rayZ);
        return K_block;
    }

    double get_K_scalar(double rotate, double rayY, double rayZ)
    {
        return ARD(rotate, rayY, rayZ, rotate, rayY, rayZ);
    }

};

/*
Cb ------ Kb
|		  |
|		  |  <----- Cov_matrix
Kb.T ---- Ks
*/

class BIRegressor
{
    MVN& prior;
    double sigma;
    double mu;
    double posterior_sigma;
    double posterior_mu;
    Matrix memory;
    Matrix prev_C_block;
    Matrix mu_x;
    Matrix Y;
    int aqcR;
    int aqcY;
    int aqcZ;

public:
    Process model;

    BIRegressor(MVN& prior);

    double UCB(Matrix& data_point, double beta = 0.001);

    double predict(double rotate, double rayY, double rayZ);

    vector acquisition(/*int step_rotate = 1, int step_rayY = 1, int step_rayZ = 1*/);

    vector step();

};

#endif // AL_H
