#include "AL.h"

double BIRegressor::UCB(Matrix &data_point, double beta)
{   /* upper confidence bound */
    Matrix C_block = prior.get_C_block(memory);
    Matrix K_block = prior.get_K_block(memory,
                                       data_point(0, 0), data_point(0, 1), data_point(0, 2));
    double K_scalar = prior.get_K_scalar(data_point(0, 0), data_point(0, 1), data_point(0, 2));

    Matrix C_inverse = gjinverse(C_block, true);
    Matrix prod1 = prod(K_block, C_inverse);
    Matrix trans = transpose(K_block);
    Matrix prod2 = prod(prod1, trans);
    double variance = K_scalar - prod2(0, 0);

    return prior.const_mean(data_point(0, 0), data_point(0, 1), data_point(0, 2))
            + std::pow(beta, 0.5) * std::pow(variance, 0.5);
};

double BIRegressor::predict(double rotate, double rayY, double rayZ)
{
    int rows = memory.size1();
    Matrix C_block = prior.get_C_block(memory);
    Matrix K_block = prior.get_K_block(memory, rotate, rayY, rayZ);

    Matrix Y(rows, 1);
    for (int i = 0; i < rows; i++)
        Y(i, 0) = memory(i, 3);

    Matrix mu_x(rows, 1);
    for (int i = 0; i < rows; i++)
        mu_x(i, 0) = prior.const_mean(memory(i, 0), memory(i, 1), memory(i, 2));
    Matrix C_inverse = gjinverse(C_block, true);
    Matrix prod1 = prod(K_block, C_inverse);
    Matrix sub1 = Y - mu_x;
    Matrix ans = prod(prod1, sub1);
    return ans(0, 0) + prior.const_mean(rotate, rayY, rayZ); // fix
}

vector BIRegressor::acquisition(/*int step_rotate, int step_rayY, int step_rayZ*/)
{
//    vector query_point = {-1};
    vector query_point;
    RNGType rng( time(0) );
    boost::uniform_int<> one_to_turn( 1, 360 );
    boost::uniform_int<> one_to_five( -10, 10 );
    boost::variate_generator< RNGType, boost::uniform_int<> >
                    get_angle(rng, one_to_turn);
    boost::variate_generator< RNGType, boost::uniform_int<> >
                    get_ray_angle(rng, one_to_five);
    int rotate, rayY, rayZ;
    bool flag = true;
    while (flag)
    {
        rotate = get_angle();
        rayY = get_ray_angle();
        rayZ = get_ray_angle();
        flag = false;
        for (int i = 0; i < memory.size1(); i++)
            if (memory(i, 0) == rotate &&
                    memory(i, 1) == rayY &&
                    memory(i, 2) == rayZ)
                flag = true;
    }
//    qDebug("query point: %d %d %d", rotate, rayY, rayZ);
    query_point.push_back(rotate);
    query_point.push_back(rayY);
    query_point.push_back(rayZ);

    {
//    double score = -1;
//    for (int r = 0; r < 180; r += step_rotate)
//    {
//        for (int y = -3; y < 3; y += step_rayY)
//        {
//            for (int z = -3; z < 3; z += step_rayZ)
//            {
//                Matrix query(1, 3);
//                query(0, 0) = r; query(0, 1) = y; query(0, 2) = z;
//                double candidate = UCB(query);
//                if (candidate > score || query_point.size() == 1)
//                {
//                   query_point.clear();
//                   query_point.push_back(r);
//                   query_point.push_back(y);
//                   query_point.push_back(z);
//                }
//            }
//        }
//    }
    }
    return query_point;
}

vector BIRegressor::step()
{
//    if (memory.size1() != 0)
//    {
//        qDebug("memory matrix:");
//        for (int i = 0 ; i < memory.size1(); i++)
//        {
//            QString str = QString::number(i) + "(";
//            for (int j = 0; j < memory.size2(); j++)
//            {
//                str.append(QString::number(memory(i, j)) + " ");
//            }
//            str.append(')');
//            qDebug(str.toLatin1());
//        }
//        qDebug("C_block matrix:");
//        for (int i = 0 ; i < prev_C_block.size1(); i++)
//        {
//            QString str = QString::number(i) + "(";
//            for (int j = 0; j < prev_C_block.size2(); j++)
//            {
//                str.append(QString::number(prev_C_block(i, j)) + " ");
//            }
//            str.append(')');
//            qDebug(str.toLatin1());
//        }
//    }

    // acqusition
    vector answer;
    int rotate;
    int rayY;
    int rayZ;
    double y = 0;
//    int rotate = aqcR;
//    int rayY = aqcY;
//    int rayZ = aqcZ++;
//    double y = model.f(rotate, rayY, rayZ, true);//0;
    // Dumb acquisition
//    if (aqcZ == 10)
//    {
//        aqcY++;
//        aqcZ = 1;
//    }
//    if (aqcY == 10)
//    {
//        aqcR++;
//        aqcY = 1;
//    }
//    qDebug("R: %d Y: %d Z: %d", aqcR, aqcY, aqcZ);
    QDateTime dateTime1 = QDateTime::currentDateTime();
    while (y == 0)
    {
        answer = acquisition();
        rotate = answer[0];
        rayY = answer[1];
        rayZ = answer[2];
        y = model.f(rotate, rayY, rayZ, true);
    }
    qDebug("rot %d Y %d Z %d dist %f", rotate, rayY, rayZ, y);

    if (memory.size1() != 0)
    {
        int rows = memory.size1();
        Matrix K_block = prior.get_K_block(memory, rotate, rayY, rayZ);
        double K_scalar = prior.get_K_scalar(rotate, rayY, rayZ);

        mu_x(rows - 1, 0) = prior.const_mean(memory(rows - 1, 0),
                                             memory(rows - 1, 1),
                                             memory(rows - 1, 2));
        Y(rows - 1, 0) = memory(rows - 1, 3);

        QDateTime dateTime3 = QDateTime::currentDateTime();
        Matrix C_inverse = gjinverse(prev_C_block, true);
        QDateTime dateTime4 = QDateTime::currentDateTime();
        qint64 millisecondsDiffinv = dateTime3.msecsTo(dateTime4);
        qDebug("Inverse time: %lld", millisecondsDiffinv);

        QDateTime dateTime5 = QDateTime::currentDateTime();
        Matrix prod1 = prod(K_block, C_inverse);
        Matrix sub1 = Y - mu_x;
        Matrix ans1 = prod(prod1, sub1);
        posterior_mu = prior.const_mean(rotate, rayY, rayZ) + ans1(0, 0);
        Matrix trans = transpose(K_block);
        Matrix ans2 = prod(prod1, trans);
        posterior_sigma = K_scalar - ans2(0, 0);
        QDateTime dateTime6 = QDateTime::currentDateTime();
        qint64 millisecondsDiffprod = dateTime5.msecsTo(dateTime6);
        qDebug("Prod time: %lld", millisecondsDiffprod);

        // update to inferred parameters
        int C_size1 = prev_C_block.size1();
        int C_size2 = prev_C_block.size2();

        K_block.resize(1, K_block.size2() + 1, true);
        K_block(0, K_block.size2() - 1) = K_scalar;

        prev_C_block.resize(C_size1 + 1, C_size2 + 1, true);

        for (int i = 0; i < C_size1; i++)
            prev_C_block(i, C_size2) = K_block(0, i);
        for (int i = 0; i < C_size2 + 1; i++)
            prev_C_block(C_size1, i) = K_block(0, i);

        sigma = posterior_sigma;
        mu = posterior_mu;
        // Update memory
        memory.resize(rows + 1, 4, true);
        Y.resize(rows + 1, 1, true);
        mu_x.resize(rows + 1, 1, true);
        memory(rows, 0) = rotate;
        memory(rows, 1) = rayY;
        memory(rows, 2) = rayZ;
        memory(rows, 3) = y;
    }
    else
    {
        memory.resize(1, 4);
        memory(0, 0) = rotate;
        memory(0, 1) = rayY;
        memory(0, 2) = rayZ;
        memory(0, 3) = y;
        sigma = prior.get_K_scalar(rotate, rayY, rayZ);
        posterior_mu = 0;
        posterior_sigma = sigma;

        prev_C_block(0, 0) = sigma;

    }
    answer.push_back(y);
    QDateTime dateTime2 = QDateTime::currentDateTime();
    qint64 millisecondsDiff = dateTime1.msecsTo(dateTime2);
    qDebug("Step time: %lld", millisecondsDiff);
    return answer;
};

BIRegressor::BIRegressor(MVN& prior): prior(prior), sigma(0), mu(0)
{
    mu_x.resize(1, 1);
    Y.resize(1, 1);
    prev_C_block.resize(1, 1);
    aqcR = 1;
    aqcY = 1;
    aqcZ = 1;
};
