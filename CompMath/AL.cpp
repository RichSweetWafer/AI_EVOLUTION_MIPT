#include "AL.h"

double BIRegressor::UCB(Matrix &data_point, double beta)
{   /* upper confidence bound */
    Matrix K_block = prior.get_K_block(memory,
                                       data_point(0, 0), data_point(0, 1), data_point(0, 2));
    double K_scalar = prior.get_K_scalar(data_point(0, 0), data_point(0, 1), data_point(0, 2));

    Matrix prod1 = prod(K_block, C_inverse);
    Matrix trans = transpose(K_block);
    Matrix prod2 = prod(prod1, trans);
    double variance = K_scalar - prod2(0, 0);

    return std::pow(beta, 0.5) * std::pow(variance, 0.5);
};

double BIRegressor::predict(double rotate, double rayY, double rayZ)
{
    Matrix K_block = prior.get_K_block(memory, rotate, rayY, rayZ);
    Matrix prod1 = prod(K_block, C_inverse);
    Matrix sub1 = Y - mu_x;
    Matrix ans = prod(prod1, sub1);
    return ans(0, 0); // fix
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
        int mem_size = memory.size1();
        for (int i = 0; i < mem_size; i++)
            if (memory(i, 0) == rotate &&
                    memory(i, 1) == rayY &&
                    memory(i, 2) == rayZ)
                flag = true;
    }
//    qDebug("query point: %d %d %d", rotate, rayY, rayZ);
    query_point.push_back(rotate);
    query_point.push_back(rayY);
    query_point.push_back(rayZ);

    return query_point;
}

vector BIRegressor::step()
{
    // Debug for matrices
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

    vector answer;
    int rotate = aqcR;
    int rayY = aqcY;
    int rayZ = aqcZ++;
    double y = model.f(rotate, rayY, rayZ, true);//0;
    answer.push_back(rotate);
    answer.push_back(rayY);
    answer.push_back(rayZ);
    answer.push_back(y);
    // Dumb acquisition
    if (aqcZ == 10)
    {
        aqcY++;
        aqcZ = 1;
    }
    if (aqcY == 10)
    {
        aqcR++;
        aqcY = 1;
    }
    // Normal acq
//    while (y == 0)
//    {
//        answer = acquisition();
//        rotate = answer[0];
//        rayY = answer[1];
//        rayZ = answer[2];
//        y = model.f(rotate, rayY, rayZ, true);
//    }
//    qDebug("rot %d Y %d Z %d dist %f", rotate, rayY, rayZ, y);

//    QDateTime dateTime1 = QDateTime::currentDateTime();
    if (memory.size1() != 0)
    {
        int rows = memory.size1();
        Matrix K_block = prior.get_K_block(memory, rotate, rayY, rayZ);
        double K_scalar = prior.get_K_scalar(rotate, rayY, rayZ);

        mu_x(rows - 1, 0) = 0;
        Y(rows - 1, 0) = memory(rows - 1, 3);

//        QDateTime dateTime3 = QDateTime::currentDateTime();
        int size = prev_C_block.size1();
        C_inverse.resize(size, size, true);
        for (int j = 0; j < size; j++)
            C_inverse(size - 1, j) = 0;
        for (int i = 0; i < size; i++)
            for (int j = 0; j <= i; j++)
                C_inverse(i, j) += invL(size -1 , i) * invL(size - 1, j);

//        qDebug("C_inverse matrix:");
//        for (int i = 0 ; i < C_inverse.size1(); i++)
//        {
//            QString str = QString::number(i) + "(";
//            for (int j = 0; j < C_inverse.size2(); j++)
//            {
//                str.append(QString::number(C_inverse(i, j)) + " ");
//            }
//            str.append(')');
//            qDebug(str.toLatin1());
//        }
//        QDateTime dateTime4 = QDateTime::currentDateTime();
//        qint64 millisecondsDiffinv = dateTime3.msecsTo(dateTime4);
//        qDebug("Inverse time: %lld", millisecondsDiffinv);

//        QDateTime dateTime5 = QDateTime::currentDateTime();
        Matrix prod1 = prod(K_block, C_inverse);
        Matrix sub1 = Y - mu_x;
        Matrix ans1 = prod(prod1, sub1);
        posterior_mu = ans1(0, 0);
        Matrix trans = transpose(K_block);
        Matrix ans2 = prod(prod1, trans);
        posterior_sigma = K_scalar - ans2(0, 0);
//        QDateTime dateTime6 = QDateTime::currentDateTime();
//        qint64 millisecondsDiffprod = dateTime5.msecsTo(dateTime6);
//        qDebug("Prod time: %lld", millisecondsDiffprod);

//        QDateTime dateTime7 = QDateTime::currentDateTime();
        // update to inferred parameters
        int C_size1 = prev_C_block.size1() + 1;

        K_block.resize(1, K_block.size2() + 1, true);
        K_block(0, K_block.size2() - 1) = K_scalar;

        prev_C_block.resize(C_size1, C_size1, true);
        for (int i = 0; i < C_size1; i++)
            prev_C_block(C_size1 - 1, i) = K_block(0, i);


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
//        QDateTime dateTime8 = QDateTime::currentDateTime();
//        qint64 millisecondsDiffinvasd = dateTime7.msecsTo(dateTime8);
//        qDebug("Resize time: %lld", millisecondsDiffinvasd);

        update_L();
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
        L(0, 0) = std::sqrt(sigma);
        invL(0, 0) = 1 / L(0, 0);
        C_inverse(0, 0) = std::pow(invL(0, 0), 2);

    }
//    QDateTime dateTime2 = QDateTime::currentDateTime();
//    qint64 millisecondsDiff = dateTime1.msecsTo(dateTime2);
//    qDebug("Step time: %lld", millisecondsDiff);
    return answer;
};

BIRegressor::BIRegressor(MVN& prior): prior(prior), sigma(0), mu(0)
{
    mu_x.resize(1, 1);
    Y.resize(1, 1);
    prev_C_block.resize(1, 1, false);
    C_inverse.resize(1, 1, false);
    L.resize(1, 1, false);
    invL.resize(1, 1, false);
    aqcR = 1;
    aqcY = 1;
    aqcZ = 1;
};

void BIRegressor::update_L() {
    double S1, S2 = 0;
    int size = prev_C_block.size1();
    L.resize(size, size, true);

    // first element in a row
    L(size - 1, 0) = prev_C_block(size - 1, 0) / L(0, 0);

    // elements in between
    for (int i = 1; i < size - 1; i++) {
        S1 = 0;
        for (int p = 0; p < i; p++)
            S1 += L(i,p) * L(size - 1,p);
        L(size - 1, i) = (prev_C_block(size - 1,i) - S1) / L(i,i);
    }

    // last element in a row
    for (int p = 0; p < size - 1; p++)
        S2 += L(size - 1, p) * L(size - 1, p);

    L(size - 1, size - 1) = (std::sqrt(prev_C_block(size - 1,size - 1) - S2));

    // For L^-1 this would be
    invL.resize(size, size, true);
    invL(size - 1, size - 1) = 1 / L(size - 1, size - 1);
    size--;
    for (int j = 0; j < size; j++)
    {
        invL(size, j) = - double(1 / L(size, size));
//        std::cout << "beforesum: " << size  << ' ' << j << "\n";
        double sum = 0;
        for (int k = 0; k < size; k++)
        {
        if (k < j || k > size)
            continue;
//        std::cout << "k: " << k << "\n";
        sum += invL(k, j)*L(size, k);
        }
        invL(size, j) *= sum;
    }
}
