#include <iostream>
#include <random>
#include <cmath>
#include <array>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/multi_array.hpp>

namespace ublas = boost::numeric::ublas;
typedef boost::multi_array<double, 3> ThreeD_Matrix; // fix
typedef ublas::matrix<double> Matrix;
typedef ublas::vector<double> Vector;

struct MDP
{
    int T;                // Time points per trial
    std::array<Matrix, 2> V;     // Policies  
    std::array<std::array<Matrix, 4>, 3> A;// Likelihood
    std::array<std::vector<Matrix>, 2> B;// Transition
    std::array<Matrix, 3> C; // Preferences over outcomes
    std::array<std::vector<double>, 2> D;// Priors over initial states
    std::array<std::vector<double>, 2> d;// Learning priors over initial states
    
    std::array<double, 5> E;                // Habits (prior probability of choosing i-th policy)            
    
    std::array<std::array<Matrix, 4>, 3> a;
    std::array<double, 5> e;
          
    double eta;
    double omega;            // Action precision
    double alpha;             // EFE precision
    // Forgetting rate
    double beta;              // Learning rate
    int NumPolicies;
    int NumFactors;
    /*......*/

    /* LABELS field for visualization */
}

// CREATES AND RETURNS MODEL
MDP explore_exploit_model(int gen_model)
{
    int T = 3;

    // Specify prior probabilities about initial states in the generative process (D)
    std::array<std::vector<double>, 2> D;
    D[0] = { 1, 0 };        /*left better, right better*/
    D[1] = { 1, 0, 0, 0 }; /*start, hint, choose left, choose right*/
    // transpose ?

    // Specify prior beliefs about initial states in the generative model (d)
    std::array<std::vector<double>, 2> d;
    d[0] = { 0.25, 0.25 };        /*left better, right better*/
    d[1] = { 1, 0, 0, 0 };          /*start, hint, choose left, choose right*/

    // Specify the probabilities of outcomes given each state in the generative process (A) 
    unsigned int Ns[2] = { D[0].size(), D[1].size() }; /*Number of states in each factors: 2 4 */

    std::array < std::array < Matrix, 4>, 3> A; /*3 time points, Ns[1] = 4 */
    for (size_t i = 0; i < Ns[1]; i++)
    {
        Matrix tmp(3, 2);
        tmp(0, 0) = 1; tmp(0, 1) = 1; /*No hint*/
        tmp(1, 0) = 0; tmp(1, 1) = 0; /*Left hint*/
        tmp(2, 0) = 0; tmp(2, 1) = 0; /*Right hit*/
        A[0][i] = tmp;
    }

    double pHA = 1;
    A[0][1](0, 0) = 1;       A[0][1](0, 1) = 1;
    A[0][1](1, 0) = pHA;     A[0][1](1, 1) = 1 - pHA;
    A[0][1](2, 0) = 1 - pHA; A[0][1](2, 1) = pHA;

    for (size_t i = 0; i < 2; i++)
    {
        Matrix tmp(3, 2);
        tmp(0, 0) = 1; tmp(0, 1) = 1; /*No hint*/
        tmp(1, 0) = 0; tmp(1, 1) = 0; /*Left hint*/
        tmp(2, 0) = 0; tmp(2, 1) = 0; /*Right hit*/
        A[1][i] = tmp;
    }

    double pWin = 0.8;
    A[1][2](0, 0) = 1;        A[1][2](0, 1) = 1;
    A[1][2](1, 0) = 1 - pWin; A[1][2](1, 1) = pWin;
    A[1][2](2, 0) = pWin;     A[1][2](2, 1) = 1 - pWin;

    A[1][3](0, 0) = 1;        A[1][3](0, 1) = 1;
    A[1][3](1, 0) = pWin;     A[1][3](1, 1) = 1 - pWin;
    A[1][3](2, 0) = 1 - pWin; A[1][3](2, 1) = pWin;

    for (size_t i = 0; i < Ns[1]; i++)
    {
        Matrix tmp(4, 2);
        tmp(0, 0) = i == 0 ? 1 : 0; tmp(0, 1) = i == 0 ? 1 : 0; /*No hint*/
        tmp(1, 0) = i == 1 ? 1 : 0; tmp(1, 1) = i == 1 ? 1 : 0; /*Left hint*/
        tmp(2, 0) = i == 2 ? 1 : 0; tmp(2, 1) = i == 2 ? 1 : 0; /*Right hit*/
        tmp(3, 0) = i == 3 ? 1 : 0; tmp(3, 1) = i == 3 ? 1 : 0; /*Right hit*/
        A[2][i] = tmp;
    }

    std::array < std::array < Matrix, 4>, 3> a;

    for (size_t t = 0; t < 3; t++)
    {
        for (size_t i = 0; i < 4; i++)
        {
            a[t][i] = A[t][i];
            size_t border = 3;
            if (t == 2)
                border++;
            for (size_t j = 0; j < border; j++)
            {
                a[t][i](j, 0) *= 200;
                a[t][i](j, 1) *= 200;
            }
        }
    }

    a[0][1](0, 0) = 0;    a[0][1](0, 1) = 0;
    a[0][1](1, 0) = 0.25; a[0][1](1, 1) = 0.25;
    a[0][1](2, 0) = 0.25; a[0][1](2, 1) = 0.25;


    std::array<std::vector<Matrix>, 2> B;
    for (size_t i = 0; i < 2; i++)
    {
        if (i == 0)
        {
            Matrix tmp(2, 2);
            tmp(0, 0) = 1; tmp(0, 1) = 0;
            tmp(1, 0) = 0; tmp(1, 1) = 1;

            B[i].push_back(tmp);
        }
        else
        {
            for (size_t j = 0; j < 4; j++)
            {
                Matrix tmp(4, 4);
                tmp(0, 0) = i == 0 ? 1 : 0; tmp(0, 1) = i == 0 ? 1 : 0;
                tmp(1, 0) = i == 1 ? 1 : 0; tmp(1, 1) = i == 1 ? 1 : 0;
                tmp(2, 0) = i == 2 ? 1 : 0; tmp(2, 1) = i == 2 ? 1 : 0;
                tmp(2, 0) = i == 3 ? 1 : 0; tmp(2, 1) = i == 3 ? 1 : 0;
                B[i].push_back(tmp);
            }

        }
    }

    int No[] = { A[0][0].size1(), A[1][0].size1(), A[2][0].size1() };

    std::array<Matrix, 3> C; /*Hints, Win/Losses, Observed behaviors*/
    for (size_t i = 0; i < 3; i++)
    {
        Matrix tmp(No[i], 3);
        for (size_t j = 0; j < No[i]; j++)
        {
            tmp(j, 0) = 0; tmp(j, 1) = 0;
        }
        C[i] = tmp;
    }

    double la = 1;
    double rs = 4;

    C[1](0, 0) = 0; C[1](0, 1) = 0; C[1](0, 2) = 0;
    C[1](1, 0) = 0; C[1](1, 1) = -la; C[1](1, 2) = -la;
    C[1](2, 0) = 0; C[1](2, 1) = rs; C[1](2, 2) = rs / 2;

    int NumPolicies = 5;
    int NumFactors = 2;

    std::array<Matrix, 2> V;
    for (size_t i = 0; i < 2; i++)
    {
        Matrix tmp(2, 5);
        tmp(0, 0) = 1; tmp(0, 1) = i == 1 ? 2 : 1; tmp(0, 2) = i == 1 ? 2 : 1; tmp(0, 3) = i == 1 ? 3 : 1; tmp(0, 4) = i == 1 ? 4 : 1;
        tmp(1, 0) = 1; tmp(1, 1) = i == 1 ? 3 : 1; tmp(1, 2) = i == 1 ? 4 : 1; tmp(0, 3) = 1;              tmp(1, 4) = 1;
        V[i] = tmp;
    }

    std::array<double, 5> E = { 1, 1, 1, 1, 1 };
    std::array<double, 5> e = { 1, 1, 1, 1, 1 };

    double eta = 1, omega = 1, beta = 1, alpha = 32, erp = 1, tau = 12;
    MDP model;
    model.T = T;// Number of time steps
    model.V = V;// allowable(deep) policies
    
    model.A = A; // state - outcome mapping
    /*for (size_t i = 0; i < 3; i++)
        for (size_t j = 0; j < 4; j++)
            model.A[i][j] = A[i][j];*/ 
    model.B = B;// transition probabilities
    /*model.B[0] = B[0];
    model.B[1] = B[1];*/
    model.C = C; // preferred states
    /*for (size_t i = 0; i < 3; i++)
        model.C[i] = C[i];*/
    model.D = D; // priors over initial states
    /*model.D[0] = D[0]; 
    model.D[1] = D[1];*/
    model.d = d; // enable learning priors over initial states
    /*model.d[0] = d[0];
    model.d[1] = d[1];*/

    if (gen_model == 1)
        model.E = E; // prior over policies
        /*for (size_t i = 0; i < 5; i++)
            model.E[i] = E[i];*/
    else if (gen_model == 2)
    {
        model.a = a; // enable learning state - outcome mappings
        /*for (size_t i = 0; i < 3; i++)
            for (size_t j = 0; j < 4; j++)
                model.a[i][j] = a[i][j];*/ 
        model.e = e; // enable learning of prior over policies
                /*for (size_t i = 0; i < 5; i++)
            model.e[i] = e[i];*/ 
    }

    model.eta = eta;// learning rate
    model.omega = omega;// forgetting rate
    model.alpha = alpha;// action precision
    model.beta = beta;// expected free energy precision

    // respecify for use in inversion script(specific to this tutorial example)
    model.NumPolicies = NumPolicies;// Number of policies
    model.NumFactors = NumFactors;// Number of state factors
    return model;
}

int main()
{
    MDP model = explore_exploit_model(1);
    auto A = model.A;
    auto B = model.B;
    auto C = model.C;
    auto D = model.D;
    auto T = model.T;
    auto V = model.V;
    auto beta = model.beta;
    auto alpha = model.alpha;
    auto eta = model.eta;
    auto omega = model.omega;

    /* No need, already normalized
    A = col_norm(A);
    B = col_norm(B);
    D = col_norm(D);
    */
    int NumPolicies = model.NumPolicies, NumFactors = model.NumFactors;

    // if we use learning priors over initial states (we do)
    std::array<std::vector<double>, 2> d_prior;
    std::array<std::vector<double>, 2> d_complexity;

    for (size_t factor = 0; factor < 2; factor++)
    {
        d_prior[factor] = model.d[factor];
        //d_complexity[factor] = spm_wnorm(d_prior[factor]);
        double sum = 0;
        std::vector<double> tmp = d_prior[factor];
        for (size_t i = 0; i < d_prior[factor].size(); i++)
        {
            sum += d_prior[factor][i];
            tmp[i] = 1 / tmp[i];
        }
        sum = 1 / sum;
        for (size_t i = 0; i < tmp.size(); i++)
        {
            tmp[i] = sum - tmp[i];
        }
        d_complexity[factor] = tmp;
    }

    /* if isfield(MDP,'a')
    // complexity of a maxtrix concentration parameters
    for modality = 1:numel(MDP.a)
        a_prior{modality} = MDP.a{modality};
        a_complexity{modality} = spm_wnorm(a_prior{modality}).*(a_prior{modality} > 0);
    end
    end*/

    // Normalise matrices before model inversion / inference
    auto a = A;
    auto b = B;
    
    /*
    % normalize C and transform into log probability
    for ii = 1:numel(C)
        C{ii} = MDP.C{ii} + 1/32;
        for t = 1:T
            C{ii}(:,t) = nat_log(exp(C{ii}(:,t))/sum(exp(C{ii}(:,t)))); ????
        end 
    end 

    % normalize D vector
    if isfield(MDP,'d')
        d = col_norm(MDP.d);
    else 
        d = col_norm(MDP.D);
    end 

    % normalize E vector
    if isfield(MDP,'e')
        E = MDP.e;
        E = E./sum(E);
    elseif isfield(MDP,'E')
        E = MDP.E;
        E = E./sum(E);
    else
        E = col_norm(ones(NumPolicies,1));
        E = E./sum(E);
    end

    % Initialize variables
    %--------------------------------------------------------------------------

    % numbers of transitions, policies and states
    NumModalities = numel(a);                    % number of outcome factors
    NumFactors = numel(d);                       % number of hidden state factors
    NumPolicies = size(V,2);                     % number of allowable policies
    for factor = 1:NumFactors
        NumStates(factor) = size(b{factor},1);   % number of hidden states
        NumControllable_transitions(factor) = size(b{factor},3); % number of hidden controllable hidden states for each factor (number of B matrices)
    end

    % initialize the approximate posterior over states conditioned on policies
    % for each factor as a flat distribution over states at each time point
    for policy = 1:NumPolicies
        for factor = 1:NumFactors
            NumStates(factor) = length(D{factor}); % number of states in each hidden state factor
            state_posterior{factor} = ones(NumStates(factor),T,policy)/NumStates(factor); 
        end  
    end 

    % initialize the approximate posterior over policies as a flat distribution 
    % over policies at each time point
    policy_posteriors = ones(NumPolicies,T)/NumPolicies; 

    % initialize posterior over actions
    chosen_action = zeros(ndims(B),T-1);
    
    % if there is only one policy
    for factors = 1:NumFactors 
        if NumControllable_transitions(factors) == 1
            chosen_action(factors,:) = ones(1,T-1);
        end
    end
    MDP.chosen_action = chosen_action;

    % initialize expected free energy precision (beta)
    posterior_beta = 1;
    gamma(1) = 1/posterior_beta; % expected free energy precision
    
    % message passing variables
    TimeConst = 4; % time constant for gradient descent
    NumIterations  = 16; % number of message passing iterations
*/

    return 0;
}


//Matrix col_norm(Matrix matrix)
//{
//    // 
//    return matrix;
//};

double nat_log(double x)
{
    return std::log(x + std::exp(-16));
}

Matrix spm_wnorm(Matrix A)
{
    for (size_t i = 0; i < A.size1(); i++)
    {
        for (size_t j = 0; j < A.size2(); j++)
        {
            A(i, j) += std::exp(-16);
        }
    }
    return A;
}