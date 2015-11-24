/* This code converts a simple 1D Gaussian Process example in Matlab to C++.
    It is a single source code file that runs using an external GP library called libgp
    It does not show any visualization, works with fake training data and assumes known hyperparameters.
    (Literally the simplest version possible!)

    Author: Jaineel Dalal <jdalal@andrew.cmu.edu>
    Date: 11/22/2015

*/
#include <iostream>
#include <Eigen/Dense>
#include "libgp/include/gp.h"
#include "libgp/include/gp_utils.h"
#include <math.h>

int main()
{
    std::cout << "Hello World!" <<'\n';

    Eigen::VectorXd TestData_X(100);
    for(int i = 0; i < 100; i++)
    {
        TestData_X(i) = (-5 + 5*i)/100;
    }

    Eigen::MatrixXd TrainData(5, 2);
    Eigen::VectorXd TrainData_X(5, 1);
    Eigen::VectorXd TrainData_Y(5, 1);
    TrainData << -4, -2,
                 -3, 0,
                 -1, 1,
                  0, 2,
                  3, -1;
    TrainData_X = TrainData.block(0, 0, 5, 1);
    TrainData_Y = TrainData.block(0, 1, 5, 1);

    //int n = 4000, m = 1000;
    double tss = 0, error, f;

    // initialize Gaussian process for 1-D input using the squared exponential
    // covariance function with additive white noise.
    libgp::GaussianProcess gp(1, "CovSum ( CovSEiso, CovNoise)");

    // initialize hyper parameter vector
    std::cout<<gp.covf().get_param_dim()<<'\n';
    Eigen::VectorXd params(gp.covf().get_param_dim());
    params << 0.0, 0.0, -2.0;
    // set parameters of covariance function
    gp.covf().set_loghyper(params);

    // add training patterns
    for(int i = 0; i < 5; i++)
    {
        double x[] = {TrainData_X(i)};
        double y = TrainData_Y(i);
        gp.add_pattern(x, y);
    }

    // Total squared error
    for(int i = 0; i < 100; i++)
    {
        double x[] = {TestData_X(i)};
        f = gp.f(x);
        double y = 1;
        error = f - y;
        tss += error*error;
    }


    /*for(int i = 0; i < n; ++i)
    {
        double x[] = {drand48()*4-2, drand48()*4-2};
        y = libgp::Utils::hill(x[0], x[1]) + libgp::Utils::randn() * 0.1;
        gp.add_pattern(x, y);
    }
    // total squared error
    for(int i = 0; i < m; ++i)
    {
        double x[] = {drand48()*4-2, drand48()*4-2};
        f = gp.f(x);
        y = libgp::Utils::hill(x[0], x[1]);
        error = f - y;
        tss += error*error;
    }
    std::cout << "mse = " << tss/m << std::endl;
*/
    std::cout<<"mse = "<<tss/100<<'\n';

    return 0;
}

