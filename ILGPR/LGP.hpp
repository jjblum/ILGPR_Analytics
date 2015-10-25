#ifndef LGP_H
#define LGP_H

#include <Eigen/Dense>
#include "gp/gp.h"
#include "gp/gp_utils.h"

const int N_MAX = 100; // maximum allowable training points

class LGP {

public:

    LGP();
    ~LGP();
    void updateCenter();
    void insertNewData(); // check that N <= N_MAX, concatenate X and y, perform ILGP updates
    void choleskyUpdate();
    double newWeight(Eigen::VectorXd);
    int size();




private:

    Eigen::MatrixXd W; // kernel width matrix
    Eigen::MatrixXd X; // locations of training data
    Eigen::VectorXd y; // values of training data
    Eigen::MatrixXd L; // lower triangular, Cholesky decomposition
    Eigen::VectorXd center; // location of the LGP
    int N; // current training points

};


#endif //LGP_H
