#ifndef ILGPR_H
#define ILGPR_H

#include <vector>
#include <Eigen/Dense>
#include <mutex>

#include "LGP.hpp"
#include "Datum.hpp"

// there is one of these for each type of scalar sensor data



class ILGPR {

public:
    ILGPR();
    ~ILGPR();
    double predict(Eigen::VectorXd x);
    void newDatum(Datum datum);
    SENSOR_TYPE type;

private:

    std::mutex queueMutex; // mutex for the data queue
    std::vector<Datum> dataQueue; // queue of points that must be processed
    std::vector<LGP> LGPs;
    Eigen::MatrixXd predictionGrid;



};

#endif // ILGPR_H
