#ifndef DATUM_H
#define DATUM_H

#include "Eigen/Dense"

const int SENSOR_TYPE_COUNT = 5;
enum SENSOR_TYPE {
    EC,
    TEMP,
    DO,
    DEPTH,
    FLOW
};

class Datum {

public:
    Datum(SENSOR_TYPE type, double valueIn, Eigen::VectorXd locationIn, long uniqueID);
    ~Datum();

private:
    SENSOR_TYPE type;
    Eigen::VectorXd location;
    double value;
    double weight; // "distance" from the LGP that this datum is assigned to. Largest distance datum is removed first
    long uniqueID;
};

#endif //DATUM_H
