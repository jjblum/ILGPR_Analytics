#ifndef DATUM_RECEIVED_FILTER_H
#define DATUM_RECEIVED_FILTER_H

#include "madara/filters/AggregateFilter.h"

#include "ILGPR.hpp"
#include "Datum.hpp"
#include <unordered_map>

class DatumReceivedFilter: public madara::filters::AggregateFilter {

public:
    DatumReceivedFilter(std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs);

private:

    void filter(madara::knowledge::KnowledgeMap &records,
                const madara::transport::TransportContext &transport_context,
                madara::knowledge::Variables &vars);
    ///////////////////////



};


#endif //DATUM_RECEIVED_FILTER_H
