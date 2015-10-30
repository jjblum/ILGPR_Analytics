#ifndef DATUM_RECEIVED_FILTER_H
#define DATUM_RECEIVED_FILTER_H

#include "madara/filters/AggregateFilter.h"

#include "ILGPR.hpp"
#include "Datum.hpp"
#include <unordered_map>
#include <mutex>
#include <regex>

#include "madara/knowledge/KnowledgeRecord.h"
#include "madara/knowledge/KnowledgeBase.h"


//typedef madara::knowledge::KnowledgeMap::iterator recordIterator;

class DatumReceivedFilter: public madara::filters::AggregateFilter {

public:
    DatumReceivedFilter(std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs_In);
    void setKB(madara::knowledge::KnowledgeBase *knowledgeIn);

private:

    madara::knowledge::KnowledgeBase *knowledge;
    std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs;
    SENSOR_TYPE getDatumType(std::string key);

    void filter(madara::knowledge::KnowledgeMap &records,
                const madara::transport::TransportContext &transport_context,
                madara::knowledge::Variables &vars);

    std::regex *datumTypeRegex;
    std::smatch match;


};


#endif //DATUM_RECEIVED_FILTER_H
