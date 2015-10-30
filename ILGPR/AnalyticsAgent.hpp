#ifndef ANALYTICS_AGENT_H
#define ANALYTICS_AGENT_H

#include <vector>
#include <string>
#include <unordered_map>
#include <memory>

#include "madara/threads/Threader.h"
#include "madara/transport/QoSTransportSettings.h"
#include "madara/knowledge/KnowledgeBase.h"
#include "madara/knowledge/containers/NativeDoubleVector.h"


namespace containers = madara::knowledge::containers;


#include "LGP.hpp"
#include "ILGPR.hpp"
#include "Datum.hpp"
#include "DatumReceivedFilter.hpp"

class AnalyticsAgent {

public:
    AnalyticsAgent(std::vector<std::string> IPsAndPorts);
    ~AnalyticsAgent();

private:

    madara::transport::QoSTransportSettings settings;
    DatumReceivedFilter *dataFilter; // is it possible to use C++11 pointers?
    madara::knowledge::KnowledgeBase knowledge;
    const std::string KNOWLEDGEBASE_NAME = std::string("analytics");

    std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> ILGPRs; // unordered map of ILGPR's for each SENSOR_TYPE that is received

};

#endif //ANALYTICS_AGENT_H
