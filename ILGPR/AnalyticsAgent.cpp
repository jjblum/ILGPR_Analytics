#include "AnalyticsAgent.hpp"

AnalyticsAgent::AnalyticsAgent(std::vector<std::string> IPsAndPorts) {

    for (int i = 0; i < IPsAndPorts.size(); i++) {
        settings.hosts.push_back(IPsAndPorts.at(i));
    }
    settings.type = madara::transport::UDP;
    dataFilter = new DatumReceivedFilter(&ILGPRs);
    settings.add_receive_filter(dataFilter);
    madara::transport::QoSTransportSettings::domains = std::string("environmentalData"); // --> ???

    knowledge = madara::knowledge::KnowledgeBase(KNOWLEDGEBASE_NAME,settings);
    knowledge.set(".id",madara::knowledge::KnowledgeRecord::Integer(999));

    dataFilter->setKB(&knowledge);
}

AnalyticsAgent::~AnalyticsAgent() {
    delete dataFilter;
    dataFilter = nullptr;
}
