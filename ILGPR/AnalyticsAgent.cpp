#include "AnalyticsAgent.hpp"

AnalyticsAgent::AnalyticsAgent(std::vector<std::string> IPsAndPorts) {

    for (int i = 0; i < IPsAndPorts.size(); i++) {
        settings.hosts.push_back(IPsAndPorts.at(i));
    }
    settings.type = madara::transport::UDP;
    //settings.add_receive_filter(&dataFilter);
    dataFilter = new DatumReceivedFilter(&ILGPRs);
    //dataFilter = std::unique_ptr<DatumReceivedFilter>(new DatumReceivedFilter(&ILGPRs));
    settings.add_receive_filter(dataFilter);
    knowledge = madara::knowledge::KnowledgeBase(KNOWLEDGEBASE_NAME,settings);
    knowledge.set(".id",madara::knowledge::KnowledgeRecord::Integer(999));
}

AnalyticsAgent::~AnalyticsAgent() {

}
