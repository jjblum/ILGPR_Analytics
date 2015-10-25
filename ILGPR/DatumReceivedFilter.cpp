#include "DatumReceivedFilter.hpp"

DatumReceivedFilter::DatumReceivedFilter(std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs) {

}

void DatumReceivedFilter::filter(madara::knowledge::KnowledgeMap &records,
                                 const madara::transport::TransportContext &transport_context,
                                 madara::knowledge::Variables &vars) {
    ////////////////////////////////////////



}
