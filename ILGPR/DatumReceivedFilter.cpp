#include "DatumReceivedFilter.hpp"

DatumReceivedFilter::DatumReceivedFilter(std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs_In) {
    ILGPRs = ILGPRs_In;

    datumTypeRegex = new std::regex("");

}

void DatumReceivedFilter::filter(madara::knowledge::KnowledgeMap &records,
                                 const madara::transport::TransportContext &transport_context,
                                 madara::knowledge::Variables &vars) {
    ////////////////////////////////////////
    //identify which ILGPR based on the sensor type received
    //for (recordIterator iter = records.begin(); iter != records.end(); iter++) {
    std::string typeString;
    SENSOR_TYPE type;
    double value;
    for (auto const &iter : records) {
        std::string key = iter.first;
        if (key.find("environmentalData")) {
            if (key.find("type")) {
                typeString = iter.second.to_string();
                type = getDatumType(typeString);
            }
            if (key.find("value")) {
                value = iter.second.to_double();
            }
        }
    }

    // create new Datum object
    //Datum datum;

    // lock the ILGPR's queue
    //std::lock_guard<std::mutex> onReceiveLock(ILGPRs->find(/*[SENSOR_TYPE from Madara packet]*/).queueMutex);

    // push new Datum into ILGPR's queue

}

SENSOR_TYPE DatumReceivedFilter::getDatumType(std::string typeString) {
    if (typeString.compare(std::string("DO")) == 0) {
        return DO;
    }
    return UNKNOWN;
}
