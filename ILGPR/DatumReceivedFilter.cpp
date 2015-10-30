#include "DatumReceivedFilter.hpp"

DatumReceivedFilter::DatumReceivedFilter(std::unordered_map<SENSOR_TYPE, ILGPR, std::hash<int>> *ILGPRs_In) {
    ILGPRs = ILGPRs_In;

    datumTypeRegex = new std::regex("");

}

void DatumReceivedFilter::setKB(madara::knowledge::KnowledgeBase *knowledgeIn) {
    knowledge = knowledgeIn;
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

    /*
     * STRATEGY: IF WE CAN GET THE PREFIX FOR THE SENDING DEVICE, YOU CAN USE CONTAINERS TO EASILY HANDLE THE REST OF THE DATA
     *           (assuming that all the info from the packet made it...
    */


    // change to madara::utility::begins_with
    for (auto const &iter : records) {
        std::string key = iter.first;
        if (key.find("environmentalData")) { // we at least know it is some environmentalData packet
            if (key.find("type")) { // type of sensor
                typeString = iter.second.to_string();
                type = getDatumType(typeString);
            }
            if (key.find("value")) { // value of the sensor --> this may or may not be a scalar, need to check for "size"
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
