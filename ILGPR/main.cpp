#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <stdexcept>
#include <vector>
#include <unordered_map>
#include <mutex>

#include "AnalyticsAgent.hpp"


const std::string ANALYTICS_IP_ADDRESS_AND_PORT = std::string("192.168.1.8:45000");
const std::string BOAT_IP_ADDRESS_AND_PORT = std::string("192.168.1.20:45000");


int main(int argc, char** argv)
{

    // gather Madara connection details from command line
    std::string ipAddress;
    std::string portString;
    int port;
    try {
        if (argc > 1) {
            ipAddress = std::string(argv[1]);
            std::cout << "IP Address = " << ipAddress.c_str() << std::endl;
        }
        else {
            throw std::invalid_argument("Must include IP and port for Madara");
        }

        if (argc > 2) {
            portString = std::string(argv[2]);
            port = boost::lexical_cast<int>(argv[2]);
            std::cout << "Port = " << port << std::endl;
        }
        else {
            throw std::invalid_argument("Must include IP and port for Madara");
        }
    }
    catch(std::invalid_argument &e) {
        std::cerr << e.what() << std::endl;
        return -1;
    }

    //////////////////////////////////////////////////////////////////////////////////
    std::vector<std::string> IPsAndPorts;
    IPsAndPorts.push_back(ANALYTICS_IP_ADDRESS_AND_PORT);
    IPsAndPorts.push_back(BOAT_IP_ADDRESS_AND_PORT);
    AnalyticsAgent analyticsAgent(IPsAndPorts);

    return 0;
}

