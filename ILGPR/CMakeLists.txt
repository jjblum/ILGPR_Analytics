project(ILGPR)
cmake_minimum_required(VERSION 2.8)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})

set(SOURCE
    ${SOURCE}
    ${CMAKE_CURRENT_SOURCE_DIR}/LGP.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Datum.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ILGPR.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticsAgent.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DatumReceivedFilter.cpp
)

set(HEADERS
    ${HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/LGP.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Datum.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/ILGPR.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/AnalyticsAgent.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DatumReceivedFilter.hpp
)

INCLUDE_DIRECTORIES(/usr/local/include)
LINK_DIRECTORIES(/usr/local/lib)
TARGET_LINK_LIBRARIES(ILGPR ACE)
TARGET_LINK_LIBRARIES(ILGPR MADARA)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
if(COMPILER_SUPPORTS_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
