cmake_minimum_required(VERSION 2.6)
 
Project(dc)

# command-line parsing using the boost program-options library
find_package(Boost COMPONENTS program_options REQUIRED) 

 
set(DC_SRC_FILES
    dc.cpp
    eigen.cpp
    octree.cpp
    # SOGReader.cpp
)

set(DC_INCLUDE_FILES
    eigen.hpp
    GeoCommon.hpp
    HashMap.hpp
    intersection.hpp
    ModelReader.hpp
    octree.hpp
    PLYReader.hpp
    PLYWriter.hpp
    # SOGReader.hpp
)

set(CMAKE_CXX_FLAGS "-fpermissive") 

ADD_EXECUTABLE(dualcontour ${DC_SRC_FILES})

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
    target_link_libraries(dualcontour ${Boost_LIBRARIES})                                                                                                                                                                                                                            
endif()

