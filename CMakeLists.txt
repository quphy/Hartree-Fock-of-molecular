cmake_minimum_required(VERSION 3.11)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

project(hf)

include(GNUInstallDirs)

include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/data)

set(CMAKE_CXX_STANDARD 17)


add_executable(hf 
                src/cal.cpp src/basis/basis.cpp src/geometry/geometry.cpp 
                src/geometry/getinfo.cpp src/integral/integral.cpp
                src/integral/kind1.cpp src/integral/kind2.cpp src/scf/scf.cpp 
                src/scf/mixing.cpp src/scf/occ.cpp src/scf/scfset.cpp
                src/hf/rhf.cpp src/hf/uhf.cpp src/hf/rohf.cpp src/main.cpp)


add_subdirectory(external/fmt)
include_directories("external/fmt/include")
target_link_libraries(hf fmt)
add_definitions(-DHAVE_FMT)


add_subdirectory(external/args)
include_directories("external/args")
target_link_libraries(hf args)
add_definitions(-DHAVE_ARGS)


add_library(json SHARED external/nlohmann/json.cpp)
include_directories("external/nlohmann")
target_include_directories(hf PUBLIC external/nlohmann/)

add_subdirectory(external/Catch2)
include_directories("external/Catch2/src")
target_link_libraries(hf Catch2WithMain)

add_subdirectory(external/libcint)
include_directories("external/libcint/src")
include_directories(${CMAKE_CURRENT_BINARY_DIR}/external/libcint/src)
target_link_libraries(hf cint)
add_definitions(-DHAVE_LIBCINT)

find_package(MPI REQUIRED)

target_link_libraries(hf MPI::MPI_CXX)


find_package(
        Armadillo
        8.500.0
        REQUIRED
)
if (${ARMADILLO_FOUND})
    message("Found Armadillo Library,")
    message("     ARMADILLO_VERSION: " ${ARMADILLO_VERSION_STRING})
    message("     ARMADILLO_LIBRARIES   : " ${ARMADILLO_LIBRARIES})
    message("     ARMADILLO_INCLUDE_DIRS: " ${ARMADILLO_INCLUDE_DIRS})

    include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS})
    target_link_libraries(hf ${ARMADILLO_LIBRARIES})
endif ()

find_package(
        HDF5
)

if (${HDF5_FOUND})
    message("Found HDF5 Library,")
    message("     HDF5_VERSION: " ${HDF5_VERSION})
    message("     HDF5_LIBRARIES   : " ${HDF5_LIBRARIES})
    message("     HDF5_INCLUDE_DIRS: " ${HDF5_INCLUDE_DIRS})

    include_directories(SYSTEM ${HDF5_INCLUDE_DIRS})
    target_link_libraries(hf ${HDF5_LIBRARIES})
endif ()

 find_package(
        GSL
        2.4
        REQUIRED
)

include_directories(SYSTEM ${GSL_INCLUDE_DIRS})
target_link_libraries(hf ${GSL_LIBRARIES})



find_package(
        Boost
        1.45.0
        REQUIRED
)

include_directories(SYSTEM ${Boost_INCLUDE_DIRS})
target_link_libraries(hf ${Boost_LIBRARIES})
