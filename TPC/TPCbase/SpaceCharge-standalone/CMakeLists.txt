cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)

# Define variables
set(MODULE O2TPCSpaceCharge)
set(sources
    ../AliTPC3DCylindricalInterpolator.cxx
    ../AliTPC3DCylindricalInterpolatorIrregular.cxx
    ../AliTPCLookUpTable3DInterpolatorD.cxx
    ../AliTPCLookUpTable3DInterpolatorIrregularD.cxx
    ../AliTPCPoissonSolver.cxx
    ../AliTPCSpaceCharge3DCalc.cxx
    ../../../STEER/STEERBase/AliSysInfo.cxx)
set(headers
    ../AliTPC3DCylindricalInterpolator.h
    ../AliTPC3DCylindricalInterpolatorIrregular.h
    ../AliTPCLookUpTable3DInterpolatorD.h
    ../AliTPCLookUpTable3DInterpolatorIrregularD.h
    ../AliTPCPoissonSolver.h
    ../AliTPCSpaceCharge3DCalc.h
    ../../../STEER/STEERBase/AliSysInfo.h)

# Configure ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
find_package(ROOT REQUIRED)
include(${ROOT_USE_FILE})

# Define include directories
include_directories(.
                    ..
                    ../../../STEER/STEERBase
                    ${ROOT_INCLUDE_DIRS})

# Generate dictionary
ROOT_GENERATE_DICTIONARY(G__${MODULE} ${headers} LINKDEF ${MODULE}LinkDef.h)

# Make shared library
add_library(${MODULE} SHARED
            ${sources}
            G__${MODULE}.cxx)
target_link_libraries(${MODULE} ${ROOT_LIBRARIES})

install(TARGETS ${MODULE}
        ARCHIVE DESTINATION lib
        LIBRARY DESTINATION lib)

install(FILES
        ${headers}
        ${MODULE}LinkDef.h
        DESTINATION include)