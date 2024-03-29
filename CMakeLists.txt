cmake_minimum_required(VERSION 3.6)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_C_COMPILER gcc)

project(SpecHap)

set(CMAKE_CXX_STANDARD 11)
add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
set(CMAKE_BUILD_TYPE DEBUG)
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXX_FLAGS} -O0 -Wall -ggdb -fkeep-inline-functions")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 ")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
set(EXTRACTOR_SRC "${CMAKE_SOURCE_DIR}/barcode_extractor")
set(HAIRS_SRC "${CMAKE_SOURCE_DIR}/hairs-src")

find_package(htslib REQUIRED)
find_package(arpack REQUIRED)
include_directories(${HTSlib_INCLUDE_DIRS} ${ARPACK_INCLUDE_DIRS})

aux_source_directory(. SPEC_SRC)
aux_source_directory(./barcode_extractor BARCODE_SRC)

add_executable(SpecHap ${SPEC_SRC})
add_executable(BarcodeExtract ${BARCODE_SRC})
target_link_libraries(SpecHap -lhts -larpack)
target_link_libraries(BarcodeExtract -lhts)

#add_custom_target(extractHAIRS
#        ALL COMMAND make
#        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/hairs-src)

install (TARGETS SpecHap DESTINATION bin)
install (TARGETS BarcodeExtract DESTINATION bin)
