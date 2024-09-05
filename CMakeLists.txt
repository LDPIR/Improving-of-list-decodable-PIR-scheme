cmake_minimum_required(VERSION 3.0.0)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(MyProject VERSION 0.1.0)

include(CTest)
enable_testing()

include_directories(/usr/local/include)
link_directories(/usr/local/lib)

set(DATA_DIR "${CMAKE_BINARY_DIR}/data/")

add_subdirectory(20xx_ours)
# add_subdirectory(2007_improving)


set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
