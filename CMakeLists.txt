cmake_minimum_required(VERSION 3.10)

project(inference CXX C)

# set(PROJECT_URL "TBD")
set(PROJECT_DESCRIPTION "workflow for generating a personal genome through variation-aware alignment and imputation")
set(CMAKE_CXX_STANDARD 11)

# add_subdirectory(varcount "${PROJECT_SOURCE_DIR}/bin")
# add_subdirectory(liftover "${PROJECT_SOURCE_DIR}/bin")
add_subdirectory(varcount)
add_subdirectory(liftover)
