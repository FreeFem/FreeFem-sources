CMAKE_MINIMUM_REQUIRED(VERSION 2.7)

PROJECT(suitesparse C)

ADD_DEFINITIONS(-O3)

FILE(GLOB UMFPACK_FILES UMFPACK/Source/*.c)
ADD_LIBRARY(umfpack SHARED ${UMFPACK_FILES})
TARGET_INCLUDE_DIRECTORIES(umfpack PRIVATE UMFPACK/Source 
                                           UMFPACK/Include
                                           AMD/Include
                                           SuiteSparse_config
                                           CHOLMOD/Include)
