cmake_minimum_required(VERSION 2.8.9)

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

SET(SRCS
    BaseModel.cc
    Debug.cc
    fusion.cc
    fusion_p.h
    IntMap.cc
    mosektask.cc
    mosektask_p.h
    SolverInfo.cc
    SolverInfo.h
    StringBuffer.cc
  
)
SET(INCLUDE_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../h/")
SET(HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/SolverInfo.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fusion_p.h
  ${INCLUDE_ROOT}/fusion.h
  ${INCLUDE_ROOT}/fusion_fwd.h
  ${INCLUDE_ROOT}/monty.h
  ${INCLUDE_ROOT}/monty_base.h
  ${INCLUDE_ROOT}/monty_iterator.h
  ${INCLUDE_ROOT}/monty_ndarray.h
  ${INCLUDE_ROOT}/monty_rc.h
  ${INCLUDE_ROOT}/monty_shape.h
  ${INCLUDE_ROOT}/mosek.h
  ${INCLUDE_ROOT}/mosektask.h
  ${INCLUDE_ROOT}/mosektask.h

)

# Create a single library for the project
add_library(mosek STATIC ${SRCS} ${HEADERS})
target_link_libraries(mosek PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libfusion64.so.9.2)
target_link_libraries(mosek PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/libmosek64.so.9.2)
add_definitions(-DNOMINMAX)


# Include settings

target_include_directories(mosek PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/../h")
