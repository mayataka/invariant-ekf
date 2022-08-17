cmake_minimum_required(VERSION 2.8.3)
project(legged_state_estimator)

if (NOT "${CMAKE_CXX_STANDARD}")
  set(CMAKE_CXX_STANDARD 11)
else()
  if ("${CMAKE_CXX_STANDARD}" LESS 11)
    set(CMAKE_CXX_STANDARD 11)
  endif()
endif()


################
## Build type ##
################
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()
if(NOT CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "WARNING: CMAKE_BUILD_TYPE is NOT set to Release, which can decrease performance significantly.")
endif()


#############
## Options ##
#############
option(OPTIMIZE_FOR_NATIVE "Enable -march=native" ON)
option(BUILD_PYTHON_INTERFACE "Build Python bindings" ON)
option(BUILD_EXAMPLES "Build examples and tests" ON)


###########
## Build ##
###########
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
find_package(Eigen3 REQUIRED)
find_package(pinocchio REQUIRED)
file(GLOB_RECURSE ${PROJECT_NAME}_SOURCES src/*.cpp)
file(GLOB_RECURSE ${PROJECT_NAME}_HEADERS include/${PROJECT_NAME}/*.h*)
add_library(
  ${PROJECT_NAME} 
  SHARED
  ${${PROJECT_NAME}_SOURCES} 
  ${${PROJECT_NAME}_HEADERS}
)
target_link_libraries(
  ${PROJECT_NAME} 
  PUBLIC
  ${PINOCCHIO_LIBRARIES}
)
target_include_directories(
  ${PROJECT_NAME} 
  PUBLIC
  ${EIGEN3_INCLUDE_DIR}
  ${PINOCCHIO_INCLUDE_DIRS}
  $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)
if (OPTIMIZE_FOR_NATIVE)
  target_compile_options(
    ${PROJECT_NAME} 
    PUBLIC
    -march=native
  )
endif()


##############
## Bindings ##
##############
if (BUILD_PYTHON_INTERFACE)
  add_subdirectory(bindings/python)
endif()


######################
#  Add Execuatables  #
######################
find_package(Threads)
find_package(Boost REQUIRED COMPONENTS system)

macro(legged_state_estimator_add_test EXACUTABLE)
  add_executable(
    ${EXACUTABLE} 
    ${PROJECT_SOURCE_DIR}/tests/${EXACUTABLE}.cpp
  )
  target_link_libraries(
    ${EXACUTABLE} 
    PRIVATE
    ${PROJECT_NAME} 
    ${Boost_LIBRARIES}
  )
  target_include_directories(
    ${EXACUTABLE} 
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
  )
endmacro()

if (BUILD_EXAMPLES)
  legged_state_estimator_add_test(propagation_speed)
  legged_state_estimator_add_test(correction_speed)
  legged_state_estimator_add_test(large_state_speed)
  legged_state_estimator_add_test(left_vs_right_error_dynamics)
endif()

macro(legged_state_estimator_add_example EXACUTABLE)
  add_executable(
    ${EXACUTABLE} 
    ${PROJECT_SOURCE_DIR}/examples/${EXACUTABLE}.cpp
  )
  target_link_libraries(
    ${EXACUTABLE} 
    PRIVATE
    ${PROJECT_NAME} 
    ${Boost_LIBRARIES}
  )
  target_include_directories(
    ${EXACUTABLE} 
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${Boost_INCLUDE_DIRS}
  )
endmacro()

if (BUILD_EXAMPLES)
  legged_state_estimator_add_example(landmarks)
  legged_state_estimator_add_example(kinematics)
endif()


#############
## Install ##
#############
include(GNUInstallDirs)
# Install lib files
install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}-config
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}
)
# Install header files
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/${PROJECT_NAME}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME}
)
# Install config file 
set(CONFIG_PACKAGE_LOCATION "share/${PROJECT_NAME}/cmake")
install(
  EXPORT ${PROJECT_NAME}-config
  FILE ${PROJECT_NAME}-config.cmake
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CONFIG_PACKAGE_LOCATION}
)
# Install alias
add_library(
  ${PROJECT_NAME}::${PROJECT_NAME} 
  ALIAS ${PROJECT_NAME}
)