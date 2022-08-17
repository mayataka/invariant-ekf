set(LEGGED_STATE_ESTIMATOR_PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
set(LEGGED_STATE_ESTIMATOR_PYTHON_BINDINGS_LIBDIR ${CMAKE_INSTALL_LIBDIR}/python${LEGGED_STATE_ESTIMATOR_PYTHON_VERSION}/site-packages/legged_state_estimator)

macro(pybind11_add_legged_state_estimator_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_link_libraries(
    ${MODULE} 
    PRIVATE
    legged_state_estimator::legged_state_estimator
  )
endmacro()

pybind11_add_legged_state_estimator_module(robot_model)
pybind11_add_legged_state_estimator_module(contact_estimator)
pybind11_add_legged_state_estimator_module(slip_estimator)
pybind11_add_legged_state_estimator_module(state_estimator_settings)
pybind11_add_legged_state_estimator_module(state_estimator)
pybind11_add_legged_state_estimator_module(noise_params)

macro(install_legged_state_estimator_pybind_module CURRENT_MODULE_DIR)
  file(GLOB PYTHON_BINDINGS_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/*.cpython*)
  file(GLOB PYTHON_FILES_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
  install(
    FILES ${PYTHON_FILES_${CURRENT_MODULE_DIR}} ${PYTHON_BINDINGS_${CURRENT_MODULE_DIR}} 
    DESTINATION ${LEGGED_STATE_ESTIMATOR_PYTHON_BINDINGS_LIBDIR}
  )
endmacro()

install_legged_state_estimator_pybind_module(legged_state_estimator)