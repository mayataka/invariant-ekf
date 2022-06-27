set(INEKF_PYTHON_VERSION ${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR})
set(INEKF_PYTHON_BINDINGS_LIBDIR ${CMAKE_INSTALL_LIBDIR}/python${INEKF_PYTHON_VERSION}/site-packages/inekf)

macro(pybind11_add_inekf_module MODULE)
  pybind11_add_module(
    ${MODULE} 
    SHARED 
    ${MODULE}.cpp
  )
  target_link_libraries(
    ${MODULE} 
    PRIVATE
    inekf::inekf
    inekf::inekf
  )
endmacro()

pybind11_add_inekf_module(robot_model)
pybind11_add_inekf_module(schmitt_trigger)
pybind11_add_inekf_module(contact_estimator)
pybind11_add_inekf_module(slip_estimator)
pybind11_add_inekf_module(state_estimator_settings)
pybind11_add_inekf_module(state_estimator)
pybind11_add_inekf_module(noise_params)

macro(install_inekf_pybind_module CURRENT_MODULE_DIR)
  file(GLOB PYTHON_BINDINGS_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_BINARY_DIR}/*.cpython*)
  file(GLOB PYTHON_FILES_${CURRENT_MODULE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
  install(
    FILES ${PYTHON_FILES_${CURRENT_MODULE_DIR}} ${PYTHON_BINDINGS_${CURRENT_MODULE_DIR}} 
    DESTINATION ${INEKF_PYTHON_BINDINGS_LIBDIR}
  )
endmacro()

install_inekf_pybind_module(inekf)