# CMake script for DiffRate
# Author: Julien Dutheil
# Created: 27/02/2020


add_executable(bpprateshift RateShift.cpp)
set (bpprateshift-targets bpprateshift)

foreach (target ${bpprateshift-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    target_link_libraries (${target} ${LIBS})
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
    target_link_libraries (${target} ${LIBS})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${bpprateshift-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})
