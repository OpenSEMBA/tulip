include("${VCPKG_ROOT_DIR}/triplets/x64-linux.cmake")

if(PORT STREQUAL "gmsh")
    list(APPEND VCPKG_CMAKE_CONFIGURE_OPTIONS "-DENABLE_MESH=ON")
endif()
