cmake_minimum_required(VERSION 3.20)

project(
    rebound_bridge
    VERSION 0.1.0
    DESCRIPTION "Standalone symplectic bridge scheduler for REBOUND"
    LANGUAGES C
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(REBOUND_BRIDGE_BUILD_EXAMPLES "Build rebound_bridge example programs" ON)
option(REBOUND_BRIDGE_BUILD_VALIDATION "Build validation comparison tools" ON)
option(REBOUND_BRIDGE_INSTALL "Install rebound_bridge targets and headers" ON)
set(REBOUND_BRIDGE_VENDORED_REBOUND_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/rebound/src")

if(NOT EXISTS "${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}/rebound.h")
    message(FATAL_ERROR "Vendored REBOUND source tree is missing rebound.h: ${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}")
endif()

file(GLOB REBOUND_BRIDGE_REBOUND_SOURCES CONFIGURE_DEPENDS "${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}/*.c")
if(NOT REBOUND_BRIDGE_REBOUND_SOURCES)
    message(FATAL_ERROR "No vendored REBOUND C sources found in ${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}")
endif()

add_library(rebound_bridge_rebound STATIC ${REBOUND_BRIDGE_REBOUND_SOURCES})
add_library(rebound_bridge::rebound ALIAS rebound_bridge_rebound)
target_include_directories(
    rebound_bridge_rebound
    PUBLIC
        $<BUILD_INTERFACE:${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(rebound_bridge_rebound PUBLIC c_std_99)
target_compile_definitions(rebound_bridge_rebound PRIVATE BUILDINGLIBREBOUND)
target_compile_definitions(rebound_bridge_rebound INTERFACE BUILDINGLIBREBOUND)
set_target_properties(
    rebound_bridge_rebound
    PROPERTIES
        EXPORT_NAME rebound
        OUTPUT_NAME rebound
        POSITION_INDEPENDENT_CODE ON
)

if(UNIX)
    target_link_libraries(rebound_bridge_rebound PUBLIC m)
endif()

add_library(rebound_bridge STATIC src/rebound_bridge.c)
add_library(rebound_bridge::rebound_bridge ALIAS rebound_bridge)
target_include_directories(
    rebound_bridge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_compile_features(rebound_bridge PUBLIC c_std_99)
target_link_libraries(rebound_bridge PUBLIC rebound_bridge::rebound)
set_target_properties(
    rebound_bridge
    PROPERTIES
        OUTPUT_NAME rebound_bridge
        POSITION_INDEPENDENT_CODE ON
)

if(REBOUND_BRIDGE_BUILD_EXAMPLES)
    set(_rebound_bridge_examples
        earth_moon_bridge
        jovian_laplace_bridge
        solar_system_bridge
    )

    foreach(_example IN LISTS _rebound_bridge_examples)
        add_executable(${_example} "examples/${_example}.c")
        target_link_libraries(${_example} PRIVATE rebound_bridge::rebound_bridge)
    endforeach()
endif()

if(REBOUND_BRIDGE_BUILD_VALIDATION)
    add_executable(
        solar_system_ias15_compare
        validation/solar_system_ias15/compare_solar_system_ias15.c
    )
    target_link_libraries(solar_system_ias15_compare PRIVATE rebound_bridge::rebound_bridge)

    add_executable(
        solar_system_whfast_bridge_bench
        validation/solar_system_ias15/bench_whfast_bridge.c
    )
    target_link_libraries(solar_system_whfast_bridge_bench PRIVATE rebound_bridge::rebound_bridge)
endif()

if(REBOUND_BRIDGE_INSTALL)
    install(
        TARGETS rebound_bridge_rebound rebound_bridge
        EXPORT rebound_bridgeTargets
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )

    install(
        DIRECTORY include/
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    )

    install(
        DIRECTORY "${REBOUND_BRIDGE_VENDORED_REBOUND_DIR}/"
        DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
        FILES_MATCHING PATTERN "*.h"
    )

    install(
        EXPORT rebound_bridgeTargets
        NAMESPACE rebound_bridge::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rebound_bridge
    )

    write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/rebound_bridgeConfigVersion.cmake"
        VERSION ${PROJECT_VERSION}
        COMPATIBILITY SameMajorVersion
    )

    configure_package_config_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/rebound_bridgeConfig.cmake.in"
        "${CMAKE_CURRENT_BINARY_DIR}/rebound_bridgeConfig.cmake"
        INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rebound_bridge
    )

    install(
        FILES
            "${CMAKE_CURRENT_BINARY_DIR}/rebound_bridgeConfig.cmake"
            "${CMAKE_CURRENT_BINARY_DIR}/rebound_bridgeConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/rebound_bridge
    )
endif()
