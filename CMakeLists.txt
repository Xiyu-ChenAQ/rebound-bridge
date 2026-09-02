cmake_minimum_required(VERSION 3.20)

project(
    rebound_bridge
    VERSION 0.2.0
    DESCRIPTION "Standalone symplectic bridge scheduler for REBOUND"
    LANGUAGES C
)

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(REBOUND_BRIDGE_BUILD_EXAMPLES   "Build rebound_bridge example programs"              OFF)
option(REBOUND_BRIDGE_BUILD_VALIDATION "Build validation comparison tools"                  OFF)
option(REBOUND_BRIDGE_BUILD_TESTS      "Build rebound_bridge tests"                         OFF)
option(REBOUND_BRIDGE_INSTALL          "Install rebound_bridge targets and headers"         ON)
set(REBOUND_BRIDGE_VENDORED_REBOUND_DIR "${CMAKE_CURRENT_SOURCE_DIR}/rebound/src")

if(MSVC)
    add_compile_options(/utf-8)
endif()

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
if(MSVC)
    target_compile_options(rebound_bridge_rebound PRIVATE /wd4477)
endif()
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

add_library(
    rebound_bridge STATIC
    # Public bridge API and lifecycle.
    src/bridge.c

    # Time stepping and operator compositions.
    src/bridge_advance.c
    src/bridge_eos.c

    # Cross-interaction operators.
    src/bridge_kick.c

)
add_library(rebound_bridge::rebound_bridge ALIAS rebound_bridge)
target_include_directories(
    rebound_bridge
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/include>
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
    # bridge_examples.c provides convenience constructors (solar-system fixtures,
    # Earth-Moon factory, etc.) that are only needed by the bundled examples.
    # It is compiled as a private object library so it never ends up in the
    # installed/exported rebound_bridge target.
    add_library(rebound_bridge_examples_obj OBJECT examples/bridge_examples.c)
    target_link_libraries(rebound_bridge_examples_obj PRIVATE rebound_bridge::rebound_bridge)
    target_include_directories(
        rebound_bridge_examples_obj
        PRIVATE
            ${CMAKE_CURRENT_SOURCE_DIR}/examples
    )

    set(_rebound_bridge_examples
        earth_moon_bridge
        jovian_laplace_bridge
        solar_system_bridge
    )

    foreach(_example IN LISTS _rebound_bridge_examples)
        add_executable(${_example} "examples/${_example}.c")
        target_link_libraries(${_example}
            PRIVATE
                rebound_bridge::rebound_bridge
                rebound_bridge_examples_obj
        )
        target_include_directories(${_example} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/examples)
    endforeach()
endif()

if(REBOUND_BRIDGE_BUILD_VALIDATION)
    # Shared utilities compiled once, linked to every validation executable.
    add_library(rebound_bridge_val_common STATIC
        validation/common/val_cli.c
        validation/common/solar_system_fixture.c
    )
    target_link_libraries(rebound_bridge_val_common PUBLIC rebound_bridge::rebound_bridge)
    target_include_directories(rebound_bridge_val_common
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/validation/common
    )

    add_executable(
        solar_system_whfast_bridge_bench
        validation/solar_system_ias15/bench_whfast_bridge.c
    )
    target_link_libraries(solar_system_whfast_bridge_bench PRIVATE rebound_bridge_val_common)

    add_executable(
        pluto_charon_whfast_bridge_bench
        validation/pluto_charon_binary/bench_pluto_charon.c
    )
    target_link_libraries(pluto_charon_whfast_bridge_bench PRIVATE rebound_bridge_val_common)

    add_executable(
        earth_moon_whfast_bridge_compare
        validation/earth_moon_whfast_bridge/compare_earth_moon_whfast_bridge.c
    )
    target_link_libraries(earth_moon_whfast_bridge_compare PRIVATE rebound_bridge_val_common)

endif()

if(REBOUND_BRIDGE_BUILD_TESTS)
    enable_testing()

    add_executable(
        rebound_bridge_step_dispatch_test
        tests/unit/bridge_step_dispatch_test.c
    )
    target_link_libraries(rebound_bridge_step_dispatch_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(NAME rebound_bridge_step_dispatch_test COMMAND rebound_bridge_step_dispatch_test)

    add_executable(
        rebound_bridge_registration_contract_test
        tests/unit/bridge_registration_contract_test.c
    )
    target_link_libraries(rebound_bridge_registration_contract_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(NAME rebound_bridge_registration_contract_test COMMAND rebound_bridge_registration_contract_test)

    add_executable(
        rebound_bridge_reversibility_test
        tests/unit/bridge_reversibility_test.c
    )
    target_link_libraries(rebound_bridge_reversibility_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(NAME rebound_bridge_reversibility_test COMMAND rebound_bridge_reversibility_test)

    add_executable(
        rebound_bridge_multibody_reversibility_test
        tests/unit/bridge_multibody_reversibility_test.c
    )
    target_link_libraries(rebound_bridge_multibody_reversibility_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(
        NAME rebound_bridge_multibody_reversibility_test
        COMMAND rebound_bridge_multibody_reversibility_test
    )

    add_executable(
        rebound_bridge_subsystem_coupling_test
        tests/unit/bridge_subsystem_coupling_test.c
    )
    target_link_libraries(rebound_bridge_subsystem_coupling_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(NAME rebound_bridge_subsystem_coupling_test COMMAND rebound_bridge_subsystem_coupling_test)

    add_executable(
        rebound_bridge_limiting_cases_test
        tests/unit/bridge_limiting_cases_test.c
    )
    target_link_libraries(rebound_bridge_limiting_cases_test PRIVATE rebound_bridge::rebound_bridge)
    add_test(NAME rebound_bridge_limiting_cases_test COMMAND rebound_bridge_limiting_cases_test)

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
        DIRECTORY src/include/
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
