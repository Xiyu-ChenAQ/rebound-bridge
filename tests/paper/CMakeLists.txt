cmake_minimum_required(VERSION 3.20)

project(rebound_bridge_paper_tests LANGUAGES C)

enable_testing()

set(REBOUND_BRIDGE_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(REBOUND_BRIDGE_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(REBOUND_BRIDGE_INSTALL OFF CACHE BOOL "" FORCE)

add_subdirectory("${CMAKE_CURRENT_LIST_DIR}/../.." rebound_bridge_root)

set(PAPER_GALILEAN_TARGETS
    galilean_horizons_compare
    galilean_whfast_compare
    galilean_long_energy
    galilean_reversibility
)

foreach(target IN LISTS PAPER_GALILEAN_TARGETS)
    add_executable(${target} "galilean/${target}.c")
    target_include_directories(
        ${target}
        PRIVATE
            "${CMAKE_CURRENT_LIST_DIR}/common"
    )
    target_link_libraries(${target} PRIVATE rebound_bridge::rebound_bridge)
endforeach()

add_executable(generalized_order_crossover generalized_order/generalized_order_crossover.c)
target_link_libraries(generalized_order_crossover PRIVATE rebound_bridge::rebound_bridge)

add_executable(bridge_limiting_cases_test checks/bridge_limiting_cases_test.c)
target_link_libraries(bridge_limiting_cases_test PRIVATE rebound_bridge::rebound_bridge)
add_test(NAME paper_limiting_cases COMMAND bridge_limiting_cases_test)

add_executable(galilean_exomoon_survival exomoon/galilean_exomoon_survival.c)
target_include_directories(
    galilean_exomoon_survival
    PRIVATE
        "${CMAKE_CURRENT_LIST_DIR}/common"
)
target_link_libraries(galilean_exomoon_survival PRIVATE rebound_bridge::rebound_bridge)

find_package(Python3 REQUIRED COMPONENTS Interpreter)
add_test(NAME paper_harness_schedule
         COMMAND "${Python3_EXECUTABLE}"
                 "${CMAKE_CURRENT_LIST_DIR}/checks/check_harness_schedule.py"
                 --horizons "$<TARGET_FILE:galilean_horizons_compare>"
                 --exomoon "$<TARGET_FILE:galilean_exomoon_survival>")

# Runtime scaling is part of the paper reproduction suite.
add_executable(complex_speed_ratio scaling/complex_speed_ratio.c)
target_include_directories(complex_speed_ratio PRIVATE "${CMAKE_CURRENT_LIST_DIR}/common")
target_link_libraries(complex_speed_ratio PRIVATE rebound_bridge::rebound_bridge)

# Two resolved subsystems: time-integration error versus spatial-model error.
add_executable(multiple_subsystems multiple_subsystems/multiple_subsystems.c)
target_link_libraries(multiple_subsystems PRIVATE rebound_bridge::rebound_bridge)
add_test(NAME paper_multiple_subsystems_kick
         COMMAND multiple_subsystems --kick-check "${CMAKE_CURRENT_BINARY_DIR}/multiple_subsystems_kick_check.csv")
