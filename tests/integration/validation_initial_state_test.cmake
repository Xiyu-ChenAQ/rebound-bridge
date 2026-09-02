if(NOT DEFINED COMPARE_EXE)
    message(FATAL_ERROR "COMPARE_EXE is required")
endif()

if(NOT DEFINED OUT_CSV)
    message(FATAL_ERROR "OUT_CSV is required")
endif()

execute_process(
    COMMAND "${COMPARE_EXE}"
        --initial-state-csv "${OUT_CSV}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)

if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "initial state check failed: ${compare_stdout}\n${compare_stderr}")
endif()

foreach(metric IN ITEMS max_abs_mass_diff max_dr_au max_dv_au_per_yr)
    if(NOT compare_stdout MATCHES "${metric}=([0-9.eE+-]+)")
        message(FATAL_ERROR "Missing ${metric} in output: ${compare_stdout}")
    endif()
    set("${metric}" "${CMAKE_MATCH_1}")
endforeach()

if(max_abs_mass_diff GREATER 1e-15)
    message(FATAL_ERROR "Initial masses differ: max_abs_mass_diff=${max_abs_mass_diff}")
endif()
if(max_dr_au GREATER 1e-13)
    message(FATAL_ERROR "Initial positions differ above machine-level tolerance: max_dr_au=${max_dr_au}")
endif()
if(max_dv_au_per_yr GREATER 1e-13)
    message(FATAL_ERROR "Initial velocities differ above machine-level tolerance: max_dv_au_per_yr=${max_dv_au_per_yr}")
endif()

file(READ "${OUT_CSV}" csv_text)
string(REGEX MATCH "^[^\r\n]*" csv_header "${csv_text}")

set(required_columns
    body
    bridge_m
    ias15_m
    mass_diff
    bridge_x
    bridge_y
    bridge_z
    bridge_vx
    bridge_vy
    bridge_vz
    ias15_x
    ias15_y
    ias15_z
    ias15_vx
    ias15_vy
    ias15_vz
    dr_au
    dv_au_per_yr
)

foreach(column IN LISTS required_columns)
    if(NOT ",${csv_header}," MATCHES ",${column},")
        message(FATAL_ERROR "CSV header is missing required column '${column}'. Header: ${csv_header}")
    endif()
endforeach()

string(REGEX MATCHALL "\n" newline_matches "${csv_text}")
list(LENGTH newline_matches newline_count)
if(NOT newline_count EQUAL 14)
    message(FATAL_ERROR "Expected header plus 13 body rows, found ${newline_count} newline(s)")
endif()
