if(NOT DEFINED COMPARE_EXE)
    message(FATAL_ERROR "COMPARE_EXE is required")
endif()

if(NOT DEFINED OUT_CSV)
    message(FATAL_ERROR "OUT_CSV is required")
endif()

execute_process(
    COMMAND "${COMPARE_EXE}"
        --years 0.01
        --samples 1
        --bridge-method eos-lf4-2
        --bridge-cross-kick force-gradient
        --force-gradient-scale 0.5
        --bridge-corrector commutator-ab
        --bridge-corrector-scale 0.25
        --inner-mu-scale 1.0000000001
        --inner-whfast-corrector 11
        --csv "${OUT_CSV}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)

if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "solar_system_ias15_compare failed: ${compare_stdout}\n${compare_stderr}")
endif()

file(READ "${OUT_CSV}" csv_text)
string(REGEX MATCH "^[^\r\n]*" csv_header "${csv_text}")

set(required_columns
    bridge_cross_kick
    force_gradient_scale
    bridge_corrector
    bridge_corrector_scale
    inner_mu_scale
    inner_whfast_corrector
    moon_abs_error_au
    earth_abs_error_au
    moon_relative_error_au
    bridge_em_a_au
    ias15_em_a_au
    em_a_diff_au
    bridge_em_e
    ias15_em_e
    em_e_diff
    bridge_p_rel
    ias15_p_rel
    p_rel_diff
)

foreach(column IN LISTS required_columns)
    if(NOT ",${csv_header}," MATCHES ",${column},")
        message(FATAL_ERROR "CSV header is missing required column '${column}'. Header: ${csv_header}")
    endif()
endforeach()

if(NOT compare_stdout MATCHES "inner_mu_scale=1\\.0000000001")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing inner_mu_scale: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "inner_whfast_corrector=11")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing inner_whfast_corrector: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "bridge_cross_kick=force-gradient")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing bridge_cross_kick: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "force_gradient_scale=0\\.5")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing force_gradient_scale: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "bridge_corrector=commutator-ab")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing bridge_corrector: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "bridge_corrector_scale=0\\.25")
    message(FATAL_ERROR "solar_system_ias15_compare output is missing bridge_corrector_scale: ${compare_stdout}")
endif()
