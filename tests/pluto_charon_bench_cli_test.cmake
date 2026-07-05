if(NOT DEFINED BENCH_EXE)
    message(FATAL_ERROR "BENCH_EXE is required")
endif()

execute_process(
    COMMAND "${BENCH_EXE}"
        --case minimal
        --years 0.01
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --repeats 1
        --warmups 0
    RESULT_VARIABLE bench_result
    OUTPUT_VARIABLE bench_stdout
    ERROR_VARIABLE bench_stderr
)

if(NOT bench_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon benchmark failed: ${bench_stdout}\n${bench_stderr}")
endif()

foreach(metric IN ITEMS case_name bridge_method whfast_dt_days bridge_mean_seconds whfast_mean_seconds mean_speedup_whfast_over_bridge)
    if(NOT bench_stdout MATCHES "${metric}=([^\\r\\n]+)")
        message(FATAL_ERROR "Missing ${metric} in benchmark output: ${bench_stdout}")
    endif()
endforeach()

if(NOT bench_stdout MATCHES "case_name=minimal")
    message(FATAL_ERROR "Expected minimal benchmark output: ${bench_stdout}")
endif()
if(NOT bench_stdout MATCHES "bridge_method=eos-lf4-2")
    message(FATAL_ERROR "Expected eos-lf4-2 benchmark output: ${bench_stdout}")
endif()
if(NOT bench_stdout MATCHES "whfast_dt_days=0\\.062500000")
    message(FATAL_ERROR "Expected fair WHFast timestep 4/64 days: ${bench_stdout}")
endif()
