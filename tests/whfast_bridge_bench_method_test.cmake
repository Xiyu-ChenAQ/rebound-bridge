if(NOT DEFINED BENCH_EXE)
    message(FATAL_ERROR "BENCH_EXE is required")
endif()

execute_process(
    COMMAND "${BENCH_EXE}"
        --years 0.001
        --dt-outer-days 0.5
        --earth-moon-ratio 20
        --jovian-ratio 50
        --bridge-method eos-lf4-2
        --repeats 1
        --warmups 0
    RESULT_VARIABLE bench_result
    OUTPUT_VARIABLE bench_stdout
    ERROR_VARIABLE bench_stderr
)

if(NOT bench_result EQUAL 0)
    message(FATAL_ERROR "solar_system_whfast_bridge_bench failed: ${bench_stdout}\n${bench_stderr}")
endif()

if(NOT bench_stdout MATCHES "bridge_method=eos-lf4-2")
    message(FATAL_ERROR "Expected bridge_method=eos-lf4-2 in output: ${bench_stdout}")
endif()
