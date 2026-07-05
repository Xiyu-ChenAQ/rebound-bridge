if(NOT DEFINED PLUTO_CHARON_EXE)
    message(FATAL_ERROR "PLUTO_CHARON_EXE is required")
endif()

if(NOT DEFINED OUT_CSV)
    message(FATAL_ERROR "OUT_CSV is required")
endif()

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --bridge-cross-kick force-gradient
        --force-gradient-scale 0.5
        --bridge-corrector commutator-ab
        --bridge-corrector-scale 0.25
        --bridge-advance-mode kdk
        --csv "${OUT_CSV}"
    RESULT_VARIABLE compare_result
    OUTPUT_VARIABLE compare_stdout
    ERROR_VARIABLE compare_stderr
)

if(NOT compare_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon comparison failed: ${compare_stdout}\n${compare_stderr}")
endif()

foreach(metric IN ITEMS bridge_method bridge_cross_kick force_gradient_scale bridge_corrector bridge_corrector_scale bridge_advance_mode ias15_sync_mode inner_mu_scale bridge_total_seconds ias15_total_seconds overall_speedup max_pc_com_abs_error_au max_pc_relative_error_au max_abs_pc_phase_arc_error_au)
    if(NOT compare_stdout MATCHES "${metric}=([^\\r\\n]+)")
        message(FATAL_ERROR "Missing ${metric} in output: ${compare_stdout}")
    endif()
    set("${metric}" "${CMAKE_MATCH_1}")
endforeach()

if(max_abs_pc_phase_arc_error_au GREATER 1e-6)
    message(FATAL_ERROR "Phase arc diagnostic is inconsistent with the position residual: max_abs_pc_phase_arc_error_au=${max_abs_pc_phase_arc_error_au}")
endif()

if(NOT force_gradient_scale MATCHES "0\\.5")
    message(FATAL_ERROR "force-gradient scale run did not report requested scale: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "bridge_corrector=commutator-ab")
    message(FATAL_ERROR "Bridge-level corrector run did not report requested corrector: ${compare_stdout}")
endif()

if(NOT compare_stdout MATCHES "bridge_corrector_scale=0\\.25")
    message(FATAL_ERROR "Bridge-level corrector scale run did not report requested scale: ${compare_stdout}")
endif()

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method forward-fg4
        --force-gradient-check
        --force-gradient-check-epsilon 1e-7
        --csv "${OUT_CSV}.force_gradient_check.csv"
    RESULT_VARIABLE force_gradient_check_result
    OUTPUT_VARIABLE force_gradient_check_stdout
    ERROR_VARIABLE force_gradient_check_stderr
)

if(NOT force_gradient_check_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon force-gradient check failed: ${force_gradient_check_stdout}\n${force_gradient_check_stderr}")
endif()

foreach(metric IN ITEMS force_gradient_max_abs_difference force_gradient_max_relative_difference force_gradient_internal_weighted_acc_residual force_gradient_main_weighted_acc_residual)
    if(NOT force_gradient_check_stdout MATCHES "${metric}=([^\\r\\n]+)")
        message(FATAL_ERROR "Missing ${metric} in force-gradient check output: ${force_gradient_check_stdout}")
    endif()
endforeach()

file(READ "${OUT_CSV}" csv_text)
string(REGEX MATCH "^[^\r\n]*" csv_header "${csv_text}")

set(required_columns
    case_name
    bridge_method
    bridge_cross_kick
    force_gradient_scale
    bridge_corrector
    bridge_corrector_scale
    bridge_advance_mode
    ias15_sync_mode
    inner_mu_scale
    dt_outer_days
    binary_ratio
    time_yr
    bridge_step_seconds
    ias15_step_seconds
    step_speedup
    bridge_cumulative_seconds
    ias15_cumulative_seconds
    cumulative_speedup
    bridge_energy_rel
    ias15_energy_rel
    energy_rel_diff
    bridge_l_rel
    ias15_l_rel
    l_rel_diff
    bridge_p_rel
    ias15_p_rel
    p_rel_diff
    pc_com_abs_error_au
    pc_com_velocity_error_au_per_yr
    pluto_abs_error_au
    charon_abs_error_au
    pc_relative_error_au
    bridge_pc_distance_au
    ias15_pc_distance_au
    pc_distance_diff_au
    bridge_pc_a_au
    ias15_pc_a_au
    pc_a_diff_au
    bridge_pc_e
    ias15_pc_e
    pc_e_diff
    bridge_pc_phase_deg
    ias15_pc_phase_deg
    pc_phase_diff_deg
    bridge_pc_true_anomaly_deg
    ias15_pc_true_anomaly_deg
    pc_true_anomaly_diff_deg
    bridge_pc_mean_anomaly_deg
    ias15_pc_mean_anomaly_deg
    pc_mean_anomaly_diff_deg
    pc_phase_arc_error_au
    bridge_com_sun_a_au
    ias15_com_sun_a_au
    com_sun_a_diff_au
    bridge_com_sun_e
    ias15_com_sun_e
    com_sun_e_diff
)

foreach(column IN LISTS required_columns)
    if(NOT ",${csv_header}," MATCHES ",${column},")
        message(FATAL_ERROR "CSV header is missing required column '${column}'. Header: ${csv_header}")
    endif()
endforeach()

string(REGEX MATCHALL "\n" newline_matches "${csv_text}")
list(LENGTH newline_matches newline_count)
if(NOT newline_count EQUAL 5)
    message(FATAL_ERROR "Expected header plus 4 sample rows, found ${newline_count} newline(s)")
endif()

set(OUTER_CSV "${OUT_CSV}.outer.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --csv "${OUTER_CSV}"
    RESULT_VARIABLE outer_result
    OUTPUT_VARIABLE outer_stdout
    ERROR_VARIABLE outer_stderr
)

if(NOT outer_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon outer comparison failed: ${outer_stdout}\n${outer_stderr}")
endif()

if(NOT outer_stdout MATCHES "case_name=outer")
    message(FATAL_ERROR "Outer run did not report case_name=outer: ${outer_stdout}")
endif()

file(READ "${OUTER_CSV}" outer_csv_text)
if(NOT outer_csv_text MATCHES "\nouter,eos-lf4-2,plain,[^,]*,none,[^,]*,")
    message(FATAL_ERROR "Outer CSV is missing outer/eos-lf4-2 sample rows: ${outer_csv_text}")
endif()

set(DKD_CSV "${OUT_CSV}.dkd.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method order2
        --bridge-advance-mode dkd
        --csv "${DKD_CSV}"
    RESULT_VARIABLE dkd_result
    OUTPUT_VARIABLE dkd_stdout
    ERROR_VARIABLE dkd_stderr
)

if(NOT dkd_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon DKD comparison failed: ${dkd_stdout}\n${dkd_stderr}")
endif()

if(NOT dkd_stdout MATCHES "bridge_advance_mode=dkd")
    message(FATAL_ERROR "DKD run did not report bridge_advance_mode=dkd: ${dkd_stdout}")
endif()

file(READ "${DKD_CSV}" dkd_csv_text)
if(NOT dkd_csv_text MATCHES "\nminimal,order2,plain,[^,]*,none,[^,]*,dkd,")
    message(FATAL_ERROR "DKD CSV is missing minimal/order2/dkd sample rows: ${dkd_csv_text}")
endif()

set(EOS_MULTIKICK_CSV "${OUT_CSV}.eos_multikick.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf4-2
        --bridge-advance-mode eos-multikick
        --csv "${EOS_MULTIKICK_CSV}"
    RESULT_VARIABLE eos_multikick_result
    OUTPUT_VARIABLE eos_multikick_stdout
    ERROR_VARIABLE eos_multikick_stderr
)

if(NOT eos_multikick_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon EOS multi-kick comparison failed: ${eos_multikick_stdout}\n${eos_multikick_stderr}")
endif()

if(NOT eos_multikick_stdout MATCHES "bridge_advance_mode=eos-multikick")
    message(FATAL_ERROR "EOS multi-kick run did not report bridge_advance_mode=eos-multikick: ${eos_multikick_stdout}")
endif()

file(READ "${EOS_MULTIKICK_CSV}" eos_multikick_csv_text)
if(NOT eos_multikick_csv_text MATCHES "\nminimal,eos-lf4-2,plain,[^,]*,none,[^,]*,eos-multikick,")
    message(FATAL_ERROR "EOS multi-kick CSV is missing minimal/eos-lf4-2/eos-multikick sample rows: ${eos_multikick_csv_text}")
endif()

set(EOS_LF8_6_4_CSV "${OUT_CSV}.eos_lf8_6_4.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method eos-lf8-6-4
        --csv "${EOS_LF8_6_4_CSV}"
    RESULT_VARIABLE eos_lf8_6_4_result
    OUTPUT_VARIABLE eos_lf8_6_4_stdout
    ERROR_VARIABLE eos_lf8_6_4_stderr
)

if(NOT eos_lf8_6_4_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon EOS LF8_6_4 comparison failed: ${eos_lf8_6_4_stdout}\n${eos_lf8_6_4_stderr}")
endif()

if(NOT eos_lf8_6_4_stdout MATCHES "bridge_method=eos-lf8-6-4")
    message(FATAL_ERROR "EOS LF8_6_4 run did not report requested method: ${eos_lf8_6_4_stdout}")
endif()

file(READ "${EOS_LF8_6_4_CSV}" eos_lf8_6_4_csv_text)
if(NOT eos_lf8_6_4_csv_text MATCHES "\nminimal,eos-lf8-6-4,plain,[^,]*,none,[^,]*,kdk,")
    message(FATAL_ERROR "EOS LF8_6_4 CSV is missing minimal/eos-lf8-6-4 sample rows: ${eos_lf8_6_4_csv_text}")
endif()

set(FORWARD_FG4_CSV "${OUT_CSV}.forward_fg4.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 4
        --binary-ratio 64
        --bridge-method forward-fg4
        --csv "${FORWARD_FG4_CSV}"
    RESULT_VARIABLE forward_fg4_result
    OUTPUT_VARIABLE forward_fg4_stdout
    ERROR_VARIABLE forward_fg4_stderr
)

if(NOT forward_fg4_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon forward-fg4 comparison failed: ${forward_fg4_stdout}\n${forward_fg4_stderr}")
endif()

if(NOT forward_fg4_stdout MATCHES "bridge_method=forward-fg4")
    message(FATAL_ERROR "forward-fg4 run did not report requested method: ${forward_fg4_stdout}")
endif()

file(READ "${FORWARD_FG4_CSV}" forward_fg4_csv_text)
if(NOT forward_fg4_csv_text MATCHES "\nminimal,forward-fg4,plain,[^,]*,none,[^,]*,kdk,")
    message(FATAL_ERROR "forward-fg4 CSV is missing minimal/forward-fg4 sample rows: ${forward_fg4_csv_text}")
endif()

set(INNER_MU_SCALE_CSV "${OUT_CSV}.inner_mu_scale.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method eos-lf4-2
        --bridge-advance-mode kdk
        --inner-mu-scale 1.000000001
        --csv "${INNER_MU_SCALE_CSV}"
    RESULT_VARIABLE inner_mu_scale_result
    OUTPUT_VARIABLE inner_mu_scale_stdout
    ERROR_VARIABLE inner_mu_scale_stderr
)

if(NOT inner_mu_scale_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon inner-mu-scale run failed: ${inner_mu_scale_stdout}\n${inner_mu_scale_stderr}")
endif()

if(NOT inner_mu_scale_stdout MATCHES "inner_mu_scale=1\\.000000001")
    message(FATAL_ERROR "inner-mu-scale run did not report requested scale: ${inner_mu_scale_stdout}")
endif()

file(READ "${INNER_MU_SCALE_CSV}" inner_mu_scale_csv_text)
if(NOT inner_mu_scale_csv_text MATCHES "\nminimal,eos-lf4-2,plain,[^,]*,none,[^,]*,kdk,outer,1\\.000000001[0-9]*,")
    message(FATAL_ERROR "inner-mu-scale CSV is missing requested scale rows: ${inner_mu_scale_csv_text}")
endif()

set(IAS15_OUTPUT_SYNC_CSV "${OUT_CSV}.ias15_output_sync.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case minimal
        --years 1
        --samples 4
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method eos-lf4-2
        --bridge-advance-mode kdk
        --ias15-sync-mode output
        --csv "${IAS15_OUTPUT_SYNC_CSV}"
    RESULT_VARIABLE ias15_output_sync_result
    OUTPUT_VARIABLE ias15_output_sync_stdout
    ERROR_VARIABLE ias15_output_sync_stderr
)

if(NOT ias15_output_sync_result EQUAL 0)
    message(FATAL_ERROR "Pluto-Charon IAS15 output-sync run failed: ${ias15_output_sync_stdout}\n${ias15_output_sync_stderr}")
endif()

if(NOT ias15_output_sync_stdout MATCHES "ias15_sync_mode=output")
    message(FATAL_ERROR "IAS15 output-sync run did not report requested mode: ${ias15_output_sync_stdout}")
endif()

file(READ "${IAS15_OUTPUT_SYNC_CSV}" ias15_output_sync_csv_text)
if(NOT ias15_output_sync_csv_text MATCHES "\nminimal,eos-lf4-2,plain,[^,]*,none,[^,]*,kdk,output,")
    message(FATAL_ERROR "IAS15 output-sync CSV is missing requested mode rows: ${ias15_output_sync_csv_text}")
endif()

set(SAMPLE_INVARIANCE_CSV_4 "${OUT_CSV}.samples4.csv")
set(SAMPLE_INVARIANCE_CSV_8 "${OUT_CSV}.samples8.csv")

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --years 1
        --samples 4
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method eos-lf4-2
        --bridge-advance-mode kdk
        --csv "${SAMPLE_INVARIANCE_CSV_4}"
    RESULT_VARIABLE sample4_result
    OUTPUT_VARIABLE sample4_stdout
    ERROR_VARIABLE sample4_stderr
)

if(NOT sample4_result EQUAL 0)
    message(FATAL_ERROR "Sample-invariance 4-sample run failed: ${sample4_stdout}\n${sample4_stderr}")
endif()

execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --years 1
        --samples 8
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method eos-lf4-2
        --bridge-advance-mode kdk
        --csv "${SAMPLE_INVARIANCE_CSV_8}"
    RESULT_VARIABLE sample8_result
    OUTPUT_VARIABLE sample8_stdout
    ERROR_VARIABLE sample8_stderr
)

if(NOT sample8_result EQUAL 0)
    message(FATAL_ERROR "Sample-invariance 8-sample run failed: ${sample8_stdout}\n${sample8_stderr}")
endif()

execute_process(
    COMMAND python -c
        "import csv,sys; p4,p8=sys.argv[1:3]; rows4=list(csv.DictReader(open(p4,newline=''))); rows8=list(csv.DictReader(open(p8,newline=''))); v4=float(rows4[-1]['pc_phase_arc_error_au']); v8=float(rows8[-1]['pc_phase_arc_error_au']); diff=abs(v4-v8); print(f'final_phase_samples4={v4:.17g}'); print(f'final_phase_samples8={v8:.17g}'); print(f'final_phase_sample_diff={diff:.17g}'); sys.exit(0 if diff <= 1e-13 else 1)"
        "${SAMPLE_INVARIANCE_CSV_4}"
        "${SAMPLE_INVARIANCE_CSV_8}"
    RESULT_VARIABLE sample_invariance_result
    OUTPUT_VARIABLE sample_invariance_stdout
    ERROR_VARIABLE sample_invariance_stderr
)

if(NOT sample_invariance_result EQUAL 0)
    message(FATAL_ERROR "Bridge comparison depends on output sample count:\n${sample_invariance_stdout}\n${sample_invariance_stderr}")
endif()

set(LONG_TIME_CSV "${OUT_CSV}.long_time_130yr.csv")
execute_process(
    COMMAND "${PLUTO_CHARON_EXE}"
        --case outer
        --years 130
        --samples 13
        --dt-outer-days 2
        --binary-ratio 32
        --bridge-method eos-lf4-2
        --bridge-advance-mode kdk
        --ias15-sync-mode output
        --csv "${LONG_TIME_CSV}"
    RESULT_VARIABLE long_time_result
    OUTPUT_VARIABLE long_time_stdout
    ERROR_VARIABLE long_time_stderr
    TIMEOUT 20
)

if(NOT long_time_result EQUAL 0)
    message(FATAL_ERROR "Long-time Pluto-Charon run failed or timed out near t=128yr: ${long_time_stdout}\n${long_time_stderr}")
endif()

if(NOT long_time_stdout MATCHES "years=130")
    message(FATAL_ERROR "Long-time run did not complete expected 130yr case: ${long_time_stdout}")
endif()
