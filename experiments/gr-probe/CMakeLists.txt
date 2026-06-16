cmake_minimum_required(VERSION 3.20)
project(gr_probe LANGUAGES C)

# Pure REBOUND 5 — no REBOUNDx. Reuse the source already downloaded by the main build.
set(REBOUND_SRC "${CMAKE_CURRENT_SOURCE_DIR}/../../build-validation/_deps/rebound_upstream-src/src" CACHE PATH "REBOUND src")
get_filename_component(REBOUND_SRC "${REBOUND_SRC}" ABSOLUTE)
if(NOT EXISTS "${REBOUND_SRC}/rebound.h")
    message(FATAL_ERROR "REBOUND src not found at ${REBOUND_SRC}")
endif()

file(GLOB REBOUND_SOURCES CONFIGURE_DEPENDS "${REBOUND_SRC}/*.c")

add_executable(gr_probe probe_gr.c ${REBOUND_SOURCES})
target_include_directories(gr_probe PRIVATE "${REBOUND_SRC}")
target_compile_definitions(gr_probe PRIVATE BUILDINGLIBREBOUND)
target_compile_features(gr_probe PRIVATE c_std_99)
if(UNIX)
    target_link_libraries(gr_probe PRIVATE m)
endif()
