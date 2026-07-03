# REBOUND Source-tree Vendoring Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Vendor the complete REBOUND upstream `src/` tree into this repository and make the default build use repository-contained source only.

**Architecture:** Keep `rebound_bridge`'s existing REBOUND-facing API unchanged. Copy the full upstream REBOUND `src/` directory into `third_party/rebound/src`, then change CMake to build `rebound_bridge_rebound` from that vendored tree instead of `_deps`, `REBOUND_SRC_DIR`, or `FetchContent`.

**Tech Stack:** C99, CMake 3.20, MSVC/Visual Studio Developer PowerShell, vendored GPL REBOUND C source.

---

## File Structure

- Create: `third_party/rebound/src/`
- Modify: `CMakeLists.txt`
- Modify: `README.md`
- Optional modify: `validation/solar_system_ias15/run.ps1`
- Do not modify as part of vendoring unless required by build failures: `include/rebound_bridge.h`
- Do not modify as part of vendoring unless required by build failures: `src/rebound_bridge.c`
- Do not include unrelated working-tree changes: `examples/solar_system_bridge.c`

`third_party/rebound/src/` is a verbatim copy of the upstream REBOUND `src/` directory currently available at `_deps/rebound-upstream/src/`. Upstream source comments must remain unchanged.

## Task 1: Vendor the REBOUND Source Tree

**Files:**
- Create: `third_party/rebound/src/*`

- [ ] **Step 1: Verify the source tree exists before copying**

Run:

```powershell
Test-Path _deps\rebound-upstream\src\rebound.h
Get-ChildItem _deps\rebound-upstream\src -Filter *.c | Measure-Object
Get-ChildItem _deps\rebound-upstream\src -Filter *.h | Measure-Object
```

Expected:

```text
True
Count is 32 for .c files
Count is 34 for .h files
```

- [ ] **Step 2: Copy the complete upstream `src/` tree without editing files**

Run:

```powershell
New-Item -ItemType Directory -Force -Path third_party\rebound | Out-Null
Copy-Item -Path _deps\rebound-upstream\src -Destination third_party\rebound\src -Recurse -Force
```

Expected:

```text
third_party\rebound\src\rebound.h exists
third_party\rebound\src\integrator_whfast.c exists
third_party\rebound\src\integrator_ias15.c exists
```

- [ ] **Step 3: Verify copied upstream files are byte-identical**

Run:

```powershell
$srcFiles = Get-ChildItem _deps\rebound-upstream\src -File
foreach ($file in $srcFiles) {
    $vendored = Join-Path "third_party\rebound\src" $file.Name
    if (!(Test-Path $vendored)) { throw "Missing vendored file: $($file.Name)" }
    $a = Get-FileHash -Algorithm SHA256 $file.FullName
    $b = Get-FileHash -Algorithm SHA256 $vendored
    if ($a.Hash -ne $b.Hash) { throw "Hash mismatch: $($file.Name)" }
}
"verified $($srcFiles.Count) vendored files"
```

Expected:

```text
verified 68 vendored files
```

- [ ] **Step 4: Commit the vendored source tree**

Run:

```powershell
git add third_party\rebound\src
git commit -m "vendor: add rebound source tree"
```

Expected:

```text
Commit includes only third_party/rebound/src files.
```

## Task 2: Switch CMake to the Vendored Tree

**Files:**
- Modify: `CMakeLists.txt`

- [ ] **Step 1: Replace external REBOUND discovery with a fixed vendored source path**

Edit `CMakeLists.txt`:

```cmake
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

option(REBOUND_BRIDGE_BUILD_EXAMPLES "Build rebound_bridge example programs" ON)
option(REBOUND_BRIDGE_BUILD_VALIDATION "Build validation comparison tools" ON)
option(REBOUND_BRIDGE_INSTALL "Install rebound_bridge targets and headers" ON)

set(REBOUND_BRIDGE_REBOUND_SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/third_party/rebound/src")

if(NOT EXISTS "${REBOUND_BRIDGE_REBOUND_SRC_DIR}/rebound.h")
    message(FATAL_ERROR "Vendored REBOUND source tree is missing rebound.h: ${REBOUND_BRIDGE_REBOUND_SRC_DIR}")
endif()

file(GLOB REBOUND_BRIDGE_REBOUND_SOURCES CONFIGURE_DEPENDS "${REBOUND_BRIDGE_REBOUND_SRC_DIR}/*.c")
if(NOT REBOUND_BRIDGE_REBOUND_SOURCES)
    message(FATAL_ERROR "No vendored REBOUND C sources found in ${REBOUND_BRIDGE_REBOUND_SRC_DIR}")
endif()
```

Remove the old `include(FetchContent)`, `REBOUND_BRIDGE_FETCH_REBOUND`, `REBOUND_BRIDGE_REBOUND_TAG`, `_rebound_bridge_default_rebound_src`, `REBOUND_SRC_DIR`, `FetchContent_Declare`, external `rebound.h` validation, and references to `${REBOUND_SRC_DIR}`.

- [ ] **Step 2: Update include directories to use the vendored path**

In the `target_include_directories(rebound_bridge_rebound ...)` block, use:

```cmake
target_include_directories(
    rebound_bridge_rebound
    PUBLIC
        $<BUILD_INTERFACE:${REBOUND_BRIDGE_REBOUND_SRC_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
```

- [ ] **Step 3: Update install headers to copy from the vendored path**

Replace:

```cmake
install(
    DIRECTORY "${REBOUND_SRC_DIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.h"
)
```

with:

```cmake
install(
    DIRECTORY "${REBOUND_BRIDGE_REBOUND_SRC_DIR}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.h"
)
```

- [ ] **Step 4: Verify CMake no longer references external REBOUND configuration**

Run:

```powershell
rg -n "FetchContent|REBOUND_SRC_DIR|rebound_upstream|REBOUND_BRIDGE_FETCH_REBOUND|REBOUND_BRIDGE_REBOUND_TAG|_deps/rebound-upstream" CMakeLists.txt
```

Expected:

```text
No matches.
```

- [ ] **Step 5: Commit the CMake switch**

Run:

```powershell
git add CMakeLists.txt
git commit -m "build: use vendored rebound source"
```

Expected:

```text
Commit includes only CMakeLists.txt.
```

## Task 3: Validate Configure and Build Without External REBOUND

**Files:**
- No source changes expected.

- [ ] **Step 1: Configure a clean build without `REBOUND_SRC_DIR`**

Run:

```powershell
$env:REBOUND_SRC_DIR = $null
cmake -S . -B build-vendored -DREBOUND_BRIDGE_BUILD_EXAMPLES=ON -DREBOUND_BRIDGE_BUILD_VALIDATION=ON
```

Expected:

```text
Configure succeeds.
Output does not show FetchContent cloning REBOUND.
Output does not require REBOUND_SRC_DIR.
```

- [ ] **Step 2: Build examples and validation targets**

Run:

```powershell
cmake --build build-vendored --config Release
```

Expected:

```text
Build succeeds.
Targets include rebound_bridge, earth_moon_bridge, jovian_laplace_bridge, solar_system_bridge.
Validation targets build if the vendored upstream tree supports all required APIs.
```

- [ ] **Step 3: If validation targets fail only because IAS15 comparison assumes external paths, temporarily disable default validation**

If Step 2 fails in validation code but core library and examples compile, edit the default option:

```cmake
option(REBOUND_BRIDGE_BUILD_VALIDATION "Build validation comparison tools" OFF)
```

Then run:

```powershell
cmake -S . -B build-vendored -DREBOUND_BRIDGE_BUILD_EXAMPLES=ON
cmake --build build-vendored --config Release
```

Expected:

```text
Core library and examples build from vendored source.
```

Commit only if this fallback edit is needed:

```powershell
git add CMakeLists.txt
git commit -m "build: make validation opt-in during vendoring"
```

## Task 4: Run Runtime Smoke Tests

**Files:**
- No source changes expected unless smoke tests expose a vendoring integration bug.

- [ ] **Step 1: Run Earth-Moon example**

Run:

```powershell
$exe = "build-vendored\Release\earth_moon_bridge.exe"
if (!(Test-Path $exe)) { $exe = "build-vendored\earth_moon_bridge.exe" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "earth_moon_bridge failed with exit code $LASTEXITCODE" }
```

Expected:

```text
Program exits with code 0.
Output includes "rebound_bridge - convergence test".
```

- [ ] **Step 2: Run Jovian Laplace example**

Run:

```powershell
$exe = "build-vendored\Release\jovian_laplace_bridge.exe"
if (!(Test-Path $exe)) { $exe = "build-vendored\jovian_laplace_bridge.exe" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "jovian_laplace_bridge failed with exit code $LASTEXITCODE" }
```

Expected:

```text
Program exits with code 0.
Output includes "Jovian Laplace resonance example".
```

- [ ] **Step 3: Run solar-system example if it is not carrying unrelated local edits**

Before running, inspect whether `examples/solar_system_bridge.c` has unrelated working-tree changes:

```powershell
git status --short examples\solar_system_bridge.c
```

If it is modified, do not edit or revert it. Run the built executable as-is:

```powershell
$exe = "build-vendored\Release\solar_system_bridge.exe"
if (!(Test-Path $exe)) { $exe = "build-vendored\solar_system_bridge.exe" }
& $exe
if ($LASTEXITCODE -ne 0) { throw "solar_system_bridge failed with exit code $LASTEXITCODE" }
```

Expected:

```text
Program exits with code 0 when the current local example state is valid.
If it fails due to pre-existing local edits, record the failure and do not revert the file.
```

## Task 5: Update Documentation

**Files:**
- Modify: `README.md`
- Optional modify: `validation/solar_system_ias15/run.ps1`

- [ ] **Step 1: Replace external REBOUND setup instructions in README**

In `README.md`, replace the "Using REBOUND" section with:

````markdown
## Using REBOUND

`rebound_bridge` vendors the REBOUND C source tree under
`third_party/rebound/src` and builds it as part of this repository. A normal
build does not download REBOUND and does not require `REBOUND_SRC_DIR`.

```powershell
cmake -S . -B build
cmake --build build --config Release
```

The build creates two static targets:

- `rebound_bridge::rebound` compiles the vendored REBOUND C sources.
- `rebound_bridge::rebound_bridge` compiles this bridge library and links it
  against `rebound_bridge::rebound`.

Downstream CMake projects link against:

```cmake
target_link_libraries(my_target PRIVATE rebound_bridge::rebound_bridge)
```
````

Do not document `REBOUND_SRC_DIR` as part of the default workflow.

- [ ] **Step 2: Update validation script only if it still injects `REBOUND_SRC_DIR`**

If `validation/solar_system_ias15/run.ps1` still contains:

```powershell
$cmakeArgs += "-DREBOUND_SRC_DIR=$ReboundSrcDir"
```

remove the external REBOUND argument path and configure with repository-contained source:

```powershell
$cmakeArgs = @(
    "-S", $repoRoot,
    "-B", $BuildDir,
    "-DREBOUND_BRIDGE_BUILD_VALIDATION=ON"
)
```

- [ ] **Step 3: Commit documentation updates**

Run:

```powershell
git add README.md validation\solar_system_ias15\run.ps1
git diff --cached --name-only
git commit -m "docs: document vendored rebound build"
```

Expected:

```text
Commit includes README.md and run.ps1 only if run.ps1 changed.
```

## Task 6: Final Verification

**Files:**
- No source changes expected.

- [ ] **Step 1: Verify no default build reference to external REBOUND remains**

Run:

```powershell
rg -n "FetchContent|REBOUND_SRC_DIR|rebound_upstream|_deps/rebound-upstream|Download REBOUND|path\\to\\rebound" CMakeLists.txt README.md cmake include src examples validation
```

Expected:

```text
No matches in default build/documentation paths.
Matches in historical notes are acceptable only if clearly marked non-default.
```

- [ ] **Step 2: Verify vendored files are still identical to source cache**

Run:

```powershell
$srcFiles = Get-ChildItem _deps\rebound-upstream\src -File
foreach ($file in $srcFiles) {
    $vendored = Join-Path "third_party\rebound\src" $file.Name
    $a = Get-FileHash -Algorithm SHA256 $file.FullName
    $b = Get-FileHash -Algorithm SHA256 $vendored
    if ($a.Hash -ne $b.Hash) { throw "Vendored file changed: $($file.Name)" }
}
"vendored source unchanged"
```

Expected:

```text
vendored source unchanged
```

- [ ] **Step 3: Rebuild from a fresh directory**

Run:

```powershell
cmake -S . -B build-vendored-final -DREBOUND_BRIDGE_BUILD_EXAMPLES=ON
cmake --build build-vendored-final --config Release
```

Expected:

```text
Configure and build succeed without network or external REBOUND path.
```

- [ ] **Step 4: Check working tree for unrelated files before final report**

Run:

```powershell
git status --short
```

Expected:

```text
Only pre-existing unrelated changes, if any, remain.
examples/solar_system_bridge.c may still be modified from before this work and must not be reverted unless the user asks.
```
