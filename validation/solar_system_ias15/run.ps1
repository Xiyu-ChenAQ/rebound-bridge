param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$BuildDir = (Join-Path (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path "build-validation"),
    [string]$ReboundSrcDir = "",
    [double]$Years = 2000.0,
    [int]$Samples = 2000
)

$ErrorActionPreference = "Stop"

$outDir = Join-Path $PSScriptRoot "out"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

$cmakeArgs = @(
    "-S", $RepoRoot,
    "-B", $BuildDir,
    "-DREBOUND_BRIDGE_BUILD_VALIDATION=ON"
)
if ($ReboundSrcDir -and (Test-Path (Join-Path $ReboundSrcDir "rebound.h"))) {
    $cmakeArgs += "-DREBOUND_SRC_DIR=$ReboundSrcDir"
} else {
    $cmakeArgs += "-UREBOUND_SRC_DIR"
}

cmake @cmakeArgs
cmake --build $BuildDir --config Release --target solar_system_ias15_compare

$exeCandidates = @(
    (Join-Path $BuildDir "Release\solar_system_ias15_compare.exe"),
    (Join-Path $BuildDir "solar_system_ias15_compare.exe")
)
$exePath = $exeCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $exePath) {
    throw "Could not locate solar_system_ias15_compare.exe under $BuildDir"
}

Push-Location $RepoRoot
try {
    & $exePath --years $Years --samples $Samples
    python (Join-Path $PSScriptRoot "plot_solar_system_ias15.py") --csv (Join-Path $outDir "solar_system_ias15_compare.csv") --outdir $outDir
} finally {
    Pop-Location
}
