param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$BuildDir = (Join-Path (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path "build-validation"),
    [string]$ReboundSrcDir = (Join-Path (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path "_deps\rebound-upstream\src"),
    [double]$Years = 2000.0,
    [int]$Samples = 2000
)

$ErrorActionPreference = "Stop"

$outDir = Join-Path $PSScriptRoot "out"
New-Item -ItemType Directory -Force -Path $outDir | Out-Null

cmake -S $RepoRoot -B $BuildDir -DREBOUND_SRC_DIR="$ReboundSrcDir" -DREBOUND_BRIDGE_BUILD_VALIDATION=ON
cmake --build $BuildDir --config Release --target solar_system_ias15_compare

Push-Location $RepoRoot
try {
    & (Join-Path $BuildDir "Release\solar_system_ias15_compare.exe") --years $Years --samples $Samples
    python (Join-Path $PSScriptRoot "plot_solar_system_ias15.py") --csv (Join-Path $outDir "solar_system_ias15_compare.csv") --outdir $outDir
} finally {
    Pop-Location
}
