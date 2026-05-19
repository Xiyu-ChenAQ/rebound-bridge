param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$BuildDir = (Join-Path (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path "build-validation-bench"),
    [double[]]$YearsList = @(10, 25, 50, 100, 200, 400, 800, 1200, 1600, 2000),
    [double]$SamplesPerYear = 4.0,
    [int]$MinSamples = 800,
    [double]$DenseYears = 400.0,
    [int]$DenseSamples = 12000
)

$ErrorActionPreference = "Stop"

$exeCandidates = @(
    (Join-Path $BuildDir "Release\solar_system_ias15_compare.exe"),
    (Join-Path $BuildDir "solar_system_ias15_compare.exe")
)
$exePath = $exeCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $exePath) {
    throw "Could not locate solar_system_ias15_compare.exe under $BuildDir"
}

$rootOutDir = Join-Path $PSScriptRoot "out"
$sweepOutDir = Join-Path $rootOutDir "efficiency_sweep"
$denseOutDir = Join-Path $sweepOutDir "dense"
$caseOutDir = Join-Path $sweepOutDir "cases"
New-Item -ItemType Directory -Force -Path $denseOutDir | Out-Null
New-Item -ItemType Directory -Force -Path $caseOutDir | Out-Null
Get-ChildItem $caseOutDir -Filter *.csv -ErrorAction SilentlyContinue | Remove-Item -Force
Get-ChildItem $denseOutDir -File -ErrorAction SilentlyContinue | Remove-Item -Force

$csvPath = Join-Path $PSScriptRoot "out\solar_system_ias15_compare.csv"
$rows = @()

Push-Location $RepoRoot
try {
    foreach ($years in $YearsList) {
        $samples = [Math]::Max($MinSamples, [int][Math]::Round($years * $SamplesPerYear))
        $output = & $exePath --years $years --samples $samples
        $bridgeTotal = [double](($output | Select-String "^bridge_total_seconds=").ToString().Split("=")[1])
        $ias15Total = [double](($output | Select-String "^ias15_total_seconds=").ToString().Split("=")[1])
        $overallSpeedup = [double](($output | Select-String "^overall_speedup=").ToString().Split("=")[1])

        $caseName = "years_{0:0000}_samples_{1:00000}.csv" -f [int][Math]::Round($years), $samples
        Copy-Item $csvPath (Join-Path $caseOutDir $caseName) -Force

        $rows += [pscustomobject]@{
            years = [double]$years
            samples = $samples
            bridge_total_seconds = $bridgeTotal
            ias15_total_seconds = $ias15Total
            overall_speedup = $overallSpeedup
        }
    }

    & $exePath --years $DenseYears --samples $DenseSamples | Out-Null
    $denseCsv = Join-Path $denseOutDir "solar_system_ias15_compare_dense.csv"
    Copy-Item $csvPath $denseCsv -Force
    python (Join-Path $PSScriptRoot "plot_solar_system_ias15.py") --csv $denseCsv --outdir $denseOutDir

    $summaryCsv = Join-Path $sweepOutDir "efficiency_sweep.csv"
    $rows | Export-Csv -NoTypeInformation -Encoding UTF8 $summaryCsv
    python (Join-Path $PSScriptRoot "plot_efficiency_sweep.py") --csv $summaryCsv --outdir $sweepOutDir
} finally {
    Pop-Location
}
