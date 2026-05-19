param(
    [string]$RepoRoot = (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path,
    [string]$BuildDir = (Join-Path (Resolve-Path (Join-Path $PSScriptRoot "..\..")).Path "build-validation-bench"),
    [double]$Years = 100.0,
    [int]$Samples = 200,
    [double[]]$DtOuterDaysList = @(1.0, 2.0, 4.0),
    [int[]]$EarthMoonRatios = @(10, 20),
    [int[]]$JovianRatios = @(25, 50)
)

$ErrorActionPreference = "Stop"

cmake --build $BuildDir --config Release --target solar_system_ias15_compare

$exeCandidates = @(
    (Join-Path $BuildDir "Release\solar_system_ias15_compare.exe"),
    (Join-Path $BuildDir "solar_system_ias15_compare.exe")
)
$exePath = $exeCandidates | Where-Object { Test-Path $_ } | Select-Object -First 1
if (-not $exePath) {
    throw "Could not locate solar_system_ias15_compare.exe under $BuildDir"
}

$rootOutDir = Join-Path $PSScriptRoot "out"
$sweepOutDir = Join-Path $rootOutDir "parameter_sweep"
$caseOutDir = Join-Path $sweepOutDir "cases"
New-Item -ItemType Directory -Force -Path $caseOutDir | Out-Null
Get-ChildItem $caseOutDir -Filter *.csv -ErrorAction SilentlyContinue | Remove-Item -Force

$csvPath = Join-Path $rootOutDir "solar_system_ias15_compare.csv"
$rows = @()

Push-Location $RepoRoot
try {
    foreach ($dtOuterDays in $DtOuterDaysList) {
        foreach ($earthMoonRatio in $EarthMoonRatios) {
            foreach ($jovianRatio in $JovianRatios) {
                $output = & $exePath `
                    --years $Years `
                    --samples $Samples `
                    --dt-outer-days $dtOuterDays `
                    --earth-moon-ratio $earthMoonRatio `
                    --jovian-ratio $jovianRatio

                $bridgeTotal = [double](($output | Select-String "^bridge_total_seconds=").ToString().Split("=")[1])
                $ias15Total = [double](($output | Select-String "^ias15_total_seconds=").ToString().Split("=")[1])
                $overallSpeedup = [double](($output | Select-String "^overall_speedup=").ToString().Split("=")[1])

                $caseName = "years_{0:0000}_samples_{1:00000}_dt_{2}_em_{3}_jov_{4}.csv" -f `
                    [int][Math]::Round($Years),
                    $Samples,
                    ($dtOuterDays.ToString("0.###") -replace "\.", "p"),
                    $earthMoonRatio,
                    $jovianRatio
                Copy-Item $csvPath (Join-Path $caseOutDir $caseName) -Force

                $caseRows = Import-Csv $csvPath
                $last = $caseRows[-1]
                $maxEarthBarycenterError = (($caseRows | ForEach-Object { [double]$_.earth_barycenter_error_au }) | Measure-Object -Maximum).Maximum
                $maxJupiterBarycenterError = (($caseRows | ForEach-Object { [double]$_.jupiter_barycenter_error_au }) | Measure-Object -Maximum).Maximum
                $maxAbsLaplaceDiff = (($caseRows | ForEach-Object { [Math]::Abs([double]$_.laplace_diff_deg) }) | Measure-Object -Maximum).Maximum
                $maxAbsBridgeEnergyRel = (($caseRows | ForEach-Object { [Math]::Abs([double]$_.bridge_energy_rel) }) | Measure-Object -Maximum).Maximum

                $rows += [pscustomobject]@{
                    years = [double]$Years
                    samples = $Samples
                    dt_outer_days = [double]$dtOuterDays
                    earth_moon_ratio = $earthMoonRatio
                    jovian_ratio = $jovianRatio
                    earth_moon_dt_days = [double]$dtOuterDays / [double]$earthMoonRatio
                    jovian_dt_days = [double]$dtOuterDays / [double]$jovianRatio
                    bridge_total_seconds = $bridgeTotal
                    ias15_total_seconds = $ias15Total
                    overall_speedup = $overallSpeedup
                    final_bridge_energy_rel = [double]$last.bridge_energy_rel
                    final_earth_barycenter_error_au = [double]$last.earth_barycenter_error_au
                    final_jupiter_barycenter_error_au = [double]$last.jupiter_barycenter_error_au
                    final_laplace_diff_deg = [double]$last.laplace_diff_deg
                    max_abs_bridge_energy_rel = $maxAbsBridgeEnergyRel
                    max_earth_barycenter_error_au = $maxEarthBarycenterError
                    max_jupiter_barycenter_error_au = $maxJupiterBarycenterError
                    max_abs_laplace_diff_deg = $maxAbsLaplaceDiff
                    case_csv = $caseName
                }
            }
        }
    }
} finally {
    Pop-Location
}

$summaryCsv = Join-Path $sweepOutDir "parameter_sweep.csv"
$rows | Sort-Object bridge_total_seconds | Export-Csv -NoTypeInformation -Encoding UTF8 $summaryCsv

$bestRuntime = $rows | Sort-Object bridge_total_seconds | Select-Object -First 1
$bestSpeedup = $rows | Sort-Object overall_speedup -Descending | Select-Object -First 1
$baseline = $rows |
    Where-Object { $_.dt_outer_days -eq 1.0 -and $_.earth_moon_ratio -eq 20 -and $_.jovian_ratio -eq 50 } |
    Select-Object -First 1

$summaryLines = @(
    "Bridge parameter sweep summary",
    "years: $Years",
    "samples: $Samples",
    "cases: $($rows.Count)",
    "summary_csv: $summaryCsv",
    "best_runtime_dt_outer_days: $($bestRuntime.dt_outer_days)",
    "best_runtime_earth_moon_ratio: $($bestRuntime.earth_moon_ratio)",
    "best_runtime_jovian_ratio: $($bestRuntime.jovian_ratio)",
    "best_runtime_bridge_seconds: $($bestRuntime.bridge_total_seconds)",
    "best_runtime_speedup: $($bestRuntime.overall_speedup)",
    "best_speedup_dt_outer_days: $($bestSpeedup.dt_outer_days)",
    "best_speedup_earth_moon_ratio: $($bestSpeedup.earth_moon_ratio)",
    "best_speedup_jovian_ratio: $($bestSpeedup.jovian_ratio)",
    "best_speedup_value: $($bestSpeedup.overall_speedup)"
)
if ($baseline) {
    $summaryLines += @(
        "baseline_bridge_seconds: $($baseline.bridge_total_seconds)",
        "baseline_speedup: $($baseline.overall_speedup)"
    )
}

$summaryPath = Join-Path $sweepOutDir "parameter_sweep_summary.txt"
$summaryLines | Set-Content -Encoding UTF8 $summaryPath
$rows | Sort-Object bridge_total_seconds | Format-Table -AutoSize
