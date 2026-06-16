# GR-potential probe

Standalone experiment: validates that REBOUNDx physics can be reimplemented as
a native bridge-style operator-split kick on top of **REBOUND 5**, without
depending on REBOUNDx (which currently requires REBOUND `<5.0.0` and is
incompatible with this project's REBOUND 5 base).

`probe_gr.c` applies the `gr_potential` acceleration — formula taken from
REBOUNDx 4.6.2 `gr_potential.c`, pure arithmetic, no REBOUNDx dependency — as a
symmetric half-kick around each WHFast drift step, the same structure the bridge
uses for its cross-kicks.

## Result

Mercury perihelion precession: **43.15 arcsec/century** (GR predicts ~42.98),
~0.4% agreement. Confirms the formula, the AU/yr/Msun unit handling (c =
63239.7263 AU/yr), and the operator-split coupling are all correct.

## Build

The `CMakeLists.txt` reuses the REBOUND 5 source the main build downloads into
`build-validation/_deps/`. Configure the main project once (or run the
validation), then:

```powershell
cmake -S . -B build
cmake --build build --config Release
./build/Release/gr_probe.exe
```

Adjust `REBOUND_SRC` in `CMakeLists.txt` if your REBOUND 5 source lives
elsewhere.
