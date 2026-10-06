# Month 2 · Week 5 — Correction pass (2026-10-06)

**Why:** Week 16 Day 3 found that every SPP processor from Week 6 on computes satellite positions at the signal **reception** time instead of the **transmission** time (up to ~±60 m range error per satellite). All weeks are being reviewed in order before fixing.

**Convention:** wrong code gets a corrected copy named `*_crctd`; originals are never edited. Corrected apps write `*_crctd` outputs. Each fix is marked `CORRECTED` in the code. `diff --strip-trailing-cr original _crctd` shows exactly what changed.

**Data:** `BRDC00IGS_R_20240070000_01D_MN.rnx` (RINEX 3.04, mixed). 2024-01-07 = Sunday = start of GPS week 2296.

## Summary

| Day | File | Issues | Effect on results |
|---|---|---|---|
| 1 | `satellite_position.cpp` | none; time argument undocumented | — |
| 2 | `rinex_nav_parser.cpp` → `_crctd` | fixed 7-line skip, memory leak, sticky `setfill` | none (latent) |
| 3 | `all_visible_sats.cpp` → `_crctd` | 7-line skip, no ephemeris age limit, no health flag, "UTC" label, 6-digit ECEF print | none: identical 11 satellites |
| 4 | `all_visible_sats_Skyplot.cpp`, `plot_sky.py` → `_crctd` | Day 3 issues + CSV precision; legend ≠ colors, "24.0" labels | plot only |

## Day 1 — Satellite position (correct)
- All 9 steps match IS-GPS-200 Table 20-IV.
- **Root of the later bug:** IS-GPS-200 defines `t` as GPS time **at transmission** (reception time minus travel time). The code just says `t_gps`, so Week 6+ passed the reception time. Passing toe in `main()` is fine for a test.
- 75 ms of travel × ~3.9 km/s ≈ 290 m along-track → up to ~±60 m in range.
- The week wrap matters here: first epoch at 0 s, toe = 597600 s → tk = +7200 s, not −597600 s.

## Day 2 — Nav parser
1. **Fixed 7-line skip** after non-GPS records. GLONASS/SBAS records are 4 lines (GPS/GAL/BDS 8), so the loop lost step and swallowed following records. Fix: skip lines starting with a space and non-`G` epoch lines; never count.
2. `*(new double)` for spare fields leaked 2 doubles per record → local variables.
3. `setfill('0')` is sticky → `07200.0` in the table → reset with `setfill(' ')`.
- **No effect on results:** all 440 GPS records come first in the BRDC file (`G 440, R 1214, E 9791, C 1079, J 96, I 146, S 11243`). Both versions parse 440.
- Synthetic test with GPS interleaved: original 4/5, `_crctd` 5/5. AddressSanitizer: original leaks 64 bytes.
- `|r| = 26,244 km` is fine: r ranges over a(1 ± e) = 26,216–26,912 km.

## Day 3 — Visible satellites
- Correct: LLA→ECEF, ENU rotation, el/az formulas, week wrap.
- Fixed: 7-line skip; **age limit** (reject |t − toe| > 2 h); **health** column (G01 unhealthy this day); label 43200 s = **12:00:00 GPS = 11:59:42 UTC** (GPS − UTC = 18 s); ECEF printed with 3 decimals instead of `-5.22293e+06`.
- Real data: identical 11 satellites and angles, all healthy.
- Range 20,404 km at 74° to 25,008 km at 6° → **67–86 ms** of signal travel.

## Day 4 — Sky plot
- C++: Day 3 fixes + CSV written with `fixed, setprecision(3)` (a new `ofstream` defaults to 6 significant digits: `2.04035e+07`). Writes `visible_sats_crctd.csv` with health.
- Python: legend said 25° but code used 30° (G21 at 27° drawn orange); legend yellow ≠ plotted yellow → colors and legend now from **one** table. Labels `24.0` (pandas `iterrows()` upcasts to float) → `G24`. Title in GPS time. Unhealthy drawn gray. Writes `sky_plot_crctd.png/.pdf`.

## Lessons
1. Document **which time** a time-dependent function expects.
2. Parse by structure, not by counting lines.
3. Test with data that exercises the edge case; the parser "worked" only because GPS came first.
4. Check ephemeris health and age before use.
5. C++ stream manipulators are sticky; new streams start with defaults.
6. GPS time ≠ UTC (18 s since 2017).
7. A plot's colors and legend must come from one definition.

## Next: Week 6
First real pseudoranges, so the transmit-time bug starts changing results. Fix with `sat_geometry` (Week 16 Day 3) and compare original vs `_crctd` errors.
