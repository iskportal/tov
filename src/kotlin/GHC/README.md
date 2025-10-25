
1. **Get an EOS file** (CSV format with pressure and energy density)
   - Should have `p` and `e` columns
   - Units should be in km^-2 for both
   - Can have comment lines starting with `#`

2. **Run the solver:**
```bash
.\gradlew.bat run --args="-i path/to/eos.csv -n 200"
```

## Example: Running for CSC

To run with CSC equation of state:

```bash
.\gradlew.bat run --args="-i ..\..\..\src\python\inputCode\csc.csv -n 200 -o export/csc_output"
```

This will:
- Solve TOV equations for 200 stars
- Filter to keep stars with R < 99 km and M > 0.05 Msun
- Write results to CSV files
- Compute tidal deformability (Lambda)
- Generate plots (PNG and ASCII)

## Test Examples

**CSC EOS:**
```bash
.\gradlew.bat run --args="-i ..\..\..\src\python\inputCode\csc.csv -n 200 -o export/csc_test"
```

**HSDD2 EOS:**
```bash
.\gradlew.bat run --args="-i ..\..\..\src\python\inputCode\hsdd2.csv -n 200 -o export/hsdd2_test"
```

**DD2 EOS:**
```bash
.\gradlew.bat run --args="-i ..\..\..\src\python\inputCode\dd2.csv -n 200 -o export/dd2_test"
```

## Command Line Options

```
-i, --input        Input EOS file (required)
-o, --output       Output folder (default: export/stars)
-n, --num-stars    Number of stars to compute (default: 200)
--rmax             Maximum radius in km (default: 100.0 km)
--dr               Radial step size in km (default: 0.0005 km)
-h, --help         Show help
```

## Output Files

For input file `csc.csv`:
- `csc_stars.csv` - Mass-radius data (M in Msun, R in km, p_c in km^-2)
- `csc_tidal.csv` - Tidal deformability (M in Msun, Lambda dimensionless)
- `csc_stars.png` - Mass-radius plot (R in km, M in Msun)
- `csc_lambda.png` - Lambda vs Mass plot (log scale for Lambda)

## CSV Format

**Stars CSV headers:** `p_c(km^-2), R(km), M(Msun), [extra columns]`
- p_c(km^-2): central pressure
- R(km): radius
- M(Msun): mass in solar masses
- Additional EOS columns are included if present (e.g., mu(pc))

**Tidal CSV headers:** `M(Msun), Lambda(dimensionless)`
- M(Msun): mass in solar masses
- Lambda(dimensionless): dimensionless tidal deformability

**Example stars CSV:**
```
p_c(km^-2),R(km),M(Msun),mu(pc)
8.541e-06,3.389500000001,0.062208463058368225,344.75789771439605
...
```

**Example tidal CSV:**
```
M(Msun),Lambda(dimensionless)
0.062208463058368225,481425.7425454705
...
```

## Build

```bash
.\gradlew.bat build
```

Builds everything and creates a standalone JAR in `build/libs/`.

## Notes

- Filtering: Only stars with R < 99 km and M > 0.05 Msun are included in output
- Lambda plot uses log scale on the y-axis
- Uses Dormand-Prince 8th order integrator for ODE solving
- All calculations in geometric units (c=G=1)
