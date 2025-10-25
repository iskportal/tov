package tovextravaganza

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator
import org.apache.commons.cli.*
import java.io.File
import kotlin.math.*

/**
 * TOVSolver - Solves the Tolman-Oppenheimer-Volkoff (TOV) equations
 * 
 * The TOV equations describe the structure of spherically symmetric, static, relativistic stars.
 * They are a generalization of Newtonian stellar structure equations for general relativity.
 * 
 * The equations integrate from the center (r=0) outward until the surface is found:
 * - dM/dr = 4πr²ε (mass within radius r)
 * - dp/dr = -(ε+p)(M + 4πr³p) / [r²(1 - 2M/r)] (pressure gradient)
 * 
 * where M is the mass within radius r, p is pressure, and ε is energy density.
 * 
 * Units (geometric units, c=G=1):
 * - Radius r: km
 * - Mass M: km (code units) or Msun (solar masses, where Msun = 1.4766 km)
 * - Pressure p: km^-2
 * - Energy density ε: km^-2
 * - Central pressure: km^-2
 */
class TOVSolver(
    val eos: EOS,              // Equation of state (pressure-energy density relation)
    val rMax: Double = 100.0,  // Maximum integration radius (km)
    val dr: Double = 0.0005,   // Radial step size (km)
    val rtol: Double = 1e-12,  // Relative tolerance for ODE integration
    val atol: Double = 1e-14   // Absolute tolerance for ODE integration
) {
    companion object {
        // Mass of sun in code units (c=G=1) in km
        // MSUN = 1.989×10^30 kg = 1.4766 km
        const val MSUN_CODE = 1.4766
    }
    
    /**
     * Inner class defining the TOV differential equations
     * This implements the right-hand side of the ODE system
     */
    private inner class TOVEquations : FirstOrderDifferentialEquations {
        override fun getDimension(): Int = 2  // Two coupled equations: dM/dr and dp/dr
        
        /**
         * Compute the derivatives for the TOV equations
         * @param t: the independent variable (radius r)
         * @param y: state vector [M(r), p(r)]
         * @param yDot: output array for derivatives [dM/dr, dp/dr]
         */
        override fun computeDerivatives(t: Double, y: DoubleArray, yDot: DoubleArray) {
            val M = y[0]  // Mass within radius r
            val p = y[1]  // Pressure at radius r
            val r = t     // Current radius
            
            // Stop integration if pressure becomes negative
            if (p <= 0.0) {
                yDot[0] = 0.0
                yDot[1] = 0.0
                return
            }
            
            // Get energy density from EOS
            val e = eos.getEnergyDensity(p)
            
            // TOV Equation 1: dM/dr = 4πr²ε
            yDot[0] = 4.0 * PI * r * r * e
            
            // TOV Equation 2: dp/dr (handle singularity at r=0)
            if (r < 1e-10) {
                // At center, dp/dr ≈ 0
                yDot[1] = 0.0
            } else {
                // dp/dr = -(ε+p)(M/r² + 4πr·p) / (1 - 2M/r)
                val tmp = 1.0 - 2.0 * M / r
                yDot[1] = -((e + p) * (M / (r * r) + 4.0 * PI * r * p)) / tmp
            }
        }
    }
    
    /**
     * Solve the TOV equations for a neutron star with given central pressure
     * @param centralP: central pressure of the star
     * @return NeutronStar object with properties (mass, radius, etc.)
     */
    fun solve(centralP: Double): NeutronStar {
        val rStart = 1e-12  // Start integration slightly away from r=0 to avoid singularity
        val equations = TOVEquations()
        // Use Dormand-Prince 8(5,3) integrator (high-order adaptive RK method)
        val integrator: FirstOrderIntegrator = DormandPrince853Integrator(dr * 1e-6, rMax, atol, rtol)
        
        // Initial conditions: M(0)=0, p(0)=centralP
        val y = DoubleArray(2)
        y[0] = 0.0
        y[1] = centralP
        
        // Integrate from center outward
        val nSteps = ((rMax - rStart) / dr).toInt() + 1
        var currentR = rStart
        var iSurf = nSteps - 1
        
        // Step outward until we reach the surface
        for (i in 0 until nSteps) {
            try {
                integrator.integrate(equations, currentR, y, currentR + dr, y)
            } catch (e: Exception) {
                // Integration failed (e.g., hit event horizon or other singularity)
                break
            }
            
            // Check if we've reached the surface (pressure drops to EOS minimum)
            val p = y[1]
            if (p <= eos.pTable[0]) {
                iSurf = i
                break
            }
            currentR += dr
        }
        
        // Calculate final radius and mass
        val R = rStart + iSurf * dr
        val mCode = y[0]
        val mSolar = mCode / MSUN_CODE  // Convert to solar masses
        
        return NeutronStar(centralP, R, mCode, mSolar, eos)
    }
    
    /**
     * Solve TOV equations for a sequence of stars with different central pressures
     * This generates the complete mass-radius relationship
     * @param numStars: number of stars in the sequence
     * @param pMin: minimum central pressure (if null, use EOS minimum)
     * @param pMax: maximum central pressure (if null, use EOS maximum)
     * @return List of all neutron star solutions
     */
    fun solveSequence(numStars: Int = 500, pMin: Double? = null, pMax: Double? = null): List<NeutronStar> {
        // Determine pressure range
        val (pMinVal, pMaxVal) = if (pMin != null && pMax != null) {
            Pair(pMin, pMax)
        } else if (pMin != null) {
            val range = eos.getPressureRange()
            Pair(pMin, range.second)
        } else if (pMax != null) {
            val range = eos.getPressureRange()
            Pair(range.first, pMax)
        } else {
            val range = eos.getPressureRange()
            // Start slightly above zero to avoid numerical issues
            Pair(kotlin.math.max(range.first, 1e-15), range.second)
        }
        
        val stars = mutableListOf<NeutronStar>()
        // Use log spacing for pressure (better coverage of mass-radius curve)
        val logPMin = log10(pMinVal)
        val logPMax = log10(pMaxVal)
        
        // Generate stars with logarithmically spaced central pressures
        for (i in 0 until numStars) {
            val frac = if (numStars == 1) 0.0 else i.toDouble() / (numStars - 1)
            val logP = logPMin + frac * (logPMax - logPMin)
            val pC = 10.0.pow(logP)  // Convert back to linear scale
            
            try {
                val star = solve(pC)
                stars.add(star)
            } catch (e: Exception) {
                // Skip failed solutions (e.g., star too compact, unphysical)
            }
        }
        
        return stars
    }
}

/**
 * Data class representing a neutron star solution
 * Contains all the relevant physical properties
 */
data class NeutronStar(
    val centralPressure: Double,  // Central pressure of the star
    val radius: Double,            // Radius (km)
    val massCode: Double,          // Mass in code units (km, c=G=1)
    val massSolar: Double,         // Mass in solar masses
    val eos: EOS                   // EOS used for this calculation
) {
    /**
     * Compactness: C = M/R
     * Important for general relativistic effects
     */
    val compactness: Double
        get() = if (radius > 0) massCode / radius else 0.0
    
    /**
     * Check if this is a valid physical solution
     * @param minMass: minimum allowed mass in solar masses
     */
    fun isValid(minMass: Double = 0.01): Boolean {
        return massSolar > minMass && radius > 0
    }
    
    override fun toString(): String =
        "NeutronStar(M=${massSolar.format(4)} Msun, R=${radius.format(4)} km, p_c=${centralPressure.toScientific()})"
}

// Helper functions for formatting numbers
private fun Double.format(decimals: Int): String = "%.${decimals}f".format(this)
private fun Double.toScientific(): String = String.format("%.3e", this)

/**
 * Main entry point - Command line interface for the TOV solver
 */
fun main(args: Array<String>) {
    // Set up command line options
    val options = Options()
    options.addOption("i", "input", true, "Input EOS file (CSV format)")
    options.addOption("o", "output", true, "Output folder (default: export/stars)")
    options.addOption("n", "num-stars", true, "Number of stars to compute (default: 200)")
    options.addOption("rmax", true, "Maximum radius (default: 100.0)")
    options.addOption("dr", true, "Radial step size (default: 0.0005)")
    options.addOption("h", "help", false, "Print help message")
    
    val parser = DefaultParser()
    
    try {
        val cmd = parser.parse(options, args)
        
        // Show help if requested or if no input file specified
        if (cmd.hasOption("h") || !cmd.hasOption("i")) {
            val formatter = HelpFormatter()
            formatter.printHelp("TOVSolver", options)
            return
        }
        
        // Parse command line arguments
        val inputFile = cmd.getOptionValue("i")
        val outputFolder = cmd.getOptionValue("o", "export/stars")
        val numStars = cmd.getOptionValue("n", "200").toInt()
        val rMax = cmd.getOptionValue("rmax", "100.0").toDouble()
        val dr = cmd.getOptionValue("dr", "0.0005").toDouble()
        
        // Load equation of state
        println("Loading EOS from $inputFile...")
        val eos = EOS.fromFile(inputFile)
        println("  ${eos.nPoints} data points")
        println("  Columns: ${eos.colNames.joinToString(", ")}")
        
        // Create solver
        val solver = TOVSolver(eos, rMax = rMax, dr = dr)
        
        // Set up output directory
        val outputDir = File(outputFolder)
        outputDir.mkdirs()
        
        // Solve for sequence of stars
        println("\nSolving TOV equations for $numStars stars...")
        val stars = solver.solveSequence(numStars)
        
        // Filter stars
        val validStars = stars.filter { it.isValid() }
        // Apply physical filters: R < 99 km (not hitting integration limit) and M > 0.05 Msun
        val filteredStars = validStars.filter { it.radius < 99.0 && it.massSolar > 0.05 }
        
        println("\nFiltered: kept ${filteredStars.size}/${validStars.size} physical solutions (R < 99 km, M > 0.05 Msun)")
        println("\nValid solutions: ${validStars.size}/${stars.size}")
        
        // Print maximum mass
        if (validStars.isNotEmpty()) {
            val maxStar = validStars.maxByOrNull { it.massSolar }
            println("Maximum mass: ${maxStar!!.massSolar} Msun at R=${maxStar.radius} km")
        }
        
        // Write results to CSV
        val baseName = File(inputFile).nameWithoutExtension
        val csvFile = File(outputDir, "${baseName}_stars.csv")
        
        csvFile.bufferedWriter().use { writer ->
            // Write header with units
            val extraCols = eos.colNames.filter { it != "p" && it != "e" }
            val header = listOf("p_c(km^-2)", "R(km)", "M(Msun)") + extraCols.map { "${it}(pc)" }
            writer.write(header.joinToString(","))
            writer.newLine()
            
            // Write data for each filtered star
            for (star in filteredStars) {
                val values = mutableListOf(
                    star.centralPressure.toScientific(),
                    "${star.radius}",
                    "${star.massSolar}"
                )
                
                // Add any extra EOS columns evaluated at central pressure
                for (col in extraCols) {
                    values.add("${eos.getValue(col, star.centralPressure)}")
                }
                
                writer.write(values.joinToString(","))
                writer.newLine()
            }
        }
        
        println("\nWrote results to: ${csvFile.absolutePath}")
        
        // Compute tidal deformability
        println("\nComputing tidal deformability...")
        val tidalCalc = TidalCalculator(solver)
        val tidalResults = filteredStars.mapNotNull { star ->
            val result = tidalCalc.compute(star.centralPressure)
            // Only keep valid tidal results (and apply R < 99 km filter)
            if (result?.isValid() == true && result.radius < 99.0) {
                Pair(star.massSolar, result.lambda)
            } else null
        }
        
        // Write tidal results
        if (tidalResults.isNotEmpty()) {
            val tidalCsvFile = File(outputDir, "${baseName}_tidal.csv")
            tidalCsvFile.bufferedWriter().use { writer ->
                writer.write("M(Msun),Lambda(dimensionless)\n")
                tidalResults.forEach { (m, l) ->
                    writer.write("$m,$l\n")
                }
            }
            println("Tidal results written to: ${tidalCsvFile.absolutePath}")
            
            // Create Lambda plot (log scale)
            val lambdaPlotFile = Plotter.createLambdaPlot(tidalResults, outputDir, baseName)
            println("Tidal plot saved to: ${lambdaPlotFile.absolutePath}")
        }
        
        // Generate plots
        println("\nGenerating plots...")
        val asciiPlot = Plotter.createAsciiPlot(csvFile)
        println(asciiPlot)
        
        val pngFile = Plotter.createPNGPlot(csvFile)
        println("PNG plot saved to: ${pngFile.absolutePath}")
        
    } catch (e: ParseException) {
        println("Error parsing command line: ${e.message}")
        val formatter = HelpFormatter()
        formatter.printHelp("TOVSolver", options)
    } catch (e: Exception) {
        println("Error: ${e.message}")
        e.printStackTrace()
    }
}
