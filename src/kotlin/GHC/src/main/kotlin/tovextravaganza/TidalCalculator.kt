package tovextravaganza

import org.apache.commons.math3.ode.FirstOrderDifferentialEquations
import org.apache.commons.math3.ode.FirstOrderIntegrator
import org.apache.commons.math3.ode.nonstiff.DormandPrince853Integrator
import kotlin.math.*

/**
 * TidalCalculator - Computes tidal deformability (Love number k2 and dimensionless Lambda)
 * 
 * The tidal deformability describes how much a neutron star is deformed by an external
 * gravitational field. It's measured by the Love number k2 and the dimensionless tidal
 * deformability Lambda.
 * 
 * This is important for gravitational wave observations of binary neutron star mergers,
 * as GW170817 measured Lambda(1.4 Msun) < 800.
 * 
 * We solve the TOV equations with additional tidal field equations to compute:
 * - k2: tidal Love number (dimensionless)
 * - Lambda: dimensionless tidal deformability = (2/3)k2 * C^(-5)
 *   where C = M/R is the compactness
 * 
 * Units:
 * - Lambda: dimensionless
 * - k2: dimensionless
 * - Mass: Msun (solar masses)
 * - Radius: km
 * - All calculations in geometric units (c=G=1)
 */
class TidalCalculator(private val tovSolver: TOVSolver) {
    private val eos = tovSolver.eos
    private val rMax = tovSolver.rMax
    private val dr = 0.001  // Larger step size for tidal equations (less stringent)
    
    /**
     * Inner class for the coupled TOV-tidal equations
     * We integrate 4 coupled ODEs: M, p, H, beta
     * H and beta are related to the tidal field perturbation
     */
    private inner class TidalEquations : FirstOrderDifferentialEquations {
        override fun getDimension(): Int = 4  // Four equations: M, p, H, beta
        
        /**
         * Compute derivatives for the coupled TOV-tidal equations
         * @param t: radius r
         * @param y: state [M, p, H, beta] where H and beta are tidal perturbation variables
         * @param yDot: derivatives [dM/dr, dp/dr, dH/dr, dbeta/dr]
         */
        override fun computeDerivatives(t: Double, y: DoubleArray, yDot: DoubleArray) {
            val M = y[0]  // Mass
            val p = y[1]  // Pressure
            val H = y[2]  // Tidal field perturbation function
            val beta = y[3]  // dH/dr
            val r = t
            
            // Stop if pressure becomes negative
            if (p <= 0.0) {
                yDot[0] = 0.0
                yDot[1] = 0.0
                yDot[2] = 0.0
                yDot[3] = 0.0
                return
            }
            
            val e = eos.getEnergyDensity(p)
            
            // Standard TOV equations (same as before)
            yDot[0] = 4.0 * PI * r * r * e  // dM/dr
            
            if (r < 1e-10) {
                yDot[1] = 0.0
            } else {
                val tmp = 1.0 - 2.0 * M / r
                yDot[1] = -((e + p) * (M / (r * r) + 4.0 * PI * r * p)) / tmp  // dp/dr
            }
            
            // Tidal equations: dH/dr = beta
            yDot[2] = beta
            
            // Tidal equations: dbeta/dr (complex expression)
            val f = max(1.0, eos.getFdedp(p))  // Speed of sound squared (capped at 1)
            val F1 = 1.0 - 2.0 * M / r  // Factor from metric
            
            // First term in dbeta/dr expression
            val term1 = (2.0 / F1) * H * (
                -2.0 * PI * (5.0 * e + 9.0 * p + f * (e + p)) +
                3.0 / (r * r) +
                (2.0 / F1) * (M / (r * r) + 4.0 * PI * r * p) * (M / (r * r) + 4.0 * PI * r * p)
            )
            
            // Second term in dbeta/dr expression
            val term2 = (2.0 * beta / r / F1) * (
                -1.0 + M / r + 2.0 * PI * r * r * (e - p)
            )
            
            yDot[3] = term1 + term2  // dbeta/dr
        }
    }
    
    /**
     * Compute tidal Love number k2 and dimensionless Lambda
     * @param centralP: central pressure of the neutron star
     * @return TidalResult with radius, mass, Lambda, and k2, or null if computation fails
     */
    fun compute(centralP: Double): TidalResult? {
        val rStart = 1e-5  // Start slightly away from center
        val equations = TidalEquations()
        val integrator: FirstOrderIntegrator = DormandPrince853Integrator(dr * 1e-6, rMax, 1e-10, 1e-8)
        
        // Initial conditions for tidal equations
        val y = DoubleArray(4)
        y[0] = 0.0  // M(0) = 0
        y[1] = centralP  // p(0) = centralP
        y[2] = rStart * rStart  // H ~ r^2 near center (even parity perturbation)
        y[3] = 2.0 * rStart  // beta = dH/dr = 2r near center
        
        // Integrate outward to surface
        val nSteps = ((rMax - rStart) / dr).toInt() + 1
        var currentR = rStart
        var iSurf = nSteps - 1
        
        for (i in 0 until nSteps) {
            try {
                integrator.integrate(equations, currentR, y, currentR + dr, y)
            } catch (e: Exception) {
                break
            }
            
            val p = y[1]
            // Surface: pressure drops to EOS minimum
            if (p <= eos.pTable[0]) {
                iSurf = i
                break
            }
            
            currentR += dr
        }
        
        // Extract values at surface
        val R = rStart + iSurf * dr
        val MCode = y[0]
        val MSolar = MCode / TOVSolver.MSUN_CODE
        val HR = y[2]  // H at surface
        val betaR = y[3]  // beta at surface
        
        val C = MCode / R  // Compactness
        
        // Check for unphysical stars
        if (R < 1e-5 || MCode < 1e-10 || C >= 0.5) {
            return TidalResult(R, MSolar, 0.0, 0.0)
        }
        
        // Compute y_R = R * beta_R / H_R (boundary condition variable)
        val yR = R * betaR / (HR + 1e-30)
        
        // Compute compactness powers
        val C2 = C * C
        val C3 = C2 * C
        val C4 = C3 * C
        val C5 = C4 * C
        
        val F = 1.0 - 2.0 * C  // Factor from Schwarzschild metric
        
        // Love number k2 formula (from Hinderer 2008)
        val num = (8.0 / 5.0) * C5 * F * F * (2.0 + 2.0 * C * (yR - 1.0) - yR)
        
        val denom = (2.0 * C * (6.0 - 3.0 * yR + 3.0 * C * (5.0 * yR - 8.0)) +
                     4.0 * C3 * (13.0 - 11.0 * yR + C * (3.0 * yR - 2.0) + 2.0 * C2 * (1.0 + yR)) +
                     3.0 * F * F * (2.0 - yR + 2.0 * C * (yR - 1.0)) * ln(F + 1e-30))
        
        val k2 = if (abs(denom) < 1e-30) 0.0 else num / denom
        
        // Dimensionless tidal deformability
        val lambda = if (C > 0) (2.0 / 3.0) * k2 / C5 else 0.0
        
        return TidalResult(R, MSolar, lambda, k2)
    }
    
    /**
     * Data class to hold tidal calculation results
     */
    data class TidalResult(
        val radius: Double,      // Star radius (km)
        val massSolar: Double,   // Star mass (solar masses)
        val lambda: Double,      // Dimensionless tidal deformability
        val k2: Double           // Tidal Love number
    ) {
        fun isValid(): Boolean = massSolar > 0.01 && radius > 0 && lambda >= 0
    }
}
