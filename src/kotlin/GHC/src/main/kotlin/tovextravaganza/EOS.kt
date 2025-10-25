package tovextravaganza

import java.io.BufferedReader
import java.io.FileReader
import kotlin.math.log10
import kotlin.math.pow

/**
 * EOS (Equation of State) class - handles reading CSV EOS files and interpolation
 * This class reads pressure-energy density tables and provides interpolation for any pressure value
 * 
 * Units:
 * - Pressure (p): km^-2
 * - Energy density (e): km^-2
 * - All values in geometric units (c=G=1)
 */
class EOS(private val dataDict: Map<String, DoubleArray>, val colNames: List<String>) {
    val pTable: DoubleArray      // Pressure table values
    val nPoints: Int              // Number of data points
    private var iLast = 0         // Cache last index for faster interpolation
    val fdedp: DoubleArray        // Pre-computed de/dp values for tidal calculations
    
    init {
        // Extract pressure table from data dictionary
        pTable = dataDict["p"]!!
        nPoints = pTable.size
        if (nPoints < 2) {
            throw IllegalArgumentException("Need at least 2 data points for interpolation.")
        }
        // Pre-compute de/dp for faster tidal calculations
        fdedp = precomputeFdedp()
    }
    
    /**
     * Pre-compute de/dp (energy density derivative with respect to pressure)
     * This is needed for tidal deformability calculations
     */
    private fun precomputeFdedp(): DoubleArray {
        val fdedp = DoubleArray(nPoints)
        val eTable = dataDict["e"]!!
        
        // Use central difference formula for interior points
        for (i in 1 until nPoints - 1) {
            fdedp[i] = (eTable[i+1] - eTable[i-1]) / (pTable[i+1] - pTable[i-1])
        }
        
        // Use boundary values for endpoints
        fdedp[0] = fdedp[1]
        fdedp[nPoints-1] = fdedp[nPoints-2]
        
        return fdedp
    }
    
    /**
     * Get de/dp at a given pressure using linear interpolation
     * Used in tidal deformability calculations
     */
    fun getFdedp(p: Double): Double {
        // Return boundary values if outside range
        if (p <= pTable[0]) return fdedp[0]
        if (p >= pTable[nPoints-1]) return fdedp[nPoints-1]
        
        // Start search from last index (cache optimization)
        var i = if (iLast < nPoints - 1) iLast else 0
        while (i > 0 && p < pTable[i]) i--
        while (i < nPoints - 1 && p > pTable[i+1]) i++
        
        // Linear interpolation
        val pI = pTable[i]
        val pIp1 = pTable[i+1]
        val fI = fdedp[i]
        val fIp1 = fdedp[i+1]
        
        val frac = (p - pI) / (pIp1 - pI)
        return fI + frac * (fIp1 - fI)
    }
    
    /**
     * Get value of any column at a given pressure using linear interpolation
     * This is the main interpolation method for pressure-energy density tables
     */
    fun getValue(colName: String, p: Double): Double {
        // Return boundary values if outside range
        if (p <= pTable[0]) return dataDict[colName]!![0]
        if (p >= pTable[nPoints-1]) return dataDict[colName]!![nPoints-1]
        
        // Search for the correct interval
        var i = iLast
        while (i > 0 && p < pTable[i]) i--
        while (i < nPoints - 1 && p > pTable[i+1]) i++
        
        // Linear interpolation between adjacent points
        val pI = pTable[i]
        val pIp1 = pTable[i+1]
        val cI = dataDict[colName]!![i]
        val cIp1 = dataDict[colName]!![i+1]
        
        val frac = (p - pI) / (pIp1 - pI)
        val result = cI + frac * (cIp1 - cI)
        
        // Cache the index for next call (optimization)
        iLast = i
        return result
    }
    
    /**
     * Get energy density at a given pressure
     * Convenience method for the TOV equations
     */
    fun getEnergyDensity(p: Double): Double = getValue("e", p)
    
    /**
     * Get the pressure range covered by this EOS
     * Returns (min pressure, max pressure)
     */
    fun getPressureRange(): Pair<Double, Double> = Pair(pTable[0], pTable[nPoints-1])
    
    companion object {
        /**
         * Load EOS from a CSV file
         * The CSV file should contain pressure (p) and energy density (e) columns
         * Additional columns are allowed and will be interpolatable
         */
        fun fromFile(filename: String): EOS {
            val (dataDict, colNames) = readCsv(filename)
            return EOS(dataDict, colNames)
        }
        
        /**
         * Read CSV file and parse the data
         * Handles headers in comments (lines starting with #)
         * Automatically identifies 'p' and 'e' columns (can be named differently)
         * Data is sorted by pressure for efficient interpolation
         */
        private fun readCsv(filename: String): Pair<Map<String, DoubleArray>, List<String>> {
            val rawData = mutableListOf<List<String>>()
            var header: List<String>? = null
            var lastCommentRow: List<String>? = null
            
            // Read the file line by line
            BufferedReader(FileReader(filename)).use { reader ->
                reader.lineSequence().forEach { line ->
                    // Handle comment lines (headers with units, etc.)
                    if (line.startsWith("#")) {
                        val parts = line.substring(1).trim().split(",").map { it.trim() }
                        if (parts.size >= 2) {
                            // Clean column names (remove units in parentheses, normalize names)
                            val cleanParts = parts.map { part ->
                                var cleaned = part.trim()
                                val idx = cleaned.indexOf('(')
                                if (idx > 0) {
                                    cleaned = cleaned.substring(0, idx).trim()
                                }
                                // Map common column name variations to standard names
                                when {
                                    cleaned.equals("pressure", ignoreCase = true) || cleaned == "p" -> "p"
                                    cleaned.equals("epsilon", ignoreCase = true) || cleaned == "e" || cleaned.equals("energy", ignoreCase = true) -> "e"
                                    else -> cleaned
                                }
                            }
                            lastCommentRow = cleanParts
                        }
                    } else if (line.isNotBlank()) {
                        // Parse data rows
                        val parts = line.split(",").map { it.trim() }
                        if (header == null) {
                            // Try to determine if this is a header or first data row
                            try {
                                parts[0].toDouble()
                                parts[1].toDouble()
                                // It's numeric data, use comment header if available
                                if (lastCommentRow != null && lastCommentRow!!.size == parts.size) {
                                    header = lastCommentRow
                                } else {
                                    // Default to p, e, col2, col3, ...
                                    header = listOf("p", "e") + List(parts.size - 2) { "col${it}" }
                                }
                            } catch (e: NumberFormatException) {
                                // This is a header row
                                header = parts
                                return@forEach
                            }
                        }
                        rawData.add(parts)
                    }
                }
            }
            
            // Default header if none was found
            if (header == null) {
                val nCols = rawData.firstOrNull()?.size ?: 2
                header = listOf("p", "e") + List(nCols - 2) { "col${it}" }
            }
            
            // Parse numeric data
            val numericData = mutableMapOf<String, MutableList<Double>>()
            val headerList = header!!
            headerList.forEach { numericData[it] = mutableListOf() }
            
            rawData.forEach { row ->
                headerList.forEachIndexed { idx, col ->
                    if (idx < row.size) {
                        try {
                            numericData[col]!!.add(row[idx].toDouble())
                        } catch (e: NumberFormatException) {
                            // Skip invalid values
                        }
                    }
                }
            }
            
            // Sort data by pressure (required for interpolation)
            val sortedData = mutableMapOf<String, DoubleArray>()
            val colArray = numericData["p"]!!.toDoubleArray()
            val sortIdx = colArray.indices.sortedBy { colArray[it] }
            
            numericData.forEach { (col, vals) ->
                sortedData[col] = sortIdx.map { vals[it] }.toDoubleArray()
            }
            
            return Pair(sortedData, headerList)
        }
    }
}
