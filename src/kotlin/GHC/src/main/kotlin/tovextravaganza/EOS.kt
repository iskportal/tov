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
class EOS(
    private val dataDict: Map<String, DoubleArray>,      // Numeric columns
    val colNames: List<String>,                          // All column names
    private val stringDict: Map<String, List<String>> = emptyMap()  // String columns
) {
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
     * Note: Does NOT update iLast to avoid interfering with other interpolation calls
     */
    fun getFdedp(p: Double): Double {
        // Return boundary values if outside range
        if (p <= pTable[0]) return fdedp[0]
        if (p >= pTable[nPoints-1]) return fdedp[nPoints-1]
        
        // Start search from last index (cache optimization, but don't update iLast)
        var i = if (iLast < nPoints - 1) iLast else 0
        while (i > 0 && p < pTable[i]) i--
        while (i < nPoints - 1 && p > pTable[i+1]) i++
        
        // Linear interpolation
        val pI = pTable[i]
        val pIp1 = pTable[i+1]
        val fI = fdedp[i]
        val fIp1 = fdedp[i+1]
        
        val frac = (p - pI) / (pIp1 - pI)
        // Don't update iLast here - let getValue() manage the cache
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
     * Get string value at given pressure (finds nearest point by energy density)
     * Uses energy density instead of pressure for phase transitions
     */
    private fun getStringValue(colName: String, p: Double): String {
        if (!stringDict.containsKey(colName)) {
            return ""
        }
        
        // Get energy density at this pressure (interpolated)
        val eAtP = getEnergyDensity(p)
        
        // Find nearest energy density index
        val eTable = dataDict["e"]!!
        var minIdx = 0
        var minDist = kotlin.math.abs(eTable[0] - eAtP)
        
        for (i in 1 until eTable.size) {
            val dist = kotlin.math.abs(eTable[i] - eAtP)
            if (dist < minDist) {
                minDist = dist
                minIdx = i
            }
        }
        
        return stringDict[colName]!![minIdx]
    }
    
    /**
     * Get all column values (except 'p' and 'e') at a given pressure
     * This interpolates all additional columns at the central pressure
     * Returns Any to handle both numeric (Double) and string values
     * @param p: pressure value
     * @return Map of column name to interpolated value
     */
    fun getAllValuesAtPressure(p: Double): Map<String, Any> {
        val result = mutableMapOf<String, Any>()
        // Get all columns except p and e (those are in the main CSV header already)
        val extraCols = colNames.filter { it != "p" && it != "e" }
        for (col in extraCols) {
            try {
                if (stringDict.containsKey(col)) {
                    result[col] = getStringValue(col, p)
                } else {
                    // Check if column exists in dataDict
                    if (dataDict.containsKey(col)) {
                        result[col] = getValue(col, p)
                    }
                }
            } catch (e: Exception) {
                // Skip columns that fail to interpolate (print warning for debugging)
                // System.err.println("Warning: Failed to interpolate column '$col' at p=$p: ${e.message}")
            }
        }
        return result
    }
    
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
            val (dataDict, colNames, stringDict) = readCsv(filename)
            return EOS(dataDict, colNames, stringDict)
        }
        
        /**
         * Read CSV file and parse the data
         * Handles headers in comments (lines starting with #)
         * Automatically identifies 'p' and 'e' columns (can be named differently)
         * Data is sorted by pressure for efficient interpolation
         */
        private fun readCsv(filename: String): Triple<Map<String, DoubleArray>, List<String>, Map<String, List<String>>> {
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
            
            // Normalize header column names (epsilon -> e, etc.)
            val headerList = header!!.map { col ->
                when {
                    col.equals("epsilon", ignoreCase = true) -> "e"
                    col.equals("pressure", ignoreCase = true) || col == "p" -> "p"
                    else -> col
                }
            }
            
            // Parse data - separate numeric and string columns
            val numericData = mutableMapOf<String, MutableList<Double>>()
            val stringData = mutableMapOf<String, MutableList<String>>()
            val columnTypes = mutableMapOf<String, String>()  // "numeric" or "string"
            
            headerList.forEach { 
                numericData[it] = mutableListOf()
                stringData[it] = mutableListOf()
            }
            
            rawData.forEach { row ->
                headerList.forEachIndexed { idx, col ->
                    if (idx < row.size) {
                        // Determine column type on first row
                        if (!columnTypes.containsKey(col)) {
                            try {
                                row[idx].toDouble()
                                columnTypes[col] = "numeric"
                            } catch (e: NumberFormatException) {
                                columnTypes[col] = "string"
                            }
                        }
                        
                        // Add value based on type
                        if (columnTypes[col] == "numeric") {
                            try {
                                numericData[col]!!.add(row[idx].toDouble())
                            } catch (e: NumberFormatException) {
                                // Try to parse as zero or NaN for invalid numeric values
                                val value = if (row[idx].trim().isEmpty() || row[idx].trim().equals("nan", ignoreCase = true)) {
                                    Double.NaN
                                } else {
                                    0.0
                                }
                                numericData[col]!!.add(value)
                            }
                        } else {
                            stringData[col]!!.add(row[idx].trim())
                        }
                    }
                }
            }
            
            // Sort data by pressure (required for interpolation)
            val sortedData = mutableMapOf<String, DoubleArray>()
            
            // Check if we have valid data
            if (rawData.isEmpty() || numericData.isEmpty() || numericData["p"] == null || numericData["p"]!!.isEmpty()) {
                throw IllegalArgumentException("No valid data found in CSV file")
            }
            
            val colArray = numericData["p"]!!.toDoubleArray()
            val sortIdx = colArray.indices.sortedBy { colArray[it] }
            
            numericData.forEach { (col, vals) ->
                if (vals.size == colArray.size) {
                    sortedData[col] = sortIdx.map { vals[it] }.toDoubleArray()
                } else {
                    // Size mismatch - use unsorted data
                    sortedData[col] = vals.toDoubleArray()
                }
            }
            
            // Sort string columns (only keep columns that actually have string data)
            val sortedStringData = mutableMapOf<String, List<String>>()
            stringData.forEach { (col, vals) ->
                if (vals.isNotEmpty()) {  // Only add if column has string data
                    if (vals.size == colArray.size) {
                        sortedStringData[col] = sortIdx.map { vals[it] }
                    } else {
                        // Size mismatch - use unsorted data
                        sortedStringData[col] = vals
                    }
                }
            }
            
            return Triple(sortedData, headerList, sortedStringData)
        }
    }
}
