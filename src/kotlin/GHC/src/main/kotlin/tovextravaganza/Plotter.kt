package tovextravaganza

import java.io.File
import java.util.Scanner
import java.awt.Color
import java.awt.Graphics2D
import java.awt.image.BufferedImage
import java.awt.BasicStroke
import javax.imageio.ImageIO
import kotlin.math.log10
import kotlin.math.pow

/**
 * Plotter - Generates plots of neutron star data
 * 
 * Creates three types of plots:
 * 1. ASCII plot: Simple text-based plot for terminal output
 * 2. PNG plot: High-quality mass-radius curve
 * 3. Lambda plot: Tidal deformability with log scale on y-axis
 * 
 * Units in plots:
 * - Mass: solar masses (Msun)
 * - Radius: km
 * - Lambda: dimensionless
 */
object Plotter {
    
    /**
     * Create a simple ASCII art plot of the mass-radius curve
     * Useful for quick visualization in the terminal
     * @param csvFile: Input CSV file with mass-radius data
     * @return Formatted ASCII string
     */
    fun createAsciiPlot(csvFile: File): String {
        val data = mutableListOf<Pair<Double, Double>>()
        
        // Read data from CSV (R in column 1, M in column 2)
        Scanner(csvFile).use { scanner ->
            if (scanner.hasNextLine()) scanner.nextLine()  // Skip header
            
            while (scanner.hasNextLine()) {
                val line = scanner.nextLine().trim()
                if (line.isEmpty()) continue
                
                val parts = line.split(",")
                if (parts.size >= 3) {
                    try {
                        val r = parts[1].toDouble()  // Radius (km)
                        val m = parts[2].toDouble()  // Mass (Msun)
                        data.add(Pair(r, m))
                    } catch (e: NumberFormatException) {
                        // Skip invalid lines
                    }
                }
            }
        }
        
        if (data.isEmpty()) {
            return "No data to plot"
        }
        
        // Find data ranges
        val rMin = data.minOf { it.first }
        val rMax = data.maxOf { it.first }
        val mMin = data.minOf { it.second }
        val mMax = data.maxOf { it.second }
        
        val rRange = rMax - rMin
        val mRange = mMax - mMin
        
        // Create 2D array for the plot
        val width = 70
        val height = 20
        val plot = Array(height) { Array(width) { ' ' } }
        
        // Plot data points as asterisks
        data.forEach { (r, m) ->
            val x = ((r - rMin) / rRange * (width - 1)).toInt().coerceIn(0, width - 1)
            val y = ((m - mMin) / mRange * (height - 1)).toInt().coerceIn(0, height - 1)
            plot[height - 1 - y][x] = '*'
        }
        
        // Build output string
        val sb = StringBuilder()
        sb.append("\n")
        sb.append("Mass-Radius Curve (ASCII Plot)\n")
        sb.append("=".repeat(width + 20)).append("\n")
        
        for (row in plot) {
            sb.append("|").append(row.joinToString("")).append("|\n")
        }
        
        sb.append("=".repeat(width + 20)).append("\n")
        sb.append("Radius: ${String.format("%.2f", rMin)} - ${String.format("%.2f", rMax)} km\n")
        sb.append("Mass: ${String.format("%.2f", mMin)} - ${String.format("%.2f", mMax)} Msun\n")
        sb.append("\n")
        
        return sb.toString()
    }
    
    /**
     * Create a high-quality PNG plot of the mass-radius curve
     * @param csvFile: Input CSV file with mass-radius data
     * @return Output PNG file
     */
    fun createPNGPlot(csvFile: File): File {
        val data = mutableListOf<Pair<Double, Double>>()
        
        // Read data from CSV
        Scanner(csvFile).use { scanner ->
            if (scanner.hasNextLine()) scanner.nextLine()  // Skip header
            
            while (scanner.hasNextLine()) {
                val line = scanner.nextLine().trim()
                if (line.isEmpty()) continue
                
                val parts = line.split(",")
                if (parts.size >= 3) {
                    try {
                        val r = parts[1].toDouble()
                        val m = parts[2].toDouble()
                        data.add(Pair(r, m))
                    } catch (e: NumberFormatException) {
                    }
                }
            }
        }
        
        if (data.isEmpty()) {
            throw IllegalStateException("No data to plot")
        }
        
        // Image dimensions
        val width = 1200
        val height = 800
        val margin = 80
        val plotWidth = width - 2 * margin
        val plotHeight = height - 2 * margin
        
        // Create image
        val image = BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
        val g = image.createGraphics()
        
        // White background
        g.color = Color.WHITE
        g.fillRect(0, 0, width, height)
        
        // Enable anti-aliasing for smooth lines
        g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING, 
                          java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
        
        // Calculate data ranges
        val rMin = data.minOf { it.first }
        val rMax = data.maxOf { it.first }
        val mMin = data.minOf { it.second }
        val mMax = data.maxOf { it.second }
        val rRange = rMax - rMin
        val mRange = mMax - mMin
        
        // Draw grid lines
        g.color = Color.LIGHT_GRAY
        g.stroke = BasicStroke(0.5f)
        
        for (i in 0..5) {
            val x = margin + (i * plotWidth / 5)
            g.drawLine(x, margin, x, height - margin)
        }
        
        for (i in 0..5) {
            val y = margin + (i * plotHeight / 5)
            g.drawLine(margin, y, width - margin, y)
        }
        
        // Draw data points and line
        g.color = Color(0, 100, 200)
        g.stroke = BasicStroke(2.0f)
        
        val path = java.awt.geom.GeneralPath()
        var first = true
        
        data.forEach { (r, m) ->
            val x = margin + ((r - rMin) / rRange * plotWidth).toInt()
            val y = height - margin - ((m - mMin) / mRange * plotHeight).toInt()
            
            if (first) {
                path.moveTo(x.toFloat(), y.toFloat())
                first = false
            } else {
                path.lineTo(x.toFloat(), y.toFloat())
            }
            
            // Draw small circles at each point
            g.fillOval(x - 3, y - 3, 6, 6)
        }
        
        // Draw the connecting line
        g.stroke = BasicStroke(1.5f)
        g.draw(path)
        
        // Add title and labels
        g.color = Color.BLACK
        g.font = java.awt.Font("Arial", java.awt.Font.PLAIN, 14)
        
        val title = "Mass-Radius Curve - ${csvFile.nameWithoutExtension}"
        val titleMetrics = g.fontMetrics
        val titleWidth = titleMetrics.stringWidth(title)
        g.drawString(title, (width - titleWidth) / 2, 30)
        
        val xLabel = "Radius (km)"
        val xMetrics = g.fontMetrics
        val xWidth = xMetrics.stringWidth(xLabel)
        g.drawString(xLabel, (width - xWidth) / 2, height - 20)
        
        val yLabel = "Mass (solar masses)"
        val yMetrics = g.fontMetrics
        val yHeight = yMetrics.stringWidth(yLabel)
        
        // Rotate for y-axis label
        g.rotate(-Math.PI / 2)
        g.drawString(yLabel, -(margin + yHeight) / 2, 20)
        g.rotate(Math.PI / 2)
        
        // Add tick labels
        g.font = java.awt.Font("Arial", java.awt.Font.PLAIN, 10)
        
        for (i in 0..5) {
            val value = rMin + (rMax - rMin) * i / 5
            val label = String.format("%.1f", value)
            val x = margin + (i * plotWidth / 5) - xMetrics.stringWidth(label) / 2
            g.drawString(label, x, height - margin + 20)
        }
        
        for (i in 0..5) {
            val value = mMin + (mMax - mMin) * i / 5
            val label = String.format("%.2f", value)
            val y = height - margin - (i * plotHeight / 5) + yMetrics.height / 3
            g.drawString(label, 10, y)
        }
        
        g.dispose()
        
        // Save PNG file
        val outputFile = File(csvFile.parentFile, "${csvFile.nameWithoutExtension}.png")
        ImageIO.write(image, "png", outputFile)
        
        return outputFile
    }
    
    /**
     * Create a Lambda plot with logarithmic y-axis
     * Lambda (tidal deformability) varies over many orders of magnitude
     * @param data: List of (mass, lambda) pairs
     * @param outputDir: Output directory
     * @param baseName: Base name for output file
     * @return Output PNG file
     */
    fun createLambdaPlot(
        data: List<Pair<Double, Double>>,
        outputDir: File,
        baseName: String
    ): File {
        if (data.isEmpty()) {
            throw IllegalStateException("No data to plot")
        }
        
        // Image dimensions
        val width = 1200
        val height = 800
        val margin = 80
        val plotWidth = width - 2 * margin
        val plotHeight = height - 2 * margin
        
        // Create image
        val image = BufferedImage(width, height, BufferedImage.TYPE_INT_RGB)
        val g = image.createGraphics()
        
        // White background
        g.color = Color.WHITE
        g.fillRect(0, 0, width, height)
        
        // Enable anti-aliasing
        g.setRenderingHint(java.awt.RenderingHints.KEY_ANTIALIASING,
                          java.awt.RenderingHints.VALUE_ANTIALIAS_ON)
        
        // Calculate data ranges (log scale for lambda)
        val mMin = data.minOf { it.first }
        val mMax = data.maxOf { it.first }
        val lMin = data.minOf { it.second }
        val lMax = data.maxOf { it.second }
        val mRange = mMax - mMin
        val lRange = log10(lMax / lMin)  // Log range for lambda
        
        // Draw grid
        g.color = Color.LIGHT_GRAY
        g.stroke = BasicStroke(0.5f)
        
        for (i in 0..5) {
            val x = margin + (i * plotWidth / 5)
            g.drawLine(x, margin, x, height - margin)
        }
        
        for (i in 0..5) {
            val y = margin + (i * plotHeight / 5)
            g.drawLine(margin, y, width - margin, y)
        }
        
        // Draw data
        g.color = Color(200, 0, 100)
        g.stroke = BasicStroke(2.0f)
        
        val path = java.awt.geom.GeneralPath()
        var first = true
        
        data.forEach { (m, l) ->
            val x = margin + ((m - mMin) / mRange * plotWidth).toInt()
            // Use log scale for lambda
            val lambdaFrac = (log10(l / lMin) / lRange)
            val y = height - margin - (lambdaFrac * plotHeight).toInt()
            
            if (first) {
                path.moveTo(x.toFloat(), y.toFloat())
                first = false
            } else {
                path.lineTo(x.toFloat(), y.toFloat())
            }
            
            g.fillOval(x - 3, y - 3, 6, 6)
        }
        
        g.stroke = BasicStroke(1.5f)
        g.draw(path)
        
        // Add title and labels
        g.color = Color.BLACK
        g.font = java.awt.Font("Arial", java.awt.Font.PLAIN, 14)
        
        val title = "Tidal Deformability - $baseName"
        val titleMetrics = g.fontMetrics
        val titleWidth = titleMetrics.stringWidth(title)
        g.drawString(title, (width - titleWidth) / 2, 30)
        
        val xLabel = "Mass (solar masses)"
        val xMetrics = g.fontMetrics
        val xWidth = xMetrics.stringWidth(xLabel)
        g.drawString(xLabel, (width - xWidth) / 2, height - 20)
        
        val yLabel = "Lambda"
        val yMetrics = g.fontMetrics
        val yHeight = yMetrics.stringWidth(yLabel)
        
        // Rotate for y-axis label
        g.rotate(-Math.PI / 2)
        g.drawString(yLabel, -(margin + yHeight) / 2, 20)
        g.rotate(Math.PI / 2)
        
        // Add tick labels
        g.font = java.awt.Font("Arial", java.awt.Font.PLAIN, 10)
        
        // Mass ticks (linear scale)
        for (i in 0..5) {
            val value = mMin + (mMax - mMin) * i / 5
            val label = String.format("%.2f", value)
            val x = margin + (i * plotWidth / 5) - xMetrics.stringWidth(label) / 2
            g.drawString(label, x, height - margin + 20)
        }
        
        // Lambda ticks (log scale)
        for (i in 0..5) {
            val logValue = log10(lMin) + (log10(lMax) - log10(lMin)) * i / 5
            val value = 10.0.pow(logValue)
            val label = String.format("%.1e", value)
            val y = height - margin - (i * plotHeight / 5) + yMetrics.height / 3
            g.drawString(label, 10, y)
        }
        
        g.dispose()
        
        // Save PNG file
        val outputFile = File(outputDir, "${baseName}_lambda.png")
        ImageIO.write(image, "png", outputFile)
        
        return outputFile
    }
}
