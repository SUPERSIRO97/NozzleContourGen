#include "Nozzle_Contour.hpp"
#include <iostream>
#include <cmath>
#include <fstream>
#include <filesystem>

// Constructor
NozzleContour::NozzleContour(double contraction_ratio, double expansion_ratio, double L_chamber,
                             double l_perc, double r_chamber,
                             std::vector<double> theta_values,
                             std::vector<double> theta_n_values,
                             std::vector<double> theta_e_values,
                             std::vector<double> r1_values,
                             std::vector<double> r2_values,
                             std::vector<double> r3_values,
                             int Nr,
                             const std::string& folderPath)
    : contraction_ratio(contraction_ratio),
      expansion_ratio(expansion_ratio),
      L_chamber(L_chamber),
      l_perc(l_perc),
      r_chamber(r_chamber),
      theta_values(theta_values),
      theta_n_values(theta_n_values),
      theta_e_values(theta_e_values),
      r1_values(r1_values),
      r2_values(r2_values),
      r3_values(r3_values),
      Nr(Nr),
      folderPath(folderPath)
{
    // Throat radius
    r_throat = r_chamber / std::sqrt(contraction_ratio);
    
    // Exit radius
    r_ext = std::sqrt(expansion_ratio) * r_throat;

    // Creation folder to store countours data
    if (!std::filesystem::exists(folderPath)) {
        std::filesystem::create_directories(folderPath);
        std::cout << "Created output folder: " << folderPath << std::endl;
    }
}

// Overload
NozzleContour::NozzleContour(double contraction_ratio, double expansion_ratio, double L_chamber,
                            double l_perc, double r_chamber,
                            std::vector<double> theta_values,
                            std::vector<double> theta_n_values,
                            std::vector<double> theta_e_values,
                            std::vector<double> r1_values,
                            std::vector<double> r2_values,
                            std::vector<double> r3_values,
                            const std::string& output_folder)
                : NozzleContour(contraction_ratio, expansion_ratio, L_chamber, l_perc, r_chamber,
                                theta_values, theta_n_values, theta_e_values,
                                r1_values, r2_values, r3_values,
                                1, output_folder) {}



// Generation uniformly spaced grid of num elements within [start, end]
std::vector<double> NozzleContour::linspace(double start, double end, std::size_t num) {

    if (num > 1) {
        std::vector<double> result(num);
        double step = (end - start) / (num - 1);
        for (std::size_t i = 0; i < num; ++i) {
            result[i] = start + i * step;
        }
        return result;

    } else if (num == 1) {
        return {end};

    } else {
        throw std::invalid_argument("linspace: 'num' must be at least 1.");
    }
}

//// GENERATION CONTOUR NOZZLE ////
void NozzleContour::DrawContour() {
    
    // Determine which variable is used as the main radius
    bool isR1 = !r1_values.empty();
    bool isR2 = !r2_values.empty();
    const std::vector<double>& mainValues = isR1 ? r1_values : r2_values;

    for (double theta_deg : theta_values)
    {
        double theta = theta_deg * M_PI / 180.0;

        // Loop over the main (reference) radius
        for (double main_r : mainValues)
        {
            std::vector<double> depValues;

            if (isR1 && isR2)
            {
                // If both r1 and r2 are defined, sice r1_value is the main radius,
                // assign r2 to the dependent radius
                depValues = isR1 ? r2_values : r1_values;
            }
            else
            {
                // if only one radius is defined, create a vector of dependent radius of size Nr
                // fot the other
                // nmaximum allawable dependent radius tangent to the main radius with slope theta
                double r_dep_max = MaxRadius(r_throat, r_chamber, theta, main_r);
                depValues = linspace(0, r_dep_max, Nr);
            }

            // Loop over dependent variable
            for (double r_dep : depValues)
            {
                // Assign radii based on which one is main
                double r1_current = isR1 ? main_r : r_dep;
                double r2_current = isR1 ? r_dep : main_r;

                // Loop divergent radius
                for (double r3_current : r3_values)
                {
                    for (double theta_n_deg : theta_n_values)
                    {
                        double theta_n = theta_n_deg * M_PI / 180.0;

                        for (double theta_e_deg : theta_e_values)
                        {
                            double theta_e = theta_e_deg * M_PI / 180.0;

                            // Clear contour vectors for this combination
                            x_contour.clear();
                            y_contour.clear();

                            // --- Chamber section ---
                            DrawChamberSection();

                            // --- Convergent section ---
                            auto [x2, theta_actual_deg] = DrawConvergentSection(r1_current, r2_current, theta);

                            // --- Divergent section ---
                            DrawDivergentSection(x2, r3_current, theta_n, theta_e);

                            // Export CSV
                            ExportContourToCSV("nozzle_contour", r1_current, r2_current, r3_current,
                                            theta_actual_deg, theta_n_deg, theta_e_deg);
                        }
                    }
                }
            }
        }
    }
}








double NozzleContour::MaxRadius(double r_throat, double r_chamber, double theta, double r_input) {
    // This function Solve the following system:
    //
    // sqrt((x2 - x1)^2 + (y2 - y1)^2) == r1 + r2;
    // theta = pi/2 - atan((y2 - y1)/(x2 - x1));
    //
    // in order to find the maximum radius tangent to the other circle
    
    double A = r_throat - r_chamber;
    double B = std::tan(M_PI / 2 - theta);
    double C = 2 * A * (1 + B * B);
    double delta = C * C - 4 * A * A * (1 + B * B);

    if (delta < 0) {
        std::cerr << "No Solution found in MaxRadius" << std::endl;
        return 0.0;
    }

    double sum_radius = (-C + std::sqrt(delta)) / 2;
    return std::abs(sum_radius - r_input); // Absolute value is nedded !
}





// Draw Chamber portion
void NozzleContour::DrawChamberSection() {
    x_contour.insert(x_contour.end(), {0.0, L_chamber});
    y_contour.insert(y_contour.end(), {r_chamber, r_chamber});
}






// Draw Convergent portion
std::pair<double, double> NozzleContour::DrawConvergentSection(double r1_current, double r2_current, double theta_current) {

    // Tangent slope converging section
    double m = -std::tan(theta_current);

    // === First circle === //
    // First circle center
    double x1 = L_chamber;
    double y1 = r_chamber - r1_current;
     
    // Tangent point to the first circle
    double P1_x = x1 + r1_current * std::cos(M_PI / 2 - theta_current);
    double P1_y = y1 + r1_current * std::sin(M_PI / 2 - theta_current);

    // === Calculate L_conv (Distance between the two circles) === //
    // setting distance point center circle 2 to straigth line == r2
    // positive solution taken
    double yp = r_throat + r2_current;
    double L_conv = P1_x - L_chamber + (-r2_current * std::sqrt(1 + m * m) + (yp - P1_y)) / m;

    // === Second circle === //
    // Second circle center
    double x2 = L_chamber + L_conv;
    double y2 = yp;

    // Tangent point to the second circle
    double P2_x = x2 + r2_current * std::cos(3 * M_PI / 2 - theta_current);
    double P2_y = y2 + r2_current * std::sin(3 * M_PI / 2 - theta_current);

    // --------- Check for the possible options ---------//

    if(P2_y < P1_y) {
        //--- Draw a tangent line to both the circles of slope m ---//
        // Points of the tangent line: y = m * (x - P1_x) + P1_y

    } else if (P2_y >= P1_y) {
        //--- The curvatures of the 2 circles dont allow to trace a straigth line of slope m ---//
        // Connect the circles directly, changing the slope at the intersection point of the two
        
        // Recompute interested quantities
        double delta_y = y2 - y1;
        double delta_x = std::sqrt((r1_current + r2_current) * (r1_current + r2_current) - delta_y * delta_y);
        x2 = x1 + delta_x;
        theta_current = M_PI/2 - atan((y2 - y1) / (x2 - x1));
    }

    // Add circle 1 and circle 2 to the contour

    // === Circle 1 ===
    std::vector<double> theta_circle1 = linspace(M_PI / 2, M_PI / 2 - theta_current, 300);
    for (std::size_t i = 1; i < theta_circle1.size(); ++i) { // Skip the first point
        double t = theta_circle1[i];
        double x = x1 + r1_current * std::cos(t);
        double y = y1 + r1_current * std::sin(t);
        x_contour.push_back(x);
        y_contour.push_back(y);
    }

    // === Circle 2 ===
    std::vector<double> theta_circle2 = linspace(3 * M_PI / 2 - theta_current, 3 * M_PI / 2, 300);
    for (std::size_t i = 0; i < theta_circle2.size(); ++i) {
        double t = theta_circle2[i];
        double x = x2 + r2_current * std::cos(t);
        double y = y2 + r2_current * std::sin(t);
        x_contour.push_back(x);
        y_contour.push_back(y);
    }

    return {x2, theta_current * 180/M_PI};
}






// Draw Divergent portion
void NozzleContour::DrawDivergentSection(double x2, double r3_current, double theta_n_current, double theta_e_current) {
    
    // === Third circle === //
    // Third circle center
    double x3 = x2;
    double y3 = r_throat + r3_current;

    // === Calculate L_div (Distance between the throat and exit plene) === //
    // Convergent distance as l_perc of cone length of 15° aperture 
    double L_div = l_perc * (r_throat * (std::sqrt(expansion_ratio) - 1) + r3_current * (1 / std::cos(15 * M_PI / 180) - 1)) / std::tan(15 * M_PI / 180);


    // Add Circle 3 to the contour
    std::vector<double> theta_circle3 = linspace(3 * M_PI / 2, 3 * M_PI / 2 + theta_n_current, 300);
    for (std::size_t i = 1; i < theta_circle3.size(); ++i) {
        double t = theta_circle3[i];
        double x = x3 + r3_current * std::cos(t);
        double y = y3 + r3_current * std::sin(t);
        x_contour.push_back(x);
        y_contour.push_back(y);
    }
       
    // Trace Bezier curve
    DrawBezierCurve(x3, y3, r3_current, L_div, theta_n_current, theta_e_current);    
}






// Draw Bezier curve
void NozzleContour::DrawBezierCurve(double x3, double y3, double r3_current, double L_div, double theta_n, double theta_e) {

    // Initial end final point convergent portion
    double x_start = x3 + r3_current * std::cos(3 * M_PI / 2 + theta_n);
    double y_start = y3 + r3_current * std::sin(3 * M_PI / 2 + theta_n);
    double x_final = x3 + L_div;
    double y_final = r_ext;

    // Compute control point for Bezier curve as intersection between,
    // Line 1: y - y_start = tan(theta_n(j)) * (x - x_start);
    // Line 2: y - y_final = tan(theta_e(i)) * (x - x_final);
    double x_control = ((y_final - y_start) + std::tan(theta_n) * x_start - std::tan(theta_e) * x_final)
                     / (std::tan(theta_n) - std::tan(theta_e));
    double y_control = y_start + std::tan(theta_n) * (x_control - x_start);

    // === Bezier quadratic curve === //
    std::vector<double> t_values = linspace(0, 1, 300);
    for (std::size_t i = 1; i < t_values.size(); ++i) {
        double t = t_values[i];
        double x = (1 - t) * (1 - t) * x_start + 2 * (1 - t) * t * x_control + t * t * x_final;
        double y = (1 - t) * (1 - t) * y_start + 2 * (1 - t) * t * y_control + t * t * y_final;
        x_contour.push_back(x);
        y_contour.push_back(y);
    }
}




// Export the generated contour to a CSV file
// The first line contains all parameters used to generate this profile
// Each file gets a unique name using a counter to avoid overwriting
void NozzleContour::ExportContourToCSV(const std::string& base_filename,
                                       double r1, double r2, double r3,
                                       double theta_actual_deg, double theta_n_deg, double theta_e_deg) {
    
    // Create a unique filename using the counter
    std::string filename = folderPath + "/" + base_filename + "_" + std::to_string(contour_counter) + ".csv";
    contour_counter++; // Increment the counter for the next file

    std::ofstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Unable to open " << filename << std::endl;
        return;
    }

    // Write the header with parameters
    // Example: # r1 = 0.5, r2 = 0.3, r3 = 0.2, theta = 30°, theta_n = 10°, theta_e = 5°
    file << "# r1 = "  << r1 << ", r2 = " << r2 << ", r3 = " << r3
     << ", theta = "   << theta_actual_deg << "°"
     << ", theta_n = " << theta_n_deg << "°"
     << ", theta_e = " << theta_e_deg << "°\n";

    // Write the contour data (x, y)
    for (size_t i = 0; i < x_contour.size(); ++i) {
        file << x_contour[i] << "," << y_contour[i] << "\n";
    }

    file.close();
    std::cout << "Contour saved in: " << filename << std::endl;
}
