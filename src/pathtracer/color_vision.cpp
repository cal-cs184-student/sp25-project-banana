#include "CGL/vector3D.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <mutex>

namespace CGL {
namespace ColorSpace {

// Debug logging system for color_vision.cpp
static std::ofstream cv_debug_log_file;
static std::mutex cv_debug_log_mutex;

// Initialize debug logging
static void init_cv_debug_log() {
    static bool initialized = false;
    if (!initialized) {
        cv_debug_log_file.open("/Users/ashvinverma/Documents/Courses/cs184/sp25-project-banana/color_vision_debug.log", std::ios::out | std::ios::trunc);
        initialized = true;
    }
}

// Debug logging macro
#define CV_DEBUG_LOG(message) \
    do { \
        std::lock_guard<std::mutex> lock(cv_debug_log_mutex); \
        init_cv_debug_log(); \
        if (cv_debug_log_file.is_open()) { \
            cv_debug_log_file << message << std::endl; \
            cv_debug_log_file.flush(); \
        } \
    } while (0)

// CIE 1931 standard observer color matching functions
// These are simplified approximations of the standard observer functions
Vector3D wave2xyz_table(double wavelength) {
    // Initialize to zero
    Vector3D xyz(0, 0, 0);
    
    // Only process wavelengths in the visible spectrum (380-780nm)
    if (wavelength < 380.0 || wavelength > 780.0) {
        return xyz;
    }

    // Approximation of CIE color matching functions
    // Based on Gaussian distributions centered at peak wavelengths
    
    // X curve - peaks around 600nm (red/orange)
    double x = std::exp(-0.5 * std::pow((wavelength - 600.0) / 80.0, 2));
    
    // Y curve - peaks around 550nm (green), corresponds to luminance sensitivity
    double y = std::exp(-0.5 * std::pow((wavelength - 550.0) / 70.0, 2));
    
    // Z curve - peaks around 450nm (blue/violet)
    double z = std::exp(-0.5 * std::pow((wavelength - 450.0) / 60.0, 2));
    
    // Scale factors to approximately normalize the response
    x *= 1.056;
    y *= 0.974;
    z *= 0.908;
    
    // Gamma correction factor to compensate for the approximation
    double gamma = 1.0;
    
    // Apply scaling based on wavelength to better match the standard observer
    if (wavelength >= 380.0 && wavelength < 420.0) {
        gamma = 0.3 + 0.7 * (wavelength - 380.0) / 40.0;
    } else if (wavelength >= 420.0 && wavelength < 700.0) {
        gamma = 1.0;
    } else if (wavelength >= 700.0 && wavelength <= 780.0) {
        gamma = 0.3 + 0.7 * (780.0 - wavelength) / 80.0;
    }
    
    // Apply gamma and store in xyz Vector3D
    xyz.x = x * gamma;
    xyz.y = y * gamma;
    xyz.z = z * gamma;
    
    CV_DEBUG_LOG("Wavelength " << wavelength << "nm -> XYZ: (" 
               << xyz.x << ", " << xyz.y << ", " << xyz.z << ")");
    
    return xyz;
}

// Convert XYZ color space to RGB using sRGB conversion matrix
Vector3D xyz2rgb(const Vector3D& xyz) {
    // sRGB conversion matrix (from XYZ to RGB)
    double r = 3.2406 * xyz.x - 1.5372 * xyz.y - 0.4986 * xyz.z;
    double g = -0.9689 * xyz.x + 1.8758 * xyz.y + 0.0415 * xyz.z;
    double b = 0.0557 * xyz.x - 0.2040 * xyz.y + 1.0570 * xyz.z;
    
    // Apply gamma correction for sRGB
    auto gamma_correct = [](double c) {
        if (c <= 0.0031308) {
            return 12.92 * c;
        } else {
            return 1.055 * std::pow(c, 1.0/2.4) - 0.055;
        }
    };
    
    r = gamma_correct(r);
    g = gamma_correct(g);
    b = gamma_correct(b);
    
    // Clamp values to [0, 1]
    r = std::max(0.0, std::min(1.0, r));
    g = std::max(0.0, std::min(1.0, g));
    b = std::max(0.0, std::min(1.0, b));
    
    CV_DEBUG_LOG("XYZ (" << xyz.x << ", " << xyz.y << ", " << xyz.z 
                << ") -> RGB: (" << r << ", " << g << ", " << b << ")");
    
    return Vector3D(r, g, b);
}

// Function to get a RGB color from a wavelength (nm)
// This gives an approximate visualization of spectral colors
Vector3D wavelength_to_rgb(double wavelength) {
    // Convert wavelength to XYZ color space
    Vector3D xyz = wave2xyz_table(wavelength);
    
    // Convert XYZ to RGB
    return xyz2rgb(xyz);
}

// Debug function to print wavelength color mappings
void print_spectral_colors() {
    CV_DEBUG_LOG("Wavelength to RGB mapping:");
    for (int wavelength = 380; wavelength <= 780; wavelength += 20) {
        Vector3D rgb = wavelength_to_rgb(wavelength);
        CV_DEBUG_LOG(wavelength << "nm: RGB(" 
                  << rgb.x << ", " << rgb.y << ", " << rgb.z << ")");
    }
}

} // namespace ColorSpace
} // namespace CGL