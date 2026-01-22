# Helicopter Blade Vibration Analysis

This project implements a numerical solution using the **Rayleigh-Ritz method** (Assumed Modes) to analyze the vibration characteristics of helicopter blades. It calculates natural frequencies and mode shapes for various configurations, considering rotation effects, mass distributions, and boundary conditions.

## Project Structure

The core logic is distributed across three main MATLAB scripts:

### 1. `vibration_analysis.m` (Main Script)
This is the entry point of the simulation. It performs the following:
*   **Parameter Definition:** Sets up physical constants (Blade radius, mass, stiffness, rotation speed, etc.).
*   **Case Management:** Iterates through 8 different scenarios:
    *   **Distributions:** Uniform vs. Non-Uniform mass and stiffness.
    *   **Boundary Conditions:** Hinged vs. Cantilevered.
    *   **Dynamics:** Non-Rotating vs. Rotating.
*   **Solver:** Calls the assembler functions, solves the generalized eigenvalue problem (`K * v = lambda * M * v`), and sorts the modes.
*   **Visualization:** Generates comparative plots for mode shapes and frequencies to analyze the effects of rotation and boundary conditions.
*   **Export:** Saves the calculated natural frequencies to `frecuencias_naturales.csv`.

### 2. `buildGamma.m`
Generates the **Assumed Mode Functions** (polynomials) required for the Rayleigh-Ritz method.
*   Constructs polynomial shape functions ($\gamma_i(x)$) that strictly satisfy the geometric boundary conditions.
*   **Hinged:** Enforces zero displacement at the root.
*   **Cantilevered:** Enforces zero displacement and zero slope at the root.
*   Uses symbolic math to ensure the polynomials also satisfy natural boundary conditions at the tip for higher accuracy.

### 3. `assembleMK.m`
Responsible for assembling the generalized **Mass (M)** and **Stiffness (K)** matrices.
*   Uses symbolic integration (`int`) to compute matrix elements based on the energy method.
*   **Mass Matrix:** Derived from the kinetic energy of the beam.
*   **Stiffness Matrix:** Derived from the potential energy (strain energy) and the work done by centrifugal forces (centrifugal stiffening).

## Requirements

*   MATLAB (R2020a or newer recommended)
*   **Symbolic Math Toolbox** (Required for symbolic integration and differentiation used in `assembleMK` and `buildGamma`)

## Usage

1.  Open `vibration_analysis.m`.
2.  (Optional) Adjust the `chi` parameter or other physical constants in the "Parameters Definition" section at the top of the file.
3.  Run the script:
    ```matlab
    vibration_analysis
    ```
4.  The script will:
    *   Display progress in the Command Window.
    *   Generate multiple figures comparing mode shapes.
    *   Save frequency results to `frecuencias_naturales.csv`.

## Author
Christian Gonzalez Perez