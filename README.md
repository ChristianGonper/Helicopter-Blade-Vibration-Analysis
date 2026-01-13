# Helicopter Blade Vibration Analysis

This project provides a comprehensive MATLAB solution for calculating the natural frequencies and mode shapes of a helicopter blade under various structural and operational conditions.

## Implementation Overview

The solution is divided into a modular set of functions and a main executable script:

1.  **`vibration_analysis.m`**: Main script that iterates through 8 scenarios (Uniform/Non-Uniform x Hinged/Cantilevered x Rotating/Non-Rotating).
2.  **`buildGamma.m`**: Generates assumed test functions (polynomials) that satisfy geometric boundary conditions.
3.  **`buildG.m`**: Computes the symbolic centrifugal tension $G(x)$ based on the mass distribution.
4.  **`assembleMK.m`**: Assembles the symbolic Mass ($M$) and Stiffness ($K$) matrices using the energy method.

## How to Run

1.  Open MATLAB and navigate to the project directory.
2.  Run the main script:
    ```matlab
    vibration_analysis
    ```
3.  The script will output progress for each case in the Command Window and generate 4 comparative figures.

## Key Scenarios Covered

| Case | Distribution | Boundary Condition | Rotation |
| :--- | :--- | :--- | :--- |
| 1 | Uniform | Hinged | 0 & $\Omega$ |
| 2 | Uniform | Cantilevered | 0 & $\Omega$ |
| 3 | Non-Uniform | Hinged | 0 & $\Omega$ |
| 4 | Non-Uniform | Cantilevered | 0 & $\Omega$ |

## Results Visualization

The script generates plots showing the first 4 mode shapes. 
- **Blue dashed lines**: Non-rotating case ($\Omega = 0$).
- **Red solid lines**: Rotating case ($\Omega = 253$ RPM).
- **Legends**: Display the natural frequency in Hz for each case.

---

**Hinged Blade Observation**: For the rotating hinged blade, the first mode frequency is approximately equal to the rotation frequency ($\lambda \approx 1$), corresponding to the rigid-body flapping mode.

**Non-Uniform Case**: Since mass and stiffness decrease linearly to zero at the tip, the natural frequencies shift significantly compared to the uniform case.
