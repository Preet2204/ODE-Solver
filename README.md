# Differential Equation Solver

## Overview
This repository contains a Python implementation of a Differential Equation Solver, designed to numerically approximate solutions to ordinary differential equations (ODEs) using the Forward Euler and Runge-Kutta 4 (RK4) methods. The project includes visualization of numerical solutions, exact solutions, and error analysis, developed as part of an academic exercise to enhance computational skills in Applied Mathematics. The solver now supports a dynamic time-step approach (`dt`) for improved flexibility and accuracy, featuring an adaptive step size control to optimize performance for higher-order and stiff ODEs.

## Features
- Implementation of the Forward Euler method for solving first-order ODEs.
- Implementation of the Runge-Kutta-Fehlberg 45 (RKF45) method for enhanced accuracy with adaptive step size control.
- Adaptive step size control within the RK4 method, adjusting `dt` dynamically based on a user-defined error tolerance to balance accuracy and efficiency.
- Visualization of numerical approximations, exact solutions, and error plots using Matplotlib with a scientific style via the `scienceplots` package.
- Error handling for invalid input parameters (e.g., negative time steps or mismatched initial conditions).
- Modular design with reusable functions, utilizing a dynamic `dt` instead of a fixed number of points.

## Installation

### Prerequisites
- Python 3.x
- pip (Python package manager)

### Dependencies
The following Python libraries are required:
- `numpy` (for numerical computations)
- `matplotlib` (for plotting)
- `scienceplots` (for scientific plot styling)

### Setup Instructions
1. Clone the repository to your local machine:
    ```bash
    git clone https://github.com/Preet2204/ODE-Solver
    cd ODE-Solver
2. Create a virtual environment (optional but recommended):
    ```bash
    python -m venv env
    source env/bin/activate  # On Windows: env\Scripts\activate
3. Install the required Dependencies
    ```bash
    pip install -r requirements.txt

## Usage
### Running the Script

1. Ensure all dependencies are installed as per the Installation section.
2. Execute the main script:
    ```bash
    python main.py

The script will generate a plot displaying the RK4 approximation, exact solution, and error curve for the default ODE (( x'' + 1000x = 0 )) with initial condition ( [0, 1, 0] ) over the interval ([0, 2]) using an initial time step ( dt = 0.0001 ) and a tolerance of ( 1 \times 10^{-6} )..

### Customizing the ODE

- Modify the ode lambda function in the script to define a different ODE (e.g., `lambda t, x, x1: -100 * x` for a second-order ODE).
- Adjust the init array to set new initial conditions (e.g., [0, 1, 0]).
- Change the dt value to alter the time step and right_t to extend the interval.

### Example Output

The plot will show:

-Solid line: RK4 approximation
-Solid line: Exact solution (( x = cos(sqrt{1000} t) ))
-Solid line: Absolute error between RK4 and exact solution
-Legend identifying each curve

## Project Structure

`main.py`: Contains the core implementation, including the condition class, `ForwardEuler` and `RungeKutta4` functions, and plotting logic.

`README.md`: This file, providing project documentation.

`requirements.txt`: List of dependencies.

`.gitignore`: Excludes virtual environment and temporary files.

## Contributing

This project is primarily for personal academic development. However, contributions are welcome. Please fork the repository, make changes, and submit a pull request with a clear description of your modifications.

## Contact
For questions or feedback, please contact Preet Siddhapura at preetsiddhapura2204@gmail.com.
