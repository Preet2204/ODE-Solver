# Differential Equation Solver

## Overview
This repository contains a Python implementation of a Differential Equation Solver, designed to numerically approximate solutions to ordinary differential equations (ODEs) using the Forward Euler and Runge-Kutta 4 (RK4) methods. The project includes visualization of the numerical solutions, exact solutions, and error analysis, developed as part of an academic exercise to enhance computational skills in Applied Mathematics.

## Features
- Implementation of the Forward Euler method for ODE solving.
- Implementation of the Runge-Kutta 4 method for higher accuracy.
- Visualization of numerical approximations, exact solutions, and error plots using Matplotlib with a scientific style.
- Error handling for invalid input parameters.
- Modular design with reusable functions.

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

The script will generate a plot displaying the Forward Euler approximation, Runge-Kutta 4 approximation, exact solution, and error curves for the default ODE ( dy/dx = 2x ) with initial condition ( y(-10) = 100 ) over the interval ([-10, 10]).

### Customizing the ODE

- Modify the `ode` lambda function in the script to define a different ODE (e.g., `lambda x, y: x + y`).
- Adjust the `init` object to set a new initial condition (e.g., `condition(0, 1)`).
- Change the `n` value to alter the number of interpolation points and `rightx` to extend the interval.

### Example Output

The plot will show:

- Black solid line: Forward Euler approximation
- Black solid line: Runge-Kutta 4 approximation
- Line: Exact solution (( y = x^2 ))
- Dashed line: Euler error
- Dotted line: RK4 error

## Project Structure

`main.py`: Contains the core implementation, including the condition class, `ForwardEuler` and `RungeKutta4` functions, and plotting logic.
`README.md`: This file, providing project documentation.
`requirements.txt`: List of dependencies.
`.gitignore`: Excludes virtual environment and temporary files.

## Contributing

This project is primarily for personal academic development. However, contributions are welcome. Please fork the repository, make changes, and submit a pull request with a clear description of your modifications.

## Contact
For questions or feedback, please contact Preet Siddhapura at preetsiddhapura2204@gmail.com.