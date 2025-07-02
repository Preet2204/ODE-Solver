import numpy as np
import matplotlib.pyplot as plt
import matplotlib_inline.backend_inline 
import scienceplots

# Set Matplotlib formats and styles for high-quality scientific plots
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'png')
plt.rcParams['figure.dpi']= 150
plt.rcParams['savefig.dpi'] = 150
plt.style.use(['science', 'no-latex', 'grid'])
fig, ax = plt.subplots()


def ForwardEuler(ode, init, n: int, rightx: float, show: bool):
    # Implement the Forward Euler method for solving first-order ODEs
    # Parameters:
    #   ode: Function defining the ODE dy/dx = f(x, y)
    #   init: Object with x and y attributes for initial conditions
    #   n: Number of interpolation points
    #   rightx: Right boundary of the x-domain
    #   show: Boolean to toggle plot display
    leftx = init[0]                 # Set Left Boundary as initial x

    # Check if Left Boundary is Greater than Right Boundary
    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Forward Euler Method.")

    # Ensure Number of points is positive
    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Forward Euler Method")

    dx = (rightx - leftx) / (n - 1)               # Step Size
    x = np.linspace(leftx, rightx, n)   
    y = np.zeros(n)
    y[0] = init[1]                          # Set Initial Condition
    for i in np.arange(1, n):
        y[i] = y[i - 1] + dx * ode(x[i - 1], y[i - 1])      # Apply Forward Euler Step

    if(show):
        ax.plot(x, y, 'o-', label='Forward Euler')
    
    return x, y

def RungeKutta4(ode, order: int, init: np.ndarray, n: int, rightx: float):
    # Implement the Runge-Kutta 4 (RK4) method for solving higher-order ODEs
    # Parameters:
    #   ode: Function defining the ODE system (e.g., f(x, y, y', y'', ...))
    #   order: Order of the ODE (number of derivatives)
    #   init: NumPy array of initial conditions [x0, y0, y'0, y''0, ...]
    #   n: Number of interpolation points
    #   rightx: Right boundary of the x-domain
    # Returns: x and y arrays of the solution

    # Validate initial condition length
    if len(init) != order + 1:
        raise Exception("Initial condition must have order + 1 elements.")
    
    # Validate ODE argument count
    if order + 1 != ode.__code__.co_argcount:
        raise Exception("The ODE Function should have %d No. of Arguments in order (x, y, y', y'', ...)")

    leftx = init[0]
    
    # Check if Left Boundary is Greater than Right Boundary
    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Runge Kutta 4 Method.")

    # Ensure Number of points is positive
    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Runge Kutta 4 Method")

    dx = (rightx - leftx) / (n - 1)     # Step Size
    x = np.linspace(leftx, rightx, n)
    y = np.zeros(n)                     # Initialize y array for the primary variable
    z = np.zeros((order - 1, n))        # Initialize array for higher derivatives (y', y'', ...)
    k = np.zeros((4, order))            # Initialize RK4 coefficients for each state variable
    state = np.zeros(order)
    state_mid = np.zeros(order)
    state_final = np.zeros(order)
    y[0] = init[1]                      # Set initial value for y
    # Set initial values for higher derivatives
    for j in np.arange(0, order - 1):
        z[j][0] = init[j + 2]

    for i in np.arange(1, n):

        # if order > 1:
        #     k[0, 1:] = dx * 

        # First stage of RK4
        state = [y[i - 1]] + list(z[:, i - 1])
        k[0] = dx * np.array(list(z[:, i - 1]) + [ode(x[i - 1], *state)])

        # Second and third stages of RK4
        for m in np.arange(0, 2):
            state_mid = [y[i - 1] + k[m][0] / 2] + list(z[:, i - 1] + k[m][1:] / 2)
            k[m + 1] = dx * np.array(list(z[:, i - 1] + k[m][1:] / 2) + [ode(x[i - 1] + dx / 2, *state_mid)])

        # Fourth stage of RK4
        state_final = [y[i - 1] + k[2][0]] + list(z[:, i - 1] + k[2][1:])
        k[3] = dx * np.array(list(z[:, i - 1] + k[2][1:]) + [ode(x[i - 1] + dx, *state_final)])

        # Update solution using weighted average of RK4 stages
        y[i] = y[i - 1] + (k[0][0] + 2 * k[1][0] + 2 * k[2][0] + k[3][0]) / 6
        if order > 1:
            z[:, i] = z[:, i - 1] + (k[0][1:] + 2 * k[1][1:] + 2 * k[2][1:] + k[3][1:]) / 6

    return x, y

# ODE Definition
n = 10000                           # Number of Points for Interpolation
rightx = 2                          # Right Boundary
init2 = np.array([0, 1, 0])         # Initial Condition: [x0, y0, y'0]
ode2 = lambda x, y, y1: -1000 * y   # Define second-order ODE system: y'' = -1000y

ode3 = lambda x, y: x * y
init3 = np.array([0, 1])

if __name__ == '__main__':
    
    x, y = RungeKutta4(ode2, 2, init2, n, rightx)
    y_exact = np.cos(np.sqrt(1000) * x)
    # y_exact = np.exp(np.pow(x, 2) /6 2)
    plt.plot(x, y, label='RK4')
    plt.plot(x, y_exact, label='Exact')
    plt.plot(x, np.abs(y_exact - y), label='Error')
    ax.legend()
    plt.show()