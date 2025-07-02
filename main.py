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


def ForwardEuler(ode, init, dt: float, rightt: float, show: bool):
    # Implement the Forward Euler method for solving first-order ODEs
    # Parameters:
    #   ode: Function defining the ODE dx/dt = f(t, x)
    #   init: Object with t and x attributes for initial conditions
    #   dt: Size of Time step
    #   rightt: Right boundary of the t-domain
    #   show: Boolean to toggle plot display
    leftt = init[0]                 # Set Left Boundary as initial t

    # Check if Left Boundary is Greater than Right Boundary
    if leftt > rightt:
        raise Exception("rightt is Smaller than leftt. For Forward Euler Method.")

    # Ensure Number of points is positive
    
    if dt <= 0:
        raise Exception("Time Step must be positive. For Forward Euler Method")

    n = int((rightt - leftt) / dt) + 1
    t = np.linspace(leftt, rightt, n)   
    x = np.zeros(n)
    x[0] = init[1]                          # Set Initial Condition
    for i in np.arange(1, n):
        x[i] = x[i - 1] + dt * ode(t[i - 1], x[i - 1])      # Apply Forward Euler Step

    if(show):
        ax.plot(t, x, label='Forward Euler')
    
    return t, x

def RungeKutta4(ode, order: int, init: np.ndarray, dt: float, rightt: float):
    # Implement the Runge-Kutta 4 (RK4) method for solving higher-order ODEs
    # Parameters:
    #   ode: Function defining the ODE system (e.g., f(t, x, x', x'', ...))
    #   order: Order of the ODE (number of derivatives)
    #   init: NumPy array of initial conditions [t0, x0, x'0, x''0, ...]
    #   dt: Size of Time step
    #   rightt: Right boundary of the t-domain
    # Returns: t and x arrays of the solution

    # Validate initial condition length
    if len(init) != order + 1:
        raise Exception("Initial condition must have order + 1 elements.")
    
    # Validate ODE argument count
    if order + 1 != ode.__code__.co_argcount:
        raise Exception("The ODE Function should have %d No. of Arguments in order (t, x, x', x'', ...)" % (order + 1))

    leftt = init[0]

    # Check if Left Boundary is Greater than Right Boundary
    if leftt > rightt:
        raise Exception("rightt is Smaller than leftt. For Runge Kutta 4 Method.")

    # Ensure Number of points is positive
    if dt <= 0:
        raise Exception("Time Step must be positive. For Runge Kutta 4 Method")

    n = int((rightt - leftt) / dt) + 1
    t = np.linspace(leftt, rightt, n)
    x = np.zeros(n)                     # Initialize x array for the primary variable
    z = np.zeros((order - 1, n))        # Initialize array for higher derivatives (x', x'', ...)
    k = np.zeros((4, order))            # Initialize RK4 coefficients for each state variable
    state = np.zeros(order)
    state_mid = np.zeros(order)
    state_final = np.zeros(order)
    x[0] = init[1]                      # Set initial value for x
    # Set initial values for higher derivatives
    for j in np.arange(0, order - 1):
        z[j][0] = init[j + 2]

    i = 1
    for i in np.arange(1, n):

        # First stage of RK4
        state = [x[i - 1]] + list(z[:, i - 1])
        k[0] = dt * np.array(list(z[:, i - 1]) + [ode(t[i - 1], *state)])

        # Second and third stages of RK4
        for m in np.arange(0, 2):
            state_mid = [x[i - 1] + k[m][0] / 2] + list(z[:, i - 1] + k[m][1:] / 2)
            k[m + 1] = dt * np.array(list(z[:, i - 1] + k[m][1:] / 2) + [ode(t[i - 1] + dt / 2, *state_mid)])

        # Fourth stage of RK4
        state_final = [x[i - 1] + k[2][0]] + list(z[:, i - 1] + k[2][1:])
        k[3] = dt * np.array(list(z[:, i - 1] + k[2][1:]) + [ode(t[i - 1] + dt, *state_final)])

        # Update solution using weighted average of RK4 stages
        x[i] = x[i - 1] + (k[0][0] + 2 * k[1][0] + 2 * k[2][0] + k[3][0]) / 6
        if order > 1:
            z[:, i] = z[:, i - 1] + (k[0][1:] + 2 * k[1][1:] + 2 * k[2][1:] + k[3][1:]) / 6
        
        t[i] = t[i - 1] + dt
    return t, x

# ODE Definition
dt = 0.0001
rightt = 2                          # Right Boundary
init2 = np.array([0, 1, 0])         # Initial Condition: [x0, y0, y'0]
ode2 = lambda t, x, x1: -1000 * x   # Define second-order ODE system: y'' = -1000y

ode3 = lambda t, x: t * x
init3 = np.array([0, 1])

if __name__ == '__main__':
    
    # t, x = ForwardEuler(ode3, init3, 0.02, rightt, True)
    
    t, x = RungeKutta4(ode2, 2, init2, dt, rightt)
    x_exact = np.cos(np.sqrt(1000) * t)
    # y_exact = np.exp(np.pow(x, 2) /6 2)
    plt.plot(t, x, label='RK4')
    plt.plot(t, x_exact, label='Exact')
    plt.plot(t, np.abs(x_exact - x), label='Error')
    ax.legend()
    plt.show()