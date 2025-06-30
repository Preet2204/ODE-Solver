import numpy as np
import matplotlib.pyplot as plt
import matplotlib_inline.backend_inline 
import scienceplots

# Set Matplotlib formats and styles for high-quality scientific plots
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'png')
plt.rcParams['figure.dpi']= 150
plt.rcParams['savefig.dpi'] = 150
plt.style.use(['science', 'no-latex', 'grid'])

class condition:
    def __init__(self, x, y):
        self.x = x              # Initial x value
        self.y = y              # Initial y value

fig, ax = plt.subplots()

def ForwardEuler(ode, init, n, rightx, show):
    leftx = init.x              # Set Left Boundary as initial x

    # Check if Left Boundary is Greater than Right Boundary
    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Forward Euler Method.")

    # Ensure Number of points is positive
    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Forward Euler Method")

    dx = (rightx - leftx) / n           # Step Size
    x = np.linspace(leftx, rightx, n)   
    y = np.zeros(n)
    y[0] = init.y                       # Set Initial Condition
    for i in np.arange(1, n):
        y[i] = y[i - 1] + dx * ode(x[i - 1], y[i - 1])      # Apply Forward Euler Step

    if(show):
        ax.plot(x, y, label='Forward Euler')
    
    return x, y
    
def RungeKutta4(ode, init, n, rightx, show):
    leftx = init.x              # Set Left Boundary as initial x

    # Check if Left Boundary is Greater than Right Boundary
    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Forward Euler Method.")

    # Ensure Number of points is positive
    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Forward Euler Method")

    dx = (rightx - leftx) / n           # Step Size
    x = np.linspace(leftx, rightx, n)   
    y = np.zeros(n)
    y[0] = init.y                       # Set Initial Condition
    for i in np.arange(1, n):
        
        K1 = dx * ode(x[i - 1], y[i - 1])                       # First Slope at current Point
        K2 = dx * ode(x[i - 1] + dx / 2, y[i - 1] + K1 / 2)     # Second Slope at midpoint
        K3 = dx * ode(x[i - 1] + dx / 2, y[i - 1] + K2 / 2)     # Third Slope at midpoint
        K4 = dx * ode(x[i - 1] + dx, y[i - 1] + K3)             # Fourth Slope at next Point
                
        y[i] = y[i - 1] + (K1 + 2 * K2 + 2 * K3 + K4) / 6       # Update y using weighted Average
    
    if(show):
        ax.plot(x, y, label='Runge-Kutta 4')
    return x, y
    

ode = lambda x, y: 2 * x            # ODE: dy/dx = 2x
init = condition(-10, 100)          # Initial Condition: y(-10) = 100
n = 1000                            # Number of Points for Interpolation
rightx = 10                         # Right Boundary

if __name__ == '__main__':

    x, eulery = ForwardEuler(ode, init, n, rightx, True)    # Compute Euler Solution
    x, rky = RungeKutta4(ode, init, n, rightx, True)        # Compute RK4 Solution
    y = x ** 2                                              # Exact Solution

    ax.plot(x, y, label='Exact')
    ax.plot(x, np.abs(y - eulery), '--', label='Euler Error')
    ax.plot(x, np.abs(y - rky), ':', label='RK4 Error')
    ax.legend()
    plt.title("Plot of ODE with Initial Condition (%f, %f)" % (init.x, init.y))
    plt.show()