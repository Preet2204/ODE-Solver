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


def ForwardEuler(ode, init, dt: float, right_t: float, show: bool):
    # Implement the Forward Euler method for solving first-order ODEs
    # Parameters:
    #   ode: Function defining the ODE dx/dt = f(t, x)
    #   init: Object with t and x attributes for initial conditions
    #   dt: Size of Time step
    #   right_t: Right boundary of the t-domain
    #   show: Boolean to toggle plot display
    left_t = init[0]                 # Set Left Boundary as initial t

    # Check if Left Boundary is Greater than Right Boundary
    if left_t > right_t:
        raise Exception("right_t is Smaller than left_t. For Forward Euler Method.")

    # Ensure Time step is positive
    if dt <= 0:
        raise Exception("Time Step must be positive. For Forward Euler Method")

    n = int((right_t - left_t) / dt) + 1
    t = np.linspace(left_t, right_t, n)
    x = np.zeros(n)
    x[0] = init[1]                          # Set Initial Condition
    for i in np.arange(1, n):
        x[i] = x[i - 1] + dt * ode(t[i - 1], x[i - 1])      # Apply Forward Euler Step

    if(show):
        ax.plot(t, x, label='Forward Euler')

    return t, x

class fehlberg:
    a = np.array([
                [          0,            0,            0,           0,        0, 0],
                [      1 / 4,            0,            0,           0,        0, 0],
                [     3 / 32,       9 / 32,            0,           0,        0, 0],
                [1932 / 2197, -7200 / 2197,  7296 / 2197,           0,        0, 0],
                [  439 / 216,           -8,   3680 / 513, -845 / 4104,        0, 0],
                [    -8 / 27,            2, -3544 / 2565, 1859 / 4104, -11 / 40, 0]
                ])
    b5 = np.array([16 / 135, 0, 6656 / 12825, 28561 / 56430, -9 / 50, 2 / 55])
    b4 = np.array([25 / 216, 0,  1408 / 2565,   2197 / 4104,  -1 / 5,      0])
    c = np.array([0, 1 / 4, 3 / 8, 12 / 13, 1, 1 / 2])

class cashkarp:
    a = np.array([
                [           0,         0,           0,              0,          0, 0],
                [       1 / 5,         0,           0,              0,          0, 0],
                [      3 / 40,    9 / 40,           0,              0,          0, 0],
                [      3 / 10,   -9 / 10,       6 / 5,              0,          0, 0],
                [    -11 / 54,     5 / 2,    -70 / 27,        35 / 27,          0, 0],
                [1631 / 55296, 175 / 512, 575 / 13824, 44275 / 110592, 253 / 4096, 0]
                ])
    b5 = np.array([    37 / 378, 0,      250 / 621,       125 / 594,            0, 512 / 1771])
    b4 = np.array([2825 / 27648, 0,  18575 / 48384,   13525 / 55296,  277 / 14336,      1 / 4])
    c = np.array([0, 1 / 5, 3 / 10, 3 / 5, 1, 7 / 8])

def RKF45step(ode, order, dt, t, x, z):
    tableau = fehlberg()

    k = np.zeros(6)
    if order > 1: k = np.zeros((6, order))

    for i in np.arange(0, 6):

        if order > 1:
            ak = np.matmul(tableau.a[:, :], k[:, :])
            state = np.array(z[:] + ak[i, 1:])
            k[i, :-1] = dt * (state[:])
            k[i, -1] = dt * ode(t + tableau.c[i] * dt, x + ak[i, 0], *state)
        else:
            ak = k[:] * tableau.a[i, :]
            k[i] = dt * ode(t + tableau.c[i] * dt, x + ak[i])

    if order > 1:
        rk4 = np.matmul(tableau.b4[:], k[:, :])
        rk5 = np.matmul(tableau.b5[:], k[:, :])
    else:
        rk4 = np.matmul(tableau.b4[:], k[:])
        rk5 = np.matmul(tableau.b5[:], k[:])
    return rk4, rk5

def RKF45(ode, order: int, init, dt: float, right_t, tol=1e-10):
    # Implement the Runge-Kutta Fehlberg 45 (RKF45) method for solving higher-order ODEs
    # Parameters:
    #   ode: Function defining the ODE system (e.g., f(t, x, x', x'', ...))
    #   order: Order of the ODE (number of derivatives)
    #   init: NumPy array of initial conditions [t0, x0, x'0, x''0, ...]
    #   dt: Size of Time step
    #   right_t: Right boundary of the t-domain
    # Returns: t and x arrays of the solution

    # Validate initial condition length
    if len(init) != order + 1:
        raise Exception("Initial condition must have order + 1 elements.")

    # Validate ODE argument count
    if order + 1 != ode.__code__.co_argcount:
        raise Exception("The ODE Function should have %d No. of Arguments in order (t, x, x', x'', ...)" % (order + 1))

    left_t = init[0]

    # Check if Left Boundary is Greater than Right Boundary
    if left_t > right_t:
        raise Exception("right_t is Smaller than left_t. For Runge Kutta 4 Method.")

    # Ensure Time step is positive
    if dt <= 0:
        raise Exception("Time Step must be positive. For Runge Kutta 4 Method")

    t = np.array([left_t])
    x = np.array([init[1]])
    z = np.zeros((1, order - 1)) if order > 1 else np.array([])
    if order > 1: z[0, :] = init[2:]

    while t[-1] < right_t:
        dt = min(dt, right_t - t[-1])

        rk4, rk5 = RKF45step(ode, order, dt, t[-1], x[-1], z[-1] if order > 1 else np.array([]))

        rk5x = 0
        if order > 1: rk5x = x[-1] + rk5[0]
        else: rk5x = x[-1] + rk5

        rk5z = 0
        if order > 1: rk5z = z[-1] + rk5[1:]

        if order > 1: error = np.max(np.abs(rk4 - rk5))
        else: error = np.abs((rk4 - rk5) / (rk5 + 1e-10))

        if error < tol:
            x = np.concatenate((x, np.array([rk5x])))
            t = np.concatenate((t, np.array([t[-1] + dt])))
            if order > 1: z = np.concatenate((z, np.array([rk5z])))
        if error > 0:
            dt *= max(0.5, min(2, 0.9 * ((tol / error) ** (1 / 5))))
            dt = max(0.000001, min(dt, 0.001))
        print(t[-1], dt)

    return t, x

# ODE Definition
dt1 = 0.0001
order1 = 2
right_t1 = 2                          # Right Boundary
init1 = np.array([0, 1, 0])         # Initial Condition: [x0, y0, y'0]
ode1 = lambda t, x, x1: -1000 * x   # Define second-order ODE system: y'' = -1000y

dt2 = 0.0001
order2 = 1
right_t2 = 2
init2 = np.array([0, 1  ])
ode2 = lambda t, x: t * x

init3 = np.array([0, 1, 0, 0])
ode3 = lambda t, x, x1, x2: x2 - x

if __name__ == '__main__':

    t, x = RKF45(ode2, order2, init2, dt2, right_t2)
    # t, x = RKF45(ode1, order1, init1, dt1, right_t1)
    # print(t, x)
    plt.plot(t, x, label='RKF45')
    x_exact = np.exp(np.pow(t, 2) / 2)
    # x_exact = np.cos(np.sqrt(1000) * t)
    plt.plot(t, x_exact, label='Exact')
    plt.plot(t, np.abs(x_exact - x), label='Error')
    print(np.max(np.abs(x_exact - x)))
    ax.legend()
    plt.show()
