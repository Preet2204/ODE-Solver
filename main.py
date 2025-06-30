import numpy as np
import matplotlib.pyplot as plt
import matplotlib_inline.backend_inline 
import scienceplots
matplotlib_inline.backend_inline.set_matplotlib_formats('pdf', 'png')
plt.rcParams['figure.dpi']= 150
plt.rcParams['savefig.dpi'] = 150
plt.style.use(['science', 'no-latex', 'grid'])

class condition:
    def __init__(self, x, y):
        self.x = x
        self.y = y

def ForwardEuler(ode, init, n, rightx, show):
    leftx = init.x

    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Forward Euler Method.")

    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Forward Euler Method")

    dx = (rightx - leftx) / n
    x = np.linspace(leftx, rightx, n)
    y = np.zeros(n)
    y[0] = init.y
    for i in np.arange(1, n):
        y[i] = y[i - 1] + dx * ode(x[i - 1], y[i - 1])

    if(show):
        plt.plot(x, y, label='Forward Euler')
        plt.title("Plot of ODE with Initial Condition (%f, %f)" % (init.x, init.y))
        plt.legend()
    return x, y
    
def RungeKutta4(ode, init, n, rightx, show):
    leftx = init.x

    if leftx > rightx:
        raise Exception("rightx is Smaller than leftx. For Forward Euler Method.")

    if n <= 0:
        raise Exception("Number of Points(n) must be positive. For Forward Euler Method")

    dx = (rightx - leftx) / n
    x = np.linspace(leftx, rightx, n)
    y = np.zeros(n)
    y[0] = init.y
    for i in np.arange(1, n):
        
        K1 = dx * ode(x[i - 1], y[i - 1])
        K2 = dx * ode(x[i - 1] + dx / 2, y[i - 1] + K1 / 2)
        K3 = dx * ode(x[i - 1] + dx / 2, y[i - 1] + K2 / 2)
        K4 = dx * ode(x[i - 1] + dx, y[i - 1] + K3)
                
        y[i] = y[i - 1] + (K1 + 2 * K2 + 2 * K3 + K4) / 6
    
    if(show):
        plt.plot(x, y, label='Runge-Kutta 4')
        plt.title("Plot of ODE with Initial Condition (%f, %f)" % (init.x, init.y))
        plt.legend()
    return x, y
    

ode = lambda x, y: 2 * x            # ODE
init = condition(-10, 100)          # Initial Condition
n = 1000                            # Number of Points for Interpolation
rightx = 10

if __name__ == '__main__':
    
    x = np.linspace(init.x, rightx, n)
    x, eulery = ForwardEuler(ode, init, n, rightx, True)
    x, rky = RungeKutta4(ode, init, n, rightx, True)
    y = x ** 2
    
    plt.plot(x, y, label='Exact')
    plt.plot(x, np.abs(y - eulery), '--', label='Euler Error')
    plt.plot(x, np.abs(y - rky), ':', label='RK4 Error')
    plt.legend()
    plt.show()