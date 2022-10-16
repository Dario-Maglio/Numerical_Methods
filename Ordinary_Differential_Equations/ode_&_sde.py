"""Basic ODE and SDE integrator script.

Before starting:
-) A and B define the temporal observation window.
-) xgrid is a vector with the time values on the x axis and the script allows to
use variable steps.
-) alpha and beta are the parameter of the desired Levy process for the SDE.

The SDE integrator starts from Ito equations in the form dx = f(x)*dt + g(x)*dw
- in which dw is a generic Levy process - and integrate it with the Heun's
algorithm (Stratonovich) and with the standard Ito formulation (better). If the
analitical solution is known, you can plot it too.

The odeint method uses the lsoda integrator from FORTRAN. lsoda differs from the
other integrators (except lsodar) in that it switches automatically between
stiff and nonstiff methods. This means that the user does not have to determine
whether the problem is stiff or not, and the solver will automatically choose
the appropriate method. It always starts with the nonstiff method.

Note:
An ordinary differential equation problem is stiff if the solution being sought
is varying slowly, but there are nearby solutions that vary rapidly, so the
numerical method must take small steps to obtain satisfactory results.

Ref:
https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.odeint.html
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.stats import levy_stable



# Parameters
A = 0.
B = 1.
STEPS = 100

# Filling the x axis
alpha = 2
beta = 0
xgrid = np.linspace(A, B, STEPS + 1)

#-------------------------------------------------------------------------------

def odesys(y, t, k1, k2, k3):
    """Define the system of equations in the form of dydt = M * y ."""
    y1, y2 = y
    dydt = [k1*y1 - k2*y1*y2,
            k2*y1*y2 - k3*y2]
    return dydt

def fun(x):
    """Define the stochastic Ito eq in the form dx = f(x)*dt + g(x)*dw.
    A drift toward the origin with decreasing fluctuations. If g(border) = 0 we
    talk of entry/exit  border depending on f(border). If f(border) = 0 we talk
    of natural border conditions.
    """
    f = -x
    g = x
    return [f, g]

def dfun(x):
    """Define the derivatives of the stochastic functions."""
    df = np.full_like(x, -1.)
    dg = np.full_like(x, 1.)
    return [df, dg]

def sfun(x, dt, z):
    """If an analitical solution to fun exists, we define it here."""
    return x * np.exp(z - 1.5 * dt)

def ito(x0, dt, z, fun, dfun):
    """Ito integration of SDE."""
    f0, g0 = fun(x0)
    dfdt0, dgdt0 = dfun(x0)
    x = x0 + (f0 - 0.5*g0*dgdt0)*dt + g0*z + 0.5*g0*dgdt0*(z**2)

    return x

def heun(x0, dt, z, fun, dfun):
    """Note: predictor-corrector uses by default the Strat integration."""
    f0, g0 = fun(x0)
    dfdt0, dgdt0 = dfun(x0)
    x1 = (f0 - 0.5*g0*dgdt0)*dt + g0*z

    f1, g1 = fun(x1)
    dfdt1, dgdt1 = dfun(x1)
    x2 = (f1 - 0.5*g1*dgdt1)*dt + g1*z

    return x0 + 0.5*(x1 + x2)

#-------------------------------------------------------------------------------

def ode_solver():
    plt.figure("ODE soultion")
    y0 = [0.5, 0.5]
    k1, k2, k3 = [10, 10, 10]
    sol = odeint(odesys, y0, xgrid, args=(k1, k2, k3))

    plt.plot(xgrid, sol[:, 0], 'b', label='prede')
    plt.plot(xgrid, sol[:, 1], 'g', label='predatori')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()

def sde_solver():
    plt.figure("SDE soultion")
    y0 = 0.5
    dt = [xgrid[i+1]-xgrid[i] for i in range(STEPS)]
    z = levy_stable.rvs(alpha, beta, size=(STEPS))
    z = np.power(dt, 1/alpha) * z
    # z = np.random.normal(size=(STEPS))
    # z = np.sqrt(dt) * z
    # assert (np.power(dt, 1/alpha) == np.sqrt(dt)).all()
    yito = np.full(STEPS + 1, y0)
    yheun = np.full(STEPS + 1, y0)
    yexact = np.full(STEPS + 1, y0)
    for index in range(STEPS):
        yexact[index + 1] = sfun(yexact[index], dt[index], z[index])
        yheun[index + 1] = heun(yheun[index], dt[index], z[index], fun, dfun)
        yito[index + 1] = ito(yito[index], dt[index], z[index], fun, dfun)

    plt.plot(xgrid, yexact, 'k', label='exact')
    plt.plot(xgrid, yheun, 'b', label='stoc heun')
    plt.plot(xgrid, yito, 'g', label='stoc ito')
    #plt.yscale('log')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()

#-------------------------------------------------------------------------------

if __name__=="__main__":

    ode_solver()

    sde_solver()
