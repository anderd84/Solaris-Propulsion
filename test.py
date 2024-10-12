import numpy as np
from scipy import linalg
rg = np.random.default_rng()
import matplotlib.pyplot as plt
from scipy import integrate
#a = np.array([[1., 2.], [3., 4.]])
#linalg.inv(a)
#A = np.array([[1,2,3],[1212122.,5.,6.],[7.,10000.,9.]])
#B = np.array([[1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,9],[1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,9],[1,2,3,4,5,6,7,8,9,1,2,3,4,5,6,7,9]])
#C = linalg.inv(A)
#D = A@B  #Matrix Multiplication
#print(f"{np.zeros([2,10])}")
#E  = np.empty([3,4])
#F = [[[11,12,13],[21,22,23],[31,32,33]],[[11,12,13],[21,22,23],[31,32,33]],[[11,12,13],[21,22,23],[31,32,33]]]
#G = np.empty([1,2])
#H = np.arange(0, 2*np.pi, 0.1)
#I = np.linspace(0,10, num=4)
#J = np.array([[1, 2+8j], [3, 4]], dtype=complex)
#K = np.array([1,2,3])
#K = K**2 #element wise Squaring
#L = 10 * np.sin(np.empty([1,2]))
#M = np.empty([3,3])
#M = np.dot(M, M)
#Ninit = 5*rg.random([10,3])  #random does from 0,1 by default
#N = np.exp(Ninit)
#O = np. arange(10)
#O = O[::-1]  # reversed a
#def f(x, y):
#    return 10 * (x+1) + y+1
#P = np.fromfunction(f, (4, 4))
#Q = rg.random(24).reshape(12,2) #If a dimension is given as -1 in a reshaping operation, the other dimensions are automatically calculated:

#zeta = 0.05
#omega_n = 2.5 #rad/s
#x_0 = 0
#v_0 = 1
#omega_d = omega_n * np.sqrt(1 - zeta**2)
#t = np.linspace(0,4/(zeta*omega_n),num = 10000)
#A_1 = x_0
#A_2 = (v_0 + x_0 * omega_n * zeta)/omega_d
#Funct = np.exp(-zeta*omega_n*t) * (A_1 * np.cos(omega_d* t) + A_2 * np.sin(omega_d* t))
#plt.plot(t, Funct)

#zeta = 0.3
#omega_n = 2.5 #rad/s
#omega_d = omega_n * np.sqrt(1 - zeta**2)
#t = np.linspace(0,4/(zeta*omega_n),num = 10000)
#phi = np.arctan(np.sqrt(1 - zeta**2)/zeta)
#Funct = 1 - np.exp(-zeta*omega_n*t) * np.sin(omega_d* t + phi)
#plt.plot(t, Funct)


#m = 100  # kg
#k = 5000  # N/m
#c = 50  # Ns/m
#x_0 = -1  # Initial position
#v_0 = 2  # Initial velocity\
#omega_n = np.sqrt(k/m)
#zeta = c / (2 * np.sqrt(m*k))
#omega_d = omega_n* np.sqrt(1 - zeta**2)
#def damped_vibrations(t, y):
#    x, v = y
#    omega = omega_n
#    F = 24 * np.sin(omega*t)
#    dxdt = v
#    dvdt = -(c/m)*v - (k/m)*x + F
#    return [dxdt, dvdt]
#
#t_span = (0, 6/(zeta*omega_n))  # Start and end times
#print(t_span[-1], zeta, omega_n)
#t_eval = np.linspace(*t_span, 1000)  # Time points at which to store the solution
#
## Solve the ODE
#sol = integrate.solve_ivp(damped_vibrations, t_span, [x_0, v_0], method='RK45', t_eval=t_eval)
#
## Plotting the results
#plt.plot(sol.t, sol.y[0], label='Displacement (x)')
#plt.plot(sol.t, sol.y[1], label='Velocity (v)')
#plt.xlabel('Time (s)')
#plt.ylabel('Displacement and Velocity')
#plt.legend()
#plt.title('Damped Vibrations')
#plt.grid(True)
#plt.show()
