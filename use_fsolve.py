from scipy.optimize import fsolve
import numpy as np
import math
import matplotlib.pyplot as plt

def target_equations(vars, *data):

    # Set knowns
    sig_1, sig_2, sig_3, alpha_deg, beta_deg = data

    # Change degree to radian
    alpha = math.radians(alpha_deg)
    beta = math.radians(beta_deg)

    # Get unknowns
    tau_1, tau_2, tau_3, sig_c = vars

    # Construct equations
    eq1 = (sig_2 - sig_c) - (sig_1 - sig_c) * np.cos(2 * alpha) + tau_1 * np.sin(2 * alpha)
    eq2 = tau_2 - (sig_1 - sig_c) * np.sin(2 * alpha) - tau_1 * np.cos(2 * alpha)
    eq3 = (sig_3 - sig_c) - (sig_2 - sig_c) * np.cos(2 * beta) + tau_2 * np.sin(2 * beta)
    eq4 = tau_3 - (sig_2 - sig_c) * np.sin(2 * beta) - tau_2 * np.cos(2 * beta)

    # Return equation values
    return [eq1, eq2, eq3, eq4]


# Knowns 1: normal stress
sig_1 = 45.
sig_2 = 70.
sig_3 = 100.

# Knowns 2: angles between sensors
alpha_deg = 45.  # between 1 & 2
beta_deg = 45.   # between 2 & 3

# Set data
data = (sig_1, sig_2, sig_3, alpha_deg, beta_deg)

# Set initial guess: tau_1, tau_2, tau_3, sig_c
ini_guess = (5, 5, 5, 35)

# Solve nonlinear system equations with initial guess
tau_1, tau_2, tau_3, sig_c = fsolve(target_equations, ini_guess, args=data)

# Calculate
radius = math.sqrt((sig_1 - sig_c)**2 + tau_1**2)

# Validate the solution
eq_1, eq_2, eq_3, eq_4 = target_equations((tau_1, tau_2, tau_3, sig_c), *data)

# print solutions
print("Results:")
print("tau_1 = %0.5f" %tau_1)
print("tau_2 = %0.5f" %tau_2)
print("tau_3 = %0.5f" %tau_3)
print("sig_c = %0.5f" %sig_c)
print("R = %0.5f" %radius)

# print validation
print("\nValidation:")
print("eq_1 = %0.5f" %eq_1)
print("eq_2 = %0.5f" %eq_2)
print("eq_3 = %0.5f" %eq_3)
print("eq_4 = %0.5f" %eq_4)

# Set a list of angle of a circle
angle = np.linspace(0, 2 * np.pi, 150)

# Set (x, y) to construct a circle
x = sig_c + radius * np.cos(angle)
y = radius * np.sin(angle)



# Set limit values for the plot
xmin = 0
xmax = 1.5 * (sig_c + radius)
ymin = np.min([-1.5 * radius, -0.5 * (sig_c + radius)])
ymax = np.max([1.5 * radius, 0.5 * (sig_c + radius)])

# Set the plot
figure, axes = plt.subplots(1)
axes.set_aspect(1)
axes.set_xlim(left=0)
axes.set_xlim(right=xmax)
axes.set_ylim(bottom=ymin)
axes.set_ylim(top=ymax)
axes.set_title('Mohr Circle')
axes.set_xlabel('sig')
axes.set_ylabel('tau')
axes.axhline(0, color='black', lw=1)    # Draw a horizontal line at tau = 0
axes.plot(x, y)

# mark reading from sensor 1, 2, and 3 and center of Mohr circle
axes.plot(sig_1, tau_1, marker="o", color="b")
axes.plot(sig_2, tau_2, marker="o", color="b")
axes.plot(sig_3, tau_3, marker="o", color="b")
axes.plot(sig_c, 0, marker="o", color="r")

space1 = 0.01 * (sig_c + radius)
space2 = 0.05 * (sig_c + radius)

plt.annotate('sig_1 = %0.2f \ntau_1 = %0.2f' %(sig_1, tau_1),
             xy=(sig_1 + space1, tau_1 + space1),
             xytext=(sig_1 + space2, tau_1 + space2),
             arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=5),
             horizontalalignment='left', verticalalignment='bottom')

plt.annotate('sig_2 = %0.2f \ntau_2 = %0.2f' %(sig_2, tau_2),
             xy=(sig_2 + space1, tau_2 + space1),
             xytext=(sig_2 + space2, tau_2 + space2),
             arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=5),
             horizontalalignment='left', verticalalignment='bottom')

plt.annotate('sig_3 = %0.2f \ntau_3 = %0.2f' %(sig_3, tau_3),
             xy=(sig_3 + space1, tau_3 + space1),
             xytext=(sig_3 + space2, tau_3 + space2),
             arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=5),
             horizontalalignment='left', verticalalignment='bottom')

plt.annotate('sig_c = %0.2f' %sig_c,
             xy=(sig_c + space1, space1),
             xytext=(sig_c + space2, space2),
             arrowprops=dict(facecolor='black', width=1, headwidth=5, headlength=5),
             horizontalalignment='left', verticalalignment='bottom')


# Show the plot
plt.show()