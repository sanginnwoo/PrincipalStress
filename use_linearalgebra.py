# import used libraries
import numpy as np
import math
import matplotlib.pyplot as plt

# Knowns 1: normal stress
sig_1 = 80.
sig_2 = 60.
sig_3 = 40.

# Knowns 2: angles between sensors
alpha_deg = 90.  # between 1 & 2
beta_deg = 45.   # between 2 & 3

# Change degree to radian
alpha = math.radians(alpha_deg)
beta = math.radians(beta_deg)

# Construct a matrix
M = np.array([[np.cos(2*alpha)-1, np.sin(2*alpha), 0, 0],
              [np.sin(2*alpha), -np.cos(2*alpha), 1, 0],
              [np.cos(2*beta)-1, 0, np.sin(2*beta), 0],
              [np.sin(2*beta), 0, -np.cos(2*beta), 1]])

# Construct the RHS vector
U = np.array([[np.cos(2*alpha)*sig_1 - sig_2],
              [np.sin(2*alpha)*sig_1],
              [np.cos(2*beta)*sig_2 - sig_3],
              [np.sin(2*beta)*sig_2]])

# Get inverse matrix
IM = np.linalg.inv(M)

# Get solution
F = np.matmul(IM, U)
sig_c = F[0][0]
tau_1 = F[1][0]
tau_2 = F[2][0]
tau_3 = F[3][0]

# Calculate the radius of Mohr circle
radius = math.sqrt((sig_1 - sig_c)**2 + tau_1**2)

# print solutions
print("Results:")
print("tau_1 = %0.5f" %tau_1)
print("tau_2 = %0.5f" %tau_2)
print("tau_3 = %0.5f" %tau_3)
print("sig_c = %0.5f" %sig_c)
print("R = %0.5f" %radius)
print("sig_p_max = %0.5f" %(sig_c + radius))
print("sig_p_min = %0.5f" %(sig_c - radius))

# Set a list of angle of a circle
angle = np.linspace(0, 2 * np.pi, 150)

# Set (x, y) to construct a circle
x = sig_c + radius * np.cos(angle)
y = radius * np.sin(angle)

# Draw the Mohr circle
figure, axes = plt.subplots(1)
axes.plot(x, y)

# Set limit values for the plot
xmax = np.ceil(0.15 * (sig_c + radius))*10
xmin = 0
ymax = np.ceil(np.max([1.5 * radius, 0.5 * (sig_c + radius)])/10.)*10
ymin = -ymax

# Set aspect ratio as 1 to show the circle
axes.set_aspect(1)

# Apply limit values of the plot
axes.set_xlim(left=0)
axes.set_xlim(right=xmax)
axes.set_ylim(bottom=ymin)
axes.set_ylim(top=ymax)

# Set title of the plot
axes.set_title('Mohr Circle')
axes.set_xlabel(r'$\sigma$')
axes.set_ylabel(r'$\tau$')

# Draw a horizontal line at tau = 0
axes.axhline(0, color='black', lw=1)

# Draw grid lines
axes.grid(color="gray", alpha=.5, linestyle='--')

# Set tick direction
axes.tick_params(direction='in')

# mark reading from sensor 1, 2, and 3 and center of Mohr circle
axes.plot(sig_1, tau_1, marker="o", color="b")
axes.plot(sig_2, tau_2, marker="o", color="b")
axes.plot(sig_3, tau_3, marker="o", color="b")
axes.plot(sig_c, 0, marker="o", color="r")

# set space between a symbol and annotation
space = 0.02 * (sig_c + radius)

# draw annotations
plt.annotate(r'$\sigma_1 = %0.1f, \tau_1 = %0.1f$' % (sig_1, tau_1),
             xy=(sig_1 + space, tau_1 + space))
plt.annotate(r'$\sigma_2 = %0.1f, \tau_2 = %0.1f$' % (sig_2, tau_2),
             xy=(sig_2 + space, tau_2 + space))
plt.annotate(r'$\sigma_3 = %0.1f, \tau_3 = %0.1f$' % (sig_3, tau_3),
             xy=(sig_3 + space, tau_3 + space))
plt.annotate(r'$\sigma_c = %0.1f$' % sig_c,
             xy=(sig_c + space, space))

# draw a box for the annotation
mybox = {'facecolor': 'white', 'edgecolor': 'black'}

# print max and min principal stresses in the plot
plt.text(xmin + (xmax-xmin) * 0.95, ymin + (ymax-ymin) * 0.05,
         r'$\sigma_{pmax}  = %0.2f, \sigma_{pmin}  = %0.2f$'
         %(sig_c + radius, sig_c - radius),
         fontsize='12', bbox=mybox, color='black',
         horizontalalignment='right', verticalalignment='bottom')

# save the plot
plt.savefig("test.svg")

# show the plot
plt.show()