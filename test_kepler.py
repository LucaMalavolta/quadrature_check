import kepler_exo as kp
import numpy as np
import matplotlib.pyplot as plt

xx = np.arange(0,1.000, 0.001)

single_orbit_RV = kp.kepler_RV_T0P(xx, np.pi/2, 1.00000, 1.0000, 0.00, np.pi/2)
plt.plot(xx, single_orbit_RV)
plt.show()
