import numpy as np
from CalcISDists import CalcISDists as calcDists

# Fill data matricies
barates = np.array([[0.8, 1.0, 1.2, 2.4, 5.8], [0.2, 3.0, 0.3, 1.5, 11.8], [-1.9, 2.8, -1.5, 0.2, 10.2]])
lastdecadegt = np.array([-211, -85-29, 26])
aris2090 = np.array([[[0.12, 0.07, 0.21], [0.04, -0.06, 0.12]],\
					[[0.08, 0.04, 0.13], [0.05, -0.04, 0.13]],\
					[[0.08, 0.04, 0.13], [0.05, -0.04, 0.13]],\
					[[0.06, 0.04, 0.10], [0.05, -0.03, 0.14]]]) * 1000

# Run the distribution fitting process
(batheteais, bathetwais, bathetgis, arthetais, arthetgis, islastdecade) = calcDists(barates, lastdecadegt, aris2090)

# Print the results to the screen

print(batheteais)
print(bathetwais)
print(bathetgis)
print(arthetais[:,0])
print(arthetgis[:,0])

# Done!
exit()