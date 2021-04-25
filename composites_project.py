# Project done by: Saksham Kumar; M.Tech, IIT Kharagpur; B.Tech, IIT Patna

"""
Important: Please ensure that the file "data.txt" is present in the directory of this code.
{That file is used to input the basic properties and the volume fraction of the matrix and fibre as entering everything by hand would be tiresome!
(Please see point 2 of introduction below.)}

Introduction:- 
-> This project has been done to perform calculations for a composite material. 
-> We just have to specify the basic properties, volume fraction of the fibre and matrix used to make the composite and the thickness and layup sequence of the composite.
-> Rest calculations will be done by this Python code. 
-> Some of the key functions performed by this code are as follows:-
-> 1. It calculates the effective properties of the resulting composite. 
	  It plots the variation of the Young’s modulus, shear modulus and Q16 as a function of angle θ. 
   2. For a given layup sequence, it calculates the A, B and D matrix for the given thickness of layers and lamina.
   3. For some applied load, it calculates the margin of safety for each layer.
 	  It also predicts whether any of the layers would fail or not.
 	  (Hashin's failure criterion has been used for doing this.)

"""

# First, importing some necessary packages
import matplotlib.pyplot as plt
import math
import numpy as np
import os

# Printing welcome message for the user
print("Welcome to the program! This program lets you calculate different properties of your composite material. :)")

# The material properties and the values of some important parameters are stored in the file <data.txt>.
# The above mentioned file should be present in the same directory as that code for its proper functioning.
# The following command checks whether the file is present in the directory or not.
exists = os.path.isfile("data.txt")

# If the file exists, we open the file and extract the required data from it..
if exists:
	for line in open("data.txt"):
		if '#' not in line:
			# Ef1, Ef2, Gf, Em, Gm are in GPa. 
			Ef1 = (float(line.split("	")[1]))
			Ef2 = (float(line.split("	")[2]))
			nu_f = (float(line.split("	")[3]))
			Gf = (float(line.split("	")[4]))
			Em = (float(line.split("	")[5]))
			nu_m = (float(line.split("	")[6]))
			Vf = (float(line.split("	")[7]))
			numberOfLayers = (int(line.split("	")[8]))
			# total thickness in mm
			totalThickness = (float(line.split("	")[9]))
			# orientation sequence
			orientationSequence = []
			sequence = line.split("	")[10]
			for element in sequence.split(","):
				angle = float(element)
				angle_radians = (angle*math.pi)/180
				orientationSequence.append(angle_radians)

			Gm = (0.5*Em)/(1 + nu_m)

	# Function calculating basic properties of the composite

	def basic_func(Ef1, Ef2, Em, Vf):
		# Axial stiffness of the composite (using rule of mixtures)
		E_parallel = Vf * Ef1 + (1 - Vf) * Em

		# Transverse stiffness of the composite (using inverse rule of mixtures)
		E_perpendicular = (Ef2 * Em) / ((Ef2 * (1 - Vf)) + (Em * Vf))

		# Poission's raio of the composite (using rule of mixtures)
		nu12 = Vf * nu_f + (1 - Vf) * nu_m

		# Shear modulus of the composite (using inverse rule of mixtures)
		G12 = (Gf * Gm) / ((Gm * Vf) + (Gf * (1 - Vf)))

		print("Axial stiffness of the composite, E_parallel (in GPa):", E_parallel)
		print("Transverse stiffness of the composite, E_perpendicular (in GPa):", E_perpendicular)
		print("Poission's raio of the composite, nu_12:", nu12)
		print("Shear modulus of the composite, G12 (in GPa):", G12)

		return (E_parallel, E_perpendicular, nu12, G12)

	E_parallel, E_perpendicular, nu12, G12 = basic_func(Ef1, Ef2, Em, Vf)

	# Lets convert everything into SI units now.
	# From GPa to Pa
	Ef1 = Ef1 * 10**9
	Ef2 = Ef2 * 10**9
	Gf = Gf * 10**9
	Em = 3.45 * 10**9
	Gm = Gm * 10**9

	E_parallel = E_parallel * 10**9
	E_perpendicular = E_perpendicular * 10**9
	G12 = G12 * 10**9

	# converting thickness from mm to meter
	totalThickness = totalThickness * 0.001

	# print(E_parallel, E_perpendicular, nu12, G12)

	# ---------- 1st Assignment Completed ----------

	# To calculate the compliance and stiffness matrix for orthotropic materials under plane stress
	E1 = E_parallel
	E2 = E_perpendicular

	# Function calculating the compliance matrix 
	def compliance(E1, E2, nu12, G12):
		# First, let's calculate the different elements of the compliance matrix.
		S11 = 1/E1
		S12 = -nu12 / E1
		S16 = 0
		S21 = -nu12/E1
		S22 = 1/E2
		S26 = 0
		S61 = 0
		S62 = 0
		S66 = 1/G12

		# Now, we assemble the elements to get the compliance matrix, S.
		S = [[S11, S12, S16], [S21, S22, S26], [S61, S62, S66]]

		print("The compliance matrix (in local (1-2) coordinate system) is:",S)

		return S

	# Now, we can define a function to calculate the stiffness matrix.
	def stiffness(E1, E2, nu12, G12):
		S_matrix = compliance(E1, E2, nu12, G12)

		# Stiffness matrix
		Q = np.linalg.inv(S_matrix)

		print("The stiffness matrix (in local (1-2) coordinate system) is:",Q)

		return Q

	S = compliance(E1, E2, nu12, G12)
	Q = stiffness(E1, E2, nu12, G12)

	# print(S)
	# print(Q)

	# ----------2nd Assignment Completed ----------

	# To rotate the compliance matrix by an angle, i.e., to transform it from local to global coordinate system
	# In other words, to calculate Sxx, Sxy .... when we're given S11, S12 and so on.
	def rotateCompliance(S, theta):
		m = math.cos(theta)
		n = math.sin(theta)

		S11 = S[0][0]
		S12 = S[0][1]
		S16 = S[0][2]
		S21 = S[1][0]
		S22 = S[1][1]
		S26 = S[1][2]
		S61 = S[2][0]
		S62 = S[2][1]
		S66 = S[2][2]
		# Different components of the rotated compliance matrix
		Sxx = pow(m, 4)*S11 + pow(m, 2)*pow(n, 2)*(2*S12 + S66) + pow(n, 4)*S22
		Sxy = pow(m, 2)*pow(n, 2)*(S11 + S22 - S66) + S12*(pow(m,4) + pow(n,4))
		Syy = pow(n, 4)*S11 + pow(m, 2)*pow(n, 2)*(2*S12 + S66) + pow(m, 4)*S22
		Sxs = 2*pow(m, 3)*n*(S11 - S12) + 2*m*pow(n, 3)*(S12 - S22) - m*n*(pow(m, 2) - pow(n, 2))*S66
		Sys = 2*m*pow(n, 3)*(S11 - S12) + 2*pow(m, 3)*n*(S12 - S22) + m*n*(pow(m, 2) - pow(n, 2))*S66
		Sss = 4*pow(m, 2)*pow(n, 2)*(S11 - S12) - 4*pow(m, 2)*pow(n, 2)*(S12 - S22) + pow((pow(m, 2) - pow(n, 2)), 2)*S66

		Syx = Sxy
		Ssx = Sxs
		Ssy = Sys

		S_rotated = [[Sxx, Sxy, Sxs], [Syx, Syy, Sys], [Ssx, Ssy, Sss]]

		# print("The rotated compliance matrix in the X-Y coordinate system is:", S_rotated)

		return S_rotated

	# Similarly, doing transformation for the stiffness matrix
	def rotateStiffness(Q, theta):
		m = math.cos(theta)
		n = math.sin(theta)

		Q11 = Q[0][0]
		Q12 = Q[0][1]
		Q16 = Q[0][2]
		Q21 = Q[1][0]
		Q22 = Q[1][1]
		Q26 = Q[1][2]
		Q61 = Q[2][0]
		Q62 = Q[2][1]
		Q66 = Q[2][2]

		Qxx = pow(m, 4)*Q11 + 2*pow(m, 2)*pow(n, 2)*(Q12 + 2*Q66) + pow(n, 4)*Q22
		Qxy = pow(m, 2)*pow(n, 2)*(Q11 + Q22 - 4*Q66) + Q12*(pow(m, 4) + pow(n, 4))
		Qyy = pow(n, 4)*Q11 + 2*pow(m, 2)*pow(n, 2)*(Q12 + 2*Q66) + pow(m, 4)*Q22
		Qxs = pow(m, 3)*n*(Q11 - Q12) + m*pow(n, 3)*(Q12 - Q22) - 2*m*n*(pow(m, 2) - pow(n, 2))*Q66
		Qys = m*pow(n, 3)*(Q11 - Q12) + pow(m, 3)*n*(Q12 - Q22) + 2*m*n*(pow(m, 2) - pow(n, 2))*Q66
		Qss = pow(m, 2)*pow(n, 2)*(Q11 + Q22 - 2*Q12) + pow((pow(m, 2) - pow(n, 2)), 2)*Q66

		Qyx = Qxy
		Qsx = Qxs
		Qsy = Qys

		Q_rotated = [[Qxx, Qxy, Qxs], [Qyx, Qyy, Qys], [Qsx, Qsy, Qss]]
		# print("The rotated stiffness matrix in the X-Y coordinate system is:", Q_rotated)

		return Q_rotated

	# Now, transforming the engineering constants from local to global coordinate system based on angle
	def rotateEnggConstants(E1, E2, G12, nu12, Q):
		
		theta = np.linspace(0, math.pi/2, num=100)
		theta_degrees = np.linspace(0, 90, num=100)
		
		# We wish to plot the variations of the Young's modulus (E1, E2), shear modulus(G12) and Q16 as a function of theta.
		# Let's also plot the normalized values for E1, E2 and G12 so that we get good values on y axis instead of very big numbers.
		# This means that apart from variations of E1, E2, G12 and Q16 with theta, we'll plot Ex/E1, Ey/E2 and Gxy/G12 too.
		# This will also give better information about how many times our engineering constants are changing with the variation in angles.
		Ex_matrix = []
		Ey_matrix = []
		Gxy_matrix = []
		Qxs_matrix = []

		Ex_by_E1_matrix = []
		Ey_by_E2_matrix = []
		Gxy_by_G12_matrix = []

		for angle in theta:
			m = math.cos(angle)
			n = math.sin(angle)
			Ex = pow((m**4/E1 + m**2 * n**2 * (1/G12 - 2*nu12/E1) + n**4 / E2), -1)
			nuxy = Ex*((nu12*(m**4 + n**4))/E1 - m**2 * n**2 *(1/E1 + 1/E2 - 1/G12))
			Ey = pow((n**4 / E1 + m**2 * n**2 * (1/G12 - 2*nu12/E1) + m**4 / E2), -1)
			eetaxy_x = Ex*(m**3 * n * (2/E1 + 2*nu12/E1 - 1/G12) - m*n**3 * (2/E2 + 2*nu12/E1 - 1/G12))
			eetaxy_y = Ey*(m*n**3 * (2/E1 + 2*nu12/E1 - 1/G12) - m**3*n*(2/E2 + 2*nu12/E1 - 1/G12))
			Gxy = pow((4 * m**2 * n**2 * (1/E1 + 1/E2 + 2*nu12/E1) + (m**2 - n**2)**2 / G12), -1)

			Qxs = pow(m, 3)*n*(Q[0][0] - Q[0][1]) + m*pow(n, 3)*(Q[0][1] - Q[1][1]) - 2*m*n*(pow(m, 2) - pow(n, 2))*Q[2][2]

			Ex_by_E1_matrix.append(Ex/E1)
			Ey_by_E2_matrix.append(Ey/E2)
			Gxy_by_G12_matrix.append(Gxy/G12)

			Ex_matrix.append(Ex)
			Ey_matrix.append(Ey)
			Gxy_matrix.append(Gxy)
			Qxs_matrix.append(Qxs)

		# Now that we have the required matrices, lets plot.
		labels1 = ['Ex/E1', 'Ey/E2', 'Gxy/G12']
		labels2 = ['Ex', 'Ey', 'Gxy', 'Qxs']

		plt.figure(1)
		plt.plot(theta_degrees, Ex_matrix, label=labels2[0])
		plt.plot(theta_degrees, Ey_matrix, label=labels2[1])
		plt.plot(theta_degrees, Gxy_matrix, label=labels2[2])
		plt.plot(theta_degrees, Qxs_matrix, label=labels2[3])
		plt.xlabel('Theta (in degrees)')
		plt.ylabel('y-axis')
		plt.title('Variations of Engineering Constants with Angle')		
		plt.legend()
		plt.savefig('Plot1.png')
		plt.show()
		
		plt.figure(2)
		plt.plot(theta_degrees, Ex_by_E1_matrix, label=labels1[0])
		plt.plot(theta_degrees, Ey_by_E2_matrix, label=labels1[1])
		plt.plot(theta_degrees, Gxy_by_G12_matrix, label=labels1[2])
		plt.xlabel('Theta (in degrees)')
		plt.ylabel('y-axis')
		plt.title('Variations of Normalized Engineering Constants with Angle')
		plt.legend()
		plt.savefig('Plot2.png')
		plt.show()

		# These two plots have been saved by the names <Plot1.png> and <Plot2.png> respectively.
		# The location of these two plots is the same directory from where this Python code is being run.

	rotateEnggConstants(E1, E2, G12, nu12, Q)

	# ---------- 3rd Assignment Completed ----------

	# Now, let's write a function for calculating A, B and D matrix.
	def calculateABD(numberOfLayers, totalThickness, orientationSequence):
		# Let's get the z-coordinates for each layer and put it inside an array named z.

		global z 
		z = []
		
		N = numberOfLayers
		t = totalThickness
		# thickness of each layer
		tk = t/N
		
		
		z0 = -N/2 * tk
		z.append(z0)


		for i in range(0, N):
			presentCoordinate = z[-1] + tk
			z.append(presentCoordinate)

		# Calculation of Aij, Bij and Dij

		# Initializing [A], [B] and [D] matrix
		Aij = np.array([[0] * 3] * 3)
		Bij = np.array([[0] * 3] * 3)
		Dij = np.array([[0] * 3] * 3)
		
		# Creating Q_layer for storing stiffness matrices (Q) for each layer
		global Q_layer
		Q_layer = []

		for k in range(1, N+1):
			Qij = rotateStiffness(Q, orientationSequence[k-1])
			Aij = Aij + np.array(Qij) * (z[k] - z[k-1])

			Bij = Bij + 0.5 * np.array(Qij) * (pow((z[k]), 2) - pow(z[k-1], 2))

			Dij = Dij + (1/3) * np.array(Qij) * (pow((z[k]), 3) - pow(z[k-1], 3))

			Q_layer.append(Qij)


		print("A matrix is (unit: Pa m):", Aij)
		print("B matrix is (unit: Pa m^2):", Bij)
		print("D matrix is (unit: Pa m^3):", Dij)
		
		return Aij, Bij, Dij

	# To calculate A, B and D matrix.
	

	A, B, D = calculateABD(numberOfLayers, totalThickness, orientationSequence)

	# ---------- 4th Assignment Completed ----------

	# Applying Hashin's failure criteria

	# First, calculating [a], [b], [c] and [d] matrices from [A], [B] and [D] matrices.
	a = np.linalg.inv(A) + np.linalg.inv(A) * B * np.linalg.inv((D - B*np.linalg.inv(A)*np.linalg.inv(B)))*B*np.linalg.inv(A)

	b = -np.linalg.inv(A)*B*np.linalg.inv(D - B*np.linalg.inv(A)*np.linalg.inv(B))

	c = np.linalg.inv(-(D - B*np.linalg.inv(A)*np.linalg.inv(B)))*B*np.linalg.inv(A)

	d = np.linalg.inv(D - B*np.linalg.inv(A)*np.linalg.inv(B))

	# print(a, b, c, d)	

	# Load: N = [[Nx], [Ny], [Nxy]] = [[100 N/mm], [0], [0]]
	# Converting into SI unit
	N = [[100*10**3], [0], [0]]

	epsilon0 = np.dot(a, N)
	kappa = np.dot(c, N)
	# epsilon matrix will contain strain in each layer, i.e., [[epsilon_x], [epsilon_y], [epsilon_xy]]
	epsilon_layer = []

	sigma_12_layer = []

	# calculating strain for each layer
	for layer in range(0, numberOfLayers):
		zavg = (z[layer] + z[layer+1])/2
		epsilon = np.array(epsilon0) + zavg*np.array(kappa)
		epsilon_layer.append(epsilon)
		sigmak =np.dot(Q_layer[layer],epsilon)

		angle = orientationSequence[layer]
		m = math.cos(angle)
		n = math.sin(angle)

		# transformion matrix for sigma
		# This will help to convert sigma back to local (1-2) coordinate system from global X-Y coordinate system
		transformation_sigma = [[m**2, n**2, 2*m*n], [n**2, m**2, -2*m*n], [-m*n, m*n, m**2 - n**2]]
		# print("transformation:",transformation_sigma)
		# print("sigmak:",sigmak)
		# Transforming sigma for applying Hashin's failure criterion
		sigma12 = np.dot(transformation_sigma, sigmak)
		sigma_12_layer.append(sigma12)

	# print(sigma_12_layer)

	# In the variable sigma_12_layer, we have got the values of stress [[sigma_11], [sigma_22], [sigma_12]] in each layer.
	# Now, we can apply Hashin's failure criteria.

	# We are doing our calculations for E-glass/Epoxy-1 composite.
	# The different strength values for this composite have been taken from internet.
	# Longitudinal tensile strength (in MPa)
	X = 783
	# Transverse tensile  strength (in MPa)
	Y = 64
	# Longitudinal compressive strength (in MPa)
	Xdash = 298
	# Transverse compressive strength (in MPa)
	Ydash = 124
	# Shear strengths S12 and S23 (both in MPa)
	S12 = 38
	S23 = 69

	# Let's convert these strength values into SI unit of Pa.
	X = X * 10**6
	Y = Y * 10**6
	Xdash = Xdash * 10**6
	Ydash = Ydash * 10**6
	S12 = S12 * 10**6
	S23 = S23 * 10**6

	# Now, let's apply Hashin's failure criterion during plane stress for each layer in the composite one by one.
	for layer in range(0, numberOfLayers):
		sigmak = sigma_12_layer[layer]
		sigma1 = sigmak[0]
		sigma2 = sigmak[1]
		sigma6 = sigmak[2]

		# Applying condition for matrix tensile failue mode
		# sigma33 = 0
		if sigma2 >= 0:
			
			value = (sigma2/Y)**2 + (sigma6/S12)**2

			if value >= 1:
				print("Layer number " + str(layer+1) + " fails due to matrix tensile failure.")
			else:
				print("Matrix failure does not occur for layer number " + str(layer+1) + ".")
				print("The margin of safety for this layer (with respect to matrix) is:", math.sqrt(1/value))

		# Condition for matrix compressive failure mode
		elif sigma2 < 0:
			
			value = (((Ydash/(2*S23))**2 - 1) * (sigma2/Ydash)) + (sigma2/(4*S23))**2 + (sigma6/S12)**2

			if value >= 1:
				print("Layer number " + str(layer+1) + " fails due to matrix compressive failure.")
			else:
				print("Matrix failure does not occur for layer number " + str(layer+1) + ".")
				# Here, there will be quadratic equation in R (margin of safety).
				# It's in the form a_m * R^2 + b_m * R - 1 = 0 such that:-
				a_m = (sigma2/(4*S23))**2 + (sigma6/S12)**2
				b_m = ((Ydash/(2*S23))**2 - 1)*sigma2/Ydash
				# print(a_m,b_m)

				# Solutions for the quadratic equation
				R1 = float((-b_m + math.sqrt(b_m**2 + 4*a_m*1))/(2*a_m))
				R2 = float((-b_m - math.sqrt(b_m**2 + 4*a_m*1))/(2*a_m))

				# We'll take only positive value as the margin of safety.
				print("The margin of safety for this layer (with respect to matrix) is:", R1)

		# Condition for tensile fibre failure mode
		if sigma1 >= 0:

			value = (sigma1/X)**2 + (sigma6/S12)**2

			if value >= 1:
				print("Layer number " + str(layer+1) + " fails due to tensile fibre failure.")
			else:
				print("Fibre failure does not occur for layer number " + str(layer+1) + ".")
				print("The margin of safety for this layer (with respect to fibre) is:", math.sqrt(1/value))

		# Condition for compressive fibre failure mode
		elif sigma1 < 0:

			value = (sigma1/Xdash)**2

			if value >= 1:
				print("Layer number " + str(layer+1) + " fails due to compressive fibre failure.")
			else:
				print("Fibre failure does not occur for layer number " + str(layer+1) + ".")
				print("The margin of safety for this layer (with respect to fibre) is:", math.sqrt(1/value))

# This sentence is for the condition if there is no data file in the directory for reading values.
# This ensures that the Python file exits the code in such case without giving errors.
else:
	print("The file is not there! Please ensure that the file with the name <data.txt> is present in the same directory as that of this code.")


