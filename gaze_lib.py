import math
import numpy as np

def toSpherical(x, y, z):
	return [math.atan2(math.sqrt(x**2 + y**2), z), math.atan2(y, x)]

def toCartisian(theta, phi, length):
	x = length*math.cos(phi)*math.sin(theta)
	y = length*math.sin(phi)*math.sin(theta)
	z = length*math.cos(theta)
	return [x, y, z]


# Calculates the spherical cooridinate angles between the head direction and the vector between the eye and the target
def calcAngle(eyePos, headDirection, target):
	diffVector = np.array(eyePos).reshape(3) - np.array(target).reshape(3)
	L = np.linalg.norm(diffVector)
	endEyeHead = eyePos - headDirection*L

	headAngles = toSpherical(headDirection[0], headDirection[1], headDirection[2])

	xdiff = target[0] - eyePos[0] 
	ydiff = target[1] - eyePos[1]
	zdiff = target[2] - eyePos[2]
	return [a_i - b_i for a_i, b_i in zip(toSpherical(xdiff, ydiff, zdiff), headAngles)]


# Calculates the eye position in the world space given 3 points (P, Q, R) on the head and an offset vector a from P 
def calcEyePosition(P, Q, R, a):
	# compute vectors PQ and PR
	P = P.reshape(3)
	pq = P-Q.reshape(3)
	pr = P-R.reshape(3)


	# create a basis for the head plane made up of vectors PQ, v2 and v1
	v1 = np.cross(pq, pr)
	v2 = np.cross(pq, v1)


	# Using the new basis, create the transformation matrix A (which transforms points from the head basis to the world basis)
	A = np.array([v1/np.linalg.norm(v1), v2/np.linalg.norm(v2), pq/np.linalg.norm(pq)])

	A = A.transpose()

	a = np.array(a)

	# return the transformed offset vector plus point P to give the location of the eye in world basis
	return A.dot(a) + P



# returns the direction the head is facing based on the 3 points on the eyetracker
def calcHeadDirection(P, Q, R):

	# the reshapes are only to convert the python numpy arrays from something like [[1],[2],[3]] to [1,2,3] : http://docs.scipy.org/doc/numpy/reference/generated/numpy.reshape.html
	P = P.reshape(3)
	pq = P-Q.reshape(3)
	pr = P-R.reshape(3)

	v1 = np.cross(pq, pr)
	v2 = np.cross(pq, v1)

	perpVector = np.cross(v1, v2)
	direction = np.cross(perpVector, v1)
	
	return direction/np.linalg.norm(direction)



def calcMidPoint(a, b, u1, u2):

	# Get the 2 equations we need to solve in order to find the 2 points of the perpendicular vector
	# https://www.youtube.com/watch?v=HC5YikQxwZA
	a = np.array(a)
	b = np.array(b)
	u1 = np.array(u1)
	u2 = np.array(u2)

	pq = np.array([u2, -u1, b - a])


	# solve the equations for t and s
	eq1 = np.dot(pq, u1)
	eq2 = np.dot(pq, u2)

	eq1A = np.dot(eq2[0], eq1)
	eq2A = np.dot(eq1[0], eq2)

	res = eq1A - eq2A

	try:
		t = float(res[2])/res[1]
		s = (eq1[2] - eq1[1]*t)/eq1[0]
	except:
		# the gaze vectors are parallel so just pick a point, we know this if res[1] will be 0
		t = 100
		s = (eq1[2] - eq1[1]*t)/eq1[0]

	# use t and s for find the points
	p = u1*t
	q = u2*s

	# find the vector PQ and find its length	
	#pqTemp = [p[i] - q[i] for i in range(0,3)]
	#lengthPQ = sqrt(pqTemp[0]*pqTemp[0] + pqTemp[1]*pqTemp[1] + pqTemp[2]*pqTemp[2])

	# use p and q to find the midpoint
	# the point seems to be refected around the vector ab, finding the midpoint between a and b reflecting the focus point about it gives us the correct focus point.
	focusPoint = (a + b - p - q)/2

	return focusPoint
