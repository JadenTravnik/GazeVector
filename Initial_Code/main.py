import numpy as np
import sys
from regressor import *
from animator import *
from gaze_lib import *
from parser import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm

if len(sys.argv) < 3:
	print('Give a Vicon filename followed by a Dlab filename')
	exit(0)

ViconDataFileName = sys.argv[1]
DlabDataFileName = sys.argv[2]

# The measured eye offset in mm
eyeOffsetLeft = [50, 60, -30]
eyeOffsetRight = [50, 60, -80]

# Get a dictionary of points from the vicon file
points = parseViconFile(ViconDataFileName)

# Put the arrays into their variables.
rForehead, rTemple, lForehead, lTemple = points['Head:RForeHead'], points['Head:RTemple'], points['Head:LForeHead'], points['Head:LTemple']

# There may be more of these points, so lets just add all of the calibration points we can
CalibrationPoints = []
calibrationNum = 0
while True:
	try:
		calibrationNum += 1 # This only works if there is a a rigid body in Vicon called the CalibrationPlane, this will need to be renamed if that changes
		CalibrationPoints.append(points['CalibrationPane:DLabCal' + str(calibrationNum)])
	except:
		calibrationNum -= 1
		break

print('Parsed Vicon Data')

#Use the vicon data size to parse the Dlab data, this also makes sure that all the data is the same size
DlabEyeDataLeft, DlabEyeDataRight = parseDLabData(DlabDataFileName, len(rForehead))

headAngles, focusPoints, eyePosLeft, eyePosRight = [], [], [], []
anglesLeft = [[] for i in range(calibrationNum)]
anglesRight = [[] for i in range(calibrationNum)]

numFrames = len(rForehead)

print('Parsed Dlab data')
print('Calculating Eye Positions and GazeVectors')
for i in range(numFrames):
	# compute the direction of the head, we don't need this at any other point in the file so only a temporary value is used
	headDirectionTemp = calcHeadDirection(lForehead[i], rForehead[i], rTemple[i])

	# Use the head direction and convert it to head angles, we will save this for later use 
	headAnglesTemp = toSpherical(headDirectionTemp[0], headDirectionTemp[1], headDirectionTemp[2])

	# Calculate the eye position by using the plane defined by the head points and the eye offset that was measured
	eyePosLeftTemp = calcEyePosition(lForehead[i], rForehead[i], rTemple[i], eyeOffsetLeft)
	eyePosRightTemp = calcEyePosition(lForehead[i], rForehead[i], rTemple[i], eyeOffsetRight)

	# calcualte the sphereical angles that would have given a perfect gaze vector
	for j in range(calibrationNum):
		anglesLeft[j].append(calcAngle(eyePosLeftTemp, headDirectionTemp, CalibrationPoints[j][i]))
		anglesRight[j].append(calcAngle(eyePosRightTemp, headDirectionTemp, CalibrationPoints[j][i]))

	headAngles.append(headAnglesTemp)
	eyePosLeft.append(eyePosLeftTemp)
	eyePosRight.append(eyePosRightTemp)


print('Parsed')
print('Finding Clusters')

# find a mean for the angles
anglesLeft = [np.mean(np.array(angles).reshape(numFrames,2), axis=0).tolist() for angles in anglesLeft]
anglesRight = [np.mean(np.array(angles).reshape(numFrames,2), axis=0).tolist() for angles in anglesRight]

# use the average to calculate the first and last clusters
anglesLeft.append(np.mean(anglesLeft, axis=0))
anglesLeft = [anglesLeft[-1]] + anglesLeft

anglesRight.append(np.mean(anglesRight, axis=0))
anglesRight = [anglesRight[-1]] + anglesRight

# If you get an error here about: ValueError: total size of new array must be unchanged,
# look at line 32
anglesLeft = np.array(anglesLeft).reshape(6, 2)
anglesRight = np.array(anglesRight).reshape(6, 2)

# calcuate the focus points from the Dlab data
leftClusters = Cluster2D(DlabEyeDataLeft, numClusters=6)
rightClusters = Cluster2D(DlabEyeDataRight, numClusters=6)

print('Found Clusters')
print('Running Regression')

# run a linear regression on the 6 clusters and predict the angles given the dlab data
newLeft = regressTest(leftClusters, anglesLeft, DlabEyeDataLeft, 2)
newRight = regressTest(rightClusters, anglesRight, DlabEyeDataRight, 2)


plt.plot(leftClusters[:,0],leftClusters[:,1])
plt.plot(anglesLeft[:,0],anglesLeft[:,1])
plt.plot(DlabEyeDataLeft[:,0],DlabEyeDataLeft[:,1])
plt.plot(newLeft[:,0],newLeft[:,1])
plt.show()


leftVector, rightVector, error = [], [], []

print('Regression Complete')
print('Computing Vectors')

for i in range(numFrames):
	
	headAnglesTemp = headAngles[i]

	# using the calculated angles and the head angles at the time, predict the new vectors and the new focus point
	leftTemp = toCartisian(newLeft[i][0] + headAnglesTemp[0], newLeft[i][1] + headAnglesTemp[1], 1)
	rightTemp = toCartisian(newRight[i][0] + headAnglesTemp[0], newRight[i][1] + headAnglesTemp[1], 1)

	focusTemp = calcMidPoint(eyePosLeft[i], eyePosRight[i], leftTemp, rightTemp)
	focusPoints.append(focusTemp)

	leftTemp = toCartisian(newLeft[i][0] + headAnglesTemp[0], newLeft[i][1] + headAnglesTemp[1], np.linalg.norm(np.array(eyePosLeft[i]) - np.array(focusTemp)))
	leftTemp = [ai + bi for ai, bi in zip(leftTemp, eyePosLeft[i])]
	leftVector.append(leftTemp)

	rightTemp = toCartisian(newRight[i][0] + headAnglesTemp[0], newRight[i][1] + headAnglesTemp[1], np.linalg.norm(np.array(eyePosRight[i]) - np.array(focusTemp)))
	rightTemp = [ai + bi for ai, bi in zip(rightTemp, eyePosRight[i])]
	rightVector.append(rightTemp)

	# calculate the distance from each point, the average of the minimum of all distances is the accuracy 
	errorArray = []
	for calibrationPoint in CalibrationPoints:
		errorArray.append(np.linalg.norm(np.array(focusTemp) - calibrationPoint[i].reshape(3)))

	# cap the error at somewhere over 250 because we dont want extream outliers to offset the average error.
	error.append(min(min(errorArray), 250))

# calculate the average error over all time steps
error = np.array(error)
avg = [np.mean(error)]*numFrames

# plot the error
if True:
	plt.plot(error, label='Error')
	plt.plot(avg, label='Avg')
	plt.ylabel('Distance (mm)')
	plt.xlabel('Time Step')
	plt.title('Min distance between Focus Point and Targets')
	plt.grid(True)
	plt.legend()
	plt.axis([0, 3500, 0, 200])

print('Vector Calculation Complete')
print('Starting Animation')

# do some reshaping for the animation
eyePosLeft = np.array(eyePosLeft).reshape(numFrames, 3, 1)
eyePosRight = np.array(eyePosRight).reshape(numFrames, 3, 1)

animationPoints = [eyePosLeft, eyePosRight, rForehead, rTemple, lForehead, lTemple] + CalibrationPoints + [focusPoints]
aniPlot(animationPoints, [[eyePosLeft, leftVector], [eyePosRight, rightVector]]) 
plt.show()
