import csv
import numpy as np

def is_number(s):
	try:
		float(s)
		return True
	except ValueError:
		return False

def parseDLabData(fileName, viconDataLength):
	eyeDataFile = csv.reader(open(fileName, 'rb'), delimiter='\t')

	eyeDataLeft, eyeDataRight = [], []

	lastLeft, lastRight = [0,0], [0,0]

	for row in eyeDataFile:
		try:	
			
			if float(row[4]) == 0 or float(row[5]) == 0:
				temp = lastLeft
			else:
				temp = [-changeRange(0, 384, -1, 1, float(row[4])), changeRange(0, 288, -1, 1, float(row[5]))]
				lastLeft = temp
			eyeDataLeft.append(temp) # left
	
		
		except:
			try:
				if float(row[6]) == 0 or float(row[7]) == 0:
					temp = lastRight
				else:
					temp = [-changeRange(0, 384, -1, 1, float(row[6])), changeRange(0, 288, -1, 1, float(row[7]))]
					lastRight = temp
				eyeDataRight.append(temp) # right
			except:
				pass

	numViconFrames = viconDataLength

	numMissingLeft = int(numViconFrames%float(len(eyeDataLeft)))/int(numViconFrames/len(eyeDataLeft))
	numMissingLeftStart = 0

	ratio = int(round(numViconFrames/float(len(eyeDataLeft) + numMissingLeftStart)))

	numMissingFromEnd = viconDataLength - len(eyeDataLeft)*ratio - numMissingLeftStart

	eyeDataLeftEx = [[0,0]]*int(numMissingLeftStart)
	eyeDataRightEx = [[0,0]]*int(numMissingLeftStart)

	for i in range(len(eyeDataLeft)):
		for j in range(ratio):
			# negative Y because its from the point of view of the camera
			eyeDataLeftEx.append([eyeDataLeft[i][0], -eyeDataLeft[i][1]])
			try:
				eyeDataRightEx.append([eyeDataRight[i][0], -eyeDataRight[i][1]])
			except:
				pass


	eyeDataRightEx = [[0,0]]*(len(eyeDataLeftEx) - len(eyeDataRightEx)) + eyeDataRightEx + [[0,0]]*numMissingFromEnd
	eyeDataLeftEx = eyeDataLeftEx + [[0,0]]*numMissingFromEnd

	return np.array(eyeDataLeftEx), np.array(eyeDataRightEx)


# Maps the value from the old to the new
def changeRange(oldMin, oldMax, newMin, newMax, value):
	oldRange = (oldMax - oldMin)  
	newRange = (newMax - newMin)  
	return (((value - oldMin) * newRange) / oldRange) + newMin	


# Parses the vicon file and returns a dictionary of point labels and data
def parseViconFile(fileName):
	viconFile = csv.reader(open(fileName, 'rb'), delimiter=',')
	i = 0

	points, lastRow, titles = [], [], []

	for row in viconFile:
		if i == 3:
			titles = filter(None, row)
			print('Parse vicon file: ' + str(titles))

		if i > 5:
			try:
				row = [float(j) for j in row[0:(2 + len(titles)*3)]]
				lastRow = row
			except:
				row = lastRow

			if row == []:
				break

			points.append([row[2+i*3:5+i*3]for i in range(0, len(titles))])	

		i += 1

	points = np.array(points)
	points = np.swapaxes(points, 0, 2)
	points = np.swapaxes(points, 0, 1)

	numFrames = points.shape[2]

	if numFrames%2 != 0:	
		points = np.delete(points, (numFrames-1), axis=2)

	numFrames = points.shape[2]

	constructedPoints = {}
	for i in range(len(points)):
		constructedPoints[titles[i]] = points[i].swapaxes(0,1).reshape(numFrames, 3, 1)
	return constructedPoints