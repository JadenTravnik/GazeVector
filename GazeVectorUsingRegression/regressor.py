import numpy as np
from sklearn.linear_model import LinearRegression as LinReg



# This function simply performs a simple linear regression http://scikit-learn.org/stable/modules/linear_model.html
def regressTest(original, target, test, numTargetCols):

	numFrames = len(test)

	linReg = LinReg(fit_intercept=True)
	
	predictedThings = linReg.fit(original, target).predict(test)
	predictedThings = predictedThings.reshape(numFrames, numTargetCols)
	
	return predictedThings




# Cluster2d takes in a 2*n array and a number of requested clusters and finds temporal clusters (see the documentation for Cluster3d)
def Cluster2D(a, numClusters=6, maxThreshold=0.05, minThreshold=0.01):
	# split the array into 2*numClusters
	a = np.array_split(a, numClusters*2, axis=0)

	clusterThreshold = len(a)/(numClusters*4)


	clusters = []
	for span in a: # for each section of time, find a group of clusters using the default values in Cluster2DwithTime
		clusters += Cluster2DwithTime(span, maxThreshold=maxThreshold, minThreshold=minThreshold, clusterThreshold=clusterThreshold)

	print('lenClusters: ' + str(len(clusters)))
	# after finding all cluster candidates, while the number of clusters is greater than the number desired
	while len(clusters) > numClusters:
		iSpl = NextToMerge(clusters) # get the index to merge, this is the cluster with the minimum distance to the next one

		# split the clusters into the first few, the 2 clusters to be merged, and the end clusters
		firstClusters, clusterA, clusterB, endClusters = clusters[:iSpl], clusters[iSpl], clusters[iSpl+1], clusters[iSpl+2:]

		# fixing a python thing that gives the entire array if given -1 as the index
		if iSpl == 0:
			clusters = [[(clusterA[i] + clusterB[i])/2 for i in range(2)]] + clusters[iSpl+2:]
		else:
			clusters = firstClusters + [[(clusterA[i] + clusterB[i])/2 for i in range(2)]] + endClusters

	return np.array(clusters)




# This is used to find the minimum distance between any 2 clusters in order to know which ones to combine
def NextToMerge(clusters):
	clusterDists = []
	for i in range(len(clusters)):
			try:
				# find the norm between the 2 clusters
				clusterDists.append(np.linalg.norm(np.array(clusters[i]) - np.array(clusters[i+1])))
			except:
				# rather that checking if there is another cluster beyond the current, it is faster to try and if at the end give a giant number
				clusterDists.append(100)
				break

	# sorting the cluster distances and keeping the index of the cluster used, we can find which indexs should be grouped
	clusterDists = [i[0] for i in sorted(enumerate(clusterDists), key=lambda x:x[1])]

	# return the minimum distance
	return clusterDists[0]




def Cluster2DwithTime(a, maxThreshold, minThreshold, clusterThreshold):
	# initialize the first cluster with the first point
	clusters = [a[0]]
	clusterNums = [1]
	c = 0

	# a function that computer the running average of the cluster
	def avgCluster(c, p, n):
		s = c*n + p
		n += 1
		return s/float(n), n
	
	for i in range(len(a)):
		# for the next point in a
		p = a[i]

		try:
			# if the point is close enough to the current cluster
			if np.linalg.norm(p-clusters[c]) <= maxThreshold:
				# merge it with the cluster and add increase the count of the number of points in the cluster
				clusters[c], clusterNums[c] = avgCluster(clusters[c], p, clusterNums[c])
				clusterNums[c] += 1
				# only if the point is in really close to the next one is it considered as the start of a new cluster, this eliminates outliers

			elif np.linalg.norm(p-a[i+1]) <= minThreshold:
				clusters.append(p)
				clusterNums.append(1)
				c += 1

		except:
			pass

	# eliminate any clusters with less than the clusterThreshold
	goodClusters = []
	for i in range(len(clusters)):
		if clusterNums[i] > clusterThreshold:
			goodClusters.append(i)

	# a simple formatting routine
	doneClusters = []
	for i in goodClusters:
		doneClusters += [clusters[i].tolist()]

	return doneClusters
