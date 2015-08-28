from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, writers
import numpy as np

def combine(i, datas):
	x, y, z = [], [], []
	j = 0
	for data in datas:
		ix, iy, iz = data[i]
		x.append(ix)
		y.append(iy)
		z.append(iz)

	return x, y, z


def aniPlot(points, lines):

	print('Making a figure')

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')

	ax.set_autoscale_on(False)
	ax.set_xlim([-500,0])
	ax.set_ylim([200,700])
	ax.set_zlim([1200,1700])



	numFrames = len(points[0])

	ix, iy, iz = combine(0, points)

	scatterPlot = ax.scatter(ix, iy, iz, color=['b', 'g','k', 'k', 'k', 'k', 'k', 'k', 'r', 'r', 'r', 'r', 'c'], s=[100, 100, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 100])

	linePlots = [ax.plot([], [], [])[0] for i in range(len(lines))]

	def init():
		for line in linePlots:
			line.set_data([],[])
			line.set_3d_properties([])
		return linePlots


	def updateLine(p1, p2, line):
		vec = []
		for i in range(3):
			vec.append([p1[i], p2[i]])
		line.set_data(vec[0:2])
		line.set_3d_properties(vec[2])
		return line

	def update_plot(num):
		num += 10
		try:
			x, y, z = combine(num, points)
			scatterPlot._offsets3d = x, y, z
		
			for i in range(len(lines)):
		
				updateLine(lines[i][0][num], lines[i][1][num], linePlots[i])
		except:
			print('Restarted animation')
			num = 0
		return scatterPlot

	print('Animating')
	ani = FuncAnimation(fig, update_plot, init_func=init, frames=numFrames, interval=1, fargs=())

	#ani.save('animation.mp4', writer = 'mencoder', fps=15)

	plt.show()