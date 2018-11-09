#!/bin/python3
# -*- coding: utf-8 -*-

"""
@author: Florent Denef
"""

from pylab import *
from mpl_toolkits.mplot3d.axes3d import Axes3D as a3d
import sys
import matplotlib.pyplot as plt


def lit_fichier_msh():
	nomf = input("Entrer le nom du fichier : ")
	try:
		f = open(nomf, 'r')
	except IOError:
		print("Erreur: problème avec l'ouverture du fichier.")
		sys.exit(0)

	line = f.readline()
	data = line.split()

	nbn = int(data[0])
	nbe = int(data[1])
	nba = int(data[2])

	coord = zeros((nbn, 2), double)
	tri = zeros((nbe, 3), int)
	ar = zeros((nba, 2), int)

	refn = zeros(nbn, int)
	reft = zeros(nbe, int)
	refa = zeros(nba, int)

	for i in range(nbn):
		line = f.readline()
		data = line.split()
		coord[i, 0] = double(data[0])
		coord[i, 1] = double(data[1])
		refn[i] = int(data[2])
	# endFor

	for i in range(nbe):
		line = f.readline()
		data = line.split()
		# Attention aux indices !
		tri[i, 0] = int(data[0]) - 1
		tri[i, 1] = int(data[1]) - 1
		tri[i, 2] = int(data[2]) - 1
		reft[i] = int(data[3])
	# endFor

	for i in range(nba):
		line = f.readline()
		data = line.split()
		# Attention aux indices !
		ar[i, 0] = int(data[0]) - 1
		ar[i, 1] = int(data[1]) - 1
		refa[i] = int(data[2])
	# endFor

	f.close()
	return [nbn, nbe, nba, coord, tri, ar, refn, reft, refa]


def trace_maillage_ind(nbn, nbe, nba, coord, tri, ar):
	triplot(coord[:, 0], coord[:, 1], tri)
	for i in range(nbn):
		x = coord[i, 0]
		y = coord[i, 1]
		plt.text(x, y, i, color = 'green')
	for i in range(nba):
		x = ((coord[ar[i, 0], 0] + coord[ar[i, 1], 0]) / 2)
		y = ((coord[ar[i, 0], 1] + coord[ar[i, 1], 1]) / 2)
		plt.text(x, y, i, color = 'orange')
	for i in range(nbe):
		x = ((coord[tri[i, 0], 0] + coord[tri[i, 1], 0] + coord[tri[i, 2], 0]) / 3)
		y = ((coord[tri[i, 0], 1] + coord[tri[i, 1], 1] + coord[tri[i, 2], 1]) / 3)
		plt.text(x, y, i, color = 'red')
	show()


def trace_maillage_ref(nbn, nbe, nba, coord, tri, ar, refn, reft, refa):
	triplot(coord[:, 0], coord[:, 1], tri)
	for i in range(nbn):
		x = coord[i, 0]
		y = coord[i, 1]
		plt.text(x, y, refn[i], color = 'green')
	for i in range(nba):
		x = ((coord[ar[i, 0], 0] + coord[ar[i, 1], 0]) / 2)
		y = ((coord[ar[i, 0], 1] + coord[ar[i, 1], 1]) / 2)
		plt.text(x, y, refa[i], color = 'orange')
	for i in range(nbe):
		x = ((coord[tri[i, 0], 0] + coord[tri[i, 1], 0] + coord[tri[i, 2], 0]) / 3)
		y = ((coord[tri[i, 0], 1] + coord[tri[i, 1], 1] + coord[tri[i, 2], 1]) / 3)
		plt.text(x, y, reft[i], color = 'red')
	show()


# def assembly(n):
# 	a = zeros((n, n))
# 	f = zeros((n, 1))
#
# 	return


if __name__ == '__main__':
	[nbn, nbe, nba, coord, tri, ar, refn, reft, refa] = lit_fichier_msh()
	x = coord[:, 0]
	y = coord[:, 1]
	psi = zeros(nbn)
	psi[23] = 1.0
	fig = figure(1)
	ax = fig.add_subplot(111, projection = '3d')
	ax.plot_trisurf(x, y, psi, triangles = tri, cmap = cm.jet, livewidth = 0.2)
	xlabel('axe des abscisses')
	ylabel('axe des ordonnees')
	ax.set_zlabel('axe des cotes')
	title('Fonction chapeau numero 23')
	ax.view_init(10, -60)
	show()
