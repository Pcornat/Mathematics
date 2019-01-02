#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

"""

import sys

import matplotlib.pyplot as plt
from pylab import *


def fonc_u(x, y):
	return 1.0 + sin(0.5 * pi * x) + x * (x - 4) * cos(0.5 * pi * y)


def fonc_uE(x):
	assert (x == 0 or x == 4)
	return


def mesT(coord: ndarray, tri: ndarray) -> double:
	"""
	Calcul d'air d'un triangle à partir de ses 3 points : 1/2 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1))

	:param coord: tableau de coordonnées des points du triangle.
	:param tri: tableau contenant le triangle dont on calcule l'air : chaque case est un sommet du triangle,
	cette case contient l'indice à mettre dans le tableau coord (cf ci-dessous)
	:return: L'air du triangle passé en paramètre.
	"""
	return 1 / 2 * builtins.abs(
		(coord[tri[1], 0] - coord[tri[0], 0]) * (coord[tri[2], 1] - coord[tri[0], 1]) -
		(coord[tri[2], 0] - coord[tri[0], 0]) * (coord[tri[1], 1] - coord[tri[0], 1])
	)


def fonc_f(x, y):
	return 0.25 * pi * pi * sin(0.5 * pi * x) + (-2.0 + 0.25 * pi * pi * x * (x - 4)) * cos(0.5 * pi * y)


def fct_kappa():
	return 1


def fct_alpha():
	return 10.0 ** 8


def coeffelem_P1_rigid(coord: ndarray, tri: ndarray) -> ndarray:
	"""
	:param coord:
	:param tri:
	:return: Matrice carrée de 3
	"""

	matk = zeros_like(tri)
	if not isinstance(coord, ndarray):
		raise TypeError
	if not isinstance(tri, ndarray):
		raise TypeError

	point1 = array([coord[tri[0], 0], coord[tri[0], 1]])  # point1[0] = x1, point1[1] = y1, déduisez la suite.
	point2 = array([coord[tri[1], 0], coord[tri[1], 1]])
	point3 = array([coord[tri[2], 0], coord[tri[2], 1]])

	coeff = 1 / (4 * mesT(coord, tri))

	matk[0, 0] = coeff * ((point2[0] - point3[0]) ** 2 + (point2[1] - point3[1]) ** 2)

	matk[1, 1] = coeff * ((point3[0] - point1[0]) ** 2 + (point3[1] - point1[1]) ** 2)

	matk[2, 2] = coeff * ((point1[0] - point2[0]) ** 2 + (point1[1] - point2[1]) ** 2)

	matk[0, 1] = matk[1, 0] = coeff * (-(point1[0] - point3[0]) * (point2[0] - point3[0]) -
									   (point1[1] - point3[1]) * (point2[1] - point3[1]))

	matk[0, 2] = matk[2, 0] = coeff * (-(point3[0] - point2[0]) * (point1[0] - point2[0]) -
									   (point3[1] - point2[1]) * (point1[1] - point2[1]))

	matk[1, 2] = matk[2, 1] = coeff * (-(point2[0] - point1[0]) * (point3[0] - point1[0]) -
									   (point2[1] - point1[1]) * (point3[1] - point1[1]))
	# matk est calculé.
	return matk


def coeffelem_P1_source(coord: ndarray, tri: ndarray) -> ndarray:
	"""

	:param coord:
	:param tri:
	:return:
	"""

	# TODO : faire la fonction
	return zeros(shape = (3, 3), dtype = double)


def lit_fichier_msh():
	"""
	Permet de lire un fichier .msh
	:return: Les différents éléments caractéristiques qui sont contenus dans le fichier.
	"""

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


def trace_maillage_ind(nbn, nbe, nba, coord: ndarray, tri: ndarray, ar):
	"""

	:param nbn:
	:param nbe:
	:param nba:
	:param coord:
	:param tri:
	:param ar:
	:return:
	"""

	triplot(coord[:, 0], coord[:, 1], tri)
	for i in range(nbn):
		x = coord[i, 0]
		y = coord[i, 1]
		plt.text(x, y, i, color = 'green')
	for i in range(nba):
		x = ((coord[ar[i, 0], 0] + coord[ar[i, 1], 0]) * 0.5)
		y = ((coord[ar[i, 0], 1] + coord[ar[i, 1], 1]) * 0.5)
		plt.text(x, y, i, color = 'orange')
	for i in range(nbe):
		x = ((coord[tri[i, 0], 0] + coord[tri[i, 1], 0] + coord[tri[i, 2], 0]) / 3)
		y = ((coord[tri[i, 0], 1] + coord[tri[i, 1], 1] + coord[tri[i, 2], 1]) / 3)
		plt.text(x, y, i, color = 'red')
	show()


def trace_maillage_ref(nbn, nbe, nba, coord, tri, ar, refn, reft, refa):
	"""

	:param nbn:
	:param nbe:
	:param nba:
	:param coord:
	:param tri:
	:param ar:
	:param refn:
	:param reft:
	:param refa:
	:return:
	"""

	triplot(coord[:, 0], coord[:, 1], tri)
	for i in range(nbn):
		x = coord[i, 0]
		y = coord[i, 1]
		plt.text(x, y, refn[i], color = 'green')
	for i in range(nba):
		x = ((coord[ar[i, 0], 0] + coord[ar[i, 1], 0]) * 0.5)
		y = ((coord[ar[i, 0], 1] + coord[ar[i, 1], 1]) * 0.5)
		plt.text(x, y, refa[i], color = 'orange')
	for i in range(nbe):
		x = ((coord[tri[i, 0], 0] + coord[tri[i, 1], 0] + coord[tri[i, 2], 0]) / 3)
		y = ((coord[tri[i, 0], 1] + coord[tri[i, 1], 1] + coord[tri[i, 2], 1]) / 3)
		plt.text(x, y, reft[i], color = 'red')
	show()


# def assembly(n, nbe):
# 	a = zeros((n, n))
# 	f = zeros((n, 1))
# 	for l in nbe:
#
# 	# endFor
# 	return


def afficherMatCreuse(nbn, coord, tri):
	"""
	Afficher la matrice creuse.
	:param nbn:
	:param coord:
	:param tri:
	:return:
	"""

	x = coord[:, 0]
	y = coord[:, 1]
	psi = zeros(nbn)
	psi[23] = 1.0
	fig = figure(1)
	ax = fig.add_subplot(111, projection = '3d')
	ax.plot_trisurf(x, y, psi, triangles = tri, cmap = cm.jet, linewidth = 0.2)
	xlabel('axe des abscisses')
	ylabel('axe des ordonnees')
	ax.set_zlabel('axe des cotes')
	title('Fonction chapeau numero 23')
	ax.view_init(10, -60)
	show()


def preTraitement():
	return


def traitement():
	return


def postTraitement():
	return


if __name__ == '__main__':
	(nbn, nbe, nba, coord, tri, ar, refn, reft, refa) = lit_fichier_msh()
	afficherMatCreuse(nbn, coord, tri)
	trace_maillage_ind(nbn, nbe, nba, coord, tri, ar)
	trace_maillage_ref(nbn, nbe, nba, coord, tri, ar, refn, reft, refa)
