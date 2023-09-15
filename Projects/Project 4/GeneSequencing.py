#!/usr/bin/python3

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import math
import time
import random

from collections import deque

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

''' Constants used to denote edit type, for readability '''
DELETE = 1
INSERT = 2
SUBMATCH = 0

class GeneSequencing:
	def __init__( self ):
		pass

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align( self, seq1, seq2, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = align_length

###################################################################################################
		''' Call the alignment algorithm '''
		align_cost, seq1_first100, seq2_first100 = self.align_sequences(seq1, seq2)
###################################################################################################

		return {'align_cost':align_cost, 'seqi_first100':seq1_first100, 'seqj_first100':seq2_first100}

	class _EditMove:
		"""
		A Helper class used to store each step within the dictionary.
		Made of a cost value, a back-pointer, and the edit type.
		Time Complexity: O(1)
		Space Complexity: O(1)
		"""
		def __init__(self, srcPos, cost, editType):
			"""
			Constructor for an EditMove, a helper class for finding edit distance
			Complexity:
				Time: O(1)
				Space: O(1)
			:param srcPos: location of back-pointing source position in the dictionary
			:param cost: cost of the edit up to now
			:param editType: the type of edit need to reach here, either INSERT, DELETE, or SUBMATCH
			"""
			self.src = srcPos
			self.cost = cost
			self.type = editType
		@classmethod
		def getDefault(cls, pos):
			"""
			A function to create the default EditMoves, or return None.
			This function handles the top and left edges of the matrix/dictionary, for both banded and unrestrained
			Complexity:
				Time: O(n)
				Space: O(1)
			:param pos: The location of the cell in the resulting dictionary
			:return: Either a default value, or None
			"""
			x, y = pos
			if pos == (0, 0): return cls((-1, -1), 0, SUBMATCH)
			if x <= 0: return cls((x, y - 1), y * INDEL, DELETE)
			elif y <= 0: return cls((x - 1, y), x * INDEL, INSERT)
			return None

	def align_sequences(self, seq1, seq2):
		"""
		Function used to get the edit distance between two strings using the Needleman/Wunsch cost function.
		This function runs two variations, one which is for a faster banded version, and one for a unrestrained variant.
		Both cases are done by creating a dictionary of values with a cost, a backtracking pointer, and the edit type.

		These algorithms were originally combined via if-statements;
		But I separated them for easier time-space complexity analysis:
		Unrestrained:
			Time: O(nm)
			Space: O(nm)
		Banded:
			Time: O(kn)
			Space: O(kn)

		:param seq1: first string to run the alligment with
		:param seq2: second string to run the alligment with
		:return: cost of alignment, first 100 characters of seq1 and seq2
		"""
		''' Trim the two sequences to the max characters allowed to be aligned
			Time-Complexity: O(n/m/MAX) further code will assume that n/m are the size of MaxCharactersToAlign '''
		seq1, seq2 = seq1[:self.MaxCharactersToAlign], seq2[:self.MaxCharactersToAlign]
		''' Terminate early if the two strings are identical
			Time-Complexity: O(n/m) '''
		if seq1 == seq2: return len(seq1)*MATCH, seq1[:100], seq2[:100]
		''' Otherwise call the appropriate algorithm '''
		if self.banded: return self._banded_algorithm(seq1, seq2)
		else: return self._unrestrained_algorithm(seq1, seq2)

	def _unrestrained_algorithm(self, seq1, seq2):
		"""
		Function used to get the edit distance between two strings using the Needleman/Wunsch cost function.
		This algorithm creates a dictionary of size nm representing a matrix.
		Each cell/value stores the cost, and the previous location.
		Each cell is found by comparing the cost of each possible edit and storing it.
		A backtracking algorithm is used to generate the edited version of the string afterwards
		Unrestrained:
			Time: O(nm)
			Space: O(nm)
		:param seq1: first string to run the alligment with
		:param seq2: second string to run the alligment with
		:return: cost of alignment, first 100 characters of seq1 and seq2
		"""
		''' Create a dictionary to represent the matrix
			Time: O(1) Space: O(1)->Expands to O(nm) '''
		matrix = {(0, 0): self._EditMove((-1, -1), 0, SUBMATCH)}
		''' Iterate over every cell of the matrix
			Time: O(nm) Space(nm) '''
		for x in range(len(seq1) + 1):
			for y in range(len(seq2) + 1):
				next_pos = (x, y)
				''' Assign all the matrix left and top edges to their default values
					Time: O(1) '''
				default = self._EditMove.getDefault(next_pos)
				if default is not None: matrix[next_pos] = default
				else:
					'''
					# Find the edit cost of a given cell, choosing the lowest value.
					# Bias towards left, top, then diagonal
					Dictionary look ups are average of O(1).
					Addition, subtraction and comparison are all O(1)
					Increases size of matrix by 1, leading to a O(nm) matrix
						Avg. Time: O(1)
						Space: O(1)
					'''
					subMatCost = (SUB if seq1[x - 1] != seq2[y - 1] else MATCH) + matrix[(x - 1, y - 1)].cost
					insCost = INDEL + matrix[(x - 1, y)].cost
					delCost = INDEL + matrix[(x, y - 1)].cost
					if delCost <= min(insCost, subMatCost):
						matrix[next_pos] = self._EditMove((x, y - 1), delCost, DELETE)
					elif insCost <= subMatCost:
						matrix[next_pos] = self._EditMove((x - 1, y), insCost, INSERT)
					else:
						matrix[next_pos] = self._EditMove((x - 1, y - 1), subMatCost, SUBMATCH)
		''' Get the final cell of the matrix '''
		last = matrix[(len(seq1), len(seq2))]
		''' Return the cost and the edited strings
			Time: O(n)	Space: O(n) '''
		return last.cost, *self._alignment_strings(seq1, seq2, self._back_track(last, matrix))

	def _banded_algorithm(self, seq1, seq2):
		"""
		Function used to get the edit distance between two strings using the Needleman/Wunsch cost function.
		This algorithm creates a dictionary of size kn representing a matrix.
		This non-square matrix is made by filling each row with k values, forming a band across the table
		Each cell/value stores the cost, and the previous location.
		Each cell is found by comparing the cost of each possible edit and storing it.
		A backtracking algorithm is used to generate the edited version of the string afterwards
		Banded:
			Time: O(kn)
			Space: O(kn)
		:param seq1: first string to run the alligment with
		:param seq2: second string to run the alligment with
		:return: cost of alignment, first 100 characters of seq1 and seq2
		"""
		''' Terminate early if the difference in sizes between the two sequences are too large.
			This is because the only valid solutions would ly outside the banded matrix '''
		if self.banded and math.fabs(len(seq1) - len(seq2)) > MAXINDELS:
			return math.inf, (alStr := "No Alignment Possible"), alStr
		''' Create a dictionary to represent the matrix
				Time: O(1) Space: O(1)->Expands to O(kn) '''
		matrix = {(0, 0): self._EditMove((-1, -1), 0, SUBMATCH)}
		''' Iterate over every cell of the matrix of the banded matrix
				Time: O(kn) Space(kn) '''
		for x in range(len(seq1)+1):
			''' Each row of the matrix is only filled with k values, offset by the row index
				These cells then form a hexagonal diamond-like shape '''
			check_y = range(max(0, x - MAXINDELS), min(len(seq2), x + MAXINDELS) + 1)
			for y in check_y:
				next_pos = (x, y)
				''' Assign all the matrix left and top edges to their default values
						Time: O(1) '''
				default = self._EditMove.getDefault(next_pos)
				if default is not None: matrix[next_pos] = default
				else:
					'''
					Find the edit cost of a given cell, choosing the lowest value.
					Bias towards left, top, then diagonal
					Distances landing outside the band are set to infinity so they will not be used.
					This causes the algorithm to produce a limited band.
					Dictionary look ups are average of O(1).
					Addition, subtraction and comparison are all O(1)
					Increases size of matrix by 1, leading to a O(nm) matrix
						Avg. Time: O(1)
						Space: O(1)
					'''
					subMatCost = (SUB if seq1[x - 1] != seq2[y - 1] else MATCH) + matrix[(x - 1, y - 1)].cost
					insCost = math.inf if y-x >= MAXINDELS else INDEL + matrix[(x - 1, y)].cost
					delCost = math.inf if x-y >= MAXINDELS else INDEL + matrix[(x, y - 1)].cost
					if delCost <= min(insCost, subMatCost):
						matrix[next_pos] = self._EditMove((x, y - 1), delCost, DELETE)
					elif insCost <= subMatCost:
						matrix[next_pos] = self._EditMove((x - 1, y), insCost, INSERT)
					else:
						matrix[next_pos] = self._EditMove((x - 1, y - 1), subMatCost, SUBMATCH)
		''' Get the final cell of the matrix '''
		last = matrix[(len(seq1), len(seq2))]
		''' Return the cost and the edited strings
				Time: O(n)	Space: O(n) '''
		return last.cost, *self._alignment_strings(seq1, seq2, self._back_track(last, matrix))

	def _back_track(self, startPt, matrix):
		"""
		Algorithm used to backtrack through the table, returning a list of edits made to achieve the final cost
		Complexity:
			Time: O(n/m)
			Space: O(n/m)
		:param startPt: starting point of the backtracking algorithm
		:param matrix: the matrix of EditMoves needed to get the final edit cost
		:return: a list of edit types used to find the final edit cost
		"""
		path = [startPt.type]
		''' while loop used to move through the matrix, visiting at most n/m cells '''
		while (startPt := matrix[startPt.src]).src != (-1, -1):
			''' back insertion Time O(1)'''
			path.append(startPt.type)
		return path

	def _alignment_strings(self, seq1, seq2, path):
		"""
		Function to produce the edited versions of the sequences.
		Complexity:
			Time: O(1) - 100 iterations, reduces O(100) to O(1)
			Space: O(1) - 200 characters, reduces O(200) to O(1)
		:param seq1: first sequence to edit
		:param seq2: second sequence to edit
		:param path: deque of edits to make
		:return: first 100 characters of the edited versions of seq1 and seq2
		"""
		str1 = str2 = ""
		seq1Idx = seq2Idx = 0
		''' Iterate through last 100 values of the path, the first 100 of the sequences 
			Time: O(1) Space: O(1) '''
		for move in path[:-101:-1]:
			''' Appends a character or '-' to each string according to the respective edit made '''
			if move == INSERT:
				str1 += seq1[seq1Idx]
				seq1Idx += 1
				str2 += '-'
			elif move == DELETE:
				str1 += '-'
				str2 += seq2[seq2Idx]
				seq2Idx += 1
			else:
				str1 += seq1[seq1Idx]
				str2 += seq2[seq2Idx]
				seq1Idx += 1
				seq2Idx += 1
		''' resulting 100 character strings '''
		return str1, str2
