# def restricted_algorithm(self, seq1, seq2, length):
# 	if math.fabs(min(length, len(seq2)) - min(length, len(seq1))) > MAXINDELS:
# 		return math.inf, "No Alignment Possible", "No Alignment Possible"
# 	swap = False
# 	if len(seq1) < len(seq2):
# 		swap = True
# 		seq1, seq2 = seq2, seq1
# 	xWidth = min(length, len(seq1)) + 1
# 	yWidth = min(length, len(seq2)) + 1
# 	table = {(0, 0): self.EditMove((-1, -1), 0, 0)}
# 	#[[[0] for _ in range(xWidth)] for _ in range(yWidth)]
# 	for i in range(yWidth):
# 		range = range(xWidth)
# 		for j in range(xWidth):
# 			# if min(i, j) <= 0: subMatCost = 0
# 			# else:
# 			# 	subMatCost = (SUB if seq2[i - 1] != seq1[j - 1] else MATCH) + table[(i - 1, j - 1)].cost
# 			# default = self.EditMove.getDefault((i, j), True, subMatCost)
# 			# if default is not None: table[(i, j)] = default
# 			# else:
# 			# 	insCost = INDEL + table[(i - 1, j)].cost
# 			# 	delCost = INDEL + table[(i, j - 1)].cost
# 			# 	if delCost <= insCost and delCost <= subMatCost:
# 			# 		table[(i, j)] = self.EditMove(j - 1, i, delCost, 2)
# 			# 	elif insCost <= subMatCost:
# 			# 		table[(i, j)] = self.EditMove(j, i - 1, insCost, 1)
# 			# 	else:
# 			# 		table[(i, j)] = self.EditMove(j - 1, i - 1, subMatCost, 0)
# 	if swap: seq1, seq2 = seq2, seq1
# 	str1, str2 = self.alignment_strings(seq1, seq2, self.trace_back(table[(yWidth - 1,MAXINDELS)], table, swap))
# 	return table[(yWidth - 1,MAXINDELS)].cost, str1[:100], str2[:100]
# 	# str1, str2 = self.alignment_strings(seq1, seq2, self.trace_back(table[(yWidth - 1, xWidth - 1)], table))
# return table[(yWidth - 1, xWidth - 1)].cost, str1[:100], str2[:100]
#
# def restricted_algorithm(self, seq1, seq2, length):
# 	if math.fabs(min(length, len(seq2)) - min(length, len(seq1))) > MAXINDELS:
# 		return math.inf, "No Alignment Possible", "No Alignment Possible"
# 	swap = False
# 	if len(seq1) < len(seq2):
# 		print("swap")
# 		swap = True
# 		seq1, seq2 = seq2, seq1
# 	xWidth = min(length, 2 * MAXINDELS)+1
# 	yWidth = min(length, len(seq2))+1
# 	table = [[[0, 0, 99, 0] for _ in range(xWidth)] for _ in range(yWidth)]
# 	for i in range(yWidth):
# 		for j in range(xWidth):
# 			if (MAXINDELS - i) > j or i + j >= yWidth + MAXINDELS: continue
# 			seq1Char = i + (j - MAXINDELS) - 1
# 			if i == 0 and seq1Char == -1: table[i][j] = [-1, -1, 0, 0]
# 			elif i == 0: table[i][j] = [j - 1, 0, (seq1Char+1) * INDEL, 1]
# 			elif seq1Char < 0: table[i][j] = [j+1, i - 1, i * INDEL, 2]
# 			else:
# 				noMatch = seq2[i-1] != seq1[seq1Char]
# 				subMatCost = table[i - 1][j][2]
# 				subMatCost += SUB if noMatch else MATCH
# 				if j + 1 >= xWidth: insCost = math.inf
# 				else: insCost = INDEL + table[i - 1][j+1][2]
# 				if j < 1: delCost = math.inf
# 				else: delCost = INDEL + table[i][j - 1][2]
# 				if delCost <= max(insCost, subMatCost): table[i][j] = [j - 1, i, delCost, 2]
# 				elif insCost <= subMatCost: table[i][j] = [j + 1, i - 1, insCost, 1]
# 				else: table[i][j] = [j, i - 1, subMatCost, 0]
# 	print(print('\n'.join(' '.join(str(x[2]) for x in row) for row in table)))
# 	if swap:
# 		seq1, seq2 = seq2, seq1
# 	str1, str2 = self.alignment_strings(seq1, seq2, self.trace_back(table[yWidth - 1][MAXINDELS], table, False))
# 	#if swap:
# 	# return table[yWidth - 1][MAXINDELS][2], str2, str1
# 	return table[yWidth - 1][MAXINDELS][2], str1[:105], str2[:105]
#
