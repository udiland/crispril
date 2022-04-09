import TreeConstruction_changed as TreeConstruction
import math

def d_f2(seq1,seq2, dicti = None):
	'''a distance function. from the article CRISPER-Cas9 knockout screening in human cells'''
	distance = 0
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	distance += len(seq2) - len(seq1)
	len1 = len(seq1) #this is also the maximum distance between two mismatches
	last_mismatch =len1
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			if last_mismatch == 0:
				Dmm = len1-1 # Dmm is the distance from the last mismatch
			else:
				Dmm = i - last_mismatch
			distance += (i)* (len1 - Dmm)/len1
			Dmm = i
	return distance/len(seq2)

def p_distance(seq1,seq2, dicti = None):
	'''the simplest distance function. will be use for tests'''
	distance = 0
	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	distance += len(seq2) - len(seq1)
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			distance += 1
	return distance/len(seq2)

def MITScore(seq1, seq2, dicti = None):
	'''frm CRISPR-MIT
	PAM cames at the end of the string'''
	distance, first_mm, last_mm = 0, -1, -1

	first_arg = 1
	M = [0, 0, 0.14, 0, 0, 0.395, 0.317, 0, 0.398, 0.079, 0.445, 0.508, 0.613, 0.851, 0.723, 0.828, 0.615, 0.804, 0.685, 0.583]

	if len(seq2) < len(seq1): #putting the longer sequence as seq2
		temp = seq1
		seq1 = seq2
		seq2 = temp
	distance += len(seq2) - len(seq1)
	for i in range(len(seq1)):
		if seq1[i] != seq2[i]:
			distance += 1
			last_mm = i
			first_arg = first_arg * (1-M[i])
			if first_mm == -1:
				first_mm = i
	if distance == 0:
		original_score = 1 ##is it right??
	else:
		d_avg = (last_mm - first_mm)/distance
		original_score = first_arg /((4*(19 - d_avg)/19 +1)*distance**2)

	return 1 - original_score

def MITScore_alternative_imp(seq1, seq2,dicti = None):
	"""
	shiran's implementation
	:param seq1, seq2: sgRNA and off target with PAM seq
	:return:
	"""
	#walk-through of the calculation: https://groups.google.com/forum/#!searchin/crispr/score/crispr/fkhX7FU3r-I/9Nc14v_j3XgJ

	mm_positions = []
	n = len(seq1) - 3
	M = [0, 0, 0, 0, 0, 0.14, 0, 0, 0.395, 0.317, 0, 0.398, 0.079, 0.445, 0.508, 0.613, 0.851, 0.723, 0.828, 0.615, 0.804, 0.685, 0.583]
	M = M[-len(seq1):] # was cheaged from M = M[-len(seq1)+3:]
	score = 1
	#first term
	for i in range(max(n-len(seq1)-1, 0), n):
		if seq1[i] != seq2[i]:
			mm_positions.append(i)
			score *= (1 - M[i])
	n_mm = len(mm_positions)

	#second term
	if n_mm <= 1:
		d_avg = n - 1
	else:
		d_diff = [mm_positions[j] - mm_positions[j-1] for j in range(1, n_mm)]
		d_avg = sum(d_diff)/len(d_diff)

	score *= (1 / ((n - 1 - d_avg) * 4 / (n-1) + 1))

	#third term
	score *= (1 / (max(1,n_mm) ** 2))

	return score

def ccTop(sgseq, target_seq, dicti = None):
	assert len(sgseq) == len(target_seq)
	max_score = sum([math.pow(1.2, i+1) for i in range(len(sgseq))])
    #max score = sum(i = 1,lengh_sg){(1.2**(i+1))}
	mm = [i+1 if sgseq[i] != target_seq[i] else 0 for i in range(len(sgseq))]
    # if mismatch-> mm[i]=i+i else->mm[i]=0
	curScore = sum(list(map(lambda x: pow(1.2, x) if x !=0 else 0, mm)))
	return curScore/max_score # a value between 0 and 1


def make_UPGMA(dm):
	'''use by the doc in http://biopython.org/DIST/docs/api/Bio.Phylo.TreeConstruction.DistanceTreeConstructor-class.html'''
	constructor = TreeConstruction.DistanceTreeConstructor()
	tree = constructor.upgma(dm)
	return tree

def make_initiale_matrix(df,seqList):
	'''input: df: distance_function. seqList: list of sequences, in the order coresponding to the names of sequences
		output: the distance according to the given distance function, arranged as list of lists: a lower triangle matrix
	'''
	res = [] # an array of arrays: a matrix
	for i in range(len(seqList)):
		j = 0
		row = []
		while j <= i:
			tempDistance = df(seqList[i],seqList[j])
			#print(i,j,tempDistance, seqList[i],seqList[j])
			row += [tempDistance]
			j += 1
		res += [row]
	return res

def make_distance_matrix(names, initiale_matrix):
	'''input: list of names of the sequences, and the output of 'make_initiale_matrix'
	output: a distance matrix, in a format adequate to the UPGMA function'''
	m = TreeConstruction._DistanceMatrix(names, initiale_matrix)
	return m


def test1():
	a = "aret"
	b = "ardw"
	c = "brdw"
	seq_list = [a,b,c]
	names = ["a", "b", "c"]
	matrix = make_initiale_matrix(d_f2,seq_list)
	#print(matrix)
	m2 = make_distance_matrix(names, matrix)
	#print(m2)
	m3 = m2.__repr__()
	#print(m3)
	upgma1 = make_UPGMA(m2)
	return upgma1

def test_s_score():
	a = "aaaaaaaaaaaaaaaaaaaa"
	b = "agaaaaaaaaaaaaaaaaaa"
	c = "aaaaaaaaaaaagggggggg"
	d = "GCCTCCCCAAAGCCTGGCCA"
	e = "ACCTCCCCATAGCCTGGCCA"
	print("real data test:", MITScore(d,e))
	print("the same:", MITScore(c,c))
	print("closer:", MITScore(a,b))
	print("farther:", MITScore(a,c))
	print("p distance")
	print("real data test:", p_distance(d,e))
	print("closer:", p_distance(a,b))
	print("farther:", p_distance(a,c))

def test_p_dist():
	s1 = "GGAATGAAAACTATGAGCAC"
	s2 = "GGTCTTAAATCTAGTATCTT"
	print(cc(s1,s2))

def test_ccTop():
	s1 = "GGAATGAAAACTATGAGCAC"
	s2 = "GGAATGAAAACTATGAGCAC"
	print(ccTop(s1,s2))

#test_ccTop()
#test_s_score()