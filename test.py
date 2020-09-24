from Bio.Seq import Seq
from Bio import motifs
import math
import numpy as np
import sys
import json

my_seq = Seq(sys.argv[3])
motifList = []

motif35 = sys.argv[1]  # user input
a_list35 = motif35.split(',')
map_object35 = map(int, a_list35)
list_of_integers35 = list(map_object35)

motif10 = sys.argv[2]  # user input
a_list10 = motif10.split(',')
map_object10 = map(int, a_list10)
list_of_integers10 = list(map_object10)

gap = [17, 18, 19, 20, 21, 22, 23, 24, 25]
for pos10 in list_of_integers10:
    possibilities = []
    for pos35 in list_of_integers35:
        for val in gap:
            if pos35 == pos10 + val:
                possibilities.append((pos35, pos10))
    motifList.append(possibilities)

# Seperating  the possible -35 & -10 pair to desire format: output - list 35,list 10
output = []


def reemovNestings(l):
    for i in l:
        if type(i) == list:
            reemovNestings(i)
        else:
            output.append(i)


reemovNestings(motifList)
list35 = []
list10 = []
for a, b in output:
    list35.append(a)
    list10.append(b)

# Slicing the Hexamer sequence based on the positions(-35,-10)
seq35 = []
seq10 = []
str35 = []
str10 = []
for i in list35:
    end = i + 6
    seq35.append(my_seq[i:end])
for i in list10:
    end = i + 6
    seq10.append(my_seq[i:end])
for i in seq35:
    str35.append(str(i))
for i in seq10:
    str10.append(str(i))

filename1 = 'fimo_35_ff_ed11.csv'
infile1 = open(filename1, 'r')
instances_regdb = list()
for line in infile1.readlines():
    if line == '\n':
        break
    lines = line.split(',')
    lines[1] = lines[1][:-1]
    instances_regdb.append(Seq(lines[1].upper()))
filename2 = 'fimo_10_ff_ed11.csv'
infile2 = open(filename2, 'r')
instances2_regdb = list()
for line in infile2.readlines():
    if line == '\n':
        break
    lines = line.split(',')
    lines[1] = lines[1][:-1]
    instances2_regdb.append(Seq(lines[1].upper()))
instances = [
    Seq("TTGACG"),
    Seq("TTTACA"),
    Seq("TTGACA"),
    Seq("TTGACA"),
    Seq("TTTACG"),
    Seq("TTTACG"),
    Seq("TTTACG"),
    Seq("CTGACA"),
    Seq("TTTACA"),
    Seq("TTTACG"),
    Seq("TTGACG"),
    Seq("CTGATA"),
    Seq("CTGATG"),
    Seq("TTTATG"),
    Seq("TTTATA"),
    Seq("TTGACA"),
    Seq("TTGACA"),
    Seq("TTGACG")
]
convert35 = []
result = []
i = 0
while i < len(instances):
    convert35.append(str(instances[i]))
    i += 1
m = motifs.create(instances_regdb[:])
pwm = m.counts.normalize(
    pseudocounts={'A': 0.49, 'C': 0.51, 'G': 0.51, 'T': 0.49})
pssm = pwm.log_odds()


def calculateX(a, b, c, d, e, f, x):
    temp1 = pssm[a, 0] + pssm[b, 1] + pssm[c, 2] + \
        pssm[d, 3] + pssm[e, 4] + pssm[f, 5]
    result.append([temp1])


i = 0
while i < len(convert35):
    calculateX(convert35[i][0], convert35[i][1], convert35[i]
               [2], convert35[i][3], convert35[i][4], convert35[i][5], i)
    i += 1
instances2 = [
    Seq("TACAGT"),
    Seq("TATTAT"),
    Seq("TACTGT"),
    Seq("TATTGT"),
    Seq("TACTAT"),
    Seq("TATAGT"),
    Seq("TATTAT"),
    Seq("TATAAT"),
    Seq("GACTGT"),
    Seq("TACAAT"),
    Seq("TATAGT"),
    Seq("GATTAT"),
    Seq("GATTAT"),
    Seq("TACAAT"),
    Seq("TACAAT"),
    Seq("GACTAT"),
    Seq("GATTGT"),
    Seq("TATTGT")
]
result2 = []
convert10 = []
i = 0
while i < len(instances2):
    convert10.append(str(instances2[i]))
    i += 1
m2 = motifs.create(instances2_regdb)
pwm2 = m2.counts.normalize(
    pseudocounts={'A': 0.49, 'C': 0.51, 'G': 0.51, 'T': 0.49})
pssm2 = pwm2.log_odds()


def calculateX2(a, b, c, d, e, f, x):
    temp1 = pssm2[a, 0] + pssm2[b, 1] + pssm2[c, 2] + \
        pssm2[d, 3] + pssm2[e, 4] + pssm2[f, 5]
    result2.append([temp1])


i = 0
while i < len(convert10):
    calculateX2(convert10[i][0], convert10[i][1], convert10[i]
                [2], convert10[i][3], convert10[i][4], convert10[i][5], i)
    i += 1
outputResult = [
    [1],
    [0.7],
    [0.86],
    [0.72],
    [0.24],
    [0.47],
    [0.36],
    [0.51],
    [0.04],
    [0.33],
    [0.58],
    [0.01],
    [0.01],
    [0.1],
    [0.15],
    [0.16],
    [0.06],
    [0.56]
]
a = []
i = 0
while i < len(outputResult):
    a.append([1, result[i][0], result2[i][0]])
    i += 1
b = []
i = 0
while i < len(outputResult):
    b += outputResult[i]
    if b[i] < 0.0001:
        b[i] = 0.0001
    b[i] = math.log(b[i])
    i += 1
x = np.asarray(a)
y = np.asarray(b)


def gradientDescent(x, y, theta, alpha, m, numIterations):
    J_history = np.zeros(shape=(numIterations, 1))
    xTrans = x.transpose()
    for i in range(0, numIterations):
        hypothesis = np.dot(x, theta)
        loss = hypothesis - y
        cost = np.sum(loss ** 2) / (2 * m)
        gradient = np.dot(xTrans, loss) / m
        theta = theta - alpha * gradient
        J_history[i][0] = cost
    return theta, J_history


m, n = np.shape(x)
numIterations = 100000  # c
alpha = 0.015  # c
theta = np.ones(n)
theta, J_history = gradientDescent(x, y, theta, alpha, m, numIterations)

# Calculating the strength for the possible pairs
fiResult35 = []
fiResult10 = []
fiStrength = []
filnStrength = []
for i in seq35:
    aa, bb, cc, dd, ee, ff = str(i)
    fi1 = pssm[aa, 0] + pssm[bb, 1] + pssm[cc, 2] + \
        pssm[dd, 3] + pssm[ee, 4] + pssm[ff, 5]
    fiResult35.append([fi1])
for i in seq10:
    aa1, bb2, cc3, dd4, ee5, ff6 = str(i)
    fi2 = pssm[aa1, 0] + pssm[bb2, 1] + pssm[cc3, 2] + \
        pssm[dd4, 3] + pssm[ee5, 4] + pssm[ff6, 5]
    fiResult10.append([fi2])
for f, b in zip(fiResult35, fiResult10):
    i = 0
    j = 0
    strength = np.array([1.0, f[i], b[j]]).dot(theta)
    fiStrength.append(strength)
    lnstrength = math.exp(strength)
    filnStrength.append(lnstrength)
    i += 1
    j += 1
finarra = [{'position35': d, 'hexamer35': b, 'position10': e, 'hexamer10': c,
            'prdStr': a} for (a, b, c, d, e) in zip(filnStrength, str35, str10, list35, list10)]


def myFunczz(e):
    return e['prdStr']


finarra.sort(reverse=True, key=myFunczz)
print(finarra)
