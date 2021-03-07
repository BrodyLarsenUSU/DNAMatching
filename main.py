import time
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np

dnaFile = "covid_DNA.fasta"
fileName = "covid_DNA.fasta"

def read_fasta(howMany = 2):
    # The first line in a FASTA record is the title line.
    # Examples:
    # >third sequence record
    # >gi|2765657|emb|Z78532.1|CCZ78532 C.californicum 5.8S rRNA gene
    # returns a list of sequences as tuples (name.)
    with open(fileName,  'r') as filePt:
        sequences = []
        fastas = filePt.read().split(">")
        fastas = fastas[1:]
        for i in range(0,howMany):
            seq = fastas[i].split("\n")
            seq_name = seq[0]
            fasta_seq = "".join(seq[1:])
            sequences.append((seq_name,fasta_seq))
        return sequences

    ##Calculates the time taken by each function call and generates graph
def timeProblems(DNAProblem, function, init=None, fit='exponential'):
        # problemList is a list of tuples [(size, arguments),...] ordered smallest to biggest
        # runs and times the function with each arguments a
        # generates a graph of run time as a function of problem size
        # fit may be 'exponential' then the time as a function of problem size is assumed
        #     to of the form time = c * a^n and the function solves for c and a
        #     where a is the base of the exponential function and c is a multiplicative factor
        # fit my be 'polynomial' then the time as a function of problem size is assumed
        #     to of the form time = c * n ^ b and the function solves for c and b
        #     where b is the power of n (the degree of the polynomial) and c is a multiplicative fac* tor
        timeLine = []
        values = []

        for (size, args) in DNAProblem:

            start_time = time.time()
            function(*args)  # use the * to unpack the tuple into arguments to the function
            elapsed = (time.time() - start_time) * 1000.0
            if elapsed > 0.0:
                timeLine.append(elapsed)
                # print(elapsed)
                values.append(size)

            print(str(size) + " Sequence Size complete")
        ##Generating the plot between time taken by each function call with n as variable and n
        plt.plot(values, timeLine, 'g')
        plt.xlabel("Problem size")
        plt.yscale('log')
        if fit == 'polynomial':
            plt.xscale('log')
        plt.ylabel("time in milliseconds")
        plt.rcParams["figure.figsize"] = [16, 9]
        plt.show()
        if fit == 'exponential':  # fit a straight line to n and log time
            slope, intercept, _, _, _ = stats.linregress([values], [np.log(t) for t in timeLine])
            print("time = %.6f * %.3f ^ n" % (np.exp(intercept), np.exp(slope)))
        elif fit == 'polynomial':  # fit a straight line to log n and log time
            slope, intercept, _, _, _ = stats.linregress([np.log(v) for v in values], [np.log(t) for t in timeLine])
            print("time = %.6f * n ^ %.3f" % (np.exp(intercept), slope))
            timems = (np.exp(intercept)) * (30000 ** slope)
            hours = (timems / 1000) / 60
            print("if n was 30,000 it would take : " + str(timems) + "milliseconds")
            print("or " + str(hours) + " hours to run")
            print("The time it takes to compute grows at a rate of n^2 because of")
            print("the nested for loops in the matchDP. The 1st for loop is getting")
            print("multiplied by itself in the 2nd for loop everytime the 2nd for loop iterates.")
            print("side note, sometimes this program works and sometimes it doesn't, sometimes ")
            print("I can run the program with the sequences being 16,000 characters long, sometimes it")
            print("crashes after 2000 character long sentences. I'm not sure if it is COLABS fault or")
            print("if it is my code. COLAB tells me i ram out of RAM")

seq = read_fasta()
(name0, seq0) = seq[0]
(name1, seq1) = seq[1]


def matchR(A, B):
    # uses recursion to solve the string matching problem
    # mod the strings to put a blank in front so that n and m represent both
    # the number of characters in the string AND the index into the string
    n = len(A)
    m = len(B)
    A = '_' + A
    B = '_' + B
    return match(n, m, A, B, )


def match(n, m, A, B, depth=0):
    # one of the strings is empty
    # print("%s A = %s B= %s" % (" " * depth, A[0:n], B[0:m]))
    if n == 0:
        sol = m
    elif m == 0:
        sol = n
    # three cases
    else:
        # print(matrixLookUp(A[n], B[m]))
        sol = max(match(n - 1, m, A, B, depth + 1) + matrixLookUp("_", A[n]),  # delete from A
                  match(n, m - 1, A, B, depth + 1) + matrixLookUp("_", B[m]),  # delete from B
                  match(n - 1, m - 1, A, B, depth + 1) + matrixLookUp(A[n], B[m]))  # substitute (if 1 != 1  )
    # print("%s solution = %d" % (" " * depth, sol))
    return sol


# prints out the solution
def printAlign(A, B, cache):
    n = len(A)
    m = len(B)
    A = '_' + A
    B = '_' + B
    # alignment = traceback(n, m, A, B)
    alignment = traceBackDP(n, m, A, B, cache)
    alignment.reverse()
    for oneAlign in alignment:
        print(oneAlign)
    # print("____________________")


def traceback(n, m, A, B):
    # n is in A, m is in B
    # base case
    # nothing to do, return empty
    if n == 0 and m == 0:
        print(cache[(n, m)])
        return []
    # no characters in A, so add a delete from A and continue
    if n == 0:
        # delete from A
        print(cache[(n, m)])
        return ["%s - %s" % ('_', B[m])] + traceback(n, m - 1, A, B)
    # no characters in B, so delete from B and continue
    if m == 0:
        # delete from B
        print(cache[(n, m)])
        return ["%s - %s" % (A[n], '_')] + traceback(n - 1, m, A, B)
    # Need to determine which alignment sub solution was used for the optimal solution
    sol = cache[(n, m)]
    # we know that sol must equal one of the sub solutions
    if sol == cache[(n - 1, m)] + matrixLookUp("_", A[n]):  # delete from B
        print(sol)
        return ["%s - %s" % (A[n], '_')] + traceback(n - 1, m, A, B)
    if sol == cache[(n, m - 1)] + matrixLookUp("_", B[m]):  # delete from A
        print(sol)
        return ["%s - %s" % ('_', B[m])] + traceback(n, m - 1, A, B)
    # must have matched the characters, check the characters
    if A[n] != B[m]:  # substitution
        print(sol)
        return ["%s x %s" % (A[n], B[m])] + traceback(n - 1, m - 1, A, B)
    # exact match
    print(sol)
    return ["%s = %s" % (A[n], B[m])] + traceback(n - 1, m - 1, A, B)


def traceBackDP(n, m, A, B, cache):
    line = []
    # print(cache[(n,m)])

    while True:
        solution = cache[(n, m)]
        if n == 0 and m == 0:
            # print("n and m == 0")
            # print(solution)
            return line

        elif n == 0:
            # line = line + ["%s - %s" % ('_', B[m])] #+ cache[(n, m-1)]
            # print("n==0")
            m = m - 1
            # print(solution)


        elif m == 0:
            # line = line + ["%s - %s" % (A[n], '_')] #+ cache[(n-1, m)]
            print("m==0")
            n = n - 1
            # print(solution)


        elif solution == cache[(n - 1, m)] + matrixLookUp("_", A[n]):  # delete from B
            # print("delete B")
            line = line + ["%s - %s" % (A[n], '_')]  # + cache[(n-1, m)]
            solution = cache[(n, m)]
            # print(solution)
            n = n - 1


        elif solution == cache[(n, m - 1)] + matrixLookUp("_", B[n]):  # delete from A
            # print("delete A")
            line = line + ["%s - %s" % ('_', B[m])]  # + cache[(n, m-1)]
            solution = cache[(n, m)]
            # print(solution)
            m = m - 1



        elif A[n] != B[m]:  # substitution
            # print("doesnt equal")
            line = line + ["%s x %s" % (A[n], B[m])]  # + cache[(n-1, m-1)]
            solution = cache[(n, m)]
            # print(solution)
            n = n - 1
            m = m - 1


        elif A[n] == B[m]:
            # print("equals")
            line = line + ["%s = %s" % (A[n], B[m])]  # + cache[(n-1, m-1)]
            solution = cache[(n, m)]
            # print(solution)
            n = n - 1
            m = m - 1


        else:
            print("error")
            n = n - 1
            m = m - 1
    print("traceback complete")
    return line


def matchDP(A, B):
    global cache
    cache = {}
    X = A
    Y = B

    n = len(A)
    m = len(B)
    # modify the strings to put a blank in front
    A = '_' + A
    B = '_' + B
    # fill in the base cases
    for i in range(0, n + 1):
        cache[(0, i)] = i  # matrixLookUp("_",B[i])
    for j in range(0, m + 1):
        cache[(j, 0)] = j  # matrixLookUp("_", A[j])
    # loop though all problems from smallest to biggest
    # what about the 0 0 case?
    cache[(0, 0)] = 0
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            # need to change the indexs into the array!
            # so from n in recursion to j and from m in recursion to i
            cache[(i, j)] = max(cache[(i - 1, j)] + matrixLookUp("_", A[i - 1]),
                                cache[(i, j - 1)] + matrixLookUp("_", B[j - 1]),
                                cache[(i - 1, j - 1)] + matrixLookUp(A[i - 1], B[j - 1]))
    return cache[(n, m)]


# looks up the number associated with th eletters
def matrixLookUp(x, y):
    # if the letter being looked up isnt A,G,C,T
    # change it to T
    if x != "A" or "G" or "C" or "T":
        x = "T"
    if y != "A" or "G" or "C" or "T":
        y = "T"

    if x == "A" and y == "A":
        return 5
    if x == "A" and y == "C":
        return -1
    if x == "A" and y == "G":
        return -2
    if x == "A" and y == "T":
        return -1

    if x == "C" and y == "A":
        return -1
    if x == "C" and y == "C":
        return 5
    if x == "C" and y == "G":
        return -3
    if x == "C" and y == "T":
        return -2

    if x == "G" and y == "A":
        return -2
    if x == "G" and y == "C":
        return -3
    if x == "G" and y == "G":
        return 5
    if x == "G" and y == "T":
        return -2

    if x == "T" and y == "A":
        return -1
    if x == "T" and y == "C":
        return -2
    if x == "T" and y == "G":
        return -2
    if x == "T" and y == "T":
        return 5

    if x == "_" and y == "A":
        return -3
    if x == "_" and y == "C":
        return -4
    if x == "_" and y == "G":
        return -2
    if x == "_" and y == "T":
        return -1

    else:
        print("ERROR: " + x + " " + y + " other that provide letters found")

#TEST correctness of code
A = 'ACGGAGAGAATC'
B = 'TGGAGAGAATTC'


A = 'GAA'
B = 'GAC'


#match to fill in the cache
#edit = matchDP(A,B)
#matchDP(A,B)
#edit = matchR(A, B)
#printAlign(A, B, cache)
#print(edit)

def generateProblemsDNA(start, end):
  return [(i, (list1[:i], list2[:i])) for i in problemSizes]

# TEST THE TIMING OF THE CODE
global cache
inputList = read_fasta(2)
# size = 10000
list1 = seq0
list2 = seq1
problemSizes = [2**i for i in range(5, 14)]
DNAProblem = [(i, (list1[:i], list2[:i])) for i in problemSizes]

#since we expect the DP algorithm to be polynomial
timeProblems(DNAProblem, matchDP, fit = 'polynomial')
printAlign(list1[:16384], list2[:16384], cache)