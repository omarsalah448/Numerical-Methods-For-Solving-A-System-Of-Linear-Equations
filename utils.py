import re
import os
import time
import copy
import matplotlib.pyplot as plt

def column(mat, i):
    return [row[i] for row in mat]

def simpleParser(eqn):
    # handle power case
    eqn=re.sub("\^","**", eqn)
    # handle digits
    eqn=re.sub(r"([0-9])([^+-.*[0-9]])", r"\1*\2", eqn)
    eqn=re.sub(r"([0-9])([a-zA-Z])", r"\1*\2", eqn)
    # handle when no coefficients
    eqn=re.sub(r"(\w*(?<![0-9]\*)[a-zA-Z])", r"1*\1", eqn)
    # handle y= case
    eqn=re.sub(".* =|.=", "", eqn)
    return eqn

def extractCoefficients(eqn):
    eqn = simpleParser(eqn)
    # only return equations
    eqn = str(re.findall("([\.*\-0-9]+(\**)[a-zA-Z][0-9]*)" ,eqn))
    # remove variables from equation
    eqn = re.sub(r"([0-9]+)(\*)([a-zA-Z])([0-9]*)", r"\1", eqn)
    # return all extracted values
    result = [float(res) for res in re.findall("[\-\.0-9]+", eqn)]
    return result

def extractResults(eqn):
    eqn = simpleParser(eqn)
    # only return numbers
    eqn = re.sub(r"(-*)([\.*0-9]+\*[a-zA-Z][0-9]*)", r"", eqn)
    # remove + sign
    eqn = re.sub("\+*", "", eqn)
    # return all extracted values
    result = [float(res)*-1 for res in re.findall("[\.\-0-9]+", eqn)]
    return result

def extractAllResults(eqns):
    b = []
    for eqn in eqns:
        b += extractResults(eqn)
    return b

def extractVariables(eqn):
    eqn = simpleParser(eqn)
    return re.findall(r"([a-zA-Z][0-9]*)", eqn)

def extractAllVariables(eqns):
    var = []
    for eqn in eqns:
        var += extractVariables(eqn)
    # remove duplicates then sort
    var = list(dict.fromkeys(var))
    var.sort()
    return var

def coefficientsCheck(eqns, variables, coefficients):
    # pass over all eqns
    for i in range(len(eqns)):
        for (j, var) in enumerate(variables):
            # if a variable is missing add a
            # 0 as its coefficient
            if var not in eqns[i]:
                coefficients[i].insert(j, 0)

def plot(x, y, xlabel, ylabel):
    plt.plot(x, y, marker='o', label=ylabel);
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show(block=False)

def outputToFile(access_method, method_name, mat, no_eqns, precision, excecution_time):
    f = open(os.path.join(os.getcwd(),"output.txt"), access_method)
    iterative = False
    no_iter = 0
    if type(mat[0]) == list:
        iterative = True
        no_iter = len(mat) - 1
    if(f):
        f.write("Method Name: {}\n".format(method_name))
        # if iterative method
        if iterative:
            f.write("Iteration\t")
        for i in reversed(range(no_eqns)):
            f.write("x{}\t\t".format(i))
        f.write("\n")
        for i in range(len(mat)):
            # if iterative method
            if iterative:
                f.write("{}\t\t".format(i))
                for j in range(no_eqns):
                    if(len(str(mat[i][j])[:precision]) > 7):
                        f.write("{}\t".format(str(mat[i][j])[:precision]))
                    else:
                        f.write("{}\t\t".format(str(mat[i][j])[:precision]))
                f.write("\n")
            # if direct method
            else:
                if(len(str(mat[i])[:precision]) > 7):
                    f.write("{}\t".format(str(mat[i])[:precision]))
                else:
                    f.write("{}\t\t".format(str(mat[i])[:precision]))
        f.write("\n")
        f.write("\n\nNumber Of Iterations = {}\n".format(no_iter))
        f.write("Precision = {}\n".format(precision))
        f.write("Excecution Time = {:.5f} secs\n".format(excecution_time))
        f.write("\n\n\n")
    f.close()


# function to print matrix for testing
def print_matrix( matrix, message ):
    print( message )
    for i in matrix:
        for j in i:
            print( "%.4f" % j, end = "\t\t" )
        print()
    print()

# function to print vector for testing
def print_vector( vector, message ):
    print( message )
    for i in vector:
        print( "%.4f" % i, end = "\t\t" )
    print()

########################
##### GAUSS-SEIDEL #####
########################

# function to check for convergence
def converge( coefficients ):
    convergence_flag = 1
    for i in range( len( coefficients ) ):
        sum = 0
        for j in range( len( coefficients[0] ) ):
            if i == j:
                continue
            sum += abs( coefficients[i][j] )
        if sum > abs( coefficients[i][i] ):
            convergence_flag = 0
    return convergence_flag

# main function to call for gauss-seidel
def gauss_seidel( eqns, tolerance, iterations, initial_values = [] ):
    coefficients = [extractCoefficients(eqn) for eqn in eqns]
    b = extractAllResults(eqns)
    variables = extractAllVariables(eqns)
    coefficientsCheck(eqns, variables, coefficients)
    # convergence check
    if not converge( coefficients ):
        return "The coefficient matrix is not diagonally dominant!"
    
    # if no  initial values entered, assume zeros
    if not len( initial_values ):
        initial_values = [ 0 for i in range( len( coefficients ) ) ]

    results = [ [0] * len( coefficients ) for i in range( iterations + 1 ) ]
    for i in range( len( initial_values ) ):
        results[0][i] = initial_values[i]

    # loop on the number of iterations
    for i in range( 1, iterations + 1 ):
        tolerance_flag = 1
        # loop on rows
        for j in range( len( coefficients ) ):
            results[i][j] = b[j]
            # loop on columns
            for k in range( len( coefficients[0] ) ):
                if k == j:
                    continue
                elif k < j:
                    results[i][j] -= coefficients[j][k] * results[i][k]
                else:
                    results[i][j] -= coefficients[j][k] * results[i-1][k]
            results[i][j] = results[i][j] / coefficients[j][j]
            
            # calculate relative error and compare with given tolerance
            relative_error = abs( ( results[i][j] - results[i-1][j] ) / results[i][j] )
            if relative_error > tolerance:
                tolerance_flag = 0
        
        if tolerance_flag:
            final_results = results[0:i] 
            break

    return final_results

##################
##### JACOBI #####
##################

# main function to call for jacobi iteration method
def jacobi( eqns, tolerance, iterations, initial_values = [] ):
    coefficients = [extractCoefficients(eqn) for eqn in eqns]
    b = extractAllResults(eqns)
    variables = extractAllVariables(eqns)
    coefficientsCheck(eqns, variables, coefficients)
    # convergence check
    if not converge( coefficients ):
        return "The coefficient matrix is not diagonally dominant!"

    # if no  initial values entered, assume zeros
    if not len( initial_values ):
        initial_values = [ 0 for i in range( len( coefficients ) ) ]

    results = [ [0] * len( coefficients ) for i in range( iterations + 1 ) ]
    for i in range( len( initial_values ) ):
        results[0][i] = initial_values[i]

    # loop on the number of iterations
    for i in range( 1, iterations + 1 ):
        tolerance_flag = 1
        # loop on rows
        for j in range( len( coefficients ) ):
            results[i][j] = b[j]
            # loop on columns
            for k in range( len( coefficients[0] ) ):
                if k == j:
                    continue
                else:
                    results[i][j] -= coefficients[j][k] * results[i-1][k]
            results[i][j] = results[i][j] / coefficients[j][j]

            # calculate relative error and compare with given tolerance
            relative_error = abs( ( results[i][j] - results[i-1][j] ) / results[i][j] )
            if relative_error > tolerance:
                tolerance_flag = 0
        
        if tolerance_flag:
            final_results = results[0:i] 
            break

    return final_results

################################
##### GAUSSIAN-ELIMINATION #####
################################

# function to carry out partial pivotting with scaling 
def pivot( augmented, index ):
    max = index
    for i in range( index, len( augmented ) ):
        if abs( augmented[i][index] ) > abs( augmented[max][index] ):
            max = i
    
    augmented[index], augmented[max] = augmented[max], augmented[index]
    return augmented

# function for the forward elimination process
def eliminate( augmented, index ):
    eliminated = copy.deepcopy( augmented )
    for i in range( index + 1, len( augmented ) ):
        for j in range( len( augmented[0] ) ):
            eliminated[i][j] = augmented[index][j]*(-1*augmented[i][index]) / augmented[index][index] + augmented[i][j]
    return eliminated 

# function to carry out the backward substitution
def substitute( coefficients, b ):
    results = []
    results = [0 for i in range( len( b ) )] 
    for i in range( len( coefficients ) - 1, -1, -1 ):
        sum = 0
        for j in range( 0, len( coefficients[0] ) ):
            sum = coefficients[i][j] * results[j] + sum
        results[i] = ( b[i] - sum ) / coefficients[i][i]
    return results

# main function to call for gaussian-elimination 
def gaussian_elimination( eqns ):
    coefficients = [extractCoefficients(eqn) for eqn in eqns]
    print(coefficients)
    b = extractAllResults(eqns)
    print(b)
    variables = extractAllVariables(eqns)
    print(variables)
    coefficientsCheck(eqns, variables, coefficients)
    # create augmented matrix
    augmented = copy.deepcopy( coefficients )
    for i in range( len( coefficients ) ):
        augmented[i].append( b[i] )

    # use partial pivotting in the eliminating process
    for i in range( len( coefficients ) - 1 ):
        augmented = pivot( augmented, i )
        augmented = eliminate( augmented, i )

    # divide the augmented matrix back into coefficients and b 
    for i in range( len( coefficients ) ):
        for j in range( len( coefficients[0] ) ):
            coefficients[i][j] = augmented[i][j]
    for i in range( len( b ) ):
        b[i] = augmented[i][len( augmented[0] ) - 1]

    # check for special cases
    if coefficients[len( coefficients ) - 1][len( coefficients[0] ) - 1] == 0:
        if b[len( coefficients ) - 1] == 0:
            return "Infinite number of solutions exist!"
        else:
            return "Solution doesn't exist!"

    results = substitute( coefficients, b )

    return results
    
########################
##### GAUSS-JORDON #####
########################

# function that carries out the scaling for gauss-jordon
def scale( augmented, index ):
    temp = augmented[index][index]
    for i in range( len( augmented[0] ) ):
        augmented[index][i] = augmented[index][i] / temp
    return augmented

def jordon_eliminate( augmented, index ):
    eliminated = copy.deepcopy( augmented )
    for i in range( len( augmented ) ):
        if i == index:
            continue
        for j in range( len( augmented[0] ) ):
            eliminated[i][j] = augmented[index][j]*(-1*augmented[i][index]) + augmented[i][j]
    return eliminated

# main function to call for gauss-jordon
def gauss_jordon( eqns ):
    coefficients = [extractCoefficients(eqn) for eqn in eqns]
    b = extractAllResults(eqns)
    variables = extractAllVariables(eqns)
    coefficientsCheck(eqns, variables, coefficients)
    # create augmented matrix
    augmented = copy.deepcopy( coefficients )
    for i in range( len( coefficients ) ):
        augmented[i].append( b[i] )

    # use scaling and partial pivotting in the eliminating process
    for i in range( len( coefficients ) ):
        augmented = pivot( augmented, i )
        augmented = scale( augmented, i )
        augmented = jordon_eliminate( augmented, i )

    # get the results from the augmented matrix
    results = []
    results = [0 for i in range( len( b ) )] 
    for i in range( len( results ) ):
        results[i] = augmented[i][len( augmented[0] ) - 1]

    return results

########################
##### LU-COMPOSITION #####
########################

def substitute_lu(a, n, b):
    y = [0]*n
    y[0] = b[0]
    # forward substitution
    for i in range(1, n):
        summ = b[i]
        for j in range(i):
            summ = summ - a[i][j] * y[j]
        y[i] = summ
    for i in range(n):
        for j in range(i):
            a[i][j] = 0
    return substitute(a, y)

# main function to call for lu-composition
def LUComp(eqns):
    coefficients = [extractCoefficients(eqn) for eqn in eqns]
    b = extractAllResults(eqns)
    variables = extractAllVariables(eqns)
    coefficientsCheck(eqns, variables, coefficients)
    n = len(coefficients)
    # decomposition
    for k in range(n-1):
        for i in range(k+1, n):
            factor = coefficients[i][k] / coefficients[k][k]
            coefficients[i][k] = factor
            for j in range(k+1, n):
                coefficients[i][j] = coefficients[i][j] - factor * coefficients[k][j]
    return substitute_lu(coefficients, n, b)