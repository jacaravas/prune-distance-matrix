#!/usr/bin/env python

import os, sys, re, argparse, csv, math
import random
import pandas as pd
import numpy as np
from tabulate import tabulate

def pruneMatrix (matrix, cutoff):
    # Note: I will be walking through and modifying matrix each pass
    #       A "one-pass" solution that ranks distances and cuts list
    #       would be possible, but you reduce your ability to modify
    #       the algorithm function 
    
    # Keep track of what was removed, when it was removed, and why 
    removed = []
    
    # I assume DISTANCE matrix of vals from 0 to 1.  Not identity matrix
    # like you might get from ANI.  Nancy Reagan taught me to "Say 'No' to ANI"
    lowestVal = 100
    
    # Save ties so that we can write tie-breaker code
    # Yes, I used a dictionary.  perl habits, but I think you'll see why
    # I prefer it to a list here (hints: uniqueness and laziness)
    lowestVals = {}
    
    # To make the loop control easily adaptable, I'm going to use a bool
    # to control it.  You can substitute new logic into the end of the
    # loop and just use it to switch the value of "done" or insert a break
    # without rewriting anything else
    done = False
    
    while not done:
        # Performance note: "iterX" functions are slow.  Performance should be
        # fine for anything not ludicrously large, though
        for rowName, row in matrix.iterrows():
            for colName, val in row.items():
                # Skip self comparisons.  They're uninformative
                if rowName == colName:
                    continue
                # Are we lower than or equal to our previous lowest value?
                if val <= lowestVal:
                    # Is it not a tie?
                    if val < lowestVal:
                        # Reset lowestVals to remove anything we found with a higher value
                        lowestVals = {}
                        
                    # Set lowestVal - in tie it just updates value
                    lowestVal = val

                    # You know how I say I'm lazy?  Here it is:
                    namesList = [rowName, colName]
                    # Sort the row/col names then make a string so that 
                    # A -> B comparison and B -> A comparison produce same value
                    namesList.sort()
                    normName = ",".join(namesList)
                    # A -> B and B -> A now have same name and value (pairwise comparison)
                    # So we can just clobber it in our dictionary and forget we were lazy
                    # in writing our loop to avoid looking at both comparisons 
                    lowestVals[normName] = [rowName, colName]
       
        # Since this code is using a hard cutoff, I'll just check that here
        # to see if anything needs to be done.
        if lowestVal > cutoff:
            done = True
            continue 
        # Now we do the work.  
        # I'll be using helper functions to keep it readable        
        
        # I'm merging the tie breaker code with the "pick which to remove"
        # code to reduce the amount of coding and debugging
        loser, loserScore = tieBreaker(lowestVals, matrix)
        # Get rid of the loser
        matrix = matrix.drop(loser, axis = 0)
        matrix = matrix.drop(loser, axis = 1)
        removed.append([loser, lowestVal, loserScore])
        # Reset loop values
        lowestVal = 100
        lowestVals = {}  
    return matrix, removed      
                

                
def tieBreaker (lowestVals, matrix):
        scores = {}
        # For each pair
        for pair in lowestVals:
            # For each sequence in pair
            for sequence in lowestVals[pair]:
                scores[sequence] = calcAverageDistance(sequence, matrix)
        # Discard the sequence with the lowest average distance
        loser = min(scores, key=lambda k: scores[k])
        return loser, scores[loser]        
                
def calcAverageDistance (sequence, matrix):
    row = matrix.loc[sequence]
    total = 0
    for colName, val in row.items():
        if colName == sequence:
            continue
        total += val
    return total/(matrix.shape[1] - 1)                        
                
def truncate(number):
    return math.floor(number * 100) / 100
    
def generateMatrix (n):
    labels = []
    for i in range(n):
        label = chr(65 + i)
        labels.append(label)
    
    matrix = {}
    # Create matrix "rows" in advance.  Will make enforcing
    # symmetry easier 
    for label1 in labels:
        row = {}
        matrix[label1] = row
    
    # Now populate matrix
    # We'll be super lazy and generate values twice [A][B] and [B][A] - previous value overwritten.
    # Note: Table is "dumb".  There is no sensible structure.  Just randomized pw distances.  
    for label1 in labels:
        for label2 in labels:
            # Self-self comparison = 0 (no differences)
            if label1 == label2:
                matrix[label1][label2] = 0.00
            else:
                # Generating range 0-0.1 for distances
                value = random.uniform(0, 0.10)
                # Truncate randomized distance to something readable
                # This also makes it much more likely to find a "tie" 
                # in the example data set we generate
                value = truncate(value)
                # Assign it to both pairs
                matrix[label1][label2] = value
                matrix[label2][label1] = value
    return matrix  
            
    
def main(arguments):
    
    # Generate a distance matrix
    dimension = 10

    distanceMatrix = generateMatrix(dimension)
    distanceDf = pd.DataFrame.from_dict(distanceMatrix, orient='index')
    print(tabulate(distanceDf, headers='keys', tablefmt='psql'))
    
    print (list(distanceDf.index.values))
    
    # Using a hard distance cutoff to control removal in this example 
    # "pruneMatrix" function.
    #
    # A slight modification to loop control in function would allow
    # you to "retain x" or "discard x" if you prefer.  When I used this,
    # I used a "retain x" style target to find the 100 most different
    # genomes in a large data set to capture maximum diversity.
    #
    # I have also used different loop controls to do a "remove until some 
    # other condition" style control.  Useful examples of this would be
    # "Keep removing sequences until your lowest pair is two different species"
    # or "Keep removing sequences until it tries to remove the last representative
    # of a species." Either approach might be a good fit for your dataset, but
    # we'll keep it simple here.
     
    cutoff = 0.05
    matrix, removed = pruneMatrix(distanceDf, cutoff)
    
    # A final report
    for sequence in removed:
        print(f'Removed {removed[0]} with pairwise distance {removed[1]} and average distance {removed[2]}')
    print()
    print ("Kept:")
    for sequence in matrix.index.values.tolist():
        print (sequence)         
    
    
if __name__=='__main__':
    main(sys.argv[1:])