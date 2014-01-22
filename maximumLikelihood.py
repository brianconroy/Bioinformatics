## MAXIMUM LIKELIHOOD PHYLOGENETIC TREE ALGORITHM ##
## Brian Conroy ##

# PART 1: Likelihood of a Tree #

import math
import random
import types

#Some test data 
tree = ("Anc1234", ("Anc2345", ("A",(),()),("B",(),())),("Anc3456",("C",(),()),("D",(),())))
mapping = {'A':'AAAA','B': 'ACGG', 'C': 'CCCC', 'D': 'CAGA', 'D': 'ACGG'}
species = ["A","B","C","D"]
branchLengths = {'A': 0.43967614295867019, 'Anc1234': 0.22935785751351712, 'B': 0.29232849303989683,
      'C': 0.62363004331713168, 'Anc3456': 0.90864948299498538,
      'D': 0.67706716537934086, 'Anc2345': 0.66529531867575342}

# Jukes-Cantor transition function
def transition(fromCharacter, toCharacter, branchLength):
    "Computes the probability of transitioning from 1 character to another"
    "over a given branch length. Assumes the Jukes-Cantor model of DNA"
    "evolution, which follows a Poisson distribution." 
    expTerm = (1.0 - math.exp((-4.0/3.0)*branchLength))
    if fromCharacter == toCharacter:
        return 1 - 0.75*expTerm
    else:
        return 0.25*expTerm

def singlePosition(tree, species, mapping, branchLengths, position, char):
    "Computes the likelihood of a tree at a single position given a list of"
    "species, dictionary of DNA, dictionary of branch lengths, "
    "and assumed ancestral character."
    if tree[0] in species:
        if char == mapping[tree[0]][position]:
            return 1
        else:
            return 0
        
    else:
        leftProbability = 0
        leftTree = tree[1]
        leftNode = leftTree[0]
        for base in ["A","T","C","G"]:
            leftSubProb = singlePosition(leftTree, species, mapping, branchLengths, position, base)
            probToSub = transition(char,base,branchLengths[leftNode])
            leftProbability += leftSubProb*probToSub

        rightProbability = 0
        rightTree = tree[2]
        rightNode = rightTree[0]
        for base in ["A","T","C","G"]:
            rightSubProb = singlePosition(rightTree, species, mapping, branchLengths, position, base)
            probToSub = transition(char,base,branchLengths[rightNode])
            rightProbability += rightSubProb*probToSub

    return leftProbability*rightProbability 

def treeLikelihood(tree, species, mapping, branchLengths):
    "Computes the likelihood of a tree given a list of species and dictionaries of"
    "DNA and branchlengths."
    numPositions = len(mapping[species[0]])
    likelihood = 1
    for p in range(0,numPositions):
        positionsProbability = 0
        for base in ["A","C","T","G"]:
            positionsProbability += 0.25*singlePosition(tree, species, mapping, branchLengths, p, base)
        likelihood = likelihood*positionsProbability

    return likelihood






# Part 2: Maximum Likelihood using NNI Heuristic #

def maximumLikelihood(species, mapping):
    "Generates the maximum likelihood tree given a list of species and"
    "dictionary of DNA character states."
    startTree = randomTree(species)
    startLengths = randomBranchLengths(startTree)
    startScore = treeLikelihood(startTree, species, mapping, startLengths)

    answerTree = startTree
    answerLengths = startLengths
    answerScore = startScore

    counter = 0 #Counts the number of iterations of the loop that don't improve score
    while counter < 4:
        NNIs = allNNIs(startTree)
        bestNNI = ()
        bestScore = 0
        bestLengths = {}
        
        for tree in NNIs:
            sLengths = {}
            for k in startLengths.keys():
                sLengths[k] = startLengths[k]
            bLengths = branchLengthChanger(tree,species,mapping,sLengths)
            score = treeLikelihood(tree,species,mapping,bLengths)
            if score > bestScore:
                bestNNI = tree
                bestScore = score
                bestLengths = bLengths
                
        if bestScore > startScore:
            startTree = bestNNI
            startLengths = bestLengths
            startScore = bestScore

            answerTree = bestNNI
            answerLengths = bestLengths
            answerScore = bestScore

            counter = 0

        else:
            startTree = bestNNI
            startLengths = bestLengths
            startScore = bestScor
            
            counter += 1

    return answerTree,answerLengths
        
        
        
    
        
                                
    

    

def randomBranchLengths(tree):
    "Generates a dictionary of random branch lengths chosen on [0,1]"
    "for the nodes and vertices of the tree."
    Nodes = leavesAndNodes(tree)
    branchLengths = {}
    for n in Nodes:
        branchLengths[n] = random.uniform(0,1)
    return branchLengths 
    
def leavesAndNodes(tree):
    "Returns a list of the leaves and nodes of a tree." 
    L = []
    if tree[1] == (): 
        return [tree[0]]
    else:
        L += [tree[0]] + leavesAndNodes(tree[1]) + leavesAndNodes(tree[2])
    return L

def randomTree(tipList):
    "Generates a random tree given a list of vertices."
    if len(tipList) == 1:
        return (tipList[0],(),()) 
    indices = range(0,len(tipList))
    groupSize = int(random.uniform(1,len(tipList)))
    leftVertices = random.sample(indices,groupSize)
    rightVertices = []
    for i in indices:
        if i not in leftVertices:
            rightVertices.append(i)
    leftTips = extract(leftVertices,tipList)
    rightTips = extract(rightVertices,tipList)
    tag = ''
    for t in range(1,5):
        tag += str(random.randint(0,9))
    
    return ('Anc'+ tag,randomTree(leftTips),randomTree(rightTips))

def extract(indices,List):
    "Helper function for randomTree that returns elements of the List"
    "that have the given indices"
    answer = []
    for i in indices:
        answer.append(List[i])
    return answer 

def allNNIs(Tree):
    "Generates all nearest neighbor interchange trees from the given Tree"
    "and returns them in a list."

    if Tree[1] == ():
        return [Tree,Tree]

    else:
        newTrees = []
        Anc = Tree[0]
        leftTree = Tree[1]
        rightTree = Tree[2]

        #If we can make swaps at the root node
        if leftTree[0][0:3] == "Anc" and rightTree[0][0:3] == "Anc":

            #make the swaps and add the new trees to the list 
            newTrees += swap(leftTree,rightTree, Anc)
     

            #save where we leave off
            rightFrame = (rightTree[0],leftTree)
            leftFrame = (leftTree[0],rightTree)

            #recurse, and attach the frame where we left off to all recursed NNI trees 
            leftNNITrees = Map(leftFrame,allNNIs(leftTree),'left')
            rightNNITrees = Map(rightFrame,allNNIs(rightTree),'right')

            #remove all trees that haven't changed
            lNNITrees = Remove(Tree,leftNNITrees)
            rNNITrees = Remove(Tree,rightNNITrees)

            #add new trees to the list
            newTrees += lNNITrees
            newTrees += rNNITrees    
            
        else:
            #save where we leave off
            rightFrame = (Anc,leftTree)
            leftFrame = (Anc,rightTree)

            #recurse, and attach the frame where we left off to all recursed NNI trees 
            leftNNITrees = Map(leftFrame,allNNIs(leftTree),'left')
            rightNNITrees = Map(rightFrame,allNNIs(rightTree),'right')

            #remove all trees that haven't changed
            lNNITrees = Remove(Tree,leftNNITrees)
            rNNITrees = Remove(Tree,rightNNITrees)

            #add new trees to the list
            newTrees += lNNITrees
            newTrees += rNNITrees 
                
        return newTrees


def Map(frame, branchList,side):
    "Helper function for allNNIs that takes list of subtrees,"
    "called branchList, and maps the remainder of the original tree"
    "onto the subtrees."
    if frame == []:
        return branchList
    completeTrees = []
    if side == 'right':
        for b in branchList:
            completeTrees.append((frame[0],frame[1],b))
    if side == 'left':
        for b in branchList:
            completeTrees.append((frame[0],b,frame[1]))
    return completeTrees

def Remove(tree,List):
    "Helper function for allNNIs that removes all trees from a list that"
    "are equal to tree."
    newList = []
    for t in List:
        if t != tree:
            newList.append(t)
    return newList 


        
def swap(leftTree, rightTree, Ancestor):
    "Helper function for allNNIs that, when given a subtree rooted in"
    "an ancestral node that has 4 grandchildren, interchanges the"
    "nearest neighbors of that subtree and returns the two resultant"
    "subtrees."
    g1 = leftTree[1]
    g2 = leftTree[2]
    g3 = rightTree[1]
    g4 = rightTree[2]

    newLeft1 = (leftTree[0],g1,g3)
    newRight1 = (rightTree[0],g2,g4)

    newLeft2 = (leftTree[0],g1,g4)
    newRight2 = (rightTree[0],g3,g2)

    return [(Ancestor,newLeft1,newRight1),(Ancestor,newLeft2,newRight2)]

def branchLengthChanger(tree,species,mapping,branchLengths):
    "Helper function for maximumLikelihood that returns the most likely"
    "dictionary of branchlengths for a tree derived from an initial"
    "random dictionary of branchLengths."
    print treeLikelihood(tree,species,mapping,branchLengths)
    iterator = 0
    speciesAndNodes=leavesAndNodes(tree)
    while iterator<6:
        lens = branchLengths
        for s in speciesAndNodes:
            singleLengthChanger(tree, species, mapping, lens, s)
        iterator += 1
    print treeLikelihood(tree,species,mapping,branchLengths)
    return branchLengths
    
def singleLengthChanger(tree, species,mapping, branchLengths, key):
    "Helper function for branchLengthChanger that adjusts a single branchlength"
    "to improve the likelihood of a tree. Returns the optimized length of the"
    "desired key." 
    originalLengths = branchLengths
    downLengths = {}
    upLengths = {}
    for k in branchLengths.keys():
        downLengths[k] = branchLengths[k]
        upLengths[k] = branchLengths[k]
    priorScore = treeLikelihood(tree, species, mapping, branchLengths)
    indicator = 0
    
    while True:
        if indicator == 0:
            oLength = originalLengths[key]
            upLengths[key] = oLength*1.05
            upScore = treeLikelihood(tree, species, mapping, upLengths)
            downLengths[key] = oLength/1.05 
            downScore = treeLikelihood(tree, species, mapping, downLengths)
            
            if upScore > priorScore and upScore >downScore:
                priorScore = upScore
                originalLengths[key] = upLengths[key]
                indicator = 1
            elif downScore > priorScore and downScore > upScore:
                priorScore = downScore
                originalLengths[key] = downLengths[key]
                indicator = 2
            else:
                break
        if indicator == 1:
            oLength = originalLengths[key]
            upLengths[key] = oLength*1.05
            upScore = treeLikelihood(tree, species, mapping, upLengths)
            if upScore > priorScore:
                priorScore = upScore
                originalLengths[key] = upLengths[key]
            else:
                break
        if indicator == 2:
            oLength = originalLengths[key]
            downLengths[key] = oLength/1.05
            downScore = treeLikelihood(tree, species, mapping, downLengths)
            if downScore > priorScore:
                priorScore = downScore
                originalLengths[key] = downLengths[key]
            else:
                break
    return originalLengths[key]


    
                                       
            
        
        
    
        

        


    

    
