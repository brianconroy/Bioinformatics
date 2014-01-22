#BRIAN CONROY
#ULTRAMETRIC TREE RENDERING ALGORITHM 

import turtle 



tipList = [ "A", "B", "C", "D", "E"]
matrix = [ [0, 8, 8, 5, 3],
       [8, 0, 3, 8, 8],
       [8, 3, 0, 8, 8],
       [5, 8, 8, 0, 5],
       [3, 8, 8, 5, 0]]



def render(tree, length, width):
    "Draws a given phylogenetic tree constrained by dimensions of" 
    "length and width."
    root = tree[0]
    leftTree = tree[1]
    rightTree = tree[2]
    if leftTree == []: 
        turtle.dot(10)
        turtle.write(root , font=("Arial", 20, "normal"))
        return
    else:
        turtle.dot(10)
        turtle.write(root, font=("Arial", 20, "normal"))
        turtle.left(90)
        turtle.forward(width)
        turtle.right(90)
        turtle.forward(length)
        render(leftTree, 0.5*length, 0.5*width) 
        turtle.back(length)
        turtle.left(90)
        turtle.backward(2*width)
        turtle.right(90)
        turtle.forward(length)
        render(rightTree, 0.5*length, 0.5*width)
        turtle.back(length)
        turtle.right(90)
        turtle.back(width)
        turtle.left(90)
        return



def buildTree(tipList, matrix):
    "Constructs a phylogenetic tree from an ultrametric matrix"
    "describing the evolutionary distances between the species"
    "listed in tipList." 
    Max = findMaxRowIndex(matrix)
    if len(tipList) == 1:
        return [tipList[0],[],[]]
    Dict = distanceGroups(tipList,matrix)
    dList = Dict.keys()
    dList.sort()
    dList.reverse()
    subTrees = []
    for d in dList:   
        indexList = Dict[d]  
        newList = subMatrix(indexList,tipList,matrix)[0]
        newMatrix = subMatrix(indexList,tipList,matrix)[1]
        subTrees.append(buildTree(newList,newMatrix))
    return assemble(dList,subTrees, tipList[Max])
    


def assemble(dList, subTrees, Max):
    "Assembles subtrees from the list subTrees according to distances, given"
    "in the list dList, along the leftmost spine of the tree from the root to"
    "the species 'Max' "
    root = dList[0]
    if len(dList) == 1:
        return [root,subTrees[0],[Max,[],[]]]
    leftTree = subTrees[0]
    newList = dList[1:]
    newTrees = subTrees[1:]
    return [root,leftTree, assemble(newList,newTrees,Max)]
    


def distanceGroups(tipList,matrix):
    "Groups indices of species based on their distances from species 'Max' "
    Max = findMaxRowIndex(matrix)
    Dict = {}
    for d in matrix[Max]:
        if d != 0:
            Dict[d] = []
    counter = 0
    for d in matrix[Max]:
        if d == 0:
            counter += 1 
        if d != 0: 
            Dict[d] += [counter]
            counter += 1
    return Dict 



def subMatrix(indexList, tipList, matrix):
    'Constructs a submatrix from the original matrix for the given'
    'list of indices, indexList, and returns that and the list'
    'of species for which the tree was constructed.'
    subMatrix = []
    for r in indexList:
        row = []
        for p in indexList:
            row += [matrix[r][p]]
        subMatrix += [row]
    subList = []
    for i in indexList:
        subList += [tipList[i]]
    return (subList, subMatrix)

            

def findMaxRowIndex(matrix):
    'Returns the index of the row in the matrix containing the '
    'max distance. '
    maxDistance = 0
    maxRowIndex = 0
    rowCounter = 0
    for r in matrix:
        if max(r) > maxDistance:
            maxRowIndex = rowCounter
            maxDistance = max(r)
        rowCounter += 1
    return maxRowIndex



        
    
    
