#BRIAN CONROY
#MAXIMUM PARSIMONY NNI HEURISTIC ALGORITHM FOR PHYLOGENETIC TREES

import turtle


def render(tree, length, width):
    "Draws a given phylogenetic tree constrained by dimensions of" 
    "length and width."
    root = tree[0]
    leftTree = tree[1]
    rightTree = tree[2]
    if leftTree == (): 
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


def allTrees(tipList):
    "Generates a list of all possible trees for a given list of tip"
    "names. Trees are represented as tuples of the form (Root, Left, Right)"
    "where Left and Right are trees or, in the case of tips, both ()." 
    if len(tipList) == 2:
        return [ ("Anc", (tipList[0], (), ()), (tipList[1], (), ())) ]
    else:
        weirdo = (tipList[0], (), ())
        allTreesNoWeirdo = allTrees(tipList[1:])
        return reduce(lambda X, Y: X+Y, map(lambda tree: addTip(weirdo, tree), allTreesNoWeirdo))



def addTip(weirdo, tree):
    "Helper function for allTrees that adds a tip to a tree and returns the"
    "list of all trees that result from adding that tip to the tree." 
    if tree[1] == (): return [ ("Anc", weirdo, tree) ]
    else:
        option1 = [ ("Anc", weirdo, tree) ]
        option2 = map(lambda newLeft: ("Anc", newLeft, tree[2]), addTip(weirdo, tree[1]))
        option3 = map(lambda newRight: ("Anc", tree[1], newRight), addTip(weirdo, tree[2]))
        return option1 + option2 + option3

tipMapping = {"Linde": "AAA" , "Atwood": "AAT", "West": "AAG", "South":"ACC", "East":"GGG"}
tree = ('Anc', ('Atwood', (), ()), ('Anc', ('Linde', (), ()), ('West', (), ())))
tree1 = ('Anc', ('Atwood', (), ()), ('Linde', (), ()))
TREE = ('Anc', ('Anc',('Anc',('Atwood', (), ()),('South',(),())),('East',(),(),())), ('Anc', ('Linde', (), ()), ('West', (), ())))


Inf = float('inf')



def bestSinglePosition(tree, character, position, tipMapping):
    'bestSinglePosition takes as input a tree, a character input we are assuming'
    'is ancestral, a number indicating the position of DNA we are exploring, and'
    'a dictionary tipmapping of the DNA of all species examined. It returns a tuple'
    'containing the parsimony score and the tree but with labeled nodes.'
    oldscore = Inf 
    root = tree[0]
    leftTree = tree[1]
    rightTree = tree[2]
    answer = ()
    Calls = {}

    for leftchar in ('A','C','T','G'): #Nested loops allow us to assign all possible characters to nodes
        for rightchar in ('A','C','T','G'):
                newscore = 0 #for each character assignment, counts number of changes
               
                if leftTree[1] == (): #if we are at a tip in the left subtree
                    if character != tipMapping[leftTree[0]][position]: #check if a change was made
                        newscore = newscore + 1 #count the change
                    newleftTree = leftTree #keep the tip   

                if leftTree[0] == 'Anc': #if we are at a node in the left subtree
                    change = 0
                    if leftchar != character:
                        change = 1
                    if (leftTree,leftchar) not in Calls.keys():
                        newscore = newscore + change + bestSinglePosition(leftTree, leftchar, position, tipMapping)[0] #CALL
                        newleftTree = bestSinglePosition(leftTree, leftchar, position, tipMapping)[1]
                        Calls[(leftTree,leftchar)] = (newscore,newleftTree)   
                    else:
                        newscore =  Calls[(leftTree,leftchar)][0]
                        newleftTree = Calls[(leftTree,leftchar)][1]
                        
                if rightTree[1] == (): #the same cases except for the rightsubtree
                    if character != tipMapping[rightTree[0]][position]:
                        newscore = newscore + 1
                    newrightTree = rightTree

                if rightTree[0] == 'Anc':
                    change = 0
                    if rightchar != character:
                        change = 1
                    if (rightTree,rightchar) not in Calls.keys():
                        newscore = change + newscore + bestSinglePosition(rightTree, rightchar, position, tipMapping)[0]
                        newrightTree = bestSinglePosition(rightTree, rightchar, position, tipMapping)[1] #CALL
                        Calls[(rightTree,rightchar)] = (newscore,newrightTree)
                    else:
                        newscore =  Calls[(rightTree,rightchar)][0]
                        newrightTree = Calls[(rightTree,rightchar)][1]

                if newscore < oldscore: #here we find the best score
                    newtree = (character, newleftTree, newrightTree)
                    oldscore = newscore
                    answer = (newscore, newtree)  
    return answer





def add(tree1, tree2):
    'Takes two identical trees except with different labellings of internal'
    'nodes, and then returns a new tree structurally identical to the originals'
    'but with internal nodes labelled as a sequence of the previous individual'
    'node labellings.'
    if tree1[1] == ():
        return tree1
    root1 = tree1[0]
    root2 = tree2[0]
    return (root1 + root2, add(tree1[1],tree2[1]), add(tree1[2],tree2[2]))



def merge(treeList):
    'Takes a list of trees that are identical except that each one is'
    'labelled for a different position, and then returns a tree with merged'
    'internal node labelling.'
    return reduce(lambda tree1, tree2: add(tree1,tree2), treeList)

Inf = float('inf')



Groodyspecies = ["Groody", "Froody", "Snoody", "Wumph", "Glumph", "Humph"]
Groodymapping = {"Groody": "AACC", "Froody": "AACG", "Snoody": "AAGC", "Wumph": "CCAA", "Glumph": "CTAA", "Humph": "AAAA"}

def maximumParsimony(tipList,tipMapping):
    'Given a dictionary of character states for the species in tipList, creates'
    'the maximum parsimony phylogenetic tree for species in the tipList,'
    'labels internal nodes, returns the max parsimony score and renders the tree.'
    allTreesL = allTrees(tipList)
    positions = len(tipMapping[tipList[0]])
    optTree = ()
    optScore = Inf
    
    for tree in allTreesL:
        treeList = []
        Score = 0
        for position in range(0,positions):
            bases = ['A','C','T','G']
            bestScore = Inf
            bestTree = ()
            for character in bases:
                possibleTree = bestSinglePosition(tree, character, position, tipMapping)[1]
                score = bestSinglePosition(tree, character, position, tipMapping)[0]
                if score <  bestScore:
                    bestScore = score
                    bestTree = possibleTree
            treeList.append(bestTree)
            Score += bestScore
        if Score < optScore:
            optTree = merge(treeList)
            optScore = Score

    render(optTree, 60, 100) 
    return optScore 



import random



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
    return ('Anc',randomTree(leftTips),randomTree(rightTips))



def extract(indices,List):
    "Helper function for randomTree that returns elements of the List"
    "that have the given indices"
    answer = []
    for i in indices:
        answer.append(List[i])
    return answer 



def NNIheuristic(tipList, tipMapping):
    "Nearest neighbor interchange heuristic to create the maximum parsimony tree given a list of species and"
    "a dictionary of character states of those species."
    #randomly generate a tree from the list of vertices
    randomStart = randomTree(tipList)

    #score and label the tree
    scoreLabel = Parsimony(randomStart,tipList, tipMapping)
    startScore = scoreLabel[0]
    startLabel = scoreLabel[1]

    #increases each iteration of the loop
    counter = 0
    #will always increase with each iteration unless no improvement was made
    threshold = 1
    startTree = randomStart

    #while our parsimony score is still improving 
    while counter < threshold:
        #generate all NNIs (unlabelled) from the start tree
        NNIs = allNNIs(startTree)
        #score and label the NNIs 
        scoredNNIs = map(lambda tree: Parsimony(tree,tipList,tipMapping), NNIs)
        #keeps track of placement within the list of unlabelled NNIs
        index = 0
        #nonzero change indicates if a better tree has been found 
        change = 0
        for t in scoredNNIs: #extract any NNIs with better parsimony scores
            if t[0] < startScore:
                threshold += 1
                counter += 1
                startTree = NNIs[index]
                startLabel = t[1]
                startScore = t[0]
                change += 1
            index += 1
        if change == 0:
            counter += 1
        
    return (startScore, startLabel)
        
                                                                                                                                  

def allNNIs(Tree):
    "Generates all nearest neighbor interchange trees from the given Tree"
    "and returns them in a list."

    if Tree[1] == ():
        return [Tree,Tree]

    else:
        newTrees = []
        leftTree = Tree[1]
        rightTree = Tree[2]

        #If we can make swaps at the root node
        if leftTree[0] == 'Anc' and rightTree[0] == 'Anc':

            #make the swaps and add the new trees to the list 
            newTree1 = swap(leftTree,rightTree)[0]
            newTree2 = swap(leftTree,rightTree)[1]
            newTrees += [newTree1,newTree2]
     

            #save where we leave off
            rightFrame = ('Anc',leftTree)
            leftFrame = ('Anc',rightTree)

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
            rightFrame = ('Anc',leftTree)
            leftFrame = ('Anc',rightTree)

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
    "Helper function for leftGeorgi that takes list of subtrees,"
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


        
def swap(leftTree, rightTree):
    "Helper function for allNNIs that, when given a subtree rooted in"
    "an ancestral node that has 4 grandchildren, interchanges the"
    "nearest neighbors of that subtree and returns the two resultant"
    "subtrees."
    g1 = leftTree[1]
    g2 = leftTree[2]
    g3 = rightTree[1]
    g4 = rightTree[2]

    newLeft1 = ('Anc',g1,g3)
    newRight1 = ('Anc',g2,g4)

    newLeft2 = ('Anc',g1,g4)
    newRight2 = ('Anc',g3,g2)

    return [('Anc',newLeft1,newRight1),('Anc',newLeft2,newRight2)]
                

                

def Parsimony(tree,tipList,tipMapping):
    'Given a dictionary of character states for the species in tipList, creates'
    'the maximum parsimony phylogenetic tree for species in the tipList,'
    'labels internal nodes, returns the max parsimony score and renders the tree.'
    allTreesL = [tree]
    positions = len(tipMapping[tipList[0]])
    optTree = ()
    optScore = Inf
    
    for tree in allTreesL:
        treeList = []
        Score = 0
        for position in range(0,positions):
            bases = ['A','C','T','G']
            bestScore = Inf
            bestTree = ()
            for character in bases:
                possibleTree = bestSinglePosition(tree, character, position, tipMapping)[1]
                score = bestSinglePosition(tree, character, position, tipMapping)[0]
                if score <  bestScore:
                    bestScore = score
                    bestTree = possibleTree
            treeList.append(bestTree)
            Score += bestScore
        if Score < optScore:
            optTree = merge(treeList)
            optScore = Score 

    return (optScore, optTree)
        
    
catList = ["Felis Cat", "Lion", "Tiger", "Wild Cat", "Cheetah", "Puma", "Sabertooth", "Homotherium", "American Cat", "Leopard", "Dog", "Wolf", "Spotted Hyena", "Striped Hyena", "Black Bear", "Brown Bear", "Cave Bear", "African Wild Cat", "Chinese Desert Cat"]     
catMapping = {"Felis Cat": "CTTACCAAAATTATTAATCACTCATTCATCGACCTACCTGCCCCATCTAACATCTCAGCATGATGAAACT\
TCGGCTCCCTTCTAGGAGTCTGCCTAATCTTACAAATCCTCACCGGCCTCTTTTTGGCCATACACTACAC\
ATCAGACACAACAACCGCCTTTTCATCAGTTACCCACATCTGTCGCGACGTTAATTATGGCTGAATCATC\
CGATATTTACACGCCAACGGAGCTTCTATATTCTTTATCTGCCTGTACATACATGTAGGACGGGGAATAT\
ACTAC",
      "Lion": "CTTGTCAAAATTATTAATCACTCATTCATTGATCTTCCCACTCCACCCAATATCTCAGCATGATGAAACT\
TTGGCTCCTTATTAGGAGTATGTTTAATCCTACAAATTCTCACCGGCCTCTTTCTAGCCATACATTACAC\
ACCAGACACAATAACCGCTTTCTCATCAGTCACCCACATTTGCCGCGATGTAAACTATGGCTGAATTATC\
CGGTACCTACACGCCAACGGAGCCTCCATATTCTTTATCTGCCTATACATGCATGTAGGACGAGGAATAT\
ACTAT",
      "Tiger": "CTTATCAAAATTATTAATCACTCATTTATTGACCTACCCGCCCCATCCAATATTTCAGCATGATGAAACT\
TTGGCTCCTTACTAGGGGTGTGCTTAATCTTACAAATCCTCACTGGCCTCTTTCTAGCCATACACTACAC\
ATCAGACACAATAACCGCATTCTCATCAGTTACCCACATTTGCCGCGACGTAAACTACGGCTGGATTATC\
CGATATCTACATGCCAACGGAGCCTCCATATTCTTTATCTGTCTATACATGCACGTAGGACGAGGAATAT\
ACTAC",
      "Wild Cat": "CTTATCAAAATTATTAATCACTCATTCATCGACCTACCCGCCCCATCTAACATCTCAGCATGATGAAACT\
TCGGCTCCCTTCTAGGAGTCTGCCTAATCTTACAAATCCTCACCGGCCTCTTTTTGGCCATACACTACAC\
ATCAGACACAATAACCGCCTTTTCATCAGTTACCCACATCTGTCGCGACGTTAATTATGGCTGAATCATC\
CGATATTTACACGCCAACGGAGCTTCTATATTCTTTATCTGCCTGTACATACATGTGGGACGGGGAATAT\
ACTAC",
      "Cheetah": "CTTATCAAAATCGTTAATCACTCATTCATCGATTTACCCACCCCACCTAACATTTCAGCATGATGAAACT\
TCGGCTCCCTACTAGGAGTCTGCCTAGTCCTACAGATCCTAACCGGCCTTTTCCTAGCCATACACTACAC\
ATCAGACACAATAACCGCCTTTTCATCAGTTACTCACATCTGCCGCGACGTCAACTACGGCTGAATTATT\
CGATACATGCACGCCAACGGAGCCTCTATATTCTTTATCTGCCTATACATGCATGTAGGACGAGGAATAT\
ACTAC",
      "Puma": "CTTATCAAAATCATTAATCACTCATTTATTGATTTACCCACCCCATCCAACATTTCAGCATGATGAAACT\
TTGGCTCCCTACTAGGGGTCTGCCTAATCCTACAAATCCTAACCGGCCTCTTCCTGGCCATACACTATAC\
ATCAGACACAATGACTGCCTTTTCATCAGTCACTCACATCTGTCGTGACGTTAACTACGGCTGAATTATT\
CGGTACATACATGCCAACGGAGCCTCCATATTCTTTATCTGCCTATACATGCACGTGGGACGAGGAATAT\
ATTAT",
      "Sabertooth": "CTAATTAAAATTATCAACCACTCATTCATTGATTTACCCACCCCATCCAACATTTCAGCATGATGAAACT\
TCGGCTCCTTATTAGGAGTGTGCTTAATCTTACAAATCCTCACTGGCTTATTTCTAGCCATACATTATAC\
ACCAGATACAACAACCGCCTTCTCATCAGTTACCCACATTTGCCGTGATGTTAATTACGGCTGAATTATC\
CGATATATACACGCCAATGGAGCTTCCATATTCTTCATCTGCCTATATATACATGTAGGTCGAGGCATAT\
ACTAC",
      "Homotherium": "CTAATTAAAATCATCAACCAATCATTCATTGACTTACCTACCCCCTCCAACATCTCAGCATGATGAAACT\
TCGGATCCCTACTAGGCATTTGCCTAATTCTTCAAATCCTCACAGGCTTATTCCTAGCCATACACTACAC\
ATCAGACACAACAACTGCTTTCTCATCAATCGCCCATATTTGCCGTGACGTAAATTATGGTTGAATTATC\
CGATATATACACGCCAATGGAGCCTCTATATTCTTCATCTGTCTATACCTACATGTAGCTCGAGGAATTT\
ATTAC",
      "American Cat": "CTTATTAAAATCATTAATCACTCATTCATTGATCTACCCACCCCATCCAACATTTCAGCATGATGAAACT\
TCGGTTCCCTACTAGGGGTCTGCCTAATCCTACAAATCCTAACCGGCCTCTTCCTGGCTATACACTACAC\
ATCAGACACAATAACCGCCTTTTCATCAGTTACTCACATCTGTCGTGACGTCAATTACGGCTGAATTATT\
CGGTATATACACGCCAACGGAGCCTCCATATTCTTTATCTGCCTATACATGCACGTAGGGCGAGGAATAT\
ATTAC",
      "Leopard": "CTCATCAAAATTATTAATCACTCATTCATTGATCTCCCCGCTCCATCCAACATCTCAACATGATGGAACT\
TTGGCTCCCTATTAGGGGTATGTTTAATCCTACAAATTCTCACCGGCCTCTTTCTAGCCATACATTATAC\
ATCAGACACAACAACCGCTTTCTCATCAGTTACCCATATCTGCCGCGATGTAAATTATGGCTGAATTATC\
CGGTATCTACACGCCAATGGAGCCTCCATATTCTTTATCTGCCTATACATACATGTAGGACGAGGGATAT\
ACTAT",
      "Dog" : "CTAGCCAAAATTGTTAATAACTCATTCATTGACCTCCCAGCGCCGTCTAACATCTCTGCTTGATGGAACT\
TCGGATCCTTACTAGGAGTATGCTTGATTCTACAGATTCTAACAGGTTTATTCTTAGCTATGCACTATAC\
ATCGGACACAGCCACAGCTTTTTCATCAGTCACCCACATCTGCCGAGACGTTAACTACGGCTGAATTATC\
CGCTATATGCACGCAAATGGCGCTTCCATATTCTTTATCTGCCTATTCCTACATGTAGGACGAGGCCTAT\
ATTAC",
      "Wolf" : "CTAGCCAAAATTGTTAATAACTCATTCATTGACCTCCCAGCGCCGTCTAACATCTCTGCTTGATGGAACT\
TCGGATCCTTACTAGGAGTATGCTTGATTCTACAGATTCTAACAGGTTTATTCTTAGCTATGCACTATAC\
ATCGGACACAGCCACAGCTTTTTCATCAGTCACCCACATCTGCCGAGACGTTAACTACGGCTGAATTATC\
CGCTATATGCACGCAAATGGCGCTTCCATATTCTTTATCTGCCTATTCCTACATGTAGGACGAGGCCTAT\
ATTAC",
      "Spotted Hyena" : "CTCATTAAAATTATCAACAAATCATTCATTGACCTCCCCACCCCATCCAACATCTCGGCATGGTGAAATT\
TCGGGTCACTATTAGGAATCTGCTTAATCTTACAAATCCTGACAGGTCTATTCCTAGCCATACACTACAC\
ATCAGACACAACAACCGCCTTCTCATCAGTGACCCACATCTGCCGAGACGTAAACTACGGCTGAATCATC\
CGATACATACACGCCAACGGAGCTTCCATATTCTTCATCTGTCTATATATACATATCGGCCGAGGAATAT\
ACTAC",
      "Striped Hyena" : "CTCATTAAAATTGTCAACGAATCATTCATCGATCTCCCCACCCCATCCAACATCTCAGCATGATGAAACT\
TCGGATCGCTATTAGGAATCTGCCTAATCTTACAGATTCTGACAGGCCTATTTCTAGCCATACACTACAC\
ATCAGACACAACAACCGCCTTTTCATCAGTAACACACATCTGCCGAGACGTCAACTATGGCTGAATTATC\
CGATATATGCACGCCAACGGGGCCTCCATGTTTTTCATCTGCCTGTTCATGCACGCCGGCCGAGGAATGT\
ACTAC",
      "Black Bear" : "TTAGCTAAAATCATCAACAACTCACTCATTGATCTCCCAGCACCATCAAATATCTCAGCATGATGAAACT\
TCGGGTCCCTCCTCGGAGTATGTTTAGTACTACAAATTCTAACGGGCCTATTCCTAGCTATACACTATAC\
ATCAGACACAACTACAGCCTTTTCATCAATCACCCATATTTGCCGAGATGTTCACTACGGATGAATTATC\
CGATACATACATGCTAACGGAGCTTCCATATTCTTTATCTGCCTGTTCATGCACGTAGGACGGGGTCTGT\
ACTAT",
      "Brown Bear" : "TTAGCTAAAATCATCAACAACTCATTTATTGACCTTCCAACACCATCAAACATCTCAGCATGATGAAACT\
TTGGATCCCTCCTTGGAGTATGTTTAATTCTACAGATTCTAACAGGCCTGTTCCTAGCCATACACTATAC\
ACCAGACACAACCACAGCTTTTTCATCGGTCACCCACATTTGCCGAGACGTTCACTACGGGTGAGTTATC\
CGATATGTACATGCAAATGGAGCCTCCATCTTCTTTATCTGCCTATTTATGCACGTAGGACGGGGCCTGT\
ACTAT",
      "Cave Bear" : "TTAGCTAAAATCATCAACAACTCATTTATTGACCTCCCAACACCATCAAACATCTCAGCATGATGAAACT\
TTGGATCCCTCCTCGGAGTATGCTTAATTCTACAGATCCTAACAGGCCTGTTTCTAGCTATACACTACAC\
ATCAGACACAACCACAGCCTTTTCATCAATCACCCATATTTGCCGAGACGTTCACTACGGTTGAGTTATC\
CGATATATACATGCAAACGGAGCCTCCATATTCTTTATCTGTCTATTCATGCACGTAGGACGGGGCCTAT\
ACTAT",
      "African Wild Cat" : "CTTATCAAAATTATTAATCACTCATTCATCGATCTACCCGCCCCATCTAACATCTCAGCATGATGAAACT\
TCGGCTCCCTTCTAGGAGTCTGCCTAATCTTACAAATCCTCACCGGCCTCTTTTTGGCCATACACTACAC\
ATCAGACACAATAACCGCCTTTTCATCAGTTACCCACATCTGTCGCGACGTTAATTATGGCTGAATCATC\
CGATATTTACACGCCAACGGAGCTTCTATATTCTTTATCTGCCTGTACATACACGTAGGACGGGGAATAT\
ACTAC",
      "Chinese Desert Cat" : "CTTATTAAAATCATCAACCATTCATTCATTGACCTACCCACCCCATCCAACATCTCAGCATGATGAAATT\
TCGGCTCCCTATTAGGAGTCTGCCTAATCTTACAAATTCTCACCGGCCTCTTTCTAGCCATACACTACAC\
ATCAGACACAGTAACCGCTTTTTCATCAGTTACTCACATCTGTCGCGATGTTAATTACGGCTGAATCATC\
CGATACATACACGCCAATGGAGCTTCCATATTCTTTATCTGCCTATACATGCACGTAGGACGAGGAATAT\
ATTAC"}    
    
        
    
            
        
        
    
        
        
    
        

            
            
    

