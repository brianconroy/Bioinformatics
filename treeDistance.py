#Aggregation clustering algorithm based on Robinson-Foulds distance metric
#Trees are binary, unrooted, and in standard Newick format


#PART 1: ROBINSON-FOULDS DISTANCE METRIC
T1 = (0, (1, (2, (3, 4))))
T2 = (0, (1, (3, (2, 4))))
T3 = (0, ((1, (2, 3)), (7, (6, (4, 5)))))
T4 = (0, ((2, (1, 3)), (6, (4, (5, 7)))))

BITS = 64

import random

import types

from Caesal import* # A data set of 450 phylogenetic trees on 51 tips
#obtained from maximum parsimony searches on regions of the
#chloroplast genome from the Caesalpinia family. 

def distance(Tree1, Tree2):
    "Computes the Robinson-Foulds distance between two binary unrooted trees."
    labels = {}
    leaves = extract(Tree1)
    for l in leaves:
        labels[l] = int(randomString(),2)
    prints1 = traverser(Tree1[1],labels,leaves)
    prints2 = traverser(Tree2[1],labels,leaves)
    distance = 0
    for p in prints1:
        if p not in prints2:
            distance += 1
    return distance
    
def traverser(Tree,labels,leaves):
    "Returns a list of fingerprints of internal edges of the Tree."
    List = []
    
    if Tree[0] in leaves and Tree[1] in leaves:
        List += [labels[Tree[0]]^labels[Tree[1]]]

    if Tree[0] not in leaves and Tree[1] not in leaves:
        leftPrints = traverser(Tree[0],labels,leaves)
        rightPrints = traverser(Tree[1],labels,leaves)
        List += leftPrints + rightPrints
        leftRecent = leftPrints[len(leftPrints)-1]
        rightRecent = rightPrints[len(rightPrints)-1]
        List += [leftRecent^rightRecent]

    if Tree[0] not in leaves and Tree[1] in leaves:
        hold = Tree[1]
        remainders = traverser(Tree[0],labels,leaves)
        closest = remainders[len(remainders)-1]
        List += remainders + [labels[hold]^closest]

    if Tree[0] in leaves and Tree[1] not in leaves:
        hold = Tree[0]
        remainders = traverser(Tree[1],labels,leaves)
        closest = remainders[len(remainders)-1]
        List += remainders + [labels[hold]^closest] 

    return List

def extract(Tree):
    "Returns a list of the leaves of a binary unrooted tree, where leaves are"
    "of type integer" 
    L = []
    if type(Tree) is types.IntType: #or types.StringType:
        return [Tree]
    else:
        if type(Tree[0]) is types.IntType: #or types.StringType:
            L += [Tree[0]] + extract(Tree[1])
        else:
            left = Tree[0][0]
            right = Tree[0][1]
            L += extract(left) + extract(right) + extract(Tree[1])
    return L
             
def randomString():
    "Returns a random string of 0s and 1s of length given by the"
    "global variable BITS."
    string = ""
    for b in range(0,BITS):
        number = random.randint(0,1)
        if number == 0:
            string += "0"
        else:
            string += "1"
    return string 


#PART 2: AGGREGATION CLUSTERING



def cluster(treeList, minClusters, maxClusters):
    "For a list of trees, this function creates clusters of the trees"
    "of all sizes between the specified minimum and maximum cluster sizes."
    "For each"
    "cluster of a given size, it returns the following statistics:"
    "The number of trees in each cluster within that cluster, the average"
    "diameter of all clusters within that cluster, and the diameter of the"
    "cluster."
    Distances = {} # A dictionary of robinson foulds distances b/ clusters

    treeCount = [] #Lists the counts of the number of trees in each cluster for given size
    clusterDiam = [] #Lists the diameters of clusters within each cluster of given size
    avgDiam= [] #Lists the average diameter over all clusters within each cluster of given size
    specificitiesL = [] 
    
    clusterList = [] #makes the trivial starting cluster 
    for tree in treeList:
        clusterList.append([tree])

    sizes = range(minClusters,len(treeList))
    sizes.reverse()
    for i in sizes: #make all clusters of desired sizes
        bestDist = float(Inf)
        BestClusters = []
        for c in clusterList: #finds the closest clusters 
            for C in clusterList:
            
                if c != C:
                    if (tuple(c),tuple(C)) in Distances.keys():
                        dist = Distances[(tuple(c),tuple(C))]
                    else:
                        dist = clusterDifference(c,C)
                        Distances[(tuple(c),tuple(C))] = dist
                    if dist < bestDist:
                        bestDist = dist
                        BestClusters = [c, C]
                
        newList = [] #Merges the selected clusters into a new cluster
        for cluster in clusterList:
            if cluster not in BestClusters:
                newList.append(cluster)
        newList.append(BestClusters[0] + BestClusters[1])
        clusterList = newList
    
        if i <= maxClusters and i >= minClusters: 

            Diameters = []
            for cluster in newList: #Computes diameters of each cluster w/in clustering
                maxDiam = 0
                for t in cluster:
                    index = 1
                    for T in cluster[1:]:
                        if t != T:
                            diam = Distances[(tuple([t]),tuple([T]))]
                            if diam > maxDiam:
                                maxDiam = diam
                    index += 1
                Diameters.append(maxDiam)
            clusterDiam.append(Diameters)

            S = 0
            for d in Diameters: #Computes the average diameter 
                S += d
            avgDiam.append(float(S)/float(len(Diameters)))

            treeNums = []
            for cluster in newList: #Counts the number of trees w/in each cluster
                treeNums.append(len(cluster))
            treeCount.append(treeNums)


            specificities = []
            for cluster in newList:
                specificities.append(specificity(cluster,1.0))
            specificitiesL.append(specificities)
            
    return treeCount,clusterDiam,avgDiam,specificitiesL
                
        
def clusterDifference(cluster1, cluster2):
    "Computes the difference between two clusters, defines as the sum of distances"
    "between each pair of trees in the clusters divided by their product."
    size1 = len(cluster1)
    size2 = len(cluster2)
    distanceSum = 0
    for t in cluster1:
        for T in cluster2:
            distanceSum += distance(t,T)
    return distanceSum/(size1*size2)

Inf = float('Inf')

def findClosestClusters(clusterList):
    "Finds the 2 clusters with minimal RF distance and returns the original"
    "list of clusters with those 2 clusters merged."
    bestDistance = Inf
    bestClusters = ()
    
    for c in clusterList:
        index = 0
        for C in clusterList[index:len(clusterList)]:
            if c != C:
                diff = clusterDifference(c,C)
                if diff<Inf:
                    bestDistance = diff
                    bestClusters = (c,C)
        index += 1

    newList = [bestClusters[0]+bestClusters[1]]
    for c in clusterList:
        if c not in bestClusters:
            newList.append(c)

    return newList


#PART 3: CONSENSUS AND SPECIFICITY

def specificity(cluster, threshold):
    "Computes the specificity of a cluster (list of trees) using a given"
    "threshold value."

    labels = {}
    leaves = extract(cluster[0])
    for l in leaves:
        labels[l] = int(randomString(),2)

    allPartitions = []
    for tree in cluster:
        partitions = traverser(tree,labels,leaves)
        allPartitions.append(partitions)
        
        
    numConPartitions = partitionCounter(allPartitions,threshold)
    specificity = float(numConPartitions)/len(allPartitions[0])
    return specificity 
    
def partitionCounter(partitionList,threshold):
    "Given a list of partitions, returns the number of"
    "partitions that occur with a frequency at least equal to the"
    "threshold value."

    consensusPartitions = 0
    alreadyCounted = []

    for partitions in partitionList: #for each list of partitions
        for p in partitions: #for each partition within that list
            occurrencesOfP = 0
            for compPartitions in partitionList:
                if p not in alreadyCounted:
                    if p in compPartitions:
                        occurrencesOfP += 1
            alreadyCounted.append(p)

            if float(occurrencesOfP)/len(partitionList) >= threshold:
                consensusPartitions += 1 
    
    return consensusPartitions       
                    
                
        
        
            

    
            
            
            
            

        
        
        
    
    

        
         
