# -*- coding: utf-8 -*- 
# creating a new network 
# drawing the created network 
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import networkx as nx 
epsilon =1e-3 
fname = "./adjacent.dat"
#fname2 = "./config_S_step0.dat"
#fname2 = "./config_S_step5000.dat"
fname2 = "./config_S_solution.dat"

E=[]

def set_Graph():
    global N,F
    global W, E
    i = 0 
    for line in open(fname,"r"):

        itemList = line.split(' ')
        del itemList[-1]
        if(i==0):
            n_V = len(itemList)
            N = np.arange(n_V)
            F = np.zeros(n_V)
            W = np.zeros((n_V,n_V))
        
        j = 0
        for x in itemList:
            W[i][j] = float(x) # I think this line is incorrect.
            W[j][i] = float(x)
            if( (i>j) and ( abs(W[i][j]) > epsilon ) ):
                E.append((i,j))
            j += 1
        i += 1
    
    k=0
    print "F=\n"
    for x in open(fname2,"r"):
        #if(i<n_V):
        F[k] = float(x) 
        k += 1

if __name__ == "__main__":
    set_Graph()
    print "N=\n", N 
    print "F=\n", F 
    print "W=\n", W
    print "E=\n", E 

    G = nx.Graph() 
    #N=[1,2,3,4,5,6,7,8,9,10]
    #E= [(1,2),(1,8),(2,3),(2,4),(4,5),(6,7),(6,8),(8,9),(8,10)]
    G.add_nodes_from(N)
    G.add_edges_from(E)
    pos=nx.circular_layout(G) #This is optional.
    nodes = G.nodes()
    color = F #[1 for n in nodes]
    print "coror=\n", color
    #nx.draw(G,with_labels=True, node_color=color, cmap=plt.cm.Blues) 
    nx.draw(G,pos,with_labels=True, node_color=color) 
    plt.axis('off')
    plt.title('an example')
    plt.savefig("image.png",format="PNG")
    plt.show()
