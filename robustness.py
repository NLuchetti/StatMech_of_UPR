# -*- coding: utf-8 -*-
"""
Created on Fri Oct 25 14:49:49 2024

@author: nicol
"""

import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import numpy
from networkx_robustness import networkx_robustness

""""""""""" Network robustness on native models """""""""""
filename = 'matrix.txt';                                                                           #Native model adjacency matrix loading
M = np.loadtxt(filename, dtype='i', delimiter=' ');
G1=nx.from_numpy_array(M);


name = open('proteins.txt', 'r');                                                                  #Protein names list
nameline = name.readlines();
nodes = [];
for i in range(0, len(nameline)):
    nodes.append(nameline[i].split()[0])
res = [];
for sub in nodes:
    res.append(re.sub('\n', '', sub));
    

" Nodes attack "
def net_attack(graph, ranked_nodes):

    fraction_removed = []

    graph1 = graph.copy()
    nnodes = len(ranked_nodes)
    n = 0

    gcc = list(nx.connected_components(graph1))[0]

    gcc_size = float(len(gcc)) / nnodes

    fraction_removed.append((float(n) / nnodes, gcc_size))

    while gcc_size > 0.02:
        
        graph1.remove_node(ranked_nodes.pop())

        gcc = list(nx.connected_components(graph1))[0]
        gcc_size = float(len(gcc)) / nnodes
        n += 1
        fraction_removed.append((float(n) / nnodes, gcc_size))
        
    return fraction_removed

" Random attack "
df = pd.DataFrame(); df["id"], df["name"] = np.arange(0, len(res)), res;
resilience_random = net_attack(G1, list(df.id));

" Betweenness attack "
betw = nx.betweenness_centrality(G1);
df2 = pd.DataFrame(); df2["id"], df2["dim"], df2["name"] = np.arange(0, len(res)), betw, res;
df2 = df2.sort_values(by=['dim'], ascending=True);
resilience_betw = net_attack(G1, list(df2.id));

" Closeness attack "
clos = nx.closeness_centrality(G1);
df3 = pd.DataFrame(); df3["id"], df3["dim"], df3["name"] = np.arange(0, len(res)), clos, res;
df3 = df3.sort_values(by=['dim'], ascending=True);
resilience_clos = net_attack(G1, list(df3.id));

" Degree attack "
G1 = nx.read_edgelist("interactions.txt");                                                                     #Interactions list as two columns
deg = [lis[-1] for lis in G1.degree()];
df1 = pd.DataFrame(); df1["id"], df1["dim"], df1["name"] = np.arange(0, len(res)), deg, res;
df1 = df1.sort_values(by=['dim'], ascending=True);
resilience_deg = net_attack(G1, list(df1.name));

" Plot "
x = [k[0] for k in resilience_random]
y = [k[1] for k in resilience_random]

x1 = [k[0] for k in resilience_deg]
y1 = [k[1] for k in resilience_deg]

x2 = [k[0] for k in resilience_betw]
y2 = [k[1] for k in resilience_betw]

x3 = [k[0] for k in resilience_clos]
y3 = [k[1] for k in resilience_clos]

plt.plot(x, y, label="Random attack", marker='o', markersize=5, color='navy', lw=0.5, markerfacecolor='None')
plt.xticks(fontname = "Times New Roman");
plt.yticks(fontname = "Times New Roman");
plt.plot(x1, y1, label="Degree-based attack", marker='p', markersize=5, color='indigo', lw=0.5, markerfacecolor='None')
plt.xticks(fontname = "Times New Roman");
plt.yticks(fontname = "Times New Roman");
plt.plot(x2, y2, label="Betw.-based attack", marker='^', markersize=5, color='darkgreen', lw=0.5, markerfacecolor='None')
plt.xticks(fontname = "Times New Roman");
plt.yticks(fontname = "Times New Roman");
plt.plot(x3, y3, label="Clos.-based attack", marker='s', markersize=5, color='darkred', lw=0.5, markerfacecolor='None', alpha=0.5)
plt.xticks(fontname = "Times New Roman");
plt.yticks(fontname = "Times New Roman");
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlim(-0.05, 1.05)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.xlabel('Fraction of removed nodes', fontsize=14, fontname = "Times New Roman")
plt.ylabel('Fraction of LCC', fontsize=14, fontname = "Times New Roman")
plt.legend(loc="best", fontsize=14, prop={'family': 'Times New Roman'})
plt.savefig('image.svg', dpi=200)                                                                                  #Save figure
plt.show();

" Jumps identification in robustness trends "
n = [1, 2, 3];
for nn in n:
    globals()["diff"+str(nn)] = [];
    globals()["names"+str(nn)] = [];
    for i in range(0, len(globals()["y"+str(nn)])-1):
            globals()["diff"+str(nn)].append(globals()["y"+str(nn)][i]*100-globals()["y"+str(nn)][i+1]*100)

    for i in range(0, len(globals()["diff"+str(nn)])):
        if globals()["diff"+str(nn)][i] >= 10:
            globals()["names"+str(nn)].append(globals()["df"+str(nn)]["name"][i]);






""""""""""" Network robustness on configuration models """""""""""
name = ['human', 'rat', 'mouse', 'monkey', 'bull', 'rabbit',
        'chicken', 'zebrafish', 'moschito', 'worm', 'yeast', 'plant'];                                            #Organism names

for nn in name:
    idx = np.arange(1, len(name)); 
    ghd = [];
    for n in idx:
        
        filename = ''+str(nn)+'_mat.txt';                                                                          #Native model adjacency matrix loading
        M = np.loadtxt(filename, dtype='i', delimiter=' ');
        
        filename1 = ''+str(nn)+'_cm_'+str(n)+'.txt';                                                               #Configuration models adjacency matrices loading
        A = np.loadtxt(filename1, dtype='i', delimiter='\t');
        
                
        " GHD measurment to evaluate nearest and farthest configuration models "            
        m_o = sum(sum(M)); m_c = sum(sum(A));
        a_o = []; a_c = [];
        for i in range(0, len(M)): 
            for j in range(0, len(M)):
                a_o.append(M[i][j]-m_o/(len(M)*(len(M)-1)));
                a_c.append(A[i][j]-m_c/(len(A)*(len(A)-1)));
         
        a = [];
        for i in range(0, len(a_o)):
            a.append((a_o[i]-a_c[i])**2);
            
        ghd.append(sum(a)/(len(M)*(len(M)-1)));
        
        globals()['num'+str(name.index(nn)+1)] = [ghd.index(min(ghd))+1, ghd.index(max(ghd))+1];

for nn in name:
    for h in globals()['num'+str(name.index(nn)+1)]:
        filename1 = ''+str(nn)+'_cm_'+str(h)+'.txt';                                                                  #Configuration models adjacency matrices loading
        A = np.loadtxt(filename1, dtype='i', delimiter='\t');
        G1=nx.from_numpy_array(A);
        
        
        " Nodes attack "
        def net_attack(graph, ranked_nodes):
        
            fraction_removed = []
        
            graph1 = graph.copy()
            nnodes = len(ranked_nodes)
            n = 0
        
            gcc = list(nx.connected_components(graph1))[0]
        
            gcc_size = float(len(gcc)) / nnodes
        
            fraction_removed.append((float(n) / nnodes, gcc_size))
        
            while gcc_size > 0.02:
                
                graph1.remove_node(ranked_nodes.pop())
        
                gcc = list(nx.connected_components(graph1))[0]
                gcc_size = float(len(gcc)) / nnodes
                n += 1
                fraction_removed.append((float(n) / nnodes, gcc_size))
                
            return fraction_removed
        
        " Random attack "
        df = pd.DataFrame(); df["id"], df["name"] = np.arange(0, len(A)), np.arange(0, len(A));
        resilience_random = net_attack(G1, list(df.id));
        
        " Betweenness attack "
        betw = nx.betweenness_centrality(G1);
        df2 = pd.DataFrame(); df2["id"], df2["dim"], df2["name"] = np.arange(0, len(A)), betw, np.arange(0, len(A));
        df2 = df2.sort_values(by=['dim'], ascending=True);
        resilience_betw = net_attack(G1, list(df2.id));
        
        " Closeness attack "
        clos = nx.closeness_centrality(G1);
        df3 = pd.DataFrame(); df3["id"], df3["dim"], df3["name"] = np.arange(0, len(A)), clos, np.arange(0, len(A));
        df3 = df3.sort_values(by=['dim'], ascending=True);
        resilience_clos = net_attack(G1, list(df3.id));
        
        " Degree attack "
        deg = [lis[-1] for lis in G1.degree()];
        df1 = pd.DataFrame(); df1["id"], df1["dim"], df1["name"] = np.arange(0, len(A)), deg, np.arange(0, len(A));
        df1 = df1.sort_values(by=['dim'], ascending=True);
        resilience_deg = net_attack(G1, list(df1.name));
    
        " Plot "
        x = [k[0] for k in resilience_random]
        y = [k[1] for k in resilience_random]
        
        x1 = [k[0] for k in resilience_deg]
        y1 = [k[1] for k in resilience_deg]
        
        x2 = [k[0] for k in resilience_betw]
        y2 = [k[1] for k in resilience_betw]
        
        x3 = [k[0] for k in resilience_clos]
        y3 = [k[1] for k in resilience_clos]
    
        plt.plot(x, y, label="Random attack", marker='o', markersize=5, color='navy', lw=0.5, markerfacecolor='None')
        plt.xticks(fontname = "Times New Roman");
        plt.yticks(fontname = "Times New Roman");
        plt.plot(x1, y1, label="Degree-based attack", marker='p', markersize=5, color='indigo', lw=0.5, markerfacecolor='None')
        plt.xticks(fontname = "Times New Roman");
        plt.yticks(fontname = "Times New Roman");
        plt.plot(x2, y2, label="Betw.-based attack", marker='^', markersize=5, color='darkgreen', lw=0.5, markerfacecolor='None')
        plt.xticks(fontname = "Times New Roman");
        plt.yticks(fontname = "Times New Roman");
        plt.plot(x3, y3, label="Clos.-based attack", marker='s', markersize=5, color='darkred', lw=0.5, markerfacecolor='None', alpha=0.5)
        plt.xticks(fontname = "Times New Roman");
        plt.yticks(fontname = "Times New Roman");
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlim(-0.05, 1.05)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('Fraction of removed nodes', fontsize=14, fontname = "Times New Roman")
        plt.ylabel('Fraction of LCC', fontsize=14, fontname = "Times New Roman")
        plt.legend(loc="best", fontsize=14, prop={'family': 'Times New Roman'})
        plt.savefig(''+str(nn)+'_'+str(h)+'.svg', dpi=200)                                                             #Save figures
        plt.show();


    