# -*- coding: utf-8 -*-
"""
Created on Mon Jun  3 13:52:50 2024

@author: nicol
"""


####### Write interactions
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

filename = "string_interactions.tsv";
df1 = pd.read_csv(filename, sep = '\t');
filename = "interactions.txt";
f = open(filename, 'w');
# f.write('#node1'+'\t'+'node2'+'\n');
for i in range(0, len(df1)):
    f.write(str(df1["node1"][i])+'\t'+str(df1["node2"][i])+'\n');
f.close();



####### Write nodes and adjacency matrix
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

filename = "interactions.txt";
G = nx.read_edgelist(filename);
nodes = list(G.nodes());
filename = "proteins.txt";
f = open(filename, 'w');
for i in range(0, len(nodes)):
    f.write(str(nodes[i])+'\n');
f.close();
M = nx.to_numpy_array(G);
filename = "metrix.txt";
np.savetxt(filename, M,fmt='%d');



####### Read and compute network metrics
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

filename = "interactions.txt";
G1 = nx.read_edgelist(filename);
filename = "metrics.txt";
df1 = pd.read_csv(filename, sep = '\t');
edges1 = np.sum(nx.to_numpy_array(G1));
density1 = edges1 / (len(df1) * (len(df1) - 1));
diam1 = nx.diameter(G1);
bar1 = nx.barycenter(G1);
mean1 = [np.mean(df1.deg), np.mean(df1.clos), np.mean(df1.bet)];
M = nx.to_numpy_array(G1);
plt.imshow(M, cmap='gray_r');
plt.ylim([len(M), 0]);
plt.xlim([0, len(M)]);
plt.xticks([0, len(M)], labels = ['1', '%d'%len(M)], fontname = 'Arial', fontsize = 13);
plt.yticks([0, len(M)], labels = ['1', '%d'%len(M)], fontname = 'Arial', fontsize = 13);
filename = 'matrix.svg';
plt.savefig(filename, dpi = 220);
plt.show();



####### Local clustering coefficient and modularity and communities computation
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

G = nx.read_edgelist("interactions.txt");
nodes = list(G.nodes);
deg = [list(G.degree())[i][1] for i in range(0, len(nodes))];
c_c = [nx.clustering(G, nodes[i]) for i in range(0, len(nodes))];
av_l_cc = np.mean(c_c);

comm_net = nx.community.louvain_communities(G, weight='weight', resolution=1, threshold=1e-07, seed=None);
community_index = {n: i for i, com in enumerate(comm_net) for n in com};
mod = nx.community.modularity(G, comm_net);

fig, ax = plt.subplots(figsize=(50, 40))
pos = nx.spring_layout(G);
node_color = [community_index[n] for n in G];
node_size = [v * 500 for v in deg];
nx.draw_networkx(
    G,
    pos=pos,
    with_labels=False,
    node_color=node_color,
    node_size=node_size,
    edge_color="grey",
    alpha=0.5,
    width=5,
    cmap=plt.cm.RdYlBu,
)

# Resize figure for label readability
ax.margins(0.1, 0.05)
fig.tight_layout()
plt.axis("off")
filename = 'graph.svg'
plt.savefig(filename, dpi=200);
plt.show();



###### Compute normalized network metrics
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

filename = "metrics.txt";
df2 = pd.read_csv(filename, sep = '\t');

df2.degree = df2.degree / len(df2);
df2.clos = (len(df2)-1)*df2.clos;
df2.bet = (df2.bet-min(df2.bet))/(max(df2.bet)-min(df2.bet));

deg = np.mean(df2.degree);
cl = np.mean(df2.clos);
bet = np.mean(df2.bet);



###### Compute highest-degree proteins
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import scipy

filename = "metrics.txt";
df1 = pd.read_csv(filename, sep = '\t'); df1 = df1.drop(columns=['clos', 'bet']);
df2 = pd.read_csv(filename, sep = '\t'); df2 = df2.drop(columns=['deg', 'bet']);
df3 = pd.read_csv(filename, sep = '\t'); df3 = df3.drop(columns=['clos', 'deg']);
G = nx.read_edgelist("interactions.txt");
nodes = list(G.nodes);
df1["nodes"], df2["nodes"], df3["nodes"] = nodes, nodes, nodes;
df1 = df1.sort_values(by = ['deg'], ascending=False);
df2 = df2.sort_values(by = ['clos'], ascending=False);
df3 = df3.sort_values(by = ['bet'], ascending=False);



###### z-score with configuration models
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import numpy
from networkx_robustness import networkx_robustness
from scipy import stats 

df = pd.read_csv('metrics.txt', sep = '\t');
M = np.loadtxt("matrix.txt", dtype='i', delimiter='\t');
G = nx.from_numpy_array(M);
nodes = list(G.nodes());

df.clos = df.clos*(len(df.clos)-1);
df.bet = (df.bet-min(df.bet))/(max(df.bet)-min(df.bet));
m_cl = [np.mean(df.clos)];
m_bet = [np.mean(df.bet)];
m_avcl = [sum(nx.clustering(G, i) for i in range(0, len(nodes)))/len(nodes)];

for i in range(1, 11):
    filename = 'cm_metrics_'+str(i)+'.txt'
    globals()["df"+str(i)] = pd.read_csv(filename, sep = '\t');
    m_cl.append(np.mean((len(globals()["df"+str(i)])-1)*globals()["df"+str(i)].closeness));
    m_bet.append(np.mean((globals()["df"+str(i)].betweenness-min(globals()["df"+str(i)].betweenness))/(max(globals()["df"+str(i)].betweenness)-min(globals()["df"+str(i)].betweenness))));
for i in range(1, 11):    
    filename1 = 'cm_'+str(i)+'.txt';
    M = np.loadtxt(filename1, dtype='i', delimiter='\t');
    G = nx.from_numpy_array(M);
    m_avcl.append(sum(nx.clustering(G, i) for i in range(0, len(nodes)))/len(nodes));
z11 = stats.zscore(m_cl);
z21 = stats.zscore(m_bet);
z31 = stats.zscore(m_avcl);



###### Shortest path length and efficiency
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt
import networkx as nx
import numpy
from networkx_robustness import networkx_robustness
from scipy import stats 

G = nx.read_edgelist("interactions.txt");
pl1 = nx.average_shortest_path_length(G);
le1 = nx.local_efficiency(G);