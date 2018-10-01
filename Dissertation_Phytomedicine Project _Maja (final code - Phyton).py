
# coding: utf-8

# # Pytomedicine Project

# Building the undirected weighted graph 

# In[40]:


# import libraries
import pandas as pd
import numpy as np
from pandas import DataFrame
import matplotlib
import matplotlib.pyplot as plt
from heapq import heappush, heappop 
from itertools import count
import networkx as nx 
import random
import warnings 
from networkx.drawing.nx_agraph import graphviz_layout
warnings.filterwarnings("ignore")
import importlib
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
import Networks_functions as nf
from networkx.algorithms import tree
from collections import Counter
from pandas import ExcelWriter
from pandas import ExcelFile
import matplotlib.patches as mpatches
import seaborn as sns; sns.set()


# In[41]:


reload(nf)


# In[42]:


# code pre each protein
# 1) NP_004574 (down regulated)
# upload the data and make a proper DataFrame
data1 = pd.read_csv('proteins_1_all.csv',delimiter=",")
genmatrix1= DataFrame(data1)
genode1=DataFrame(data1,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])


# In[43]:


# create a tulpe from a list to use it as nodes in network construction
my_list1 = list(genmatrix1["Official Symbol Interactor A"].values)
my_list2 = list(genmatrix1["Official Symbol Interactor B"].values)
my_list_1=my_list1+my_list2
my_tuple_1=tuple(my_list_1)


# In[44]:


# calculate wight based on # of interactions
weight_node_1=Counter(my_list_1)


# In[45]:


# create a list with edges
list1=list(weight_node_1.keys())
list2=list(weight_node_1.values())


# In[46]:


# create a graph
G1=nx.Graph()


# In[47]:


# add nodes
G1.add_node(my_tuple_1)


# In[48]:


# add edges
G1.add_edge('NP_004574', list1[0], weight= list2[0])

for value in range(0,len(list1)):
   variable = G1.add_edge('NP_004574', list1[value], weight= list2[value])


# In[49]:


# add additional nodes to build a whole network (each step performed as above)
# 2) EAW53700 (down regulated)
data2 = pd.read_csv('proteins_2_all.csv')
genmatrix2= DataFrame(data2)
genode2=DataFrame(data2,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list3 = list(genmatrix2["Official Symbol Interactor A"].values)
my_list4 = list(genmatrix2["Official Symbol Interactor B"].values)

my_list_2=my_list3+my_list4
my_tuple_2=tuple(my_list_2)

from collections import Counter
weight_node_2=Counter(my_list_2)

list3=list(weight_node_2.keys())
list4=list(weight_node_2.values())

import networkx as nx
G2=nx.Graph()
G2.add_node(my_tuple_2)

G2.add_edge('EAW53700', list3[0], weight= list4[0])

for value in range(0,len(list3)):
   variable = G2.add_edge('EAW53700', list3[value], weight= list4[value])


# In[50]:


# create a network netween two molecus and color the nodes with the degree >1 in blue
# create a graph
U_two_nodes=nx.Graph()
U_two_nodes.add_edges_from(list(G1.edges())+list(G2.edges()))
U_two_nodes.add_nodes_from(list(G1.nodes())+list(G2.nodes()))

# find the degree
info=dict(U_two_nodes.degree(nbunch=None, weight=None))
#DataFrame(info, index=[0])

# select overlap i.e. nodes with degree >1
overlap_nodes=list(k for k, v in info.items() if v > 1)

# build two molecules network, and select nodes with the degree <1
single_nodes=list(k for k, v in info.items() if v <= 1)

#create a graph
fig = plt.figure()
nx.draw(U_two_nodes, with_labels=True, nodelist=overlap_nodes, node_color='b', font_size=6)
plt.show()


# In[51]:


#  diffrent visualisation of the two molecule network
fig = plt.figure()
layout = nx.spring_layout(U_two_nodes,iterations=100)
nx.draw_networkx_nodes(U_two_nodes, layout, nodelist=overlap_nodes, node_color='b', node_size=len(overlap_nodes))
nx.draw_networkx_nodes(U_two_nodes, layout, nodelist=single_nodes, node_color='r', node_size=len(single_nodes))
nx.draw_networkx_edges(U_two_nodes, layout, width=10.5, dge_color="ccccc")
plt.show()


# In[52]:


#3) CAG46507 (down regulated)
data3 = pd.read_csv('proteins_3_all.csv',delimiter=",")
genmatrix3= DataFrame(data3)
genode3=DataFrame(data3,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list5 = list(genmatrix3["Official Symbol Interactor A"].values)
my_list6 = list(genmatrix3["Official Symbol Interactor B"].values)
my_list_3=my_list5+my_list6
my_tuple_3=tuple(my_list_3)

from collections import Counter
weight_node_3=Counter(my_list_3)
weight_node_3

list5=list(weight_node_3.keys())
list6=list(weight_node_3.values())

import networkx as nx
G3=nx.Graph()
G3.add_node(my_tuple_3)

G3.add_edge('CAG46507', list5[0], weight= list6[0])

for value in range(0,len(list5)):
   variable = G3.add_edge('CAG46507', list5[value], weight= list6[value])



# In[53]:


# 4) EAW98517 (down regulated)
data4 = pd.read_csv('proteins_4_all.csv',delimiter=",")
genmatrix4= DataFrame(data4)
genode4=DataFrame(data4,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list7 = list(genmatrix4["Official Symbol Interactor A"].values)
my_list8 = list(genmatrix4["Official Symbol Interactor B"].values)
my_list_4=my_list7+my_list8
my_tuple_4=tuple(my_list_4)

from collections import Counter
weight_node_4=Counter(my_list_4)
weight_node_4

list7=list(weight_node_4.keys())
list8=list(weight_node_4.values())

import networkx as nx
G4=nx.Graph()
G4.add_node(my_tuple_4)
G4.add_edge('EAW98517', list7[0], weight= list8[0])

for value in range(0,len(list7)):
   variable = G4.add_edge('EAW98517', list7[value], weight= list8[value])


# In[54]:


# 5) NP_006588 (up reguated)
data5 = pd.read_csv('proteins_5_all.csv',delimiter=",")
genmatrix5= DataFrame(data5)
genode5=DataFrame(data5,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list9 = list(genmatrix5["Official Symbol Interactor A"].values)
my_list10 = list(genmatrix5["Official Symbol Interactor B"].values)
my_list_5=my_list9+my_list10
my_tuple_5=tuple(my_list_5)

from collections import Counter
weight_node_5=Counter(my_list_5)
weight_node_5

list9=list(weight_node_5.keys())
list10=list(weight_node_5.values())

import networkx as nx
G5=nx.Graph()
G5.add_node(my_tuple_5)

G5.add_edge('NP_006588', list9[0], weight= list10[0])

for value in range(0,len(list9)):
   variable = G5.add_edge('NP_006588', list9[value], weight= list10[value])



# In[55]:


# 6) ABB01006 (up reguated)
data6 = pd.read_csv('proteins_6_all.csv',delimiter=",")
genmatrix6= DataFrame(data6)
genode6=DataFrame(data6,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list11 = list(genmatrix6["Official Symbol Interactor A"].values)
my_list12 = list(genmatrix6["Official Symbol Interactor B"].values)
my_list_6=my_list11+my_list12
my_tuple_6=tuple(my_list_6)

from collections import Counter
weight_node_6=Counter(my_list_6)

list11=list(weight_node_6.keys())
list12=list(weight_node_6.values())

import networkx as nx
G6=nx.Graph()
G6.add_node(my_tuple_6)

G6.add_edge('ABB01006', list11[0], weight= list12[0])

for value in range(0,len(list11)):
   variable = G6.add_edge('ABB01006', list11[value], weight= list12[value])


# In[56]:


# 7) CAI64497 (up reguated)
data7 = pd.read_csv('proteins_7_all.csv',delimiter=",")
genmatrix7= DataFrame(data7)
genode7=DataFrame(data7,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list13 = list(genmatrix7["Official Symbol Interactor A"].values)
my_list14 = list(genmatrix7["Official Symbol Interactor B"].values)
my_list_7=my_list13+my_list14
my_tuple_7=tuple(my_list_7)

from collections import Counter
weight_node_7=Counter(my_list_7)

list13=list(weight_node_7.keys())
list14=list(weight_node_7.values())

import networkx as nx
G7=nx.Graph()
G7.add_node(my_tuple_7)

G7.add_edge('CAI64497', list13[0], weight= list14[0])

for value in range(0,len(list13)):
   variable = G7.add_edge('CAI64497', list13[value], weight= list14[value])


# In[57]:


# 8) NP_005339 (up reguated)

data8 = pd.read_csv('proteins_8_all.csv',delimiter=",")
genmatrix8= DataFrame(data8)
genode8=DataFrame(data8,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list13 = list(genmatrix8["Official Symbol Interactor A"].values)
my_list14 = list(genmatrix8["Official Symbol Interactor B"].values)
my_list_8=my_list13+my_list14
my_tuple_8=tuple(my_list_8)

from collections import Counter
weight_node_8=Counter(my_list_8)

list13=list(weight_node_8.keys())
list14=list(weight_node_8.values())

import networkx as nx
G8=nx.Graph()
G8.add_node(my_tuple_8)

G8.add_edge('NP_005339', list13[0], weight= list14[0])

for value in range(0,len(list13)):
   variable = G8.add_edge('NP_005339', list13[value], weight= list14[value])


# In[58]:


# 9) P43358 (up reguated)

data9 = pd.read_csv('proteins_9_all.csv',delimiter=",")
genmatrix9= DataFrame(data9)
genode9=DataFrame(data9,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list15 = list(genmatrix9["Official Symbol Interactor A"].values)
my_list16 = list(genmatrix9["Official Symbol Interactor B"].values)
my_list_9=my_list15+my_list16
my_tuple_9=tuple(my_list_9)

from collections import Counter
weight_node_9=Counter(my_list_9)

list15=list(weight_node_9.keys())
list16=list(weight_node_9.values())

import networkx as nx
G9=nx.Graph()
G9.add_node(my_tuple_9)

G9.add_edge('P43358', list15[0], weight= list16[0])

for value in range(0,len(list15)):
   variable = G9.add_edge('P43358', list15[value], weight= list16[value])


# In[59]:


# 10) NP_037473 (up reguated)

data10 = pd.read_csv('proteins_10_all.csv',delimiter=",")
genmatrix10= DataFrame(data10)
genode10=DataFrame(data10,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list17 = list(genmatrix10["Official Symbol Interactor A"].values)
my_list18 = list(genmatrix10["Official Symbol Interactor B"].values)
my_list_10=my_list17+my_list18
my_tuple_10=tuple(my_list_10)

from collections import Counter
weight_node_10=Counter(my_list_10)

list17=list(weight_node_10.keys())
list18=list(weight_node_10.values())

import networkx as nx
G10=nx.Graph()
G10.add_node(my_tuple_10)

G10.add_edge('NP_037473', list17[0], weight= list18[0])

for value in range(0,len(list17)):
   variable = G10.add_edge('NP_037473', list17[value], weight= list18[value])



# In[60]:


# 11) EAW86495 (up reguated)

data11 = pd.read_csv('proteins_11_all.csv',delimiter=",")
genmatrix11= DataFrame(data11)
genode11=DataFrame(data11,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list19 = list(genmatrix11["Official Symbol Interactor A"].values)
my_list20 = list(genmatrix11["Official Symbol Interactor B"].values)
my_list_11=my_list19+my_list20
my_tuple_11=tuple(my_list_11)

from collections import Counter
weight_node_11=Counter(my_list_11)

list19=list(weight_node_11.keys())
list20=list(weight_node_11.values())

import networkx as nx
G11=nx.Graph()
G11.add_node(my_tuple_11)

G11.add_edge('EAW86495', list19[0], weight= list20[0])

for value in range(0,len(list19)):
   variable = G11.add_edge('EAW86495', list19[value], weight= list20[value])


# In[61]:


# 12) AAH51814 (up reguated)

data12 = pd.read_csv('proteins_12_all.csv',delimiter=",")
genmatrix12= DataFrame(data12)
genode12=DataFrame(data11,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list21 = list(genmatrix12["Official Symbol Interactor A"].values)
my_list22 = list(genmatrix12["Official Symbol Interactor B"].values)
my_list_12=my_list21+my_list22
my_tuple_12=tuple(my_list_12)

from collections import Counter
weight_node_12=Counter(my_list_12)

list21=list(weight_node_12.keys())
list22=list(weight_node_12.values())

import networkx as nx
G12=nx.Graph()
G12.add_node(my_tuple_12)

G12.add_edge('AAH51814', list21[0], weight= list22[0])

for value in range(0,len(list21)):
   variable = G12.add_edge('AAH51814', list21[value], weight= list22[value])


# In[62]:


# 13) NP_001737 (up reguated)

data13 = pd.read_csv('proteins_13_all.csv',delimiter=",")
genmatrix13= DataFrame(data13)
genode13=DataFrame(data13,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list23 = list(genmatrix13["Official Symbol Interactor A"].values)
my_list24 = list(genmatrix13["Official Symbol Interactor B"].values)
my_list_13=my_list23+my_list24
my_tuple_13=tuple(my_list_13)

from collections import Counter
weight_node_13=Counter(my_list_13)

list21=list(weight_node_13.keys())
list22=list(weight_node_13.values())

import networkx as nx
G13=nx.Graph()
G13.add_node(my_tuple_13)

G13.add_edge('NP_001737', list21[0], weight= list22[0])

for value in range(0,len(list21)):
   variable = G13.add_edge('NP_001737', list21[value], weight= list22[value])


# In[63]:


# 14) AAW67757 (up reguated)

data14 = pd.read_csv('proteins_14_all.csv',delimiter=",")
genmatrix14= DataFrame(data14)
genode14=DataFrame(data14,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list25 = list(genmatrix13["Official Symbol Interactor A"].values)
my_list26 = list(genmatrix13["Official Symbol Interactor B"].values)
my_list_14=my_list25+my_list26
my_tuple_14=tuple(my_list_14)

from collections import Counter
weight_node_14=Counter(my_list_14)

list25=list(weight_node_14.keys())
list26=list(weight_node_14.values())

import networkx as nx
G14=nx.Graph()
G14.add_node(my_tuple_14)

G14.add_edge('AAW67757', list25[0], weight= list26[0])

for value in range(0,len(list25)):
   variable = G14.add_edge('AAW67757', list25[value], weight= list26[value])


# In[64]:


# 15) NP_006380 (up reguated)

data15 = pd.read_csv('proteins_15_all.csv',delimiter=",")
genmatrix15= DataFrame(data15)
genode15=DataFrame(data15,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list27 = list(genmatrix15["Official Symbol Interactor A"].values)
my_list28 = list(genmatrix15["Official Symbol Interactor B"].values)
my_list_15=my_list27+my_list28
my_tuple_15=tuple(my_list_15)

from collections import Counter
weight_node_15=Counter(my_list_15)

list27=list(weight_node_15.keys())
list28=list(weight_node_15.values())

import networkx as nx
G15=nx.Graph()
G15.add_node(my_tuple_15)

G15.add_edge('NP_006380', list27[0], weight= list28[0])

for value in range(0,len(list27)):
   variable = G15.add_edge('NP_006380', list27[value], weight= list28[value])


# In[65]:


# 16) BAD93042 (up reguated)

data16 = pd.read_csv('proteins_16_all.csv',delimiter=",")
genmatrix16= DataFrame(data16)
genode16=DataFrame(data16,columns=['Official Symbol Interactor A', 'Official Symbol Interactor B'])

my_list29 = list(genmatrix16["Official Symbol Interactor A"].values)
my_list30 = list(genmatrix16["Official Symbol Interactor B"].values)
my_list_16=my_list29+my_list30
my_tuple_16=tuple(my_list_16)

from collections import Counter
weight_node_16=Counter(my_list_16)

list29=list(weight_node_16.keys())
list30=list(weight_node_16.values())

import networkx as nx
G16=nx.Graph()
G16.add_node(my_tuple_16)

G16.add_edge('BAD93042', list29[0], weight= list30[0])

for value in range(0,len(list29)):
   variable = G16.add_edge('BAD93042', list29[value], weight= list30[value])


# Network with all 16 proteins 

# In[66]:


# Final graph
# all 16 proteins together (diffrent way)
U=nx.Graph()
U.add_edges_from(list(G1.edges())+list(G2.edges())+list(G3.edges())+list(G4.edges())+list(G5.edges())+list(G6.edges())+list(G7.edges())+list(G8.edges())+list(G9.edges())+list(G10.edges())+list(G11.edges())+list(G12.edges())+list(G13.edges())+list(G14.edges())+list(G15.edges())+list(G16.edges()))
# U.add_nodes_from(U.nodes()+g.nodes())
nx.draw(U)
plt.show()


# In[67]:


# . save it as U_whole
U_whole=U.copy()


# In[68]:


U_whole.nodes()


# Descriptive Statisctics

# In[33]:


# no of nodes
list_nodes = pd.DataFrame(list(U_whole.nodes()))
#list_nodes.to_csv('list_edges.csv')
len(list_nodes)


# In[34]:


# no. of edges
list_edges = pd.DataFrame(list(U_whole.edges()))
list_edges.to_csv('list_edges.csv')
#len(list_edges)
list_edges


# In[35]:


# sort nodes for the higest degree
a1=dict(nx.degree(U_whole))
a1_sorted= sorted(a1,key=a1.get, reverse=True)

for r in a1_sorted:
     print(r,a1[r])


# In[38]:


# sort nodes for the higest betweeness centrality
a2=dict(nx.betweenness_centrality(U_whole))
a2_sorted= sorted(a2,key=a2.get, reverse=True)

for r in a2_sorted:
     print(r,a2[r])


# Degree Trimming

# Trimming with the degree higher than 6

# In[39]:


# remove nodes with the degree lower than 6
outdeg = dict(nx.degree(U))
to_remove = [n for n in outdeg if outdeg[n] <= 6]
#create a new variable
U_degree6=U.remove_nodes_from(to_remove)
# keep nodes with the degree higher than 6 ( looks similar like the trim network)
to_keep = [n for n in outdeg if outdeg[n] > 6]
U_degree6=U.subgraph(to_keep)


# In[40]:


#X odes of degree > 6 and all organismus is the same
nx.draw(U_degree6,with_labels=True,font_size=6)
plt.show()


# In[41]:


# count nb of nodes when degree higher than 6
nx.number_of_nodes(U_degree6)


# In[42]:


# what are its nodes
nx.nodes(U_degree6)


# In[43]:


# count nb of edges when degree higher than 6
nx.number_of_edges(U_degree6)


# In[46]:


# sort for the higest degree 
a3=dict(nx.degree(U_degree6))
a3_sorted= sorted(a3,key=a3.get, reverse=True)

for r in a3_sorted:
     print r,a3[r]


# In[47]:


# sort for the higest betweenness centrality
a4=nx.betweenness_centrality(U_degree6)
a4_sorted= sorted(a4,key=a4.get, reverse=True)

for r in a4_sorted:
     print r,a4[r]


# Find Minimum spanning tree

# In[48]:


#print(networkx.__version__)
X= nx.to_numpy_matrix(U_whole)
mst = minimum_spanning_tree(X)
mst_new= nx.from_scipy_sparse_matrix(mst)
mst_edges= nx.generate_edgelist(mst_new)
temp_tup=tuple(nx.nodes(U_whole))
temp_U = nx.parse_edgelist(mst_edges)
Z=nx.Graph()
Z.add_node(temp_tup)
edges = nx.edges(temp_U)
Z.add_edges_from(edges)

#print(list(mst_edges))
plt.axis('off')
nx.draw(Z, node_size=4)
plt.show()


# 1a. In a degree trimmed network look for one identifier

# In[49]:


# response to stress in the degree trimmed network color in blue
color_map = []
for node in U_degree6:
    if node == "NP_005339":
        color_map.append('blue')# response to stress
    elif node == "BAD93042":
        color_map.append('blue')
    elif node == "CAI64497":
        color_map.append('blue')
    else: color_map.append('red') # others    
nx.draw(U_degree6,node_color = color_map,with_labels = True)
plt.show()


# Shortest path analysis among these three proteins with response to stress identifier

# In[50]:


# shortest path between two nodes with response to stress
nodes= nx.nodes(U_degree6)
nx.shortest_path(U_degree6,source='BAD93042',target="NP_005339")


# In[51]:


# shortest path between two nodes with response to stress
nodes= nx.nodes(U_degree6)
nx.shortest_path(U,source='BAD93042',target="CAI64497")


# In[52]:


# shortest path between two nodes with response to stress
nodes= nx.nodes(U)
nx.shortest_path(U,source='NP_005339',target="CAI64497")


# Finding the meaning

# 1b. Creating subgraphs

# Meaning in the network

# In[67]:


# two identifiers function
nf.choose_identifier('ATP binding','response to stress',path_dic_identifier,U_degree6)


# In[68]:


# select one indentifier in the MST
nf.choose_identifier_MST('response to stress',path_dic_identifier,U)


# In[61]:


# select one identiofier in the whole function
nf.choose_identifier_whole('response to stress',path_dic_identifier,U_whole)


# The knowlege graph expansion

# In[69]:


# import data
import pandas as pd
df_data= pd.read_excel('Finalised-knowledge-graph (3).xlsx')

# delete NAN
df_data1=df_data.dropna()

# create a dictionary for the knowledge expansion via set operation
dic={}
dic2={}
for i in list(df_data1.ID.unique()):
    df_per_protein= df_data1[df_data1['ID']==i] # search for specific protein
    df_grouped=pd.DataFrame(df_per_protein.groupby("Species")['function'].apply(list)) # create a group of functions per species in a dataframe
    if 'Arabidopsis thaliana' in df_grouped.index:
        set_plant=set(df_grouped['function']['Arabidopsis thaliana'])
    else:
        set_plant=set()
    
    if 'Mus musculus' in df_grouped.index:
        set_mouse=set(df_grouped['function']['Mus musculus'])
    else:
        set_mouse=set()
        
    if 'Rattus norvegicus' in df_grouped.index:
        set_rat=set(df_grouped['function']['Rattus norvegicus'])
    else:
        set_rat=set()
        
    if 'Homo Sapiens' in df_grouped.index:
        set_human=set(df_grouped['function']['Homo Sapiens'])
    else:
        set_human=set()
    
    intersection_human_rat=set_human.intersection(set_rat)
    
    intersection_human_plant=set_human.intersection(set_plant)
    
    difference_mouse=list(set_mouse.difference(set_human))
    
    difference_rat=list(set_rat.difference(set_human))
    
    difference_plant=list(set_plant.difference(set_human))
    
    dic[i] = list(set_human),difference_mouse,difference_rat,difference_plant
    dic2[i] = len(list(set_human)),len(difference_mouse),len(difference_rat),len(difference_plant)
new_df = pd.DataFrame.from_dict(dic)
new_df = new_df.transpose()
new_df.columns=['Function_Human','Function_Mouse','Function_Rat','Function_Plant']
new_df.head()

new_df2 = pd.DataFrame.from_dict(dic2)
new_df2 = new_df2.transpose()
new_df2.columns=['Function_Human','Function_Mouse','Function_Rat','Function_Plant']
new_df2


# In[70]:


# save new_df as csv file
writer = ExcelWriter('Features.xlsx')
new_df.to_excel(writer,'protein_name')
writer.save()


# In[71]:


# heat map
sns.heatmap(new_df2, cmap="Blues")
plt.figure(figsize=(10,10))
plt.show()


# Validation

# In[72]:


# import data
import pandas as pd
df_data_VAL= pd.read_excel('functions_data.xls', sheetname='knowledge_data-2')

# create a list of unique proteins
proteins=list(df_data_VAL.protein.unique())

# group fatures per protein
df_grouped=pd.DataFrame(df_data_VAL.groupby('protein')['identifier'].apply(list))

# count the amount of features per protein
df_len=pd.DataFrame(df_data_VAL.groupby('protein')['identifier'].apply(len))

#create a heat map
sns.heatmap(df_len, cmap="Blues")
plt.figure(figsize=(100,100))
plt.show()


# In[73]:


df_len

