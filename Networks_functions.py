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
import matplotlib.patches as mpatches

# function 0: creating a dictionary of identifiers

def first_function(attribute, dictionary):
    list_identifier=[]
    for i in dictionary:
        if attribute in dictionary[i]:
            list_identifier.append(i)
    return(list_identifier)


# function 1: create the subnetwork with the one attribiute, and color the second attribute  within

def choose_identifier(attribute_1,attribute_2,path_dic_identifier,U_degree6,):

    search_attribute_1=first_function(attribute_1,path_dic_identifier)
    sub_1=U_degree6.subgraph(search_attribute_1)
    path_dic_identifier_new = {k: path_dic_identifier[k] for k in search_attribute_1}
    search_attribute_2=first_function(attribute_2,path_dic_identifier_new)
    sub_2=U_degree6.subgraph(search_attribute_2)
    sub=nx.nodes(sub_1)+nx.nodes(sub_2)
    sub_new=U_degree6.subgraph(sub)

    labels={}
    for i in nx.nodes(sub_new):
        labels[i] = i

    plt.figure()

    pos = nx.spring_layout(sub_new)

    nx.draw_networkx_nodes(sub_new, pos, nodelist=search_attribute_1, node_color='r')
    nx.draw_networkx_nodes(sub_new, pos, search_attribute_2, node_color='b')
    nx.draw_networkx_labels(sub_new, pos, labels, font_size=10)
    nx.draw_networkx_edges(sub_new, pos, edgelist=nx.edges(sub_new))
    plt.legend((str(attribute_1),str(attribute_2)))
    plt.axis('off')
    plt.show()
  
    degree_nodes=nx.degree(sub_new)

    degree_nodes_sorted= sorted(degree_nodes,key=degree_nodes.get, reverse=True)

    print ('degree:')
    
    print("\n")
    
    for r in degree_nodes_sorted:
        print(r,degree_nodes[r])
            
    centrality_nodes=nx.betweenness_centrality(sub_new)

    centrality_nodes_sorted= sorted(centrality_nodes,key=centrality_nodes.get, reverse=True)

    print("\n")
    
    print ('betweeness centrality:')
    
    print("\n")
    
    for r in centrality_nodes_sorted:
        
        print(r,centrality_nodes[r])

    print("\n")
    
    print ('graphical representation:')
            
    plt.figure()

    degree_nodes=nx.degree(sub_new)
    centrality_nodes=nx.betweenness_centrality(sub_new)
    x=list(degree_nodes.keys())
    z=list(centrality_nodes.values())
    y=list(degree_nodes.values())
    y_sort = sorted(y, key=int, reverse=False)
    fig, ax = plt.subplots(figsize=(20, 10)) 
    #ax1=plt.subplots()
    #ax2=ax1.twinx()
    ax.plot(x,y_sort, c='b')
    ax1=ax.twinx()
    #df.sort()
    ax1.plot(x,z, c='r')
    red_patch = mpatches.Patch(color='red', label='betweeness centrality')
    blue_patch = mpatches.Patch(color='blue', label='degree')
    plt.legend(handles=[red_patch, blue_patch])

    return(plt.show())

# function 2: create the subnetwork of a MST coloring 1 function

def choose_identifier_MST(attribute_1,path_dic_identifier,U):

    search_attribute_1=first_function(attribute_1,path_dic_identifier)
    sub_1=U.subgraph(search_attribute_1)

    labels={}
    for i in nx.nodes(sub_1):
        labels[i] = i

    plt.figure()

    pos = nx.spring_layout(U)

    nx.draw_networkx_nodes(U, pos, nodelist=search_attribute_1, node_color='b')
    nx.draw_networkx_nodes(U, pos, nodelist=list(set(nx.nodes(U))-set(search_attribute_1)), node_color='r')
    nx.draw_networkx_labels(U, pos, labels, font_size=10)
    nx.draw_networkx_edges(U, pos, edgelist=nx.edges(U))
    plt.axis('off')
    plt.legend((str(attribute_1),'others'))
  
    print("degree of",attribute_1,":")
    
    print("\n")
    
    degree_nodes=nx.degree(U)

    degree={}
    for i in labels.keys():
        degree[i]=degree_nodes[i]

    degree_nodes_sorted= sorted(degree,key=degree.get, reverse=True)

    for r in degree_nodes_sorted:
        print(r,degree[r])
            
    print("\n")
            
    print("betweeness centrality of",attribute_1,":")
    
    print("\n")
            
    centrality_nodes=nx.betweenness_centrality(U)
    
    centrality={}
    for i in labels.keys():
        centrality[i]=centrality_nodes[i]

    centrality_nodes_sorted= sorted(centrality,key=centrality.get, reverse=True)

    for r in centrality_nodes_sorted:
        print(r,centrality[r])

    print("\n")
    
    print ('graphical representation:')
            
    plt.figure()

    degree_nodes=nx.degree(U)
    centrality_nodes=nx.betweenness_centrality(U)
    x=list(degree_nodes.keys())
    z=list(centrality_nodes.values())
    y=list(degree_nodes.values())
    y_sort = sorted(y, key=int, reverse=False)
    fig, ax = plt.subplots(figsize=(30, 15)) 
    #ax1=plt.subplots()
    #ax2=ax1.twinx()
    ax.plot(x,y_sort, c='b')
    ax1=ax.twinx()
    #df.sort()
    ax1.plot(x,z, c='r')
    red_patch = mpatches.Patch(color='red', label='betweeness centrality')
    blue_patch = mpatches.Patch(color='blue', label='degree')
    plt.legend(handles=[red_patch, blue_patch])

    return(plt.show())

#function 3: Draw the whole graph + degree + centrality

def choose_identifier_whole(attribute_1,path_dic_identifier,U_whole):

    search_attribute_1=first_function(attribute_1,path_dic_identifier)
    sub_1=U_whole.subgraph(search_attribute_1)

    labels={}
    for i in nx.nodes(sub_1):
        labels[i] = i

    plt.figure()

    pos = nx.spring_layout(U_whole)

    nx.draw_networkx_nodes(U_whole, pos, nodelist=search_attribute_1, node_color='b')
    nx.draw_networkx_nodes(U_whole, pos, nodelist=list(set(nx.nodes(U_whole))-set(search_attribute_1)), node_color='r')
    nx.draw_networkx_edges(U_whole, pos, edgelist=nx.edges(U_whole))
    plt.axis('off')
    plt.legend((str(attribute_1),'others'))
    plt.show()
    
    print("degree of",attribute_1,":")
    
    print("\n")
    
    degree_nodes=nx.degree(U_whole)
    
    degree={}
    for i in labels.keys():
        degree[i]=degree_nodes[i]
    
    degree_nodes_sorted= sorted(degree,key=degree.get, reverse=True)

    for r in degree_nodes_sorted:
    
        print(r,degree[r])
            
    print("\n")
            
    print("betweeness centrality of",attribute_1,":")
    
    print("\n")
            
    centrality_nodes=nx.betweenness_centrality(U_whole)
    
    centrality={}
    for i in labels.keys():
        centrality[i]=centrality_nodes[i]

    centrality_nodes_sorted= sorted(centrality,key=centrality.get, reverse=True)

    for r in centrality_nodes_sorted:
      
        print(r,centrality[r])

    return(plt.show())