from sklearn import tree
from matplotlib import pyplot as plt
import graphviz
import numpy as np
from matplotlib.colors import ListedColormap, to_rgb, to_hex
import pygraphviz as pgv
import ast

def plot_clf_regions(X, y, kernelA, kernelB, clf, tick_size=18, colors = ['crimson', 'white', 'lightgreen']):
    
    fig = plt.figure(figsize=[10,8])
    ax = fig.add_subplot(111)
    #print(colors)
    
    ax.plot(X[:,0][y==-1], X[:,1][y==-1], "o", color=colors[0])
    ax.plot(X[:,0][y==0], X[:,1][y==0], "^", color='black')
    ax.plot(X[:,0][y==1], X[:,1][y==1], "x", color=colors[2])
    
    xx, yy = np.meshgrid(np.linspace(0.0, X[:, 0].max()+1, 100), np.linspace(0.0, X[:, 1].max()+1, 100))
    pred = clf.predict(np.c_[(xx.ravel(), yy.ravel())])
    #print(pred)
    
    _colors = colors.copy()
    cl_set = set(pred)
    if -1 not in cl_set:
        _colors.remove(colors[0])
    if 0 not in cl_set:
        _colors.remove(colors[1])
    if 1 not in cl_set:
        _colors.remove(colors[2])
    
    ax.contourf(xx, yy, pred.reshape(xx.shape), cmap=ListedColormap(_colors), alpha=0.25)
    
    ax.set_xlabel("{} (MFLOPs)".format(kernelA), size=18)
    ax.set_ylabel("{} (MFLOPs)".format(kernelB), size=18)
    ax.tick_params(labelsize=18)
    
    return fig

def get_tree_gviz(clf, kernelA, kernelB, colors = ['crimson', 'white', 'lightgreen']):
    
    _cls = sorted(list(clf.classes_.astype(str)))
    _colors = colors.copy()
    if '-1' not in _cls:
        _colors.remove(colors[0])
    if '0' not in _cls:
        _colors.remove(colors[1])
    if '1' not in _cls:
        _colors.remove(colors[2])
    
    dot_data = tree.export_graphviz(clf, out_file=None, 
                                feature_names=[kernelA, kernelB],
                                class_names = _cls,
                                filled=True, impurity=True)
    #g = pydotplus.graph_from_dot_data(dot_data)
    dg = pgv.AGraph(dot_data)
    #print(dg)
    for node in dg.nodes():
        
        bLeaf = False
        
        label_data = node.attr['label'].split("\\n")
        
        if len(label_data) == 4:
            bLeaf = True
        elif len(label_data) == 3:
            print("Decision tree not found. There is only one class")
            return -1
        
        impurity = float(label_data[-4].split("=")[1])
        cls = int(label_data[-1].split("=")[1])
        values = ast.literal_eval(label_data[-2].split("=")[1].strip())
        #print(impurity, cls, values)
        
        # set color
        r, g, b = to_rgb(_colors[np.argmax(values)])
        f = impurity*3/2.
        rgb = (min(f + (1-f)*r, 1.0),\
               min(f + (1-f)*g, 1.0),\
               min(f + (1-f)*b,1.0), 0.5)
        #print(rgb)
        node.attr['fillcolor'] = to_hex(rgb, keep_alpha=True)
        
        new_label = ""
        if not bLeaf:
            new_label += label_data[0] +"\n\n"
        new_label += label_data[-4] + "\n"
        new_label += "num " + label_data[-3] + "\n"
        new_label += "majority " + label_data[-1] 
        
        node.attr['label'] = new_label
        
        #print(impurity, cls, to_hex(rgb))
        
    graph = graphviz.Source(dg.string()) 
    return graph

def gini_impurity(x):
    ig = 0.0
    if len(set(x)) == 1:
        return ig
    else:
        n = float(len(x))
        for i in set(x):
            p = x.count(i)/n
            ig += p*(1-p)
        return ig