# useful libraries to import

import pandas as pd
import numpy as np
import  sklearn.decomposition
import matplotlib.pyplot as plt
import seaborn as sns
from more_itertools import locate

def plot_pca( pca , 
             bigwig_metadata=None,
             metadata_label_column=None, 
             alpha=0.5, 
             lw=0, 
             figsize=(10,8),
             label_display=False,
             filter_label_column=False,
             filterby=False):
    
    """ 
    Skeleton for plotting PCA and annotating the plot. 
    Can be modified/extended to answer various questions.
    """
    
    plt.figure(figsize=figsize)

    if metadata_label_column is not None:
        if bigwig_metadata is None: 
            raise ValueError("must provide metadata table to label by a metadata column") 
        labels = [bigwig_metadata.query(
                    "`File accession`==@ file_accession ").loc[:,metadata_label_column].values[0]
                  for file_accession in pca.feature_names_in_]
        le = sklearn.preprocessing.LabelEncoder()
        le.fit(labels)
        labels_transformed = le.transform(labels)

        if filterby is not False:
            labels = [bigwig_metadata.query(
                        "`File accession`==@ file_accession ").loc[:,filter_label_column].values[0]
                        for file_accession in pca.feature_names_in_]
            filter_idx = list(locate(labels, lambda x: x == filterby))
            plt.scatter(np.delete(pca.components_[0], filter_idx, 0),
                        np.delete(pca.components_[1], filter_idx, 0),
                        c = np.delete(labels_transformed, filter_idx, 0),
                        alpha=0.8,
                        lw=0)
        else:
            plt.scatter(pca.components_[0],
                        pca.components_[1],
                        c = labels_transformed,
                        alpha=alpha,
                        lw=lw)
    
    else: 
        labels_transformed = None
        plt.scatter(pca.components_[0],
                    pca.components_[1],
                    alpha=alpha,
                    lw=lw)

    

    plt.xlabel("PC-1 (%f)" % pca.explained_variance_ratio_[0])
    plt.ylabel("PC-2 (%f)" % pca.explained_variance_ratio_[1])


    # display labels for data points on the plot
    if label_display is True:
        for idx, (x,y) in enumerate(zip(pca.components_[0],pca.components_[1])):
            label = labels[idx]
            plt.annotate(label, # this is the text
                        (x,y), # these are the coordinates to position the label
                        textcoords="offset points", # how to position the text
                        xytext=(0,10), # distance from text to points (x,y)
                        ha='center',
                        fontsize=8) # horizontal alignment can be left, right or center
    