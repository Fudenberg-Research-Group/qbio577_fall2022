# useful libraries to import

import pandas as pd
import numpy as np
import  sklearn.decomposition
import matplotlib.pyplot as plt
import seaborn as sns 

def plot_pca( pca , 
             bigwig_metadata=None,
             metadata_label_column=None, 
             alpha=0.5, 
             lw=0, 
             figsize=(8,8)):
    
    """ 
    Skeleton for plotting PCA and annotating the plot. 
    Can be modified/extended to answer various questions.
    """
    
    plt.rcParams["font.weight"] = "bold"
    plt.rcParams["axes.labelweight"] = "bold"
    if metadata_label_column is not None:
        if bigwig_metadata is None: 
            raise ValueError("must provide metadata table to label by a metadata column") 
        labels = [bigwig_metadata.query(
                    "`File accession`==@ file_accession ").loc[:,metadata_label_column].values[0]
                  for file_accession in pca.feature_names_in_]
        le = sklearn.preprocessing.LabelEncoder()
        le.fit(labels)
        labels_t = le.transform(labels)
        #print(np.unique(labels))
    else: 
        labels_t = None
     
    fig, ax = plt.subplots(figsize=figsize)
    sns.scatterplot(pca.components_[0],
                pca.components_[1],
                hue = labels_t,
                alpha=alpha,
                lw=lw, 
                palette='colorblind')
    ax.xaxis.set_tick_params(labelsize=10)
    if metadata_label_column is not None:
        legend_labels, _= ax.get_legend_handles_labels()
        ax.legend(legend_labels, np.unique(labels), bbox_to_anchor=(1,1))
        sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1))

    #kwargs={'fontsize':'4'}
    #if labels is not None:
    #    for i in range(0, len(pca.components_[0])):
    #        ax.text(pca.components_[0][i], pca.components_[1][i], f'{labels[i]}', **kwargs)

