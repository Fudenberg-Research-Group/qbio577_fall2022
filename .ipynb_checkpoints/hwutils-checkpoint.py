# useful libraries to import

import pandas as pd
import numpy as np
import  sklearn.decomposition
import matplotlib.pyplot as plt

def plot_pca( pca , 
             bigwig_metadata=None,
             metadata_label_column=None, 
             alpha=0.5, 
             lw=0, 
             figsize=(8,8)):
    
    """ 
    Visualize PCA results in 2-D plots using PC1 and PC2.
    
    Parameters:
    -----------
    pca: pca result 
    bigwig_metadata: dataframe of bigwig metadata (default: None)
    metadata_label_column: the name of the column you are interested in (default: None)
    alpha: transparency value (default: 0.5)
    lw: line width (default: 0)
    figsize: size of the figure (default: (8,8))
    """
    
    
    if metadata_label_column is not None:
        if bigwig_metadata is None: 
            raise ValueError("must provide metadata table to label by a metadata column") 
        labels = [bigwig_metadata.query(
                    "`File accession`==@ file_accession ").loc[:,metadata_label_column].values[0]
                  for file_accession in pca.feature_names_in_]
        le = sklearn.preprocessing.LabelEncoder()
        le.fit(labels)
        labels = le.transform(labels)
        
        inverse_labels=le.inverse_transform(labels)
        plt.figure(figsize=figsize)
        scatter=plt.scatter(pca.components_[0],
                pca.components_[1],
                c = labels,
                cmap='Spectral',
                alpha=alpha,
                lw=lw)
        classes=np.unique(inverse_labels).tolist()
        plt.legend(handles=scatter.legend_elements(num=len(classes))[0], labels=classes,loc='center left', bbox_to_anchor=(1, 0.5))
        
         
    else: 
        labels = None
        plt.figure(figsize=figsize)
        plt.scatter(pca.components_[0],
                    pca.components_[1],
                    c = labels,
                    alpha=alpha,
                    lw=lw)
  

