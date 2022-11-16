import pandas as pd
import numpy as np
import  sklearn.decomposition
import matplotlib.pyplot as plt
import seaborn as sns 
import colorcet as cc
import math  

def pca_plot_celltype(pca, df, celltype):
    pca_components = pd.DataFrame(pca.components_)
    pca_components.columns = df.columns[3:]
    pca_components_t = pca_components.transpose()
    acc_id_Error = metadata_filter[metadata_filter["Biosample term name"] != celltype]["File accession"]
    df_filtered_H3K9me3 = pca_components_t.drop(pca_components_t[(pca_components_t.index.isin(acc_id_Error))].index, 0)
    labels = [metadata_filter.query(
                        "`File accession`==@ file_accession ").loc[:,"Biosample term name"].values[0] + "_" + metadata_filter.query(
                        "`File accession`==@ file_accession ").loc[:,"Experiment target"].values[0]
                      for file_accession in pca_components_t.index]

    for i in range(0, len(labels)):
        if labels[i] != celltype + "_H3K9me3-human" and labels[i] != celltype + "_H3K36me3-human":
            labels[i] = "others"
    fig, ax = plt.subplots(figsize=(5,5))


    colors = ["#FF9505",  "#40531B", "#DADADA"]
    # Set your custom color palette
    sns.set_palette(sns.color_palette(colors))
    hue_order = [celltype + "_H3K9me3-human", celltype + "_H3K36me3-human", 'others']

    sns.scatterplot(pca_components_t[0],
                    pca_components_t[1],
                    hue = labels,
                    alpha=0.8, hue_order=hue_order)

    legend_labels, _= ax.get_legend_handles_labels()
    ax.legend(legend_labels, np.unique(labels), bbox_to_anchor=(1,1))
    sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1), ncol=math.ceil(3/12))

    plt.title("PC1 vs PC2 \n H3K9me3 and H3K36me3 across " + celltype, fontsize = 20, weight="bold") # title with fontsize 20
    plt.xlabel('PC1 (41.82% varaince expained)', fontsize = 15, weight="bold") # x-axis label with fontsize 15
    plt.ylabel('PC2 (10.5% varaince expained)', fontsize = 15, weight="bold") # y-axis label with fontsize 15