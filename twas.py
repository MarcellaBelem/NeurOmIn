import pandas as pd
import streamlit as st
import os
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import plotly.express as px

@st.cache_data
def read_files(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo twas especifico
    data= aqruivo df com as variantes
    Essa função vai ler os aqruivos iniciados com TWASAtlas_ e apresenta-lo
    output: df com todas os genes identificados por TWAS para determinada condição
    '''

    for filename in os.listdir(folder_path):
        if filename.startswith("TWASAtlas"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep=',')

    return df

def show_association(data, folder_path, disease):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo twas especifico
    data= aqruivo df com as variantes e o loci
    Essa função vai ler os aqruivos iniciados com TWASAtlas_ e identificar os genes encontrados no df das varinates e busca-lo no 
    TWAS especifico, alem de apresentar um grafico que mostra a significancia do gene no tecido
    output: df com todas os genes identificados por TWAS para determinada condição, um grafico de comparação entre
    tecido e gene
    '''
    twas=read_files(folder_path)
    st.dataframe(twas)
    st.write(twas.shape)

    data=data[0]
    genes=list(set(data["SYMBOL"]))
    
   
    group_select=st.segmented_control("select:", options=["all","exclusive"], default="exclusive", selection_mode="single",key="select")
    col1,col2=st.columns(2)
    
    if group_select=="exclusive":
    #fazer plot de relação com tecidos
        
        df_long = twas[twas['Gene Symbol'].isin(genes)]
        
        
        tissue_selected = col1.multiselect("Select tissues (optional):", 
                                        sorted(df_long["Tissues"].unique()), default=df_long["Tissues"].unique())

        if tissue_selected:
            df_long = df_long[df_long["Tissues"].isin(tissue_selected)]
            df_long=df_long[["Gene Symbol","Tissues","P-value"]]
            df_long["P-value"] = df_long["P-value"].replace(["'-", "NA", "na", "None", ""], np.nan)
            df_long["P-value"] = pd.to_numeric(df_long["P-value"])

            df_long = df_long.dropna(subset=["P-value"]) 

            df_long['log_10(P)'] = -np.log10(df_long['P-value'])

            st.subheader(f"Integrating with {disease} genes")
            fig=px.scatter(
            df_long,
            x='Tissues',
            y='Gene Symbol',
            size='log_10(P)',         
            color='Tissues',        
            color_continuous_scale='Viridis',
            title=f'Gene Presence Across Tissues in {disease}')

            fig.update_traces(marker=dict(opacity=0.7, line=dict(width=0.5, color='DarkSlateGrey')))
            fig.update_layout(
            xaxis_title='Tissues',
            yaxis_title='Genes',
            plot_bgcolor='white',
            showlegend=False,
            height=500)

            st.plotly_chart(fig)
    elif group_select == "all":
        #preselected_genes = twas[:5,"Gene Symbol"].unique()
        gene_selected = col1.multiselect("Select a gene:", sorted(twas["Gene Symbol"].unique()), 
                                        default=twas["Gene Symbol"].unique()[:5])
        tissue_selected = col2.multiselect("Select tissues (optional):", 
                                        sorted(twas["Tissues"].unique()), default=twas["Tissues"].unique())

        # Filtra os dados dinamicamente
        filtered_df = twas[twas["Gene Symbol"].isin(gene_selected)]
        if tissue_selected:
            filtered_df = filtered_df[filtered_df["Tissues"].isin(tissue_selected)]
            df_long=filtered_df[["Gene Symbol","Tissues","P-value"]]
            df_long["P-value"] = df_long["P-value"].replace(["'-", "NA", "na", "None", ""], np.nan)
            df_long["P-value"] = pd.to_numeric(df_long["P-value"])

            df_long = df_long.dropna(subset=["P-value"]) 

            df_long['log_10(P)'] = -np.log10(df_long['P-value'])

            st.subheader(f"Integrating with {disease} genes")
            fig=px.scatter(
            df_long,
            x='Tissues',
            y='Gene Symbol',
            size='log_10(P)',         
            color='Tissues',        
            color_continuous_scale='Viridis',
            title=f'Gene Presence Across Tissues in {disease}')

            fig.update_traces(marker=dict(opacity=0.7, line=dict(width=0.5, color='DarkSlateGrey')))
            fig.update_layout(
            xaxis_title='Tissues',
            yaxis_title='Genes',
            plot_bgcolor='white',
            showlegend=False,
            height=500)

            st.plotly_chart(fig)

 
        