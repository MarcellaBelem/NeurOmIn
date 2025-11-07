import pandas as pd
import streamlit as st
import os
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx

@st.cache_data
def read_files(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo pgs especifico
    Essa função vai ler os aqruivos iniciados com PGS e concatena-los
    output: df com todas os PRS identificados para determinada condição
    '''
    all_data=[]

    for filename in os.listdir(folder_path):
        if filename.startswith("PGS"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep='\t', comment="#")
            
            study= os.path.splitext(filename)[0]
            df["PGS_id"]=study

            all_data.append(df)
            data=pd.concat(all_data, ignore_index=True)

    
    return data

def read_files_adicionales(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo pgs especifico
    Essa função vai ler os aqruivos iniciados com pgs e concatena-lo, objetivo obter as metricas e ids das amostras para ver a ancestralidade
    output: df com todas os id_sample identificados 
    '''
    all_data=[]
    for filename in os.listdir(folder_path):
        if filename.startswith("pgs"):
            file_path = os.path.join(folder_path, filename)
            df= pd.read_csv(file_path, sep=',')
            all_data.append(df)

    metrics=all_data[0]
    samples=all_data[1]

    metrics["PSS"]=metrics["PGS Sample Set ID\n(PSS)"].str.split("|").str.get(0)
    data=samples[samples["PGS Sample Set ID\n(PSS)"].isin(metrics["PSS"])].copy()


    data = data.merge(metrics[["PSS", "Evaluated Score"]],
    left_on="PGS Sample Set ID\n(PSS)",
    right_on="PSS", how="left")

    data=data.rename(columns={"Evaluated Score": "PGS_id", "PGS Sample Set ID\n(PSS)":"PGS Sample Set ID (PSS)"})
    data["PGS_id"]=data["PGS_id"].str.split(" ").str.get(0)
    data=data[["PGS_id","PGS Sample Set ID (PSS)","Sample Numbers","Age of Study Participants","Sample Ancestry"]]

    expandir=st.toggle("See infos about the samples", value=False)
    if expandir == True:
        st.subheader("Informações sobre as amostras do PGS")
        st.dataframe(data)
        st.write(data.shape)

    return data


#mostrar um mapa entre ancestralidade e variante
#olhar os scores para as variantes possiveis
def show_score(data, folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo PGS especifico
    data= aqruivo df com as variantes
    Essa função vai ler os aqruivos iniciados com PRS e associar a nossa variantes para ver a intersseção
    output: df com todas as varinates e seus respectivos PRS calculados anteriormente
    '''
    prs=read_files(folder_path)
    st.dataframe(prs)
    st.write(prs.shape)
    data=data[0]

    all=[]
    #preciso separar por rs
    data["rsid_only"] = data["rsid"].str.split(",").str.get(0)
    #preciso pegar 
    
    #filtro
    rs=prs[prs["rsID"].isin(data["rsid_only"])]

    all.append(rs)
     
    st.subheader("Only variants matched from PGS")
    associations= pd.concat(all, ignore_index=True)
    st.dataframe(associations)
    st.write(associations.shape)
    return associations

def make_net_ancestry(data, adicionales, folder_path,disease):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo PGS especifico
    data= aqruivo df com as variantes
    adicionales: é o arquivo com os id_smples e ancestralidade
    disease: str para colocar o nome da doença
    Essa função vai ler os aqruivos iniciados com PRS e associar a nossa variantes para ver a intersseção e 
    gerar uma rede de variantes e a ancestralidade encontrada para ela
    output: rede de interação entre variantes e ancestralidade
    '''
    prs=read_files(folder_path)

    data=data[0] #com minhas variantes

    info=adicionales

    prs=prs.merge(info[["PGS_id", "Sample Ancestry"]],
                  on="PGS_id",
                  how="left")
    
    prs["Ancestry"]=prs["Sample Ancestry"].str.split(",").explode("Ancestry")
    prs.loc[prs["Ancestry"] == " NR", "Ancestry"] = "Not reported"

    rs = list(set(data["rsid_only"]))


    df=prs[prs["rsID"].isin(rs)]
    
    st.subheader("Ancestry-Variants Network")
    col1,col2=st.columns(2)
    nos=col1.segmented_control( "Select one:",["ALL","Exclusive"],selection_mode="single", default="Exclusive")
    ancestrys=prs["Ancestry"].unique().copy()
    ancestry=col2.pills("select ancestry to see:", ancestrys, selection_mode="multi",key="ancestry")

    if nos == "ALL":
        nodes=st.checkbox("Nodes significants (|weight| 0.3)",value=True)
        if nodes:
            df = prs.dropna(subset=["rsID"])
            df=df[df["effect_weight"].abs()> 0.3]
        else:
            df = prs.dropna(subset=["rsID"])
    if ancestry:
        df=df[df["Ancestry"].isin(ancestry)]

    #rede
    net= Network(directed=False,height="600px", width="100%", notebook=False, bgcolor="#FFFFFF", font_color="black", cdn_resources='in_line')
    
    ancestry_colors = {
        "European": "#1f77b4",       # azul
        "African": "#ff7f0e",        # laranja
        "South Asian\n(Indian)": "#2ca02c",     # verde
        "Admixed": "#d62728",        # vermelho
        "Not reported": "#373737",   # gray
        "Other\n(Ashkenazi Jewish)": "#BF6192"
        # adicione outras ancestries se necessário
    }
    for _,row in df.iterrows():
        cor_peso = "blue" if row["effect_weight"] > 0 else "red"
        cor_ancestralidade=ancestry_colors.get(row["Ancestry"], "#7f7f7f")
        net.add_node(row["Ancestry"], label=row["Ancestry"], color=cor_ancestralidade, size=15, shape="square")
        net.add_node(row["rsID"], label=row["rsID"], color="orange", size=10, shape="dot")
        net.add_edge(row["Ancestry"], row["rsID"], value=row["effect_weight"], color =cor_peso, width=0.5, arrows="to")

    col1,col2=st.columns(2)
    repulsao=col1.slider("select distance between nodes?", 100 , 500 , 150, key="r1")
    edger_length=col2.slider("select edger length?", 100 , 500 , 150, key="e1")

    net.repulsion(node_distance=repulsao, spring_length=edger_length)
    net.toggle_physics(True)

    net.set_options("""
    {
    "interaction": {
        "navigationButtons": true,
        "zoomView": true,
        "dragNodes": true,
        "selectable": true,
        "selectConnectedEdges": true
    },
    "physics": {
        "stabilization": true
        }
    }
    """)


    html_content = net.generate_html()
    st.components.v1.html(html_content, height=600, scrolling=True)

    col1.download_button(
    label="⬇️ Download Network (HTML)",
    data=html_content,
    file_name=f"rede_variants_ancestry_{disease}.html",
    mime="text/html")