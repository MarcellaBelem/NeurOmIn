import pandas as pd
import streamlit as st
import os
import scipy as sp
import matplotlib.pyplot as plt
from pyvis.network import Network
import networkx as nx

@st.cache_data
def read_files(folder_path="databases/"):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo TFlink 
    Essa função vai ler o TFLink
    output: df com todas os genes alvos e seu reguladores (TF)
    '''

    for filename in os.listdir(folder_path):
        if filename.startswith("tflink_database"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep='\t')

    return df

#fazer uma rede com peso da qtd de artigos confirmando interação
def make_net_tf(data,disease):
    '''
    data: arquivo com as varintes e seu loci
    disease: o nome da doença
    Essa função vai associar os genes encontrados nas variantes e formar uma rede de interação de TF-alvo
    output: rede de interção de TF e gene-alvo
    '''
    tf_data=read_files()

    data=data[0]
    genes = list(set(data["SYMBOL"]))

    st.subheader("Transcrition Factor-Gene Network")
    col1=st.columns(1)[0]
    gene_select=col1.multiselect("select genes specifics:", options=genes, key="select1")
    
    if not gene_select:
        df=tf_data[tf_data["Name.Target"].isin(genes)]
    else:
        df=tf_data[tf_data["Name.Target"].isin(gene_select)]
    
    df["weight"]=df["PubmedID"].str.split(";").apply(len)
    csv = df.to_csv(index=False).encode('utf-8')

    nos=st.checkbox("Nodes significants (3 ou more observations)",value=True)
    if nos:
        df=df[df["weight"]>= 3]

    #rede
    net= Network(height="600px", width="100%",directed= True, notebook=False, bgcolor="#FFFFFF", font_color="black", cdn_resources='in_line')
    
    for _,row in df.iterrows():
        net.add_node(row["Name.TF"], label=row["Name.TF"], color="blue", size=15, shape="square")
        net.add_node(row["Name.Target"], label=row["Name.Target"], color="orange", size=10, shape="dot")
        net.add_edge(row["Name.TF"], row["Name.Target"], value=row["weight"], color = "#615456",arrows="to")

    col1,col2=st.columns(2)
    repulsao=col1.slider("select distance between nodes?", 100 , 200 , 150, key="r1")
    edger_length=col2.slider("select edger length?", 100 , 200 , 150, key="e1")

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
    file_name=f"rede_tf_{disease}.html",
    mime="text/html")

    col2.download_button(
    label="⬇️ Download file (csv)",
    data=csv,
    file_name=f"rede_tf_{disease}.csv",
    mime="text/plain")

    


def read_files_mirna(folder_path="databases/"):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo mirTarBase 
    Essa função vai ler o mirTarbase
    output: df com todas os microRNAs humanos ee cm alvos humanos
    '''

    for filename in os.listdir(folder_path):
        if filename.startswith("miRTarBase"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep=',')

    df=df[["miRTarBase ID","miRNA","Species (miRNA)","Target Gene","Species (Target Gene)","Experiments",
           "Support Type","References (PMID)"]]

    df=df[(df["Species (miRNA)"]== "hsa") & (df["Species (Target Gene)"]=="hsa")]


    return df

#fazer uma rede com peso da qtd de artigos confirmando interação
def make_net_mirna(data,disease):
    '''
    data: arquivo com as varintes e seu loci
    disease: o nome da doença
    Essa função vai associar os genes encontrados nas variantes e formar uma rede de interação de miRNA
    output: rede de interção de miRNA e gene-alvo
    '''
    mirna_data=read_files_mirna()

    data=data[0]
    genes = list(set(data["SYMBOL"]))

    st.subheader("miRNA-Gene Network")
    col1=st.columns(1)[0]
    gene_select=col1.multiselect("select genes specifics:", options=genes, key="select2")
    
    if not gene_select:
        df=mirna_data[mirna_data["Target Gene"].isin(genes)]
    else:
        df=mirna_data[mirna_data["Target Gene"].isin(gene_select)]

    
    df["weight"]=df["Experiments"].str.split("//").apply(len)
    df= (
    df.groupby(["miRNA", "Target Gene"], as_index=False)
      .agg({
          "weight": "sum",  # soma o número total de experimentos
          "Experiments": lambda x: "//".join(sorted(set("//".join(x).split("//"))))  # concatena únicos
          #"Target Gene": lambda x: ", ".join(sorted(set(x)))    
      }))
    #df["weight"] = df["Experiments"].str.split("//").apply(lambda x: len(set(x)))
    
    csv = df.to_csv(index=False).encode('utf-8')


    nos=st.checkbox("Nodes significants (2 ou more Experiments)",value=True,  key="nodes_mirna")
    if nos:
        df=df[df["weight"]>= 2]

    #rede
    net= Network(height="600px",directed= True, width="100%", notebook=False, bgcolor="#FFFFFF", font_color="black", cdn_resources='in_line')
    
    for _,row in df.iterrows():
        net.add_node(row["miRNA"], label=row["miRNA"], color="blue", size=15, shape="square")
        net.add_node(row["Target Gene"], label=row["Target Gene"], color="orange", size=10, shape="dot")
        net.add_edge(row["miRNA"], row["Target Gene"], value=row["weight"], color = "#615456", arrows={"to": {"enabled": True, "type": "vee", "scaleFactor": 1.5}})
    

    col1,col2=st.columns(2)
    repulsao=col1.slider("select distance between nodes?", 100 , 200 , 150, key="r2")
    edger_length=col2.slider("select edger length?", 100 , 200 , 150, key="e2")

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
    file_name=f"rede_miRNA_{disease}.html",
    mime="text/html")
   
    col2.download_button(
    label="⬇️ Download file (csv)",
    data=csv,
    file_name=f"rede_miRNA_{disease}.csv",
    mime="text/plain")


def make_graph(data):
    '''
    Esta função recebe um arqivo de texto de interações e forma uma rede
    :param txt: aruivo.txt contendo interações e valores de peso
    :return: retorna um objeto Graph (um grafo não direcionado)
    '''
    tf_data=read_files()

    data=data[0]
    genes = list(set(data["SYMBOL"]))

    df=tf_data[tf_data["Name.Target"].isin(genes)]
    df["weight"]=df["PubmedID"].str.split(";").apply(len)


    nos=st.checkbox("Nodes significants (3 ou more observations)",value=True)
    if nos:
        df=df[df["weight"]>= 3]

    g = nx.DiGraph()
    for _, row in df.iterrows():
        source=row["Name.TF"]
        target=row["Name.Target"]
        mode=row["weight"]

        g.add_node(source, type="TF", color="blue",  node_shape="s")
        g.add_node(target, type="Gene", color="orange",  node_shape="o")
        g.add_edge(source, target, weigh=mode)

    return g


def plot_network(graph):
    '''
    Esta função faz o plot do grafo, colorindo os nodes por logFC e as arestas pelo weight
    :param graph: objeto grafo
    :return: a visualização do grafo
    '''
    # torna-lo direcionado
    #graph = nx.DiGraph(graph)

    node_colors = "red"  # Obter cores dos nós

    edge_colors = "black"#[cor["color"] for edge, weight, cor in graph.edges(data=True)]  ## Obter cores das arestas

    edge_labels = nx.get_edge_attributes(graph, 'weight')

    options = {
        'node_color': node_colors,
        'edge_color': edge_colors,
        'node_size': 100,
        'width': 3,
        'font_size': 10,
        'font_color': "white"}  # config do grafo

    pos = nx.spring_layout(graph)  # layout do graph

    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.title("Rede regulatória", color='black', fontsize=14, pad=15)
    ax.set_title("Rede regulatória centrada em TF", color='black', fontsize=14, pad=15)
    nx.draw(graph, pos,  ax=ax, with_labels=True, node_color=node_colors, node_size=100, edge_color=edge_colors, width=3.0,
            font_size=10, font_color='black', font_weight='bold')  # desenhar o grafo

    nx.draw_networkx_edge_labels(graph, pos, ax=ax, edge_labels=edge_labels, font_size=8,
                                 font_color='purple')  # desenhar as arestas

    # nx.draw(g, options,pos= pos,with_labels=True, font_weight= 'bold') #desenha o graph

    st.pyplot(fig)  # plotar




def rede(graph):
    net = Network(height="600px", width="100%", directed=True, bgcolor="#FDFDFD", font_color="black", cdn_resources='in_line')

    net.from_nx(graph)

    col1,col2=st.columns(2)
    repulsao=col1.slider("select distance between nodes?", 100 , 200 , 150)
    edger_length=col2.slider("select edger length?", 100 , 200 , 150)

    net.repulsion(node_distance=repulsao, spring_length=edger_length)
    net.toggle_physics(True)

    html_content = net.generate_html()
    st.components.v1.html(html_content, height=600, scrolling=True)
        ##add botoes de download, e de zoom, select essas coisas, integrar cytoscape