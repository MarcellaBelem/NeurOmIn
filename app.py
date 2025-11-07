#entrar no ambiente virtual
#app-env\Scripts\activate

import streamlit as st
import variants
import gwas
import twas
import network
import annotate
import poligenic
import proteomic


st.set_page_config(page_title="NeurOmIn", layout="wide")

st.markdown(
    """
    <style>
        /* Fundo da página principal */
        .stApp {
            background-color: #f0f4f8; /* cinza azulado suave */
        }

        /* Fundo e texto da sidebar */
        [data-testid="stSidebar"] {
            background-color: #1e3d59; /* azul escuro */
            color: white;
        }

        /* Ajuste das labels na sidebar */
        [data-testid="stSidebar"] * {
            color: white !important;
        }

        /* Cabeçalho */
        h1, h2, h3, h4 {
            color: #1e3d59;
        }
    </style>
    """,
    unsafe_allow_html=True
)

st.sidebar.image("neuromin_logo.png", use_container_width=True)

st.title("NeurOmIn - Neurodegenerative Omics Integration ")

# Menu lateral
menu = st.sidebar.radio(
    "Menu",
    ["HOME", "Genomics","Transcriptomics","Proteomics", "ABOUT US"]
)
if menu == "HOME":
    st.header("A database for exploring SNVs associated with neurodegenerative disease with integration omics data")

elif menu == "Genomics":
    st.title("Visualização de dados")
    doenças=st.pills("Select the disease",["PD"],selection_mode="single", key="disease", default="PD",label_visibility="visible", width="stretch")
    if doenças == "PD":
        st.header("Parkinson's disease")
        titulos_abas=["Variants", "GWAS", "PRS"]
        abas = st.tabs(titulos_abas)
        # Adicionar conteúdo a cada aba
        with abas[0]:
            st.header('Visualização das variantes')
            st.write('Aqui realizamos a análise descritivas das variantes...')
            variants.read_files(folder_path="data/")
            variants.show_variants(folder_path="data/")
        

        with abas[1]:
            st.header('GWAS')
            st.write('Aqui realizamos a integração com GWAS...')
            data=variants.read_files(folder_path="data/")
            gwas.show_association(data,folder_path="databases/")

        with abas[2]:
            st.header("Poligenic Score Risk")
            st.write()
            data=variants.read_files(folder_path="data/")
            poligenic.show_score(data, folder_path="pgs/")
            adicional=poligenic.read_files_adicionales(folder_path="pgs/")
            poligenic.make_net_ancestry(data,adicional, folder_path="pgs/", disease=doenças)
    
    else:
        st.header("Alzheimer's Disease")

        titulos_abas=["Variants", "GWAS", "PRS"]
        abas = st.tabs(titulos_abas)
        
        # Adicionar conteúdo a cada aba
        with abas[0]:
            st.header('Visualização das variantes')
            st.write('Aqui realizamos a análise descritivas das variantes...')
            #variants.show_variants()
        

        with abas[1]:
            st.header('GWAS')
            st.write('Aqui realizamos a integração com GWAS...')
            #data=variants.read_files()
            #gwas.show_association(data)

        with abas[2]:
            st.header("Poligenic Score Risk")
            st.write()
            #data=variants.read_files()
            #poligenic.show_score(data)
            #adicional=poligenic.read_files_adicionales()
            #poligenic.make_net_ancestry(data,adicional)

elif menu == "Transcriptomics":
    doenças=st.pills("Select disease",["PD"],selection_mode="single", key="disease", default="PD",label_visibility="visible", width="stretch")
    if doenças == "PD":
        st.header("Parkinson's disease")
        titulos_abas=["TWAS", "Networks", "Annotation"]
        abas = st.tabs(titulos_abas)

        with abas[0]:
            st.header('TWAS')
            st.write('Aqui realizamos a integração com TWAS...')
            data=variants.read_files("data/")
            twas.show_association(data, folder_path="databases/", disease=doenças)

        
        with abas[1]:
            st.header("Network Analysis and Visualization")
            
            st.write('Aqui realizamos a visualização das redes...')
            data=variants.read_files(folder_path="data/")
            network.make_net_tf(data,disease=doenças)
            network.make_net_mirna(data,disease=doenças)
            #g= network.make_graph(data)
            #network.plot_network(g)
            #network.rede(g)

        with abas[2]:
            st.header("Annotation")
            
            st.write('Aqui realizamos a anotação dos genes identificados...')
            data=variants.read_files(folder_path="data/")
            results=annotate.enrichiment(data)
            annotate.plot_enrichiment(results)
    elif doenças == "AD":
        st.header("Alzheimer's disease")
        titulos_abas=["TWAS", "Networks", "Annotation"]
        abas = st.tabs(titulos_abas)

        with abas[0]:
            st.header('TWAS')
            st.write('Aqui realizamos a integração com TWAS...')
            data=variants.read_files("../../../../OneDrive/bio e/parkison/")
            twas.show_association(data, folder_path="../../../../OneDrive/bio e/parkison/", disease=doenças)

        
        with abas[1]:
            st.header("Network Analysis and Visualization")
            
            st.write('Aqui realizamos a visualização das redes...')
            data=variants.read_files(folder_path="../../../../OneDrive/bio e/parkison/")
            network.make_net_tf(data,disease=doenças)
            network.make_net_mirna(data,disease=doenças)
           

        with abas[2]:
            st.header("Annotation")
            
            st.write('Aqui realizamos a anotação dos genes identificados...')
            data=variants.read_files(folder_path="../../../../OneDrive/bio e/parkison/")
            results=annotate.enrichiment(data)
            annotate.plot_enrichiment(results)

elif menu == "Proteomics":
    titulos_abas=["Alphafold"]
    abas = st.tabs(titulos_abas)

    with abas[0]:
        st.header('Alphafold')
        st.write('Aqui realizamos a integração com Alphafold e fazemos a visualização da proteina e sua mutação...')
        data=variants.read_files(folder_path="data/")
        dado=data[1]
        dado["Uniprot_id"] = dado["UNIPROT_ISOFORM"].str.split("-").str.get(0)
        dado=dado[dado["SYMBOL"]!= "-"]
        dado["Feature"] = dado["Feature"].str.split(".").str.get(0)# separar a isoforma

        uniprots = {row["SYMBOL"]: [row["Uniprot_id"], row["Feature"]] for _, row in dado.iterrows()}
        uniprots['VPS35'][0] = 'Q96QK1'
        uniprots['LRRK2'][0] = 'Q55007'
        uniprots['PARK7'][0] = 'Q99497'
   
        #uniprots={"SNCA":["P37840-3","ENST00000394991"],"VPS35":["Q96QK1","ENST00000299138"]}  
        uniprot_id=st.selectbox("choose a gene:", uniprots, index=None)
        if uniprot_id:
            pdb_id=proteomic.proteina(uniprot_id, uniprots)
        # proteomic.find_aa(data)
            proteomic.show_variant(data,pdb_id, uniprot_id, uniprots)

elif menu== "DY":
    titulos_abas=["Upload", "Analysis"]
    abas = st.tabs(titulos_abas)
    with abas[0]:
        st.header("Upload")
        
        st.write('Aqui realizamos a visualização das redes...')


    with abas[1]:
        st.header("Analysis")
        
        st.write('Aqui realizamos a anotação dos genes identificados...')

elif menu == "ABOUT US":
    st.header("Laboratorio de Bioinformatica e ciência de Dados (Laboratory of Bioinformatics and Data Science - LBCD), intituto de ciencias biologicas, UFPA")
    st.image("lbcd_logo.jpeg")