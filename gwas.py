import pandas as pd
import streamlit as st
import os

@st.cache_data
def read_files(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo gwas especifico
    Essa função vai ler os aruivos iniciados com gwas_ e organiza sua apresentação
    output: df com todas as variantes cattalogadas para a doença 
    '''

    for filename in os.listdir(folder_path):
        if filename.startswith("gwas_"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep='\t')
    
    df=df[["STRONGEST SNP-RISK ALLELE","SNPS","RISK ALLELE FREQUENCY" ,
     "P-VALUE","OR or BETA" , "95% CI (TEXT)", "MAPPED_TRAIT","DISEASE/TRAIT","REGION" ,"CHR_ID","CHR_POS" 
     ,"REPORTED GENE(S)" ,"MAPPED_GENE", "CONTEXT", "INTERGENIC" ,"PUBMEDID", "FIRST AUTHOR", "LINK" , "STUDY ACCESSION" ]]
    
    data=df.rename(columns={"STRONGEST SNP-RISK ALLELE": "Variants and risk allele",  "SNPS": "RSID", "95% CI (TEXT)":"CI",
                                    "OR or BETA":"Odds_ratio", "CONTEXT":"Type_variant"})

    
    return data

def show_association(data, folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso arquivo gwas especifico
    data= aqruivo df com as variantes
    Essa função vai ler os aqruivos iniciados com gwas_ e carregar o df das varinates para achar correspondencia entre nosso db e o gwas
    identificando por rs, posição e gene, destacando as variantes identificadas por rs e posição.
    output: df com todas as variantes cattalogadas para a doença no gwas, df com a interseção das variantes e gwas
    '''
    gwas=read_files(folder_path)
    st.dataframe(gwas)
    st.write(gwas.shape)
    data=data[0] #o data é uma tupla visto no codigo de variantes a função retornar dois df

    all=[]
    #preciso separar por rs
    data["rsid_only"] = data["rsid"].str.split(",").str.get(0)
    #preciso separar por chr e posição
    data["chr_pos"] ="chr"+data["Variants"].str.split("_").str.get(0)+":"+data["Variants"].str.split("_").str.get(1)
    #preciso separar por gene e marcar a linha
    genes = set(data["SYMBOL"])

    gene = gwas[gwas["REPORTED GENE(S)"].apply(lambda x: any(g.strip() in genes for g in str(x).split(",")))]
    
    #filtro
    rs=gwas[gwas["RSID"].isin(data["rsid_only"])]
    pos=gwas[gwas["RSID"].isin(data["chr_pos"])]


    all.append(rs)
    all.append(pos)
    all.append(gene)
    
    associations= pd.concat(all, ignore_index=True)

    #marcar as linhas 

    mask_rs = associations["RSID"].isin(data["rsid_only"])
    mask_pos = associations["RSID"].isin(data["chr_pos"])
    

    def highlight_rows(row):
        if mask_rs.loc[row.name]:
            return ["background-color: lightyellow"] * len(row)
        elif mask_pos.loc[row.name]:
            return ["background-color: lightblue"] * len(row)


    associations_mask = associations.style.apply(highlight_rows, axis=1)
    st.dataframe(associations_mask)
    st.write(associations.shape)

    #st.dataframe(
    #associations.style.set_properties(
    #    subset=associations[associations["RSID"].isin(data["rsid_only"])], **{'background-color': '#e8f4f8', 'color': '#0b3d91', 'font-weight': 'bold'}
    #))

    def color_by_expression(val):
        color = '#90EE90' if val > 10 else '#ffcccb'
        return f'background-color: {color}'
    
    #st.dataframe(df.style.applymap(color_by_expression, subset=[""]))
    return associations
