import streamlit as st
import py3Dmol
import requests
import myvariant
import pandas as pd
import re
import json


#fazer um dicionario que vai servir de lista para o selectbox
def proteina(uniprot_id, uniprots): 
    '''
    Esta função faz a leitura do gene e apresenta sua estrutura proteica
    '''
   #os uniprots posso pegar de data, coluna UNIPROT_ISOFORM
    #uniprots={"SNCA":"P37840","VPS35":"Q96QK1"}  

    url = f"https://alphafold.ebi.ac.uk/files/AF-{uniprots[uniprot_id][0]}-F1-model_v6.pdb"
    #url =f"https://www.alphafold.ebi.ac.uk/entry/{uniprots[uniprot_id]}"
    r = requests.get(url)
    pdb_data = r.text

    st.write(f"Protein struture for {uniprot_id}({uniprots[uniprot_id][0]})")

    viewer = py3Dmol.view(width=800, height=500)
    viewer.addModel(pdb_data, "pdb")
    viewer.setStyle({'cartoon': {'color': 'lightblue'}})

    #for pos in mutations:
    # viewer.addStyle({'resi': pos}, {'stick': {'color': 'red'}})

    viewer.zoomTo()
    st.components.v1.html(viewer._make_html(), height=500, width=800)
    return pdb_data


def find_aa(data):
    mv=myvariant.MyVariantInfo()
    data=data[0]
    variant=pd.DataFrame()
    st.dataframe(data)
    data['variant_id']= ("chr"+data["Variants"].str.split("_").str.get(0)+":g."+data["Variants"].str.split("_").str.get(1)+data['Variants'].str.split("_").str.get(2).str.split("/").str[0]+">"+data['Variants'].str.split("/").str.get(1))

    variants = data["variant_id"].tolist()
    st.dataframe(variants)
    #preciso separar por gene e marcar a linha
    res = mv.getvariant(variants, fields="snpeff.ann.hgvs_p,snpeff.ann.hgvs_c,snpeff.ann.gene_name,snpeff.ann.effect,snpeff.ann.impact")
    print(res)
    #results_df = pd.json_normalize(res)
    #st.write(results_df)
    st.write("teste")
    vari = mv.getvariants(["chr1:g.218631822G>A", "chr9:g.107620835G>A"],  fields="snpeff")
    print(vari)
    st.write(vari)

def find_hgvsp(data):
    #data=data[0]
   
    data['variant_id']= ("chr"+data["Variants"].str.split("_").str.get(0)+":g."+data["Variants"].str.split("_").str.get(1)+data['Variants'].str.split("_").str.get(2).str.split("/").str[0]+">"+data['Variants'].str.split("/").str.get(1))

    
    data=data[data["Consequence"] == "missense_variant"]
    variant = data["variant_id"].tolist()
  



    url = "https://rest.ensembl.org/vep/human/hgvs"

    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {"hgvs_notations": variant}


    response = requests.post(url, headers=headers, json=payload, timeout=30)

    data = response.json()
 

    results = []
    for entry in data:
        for tc in entry.get("transcript_consequences", []):
            results.append({
                "gene": tc.get("gene_symbol"),
                "transcript": tc.get("transcript_id"),
                #"hgvsc": tc.get("hgvsc"),
                #"hgvsp": tc.get("hgvsp"),
                "amino_acids": tc.get("amino_acids"),
                "protein_start": tc.get("protein_start"),
                "protein_end": tc.get("protein_end"),
                "impact": tc.get("impact"),
                "consequence": ", ".join(tc.get("consequence_terms", []))
            })

    df = pd.DataFrame(results).drop_duplicates()

    return df



def show_variant(data, pdb_id, uniprot_id, uniprots):
    '''
    data: df com as variantes e o loci
    Esta função apresenta a mutação missense e su provavel impacto na estrutura proteica
    '''
    
    #converter a coordenada para ver a proteina

##selcionar so o transcrito para aparecer coluna Feature, , tenho que tirar o .1

    data=data[0]
    #data["Feature"] = data["Feature"].str.split(".").str.get(0)# separar a isoforma

    data=data[data["SYMBOL"]==uniprot_id]
    data=data[data["Consequence"] == "missense_variant"]
    variants=data["Variants"].tolist()

 
    variant_select=st.selectbox("choose a variant to see impact in protein:", variants, index=1)
    
    if variant_select:
    
    #fazer uma seleção para identificar essas posições
        variant=data[data["Variants"]== variant_select]
     
        #transcrito=uniprots[uniprot_id][1]
   
        df_hgvsp=find_hgvsp(variant)
        if not df_hgvsp.empty:
            df_hgvsp=df_hgvsp[df_hgvsp["transcript"]==transcrito]
            st.dataframe(df_hgvsp)
        #match = re.search(r"p\.[A-Za-z]+(\d+)[A-Za-z*]+", df_hgvsp["hgvsp"])
        #if match:
            #mutation_pos = int(match.group(1))
            mutation_pos=[int(df_hgvsp["protein_start"].iloc[0]), int(df_hgvsp["protein_end"].iloc[0])]
            aa=df_hgvsp["amino_acids"].iloc[0]
            st.write(f"Posição mutada: {mutation_pos}, AA mutados: {aa}")
            viewer = py3Dmol.view( width=800, height=500)
            viewer.addModel(pdb_id, "pdb")
            viewer.setStyle({'cartoon': {'color': 'lightgrey'}})
            viewer.addStyle({'resi': mutation_pos}, {'stick': {'color': 'red'}})

            viewer.zoomTo()
            st.components.v1.html(viewer._make_html(), height=500, width=800)
        else:
            st.write("we dont have information about this variant yet, try another one.")

    #fazer a predição do efeito da variante