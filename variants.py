#leitura do arquivo e interação com graficos

import streamlit as st
import pandas as pd
import os
import plotly.express as px
#folder_path="../../../../OneDrive/bio e/parkison/"
@st.cache_data
def read_files(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso de variantes
    Essa função vai ler os aruivos iniciados com variantes_ e concatena-los em um unico df
    output: df com todas as variantes[0], df com todas as variantes e o escore separado[1]
    '''
    all_data = []

    for filename in os.listdir(folder_path):
        if filename.startswith("variants_"):
            file_path = os.path.join(folder_path, filename)
            df = pd.read_csv(file_path, sep='\t')

            all_data.append(df)
            data=pd.concat(all_data, ignore_index=True)

    data=data.rename(columns={"X.Uploaded_variation": "Variants",  "Existing_variation": "rsid", "REVEL":"REVEL_score",
                                    "CADD_PHRED":"CADD_score", "MetaRNN_pred":"MetaRNN_Class"})

    data["Frequency_ALT"] = data["freqs"].str.split("/").str.get(1)
    dado=data[["Variants","SYMBOL","rsid", "Consequence","sift_class","PolyPhen_class", "am_class","CADD_class",
                    "MetaRNN_class","REVEL_class", "UNIPROT_ISOFORM","Feature"]]
    data=data[["Variants","SYMBOL","rsid","Consequence","SIFT","PolyPhen", "am_class","am_pathogenicity","CLIN_SIG","CADD_class", "CADD_score",
                    "MetaRNN_Class","MetaRNN_score","REVEL_class","REVEL_score","Frequency_ALT"]]
    return data,dado

def show_variants(folder_path):
    '''
    folder_path: especifica o caminho de leitura dos arquivos, nesse caso de variantes
    Essa função vai ler os aruivos iniciados com variantes e apresentar um df interativo e com poder de seleção do gene especifico,
    além de apresentar uma sumarização da classificação por preditor e um heatmap com a predição por variante
    output: df com todas as variantes interativo, um grafico de preditores, um heatmap por variante, um df so com as variantes novas
    '''
    data,dado=read_files(folder_path)
    genes=sorted(set(data["SYMBOL"]))

    gene_select = st.multiselect(
        "Select an specific gene...",
        genes,
        placeholder="Choose a gene...", 
         key="gene_selector")

    def color_by_expression(val):
        if pd.isna(val):
            return ''
        elif val == "'-":
            return ''
        color = '#90EE90' if val <0.5 else '#ffcccb'
        return f'background-color: {color}'
    
    
    if not gene_select:
        st.dataframe(data)
        #st.dataframe(data.style.applymap(color_by_expression, subset=["am_pathogenicity"]))
        st.write(data.shape)
    else:
        data=data[data["SYMBOL"].isin(gene_select)]
        dado=dado[dado["SYMBOL"].isin(gene_select)]
        st.dataframe(data)
        st.write(data.shape)

    col1,col2=st.columns(2)
    top = st.container()
    bottom = st.container()

    
    col1.subheader("Impact Predicts")
    predictors = ["am_class", "PolyPhen_class", "sift_class", "MetaRNN_class", "REVEL_class","CADD_class"]

    melted = dado.melt(
            id_vars=["Variants", "SYMBOL"],
            value_vars=predictors,
            var_name="Preditor",
            value_name="Classificação"
        )
    
    predictor_select = col1.selectbox("Choose a predictor:", 
                                    sorted(predictors),
                                    key="predictor_selector")
    if predictor_select:
        df_predictor = melted[melted["Preditor"] == predictor_select]

        fig_general = px.histogram(
            df_predictor,
            x="Preditor",
            color="Classificação",
            barmode="group",
            text_auto=True,
            title=f"Variants Classification by {predictor_select}",
            labels={'x':'Predictor', 'y':'count'}
        )
        col1.plotly_chart(fig_general, use_container_width=True)

        col2.subheader("Impact per Variant")

        variants_list=list(sorted(dado["Variants"].unique()))
        variants=["All"]+variants_list
        variant_select = col2.selectbox("Choose a variant:", 
                                    variants,
                                    key="variant_selector")

        if variant_select == "All":
            heatmap_df = melted.pivot_table(
                index="Variants",
                columns="Preditor",
                values="Classificação",
                aggfunc="first"
            )

            class_map = {"Deleterious":1, "D":1,"likely_deleterious":1, "T":2, "likely_benign":2,"benign":1, "-":0, "ambiguous":4}
            z_numeric = heatmap_df.map(lambda x: class_map.get(x, 0))

            fig_variant = px.imshow(
                z_numeric, 
                text_auto=False,
                aspect="auto",
                color_continuous_scale="RdYlBu", 
            )

            fig_variant.update_layout(
                title=f"Variant {variant_select} classification by predictors",
                xaxis_title="Predictor",
                yaxis_title="Variants",
                coloraxis_showscale=False
            )

            col2.plotly_chart(fig_variant, use_container_width=True)
        else:
            df_variant = melted[melted["Variants"] == variant_select]

            heatmap_df = df_variant.pivot_table(
                index="Variants",
                columns="Preditor",
                values="Classificação",
                aggfunc="first"
            )

            class_map = {"Deleterious":1, "D":1,"likely_deleterious":1, "T":2, "likely_benign":2,"benign":1, "-":0, "ambiguous":4}
            z_numeric = heatmap_df.applymap(lambda x: class_map.get(x, 0))

            fig_variant = px.imshow(
                z_numeric, 
                text_auto=False,
                aspect="auto",
                color_continuous_scale="RdYlBu", 
            )

            fig_variant.update_layout(
                title=f"Variant {variant_select} classification by predictors",
                xaxis_title="Preditor",
                yaxis_title="Variants",
                coloraxis_showscale=False
            )

            col2.plotly_chart(fig_variant, use_container_width=True)
    st.subheader("Only new variants")
    new=data[data["rsid"] == "-"]
    st.dataframe(new)
    st.write(new.shape)


 
def preditor_grafic_variants(folder_path):
    _,dado =read_files(folder_path)
    st.subheader("Impact per Variant")

    variants=list(dado["Variants"].unique()) + ["all"]
    variant_select = st.selectbox("Choose a variant:", 
                                  sorted(variants),
                                   key="variant_selector1")

    if variant_select:
        df_variant = melted[melted["Variants"] == variant_select]

        heatmap_df = df_variant.pivot_table(
            index="Variants",
            columns="Preditor",
            values="Classificação",
            aggfunc="first"
        )

        fig_variant = px.imshow(
            heatmap_df.applymap(str), 
            text_auto=True,
            aspect="auto",
            color_continuous_scale="RdYlBu", 
        )

        fig_variant.update_layout(
            title=f"Variant {variant_select} classification by predictors",
            xaxis_title="Preditor",
            yaxis_title="Variants",
            coloraxis_showscale=False
        )

        st.plotly_chart(fig_variant, use_container_width=True)