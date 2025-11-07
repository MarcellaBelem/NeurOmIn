
import streamlit as st
from gprofiler import GProfiler
import matplotlib.pyplot as plt
import seaborn as sns


def enrichiment(data):
    '''
    data= df com a coluna SYMBOL dos genes
    Esta função faz a anotação funcional dos genes 
    '''
    
    gp = GProfiler(user_agent='ExampleTool', return_dataframe=True)

    data=data[0]
    genes = list(set(data["SYMBOL"]))

    results = gp.profile(organism='hsapiens', query=genes)
    st.dataframe(results.head())

    return results


def plot_enrichiment(result):
    '''
    Esta função faz um barplot do enriquecimento funcional
    :param result: df com o resultado do enriquecimento
    :param top: qts de enriquecimentos para plotar
    :return: a visualização do barplot
    '''
    #selecionar so os enriquecimentos significativos
    significant_results = result[result['p_value'] < 0.05]

    # organizar pelo valor de p
    significant_results = significant_results.sort_values(by='p_value', ascending=True)

    # Top resultados para plotar
    #top= st.selectbox("Top results", 5, 10, 20)
    options = [5, 10, 15, 20]
    top = st.pills("Top Results", options, selection_mode="single", default=5)
    top_results = significant_results.head(top) #escolher qts plotar
    
    col1,col2=st.columns(2)
    #plotar
    fig,ax=plt.subplots(figsize=(18, 10))
    sns.barplot(top_results, x="query_size", y="name", hue="p_value", legend=True, ax=ax)
    ax.set_xlabel('Counts')
    ax.set_ylabel('Process')
    ax.set_title(f'Top {top} Enriched Terms')
    ax.invert_yaxis()
    col1.pyplot(fig)

    fig_dot,ax=plt.subplots(figsize=(12, 8))
    sns.scatterplot(data=top_results, x="query_size", y="name",
                    size="query_size", hue="p_value", palette="viridis", ax=ax, legend="full", sizes=(50, 400))
    
    ax.set_xlabel('Counts')
    ax.set_ylabel('Process')
    ax.set_title(f'Top {top} Enriched Terms')
    ax.invert_yaxis()
    col2.pyplot(fig_dot)