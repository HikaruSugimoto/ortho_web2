import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import zipfile
import base64
from PIL import Image
from io import StringIO
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import os
#plt.rcParams['font.family']= 'sans-serif'
#plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'
plt.rcParams['xtick.major.width'] = 1.2
plt.rcParams['ytick.major.width'] = 1.2
plt.rcParams['font.size'] = 12 
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['figure.dpi'] = 300

    
def set_species_name(name):
    return name.split(' (')[1].replace(' ', '_')[:-1]

def upload_demo_data():
    if os.path.isfile('demo.zip'):
        os.remove('demo.zip')
    with zipfile.ZipFile('demo.zip', 'x') as csv_zip:
        csv_zip.writestr("demo.csv", pd.read_csv("./demo/demo.csv").to_csv(index=False))   
    with open("demo.zip", "rb") as file:
        zip_data = file.read()
        b64 = base64.b64encode(zip_data).decode()
        zip_filename = 'demo.zip'
        href = f'<a href="data:application/zip;base64,{b64}" download="{zip_filename}">Download demo data (Golden hamster gene symbols)</a>'
        st.sidebar.markdown(href, unsafe_allow_html=True)
    if os.path.isfile('demo.zip'):
        os.remove('demo.zip')

def download_result(df):
    if os.path.isfile('result.zip'):
        os.remove('result.zip')
    with zipfile.ZipFile('result.zip', 'x') as csv_zip:
        csv_zip.writestr("result.csv", df.to_csv(index=False))   
    with open("result.zip", "rb") as file:
        zip_data = file.read()
        b64 = base64.b64encode(zip_data).decode()
        zip_filename = 'result.zip'
        href = f'<a href="data:application/zip;base64,{b64}" download="{zip_filename}">Download result data</a>'
        st.markdown(href, unsafe_allow_html=True)
    if os.path.isfile('result.zip'):
        os.remove('result.zip')

def display_phylogenetic_tree(Phylo_s):
    f = open('./Gene_Trees/'+Phylo_s+'_tree.txt', 'r')
    data = f.read()
    tree: Tree = Phylo.read(file=StringIO(data), format="newick")

    fig = plt.figure(figsize=(3,6), dpi=100)
    axes = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, axes=axes, do_show=False)

    axes = plt.gca()
    axes.spines['top'].set_visible(False)
    axes.spines['right'].set_visible(False)
    axes.spines['bottom'].set_visible(False)
    axes.spines['left'].set_visible(False)
    axes.set_xticks([])
    axes.set_yticks([])
    axes.set_xlabel('')
    axes.set_ylabel('')
    #axes.set_title(df_ratio_all["OG"][i])
    #axes.set_xlim(-0.5, 1)
    st.pyplot(fig)       

def extract_first(x):
    if pd.isna(x):
        return x
    parts = x.split(",")
    unique_parts = []
    for part in parts:
        part = part.strip()
        if part not in unique_parts:
            unique_parts.append(part)
    return ",".join(unique_parts)

#Main
st.set_page_config(layout="wide")
st.title("A database of mammalian orthologs and phylogenetic trees")
upload_demo_data()

#Input
options1 = ['Golden hamster (Mesocricetus auratus)','Human (Homo sapiens)',
            'Mouse (Mus musculus)','Rat (Rattus norvegicus)', 'Dog (Canis lupus familiaris)',
            'Elephant (Loxodonta africana)','Lesser hedgehog tenrec (Echinops telfairi)', 
            'Thirteen-lined ground squirrel (Ictidomys tridecemlineatus)',
            'Gray mouse lemur (Microcebus murinus)', 'Little brown bat (Myotis lucifugus)',
            'Squirrel monkey (Saimiri boliviensis boliviensis)', 
            'Tasmanian devil (Sarcophilus harrisii)',
            'Pig (Sus scrofa)', 'Arctic ground squirrel (Urocitellus parryii)',
            'American black bear (Ursus americanus)']
options2 = ['Mouse (Mus musculus)','Rat (Rattus norvegicus)','Human (Homo sapiens)',
            'Golden hamster (Mesocricetus auratus)',
            'Dog (Canis lupus familiaris)',
            'Elephant (Loxodonta africana)','Lesser hedgehog tenrec (Echinops telfairi)', 
            'Thirteen-lined ground squirrel (Ictidomys tridecemlineatus)',
            'Gray mouse lemur (Microcebus murinus)', 'Little brown bat (Myotis lucifugus)',
            'Squirrel monkey (Saimiri boliviensis boliviensis)', 
            'Tasmanian devil (Sarcophilus harrisii)',
            'Pig (Sus scrofa)', 'Arctic ground squirrel (Urocitellus parryii)',
            'American black bear (Ursus americanus)']
Input_s= st.sidebar.selectbox('Input species',options1)
Output_s= st.sidebar.selectbox('Output species',options2)      
Input_name = st.sidebar.file_uploader("Upload input species' gene symbols", type="csv")

Input_s = set_species_name(Input_s)
Output_s = set_species_name(Output_s)

#Home
if Input_name is None:
    #image = Image.open('./Fig/0.png')
    #st.image(image, caption='',use_column_width=True)
    st.subheader('Overview')
    st.write("This web application outputs orthologs of the other species shown below and phylogenetic trees of the orthologs. Upload gene symbols of the species of interest in the following format.")
    image = Image.open('./Fig/1.png')
    st.image(image, caption='',use_column_width=True)
    
    st.subheader('Methods')
    f=open('./Fig/3.txt', 'r')
    st.write(f.read())
    f.close()
        
    st.subheader('License')
    f=open('./Fig/4.txt', 'r')
    st.write(f.read())
    f.close()
       
    st.subheader('Help')
    st.write("Download codes and run this web application on a local machine: https://github.com/HikaruSugimoto/ortho_web")
    st.write("Please send e-mail to Hikaru Sugimoto and Shinya Kuroda if you have any questions or requests.")
#Output
if Input_name is not None:
    st.subheader('Uploaded data')
    Input_name = pd.read_csv(Input_name, names=[Input_s])
    st.write(Input_name)  
    
    #Ortholog conversion
    df_ortho=pd.read_csv('./Orthologues/Orthologues_'+Input_s+'/'+Input_s+'__v__'+Output_s+'.tsv', delimiter='\t')
    df_ortho=df_ortho.reset_index(drop=False)
    
    #Load ID conversion file
    ID_input=pd.read_csv('./ID_conversion/'+Input_s+'.tsv', delimiter='\t')
    ID_input=ID_input.rename(columns={"Protein accession": Input_s}, inplace=False)
    ID_input=ID_input[[Input_s,"Symbol"]].dropna().reset_index(drop=True)

    ID_output=pd.read_csv('./ID_conversion/'+Output_s+'.tsv', delimiter='\t')
    ID_output=ID_output.rename(columns={"Protein accession": Output_s}, inplace=False)
    ID_output=ID_output[[Output_s,"Symbol"]].dropna().reset_index(drop=True)
    
    #Convert input species ID to gene symbol
    Input_all=pd.DataFrame()
    All=df_ortho[Input_s].str.split(', ', expand=True).add_prefix('input_')
    All=pd.concat([df_ortho["index"],All],axis=1)
    for i in range (0,len(All.T)-1):
        conv=pd.merge(ID_input, All.rename(columns={"input_"+str(i): Input_s}, inplace=False).astype('object'),
                                on=Input_s, how='inner')
        Input_all=pd.concat([Input_all,conv[["index","Symbol"]]],axis=0)
    Input_all = Input_all.groupby("index").agg({'Symbol': lambda x: ', '.join(x)}).reset_index(drop=False)

    #Convert output species ID to gene symbol
    Output_all=pd.DataFrame()
    All=df_ortho[Output_s].str.split(', ', expand=True).add_prefix('output_')
    All=pd.concat([df_ortho["index"],All],axis=1)
    for i in range (0,len(All.T)-1):
        conv=pd.merge(ID_output, All.rename(columns={"output_"+str(i): Output_s}, inplace=False).astype('object'),
                                on=Output_s, how='inner')
        Output_all=pd.concat([Output_all,conv[["index","Symbol"]]],axis=0)
    Output_all = Output_all.groupby("index").agg({'Symbol': lambda x: ', '.join(x)}).reset_index(drop=False)
    
    #combine
    In2Out=pd.merge(Input_all, Output_all,on='index', how='outer')
    In2Out=pd.merge(df_ortho, In2Out,on='index', how='outer')
    In2Out=In2Out.rename(columns={"Symbol_x": Input_s+'_gene_symbol'}, inplace=False)
    In2Out=In2Out.rename(columns={"Symbol_y": Output_s+'_gene_symbol'}, inplace=False)
    In2Out=In2Out.drop("index", axis=1)
    In2Out=In2Out[["Orthogroup",Input_s+'_gene_symbol',Output_s+'_gene_symbol',Input_s,Output_s]]
    In2Out=In2Out.reset_index(drop=False)
    In2Out=In2Out.rename(columns={Input_s: Input_s+'_accession'}, inplace=False)
    In2Out=In2Out.rename(columns={Output_s: Output_s+'_accession'}, inplace=False)
    
    #Convert input species gene symbol to output species gene symbol
    All=In2Out[Input_s+'_gene_symbol'].str.split(', ', expand=True).add_prefix('input_')
    All=pd.concat([In2Out["index"],All],axis=1)
    In2Out1=In2Out.drop(Input_s+'_gene_symbol', axis=1)

    Inp_all=pd.DataFrame()
    for i in range (0,len(All.T)-1):
        conv=All.rename(columns={"input_"+str(i): Input_s+'_gene_symbol'}, inplace=False)[["index",Input_s+'_gene_symbol']]
        conv=pd.merge(In2Out1, conv[["index",Input_s+'_gene_symbol']].dropna(),on='index', how='inner')
        Inp_all=pd.concat([Inp_all,conv[['index',"Orthogroup",Input_s+'_gene_symbol',
                                        Output_s+'_gene_symbol',Input_s+'_accession',Output_s+'_accession']]],axis=0)
    Inp_all=Inp_all.sort_values('index')
    Inp_all=Inp_all.drop_duplicates()
    Inp_all=Inp_all.drop('index', axis=1)    
    
    Inp_all1=pd.merge(Input_name.rename(columns={Input_s: Input_s+'_gene_symbol'}, inplace=False),
                    Inp_all,on=Input_s+'_gene_symbol', how='inner')
    Inp_all1=pd.merge(Input_name.rename(columns={Input_s: Input_s+'_gene_symbol'}, inplace=False),
                    Inp_all1,on=Input_s+'_gene_symbol', how='outer')
    Inp_all1=Inp_all1[[Input_s+'_gene_symbol',Output_s+'_gene_symbol',Input_s+'_accession',Output_s+'_accession','Orthogroup']]
    Inp_all1=Inp_all1.drop_duplicates()

    Inp_all1['extracted']=Inp_all1[Output_s+'_gene_symbol'].apply(extract_first)
    Inp_all1=Inp_all1.drop(Output_s+'_gene_symbol', axis=1)
    Inp_all1=Inp_all1.rename(columns={'extracted': Output_s+'_gene_symbol'}, inplace=False)
    Inp_all1=Inp_all1[[Input_s+'_gene_symbol',Output_s+'_gene_symbol',
                    Input_s+'_accession',Output_s+'_accession','Orthogroup']]
    Inp_all1=Inp_all1.rename(columns={Input_s+'_gene_symbol':'Input'}, inplace=False)
    
    st.write(Inp_all1.set_index('Input'))
    download_result(Inp_all1)
    
    options3 = Inp_all1["Orthogroup"].dropna()
    Phylo_s= st.selectbox('Select the number of the ortholog group for which you want to display the phylogenetic tree.',options3)
    display_phylogenetic_tree(Phylo_s)
