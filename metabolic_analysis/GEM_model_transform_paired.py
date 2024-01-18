import riptide
import sys
import os
from os.path import join  
import cobra
import json
import pandas as pd


from script.ECMpy_function import *
from script.AutoPACMEN_function import *
import re
import matplotlib.pyplot as plt
import argparse


"""
python GEM_model_transform_paired.py --input_c input_file/model_input_control.tsv --input_t input_file/model_input_treatment.tsv --output output_file/model_output_difference.tsv
"""

def read_config():
    parser = argparse.ArgumentParser(description='This code is used to integrate transcriptome data')
    parser.add_argument('--input_c', '-c', required=True, help='input rpkm file / C')
    parser.add_argument('--input_t', '-t', required=True, help='input rpkm file / T')
    parser.add_argument('--output', '-o', required=False, help='ouput flux file')
    arg = parser.parse_args()
    input_rpkm_c_path = arg.input_c
    input_rpkm_t_path = arg.input_t
    output_flux_path = arg.output

    return input_rpkm_c_path,input_rpkm_t_path,output_flux_path   

def bing_transcriptome(inputfile1,inputfile2):

    # 模型背景文件路径,需要修改
    model_path="model/ecMTM_TurNup_C13.json"

    # data.to_csv(inputfile, index=False, sep='\t', encoding='utf-8')
    my_model=get_enzyme_constraint_model(model_path)
    
    my_model.reactions.get_by_id('R2399').bounds = (-1000, 0.0)#葡萄糖打开
    my_model.reactions.get_by_id('R2402').bounds = (-1000, 0.0)#果糖打开
    my_model.reactions.get_by_id('R2412').bounds = (-1000, 0.0)#蔗糖打开
    my_model.reactions.get_by_id('R2417').bounds = (-1000, 0.0)#纤维二糖打开
    my_model.reactions.get_by_id('R2411').bounds = (-1000, 0.0)#麦芽糖打开
    my_model.reactions.get_by_id('R2401').bounds = (-1000, 0.0)#半乳糖打开
    my_model.reactions.get_by_id('R2408').bounds = (-1000, 0.0)#木糖糖打开
    my_model.reactions.get_by_id('R2405').bounds = (-1000, 0.0)#阿拉伯糖打开
    transcript_abundances_1 = riptide.read_transcription_file(inputfile1)#读取转录组数据
    riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)#整合数据
    model=riptide_object_1_a.model#得到模型
    solution = model.optimize()#求解
    reac=dict(solution.fluxes[solution.fluxes> 1e-10])#选取通量大的反应


    reactionlist=[]#反应列表信息
    fluxlist=[]#通量列表信息
    reactionlist2=[]#反应式
    for i in reac:
    # print( enz_model.reactions.get_by_id(i).gene_reaction_rule,":",enz_model.reactions.get_by_id(i), ":",reac[i])
    # print(i)
    # genelist.append(enz_model.reactions.get_by_id(i).gene_reaction_rule)
    # ID=i
        reID=str(model.reactions.get_by_id(i)).split(':')[0]
        # if ID=='R38_num':
        #     print(ID)
        #     print(reac[i])
        # reID= reID.split('_')[0]
        reactionlist.append(reID)
        fluxlist.append(reac[i])
        reactionlist2.append(model.reactions.get_by_id(i))
    df1=pd.DataFrame({'reactionid':reactionlist,'equation':reactionlist2,'flux':fluxlist})

    transcript_abundances_1 = riptide.read_transcription_file(inputfile2)#读取转录组数据
    riptide_object_1_a = riptide.contextualize(model=my_model, transcriptome=transcript_abundances_1)#整合数据
    model=riptide_object_1_a.model#得到模型
    solution = model.optimize()#求解
    reac=dict(solution.fluxes[solution.fluxes> 1e-10])#选取通量大的反应


    reactionlist=[]#反应列表信息
    fluxlist=[]#通量列表信息
    reactionlist2=[]#反应式
    for i in reac:
    # print( enz_model.reactions.get_by_id(i).gene_reaction_rule,":",enz_model.reactions.get_by_id(i), ":",reac[i])
    # print(i)
    # genelist.append(enz_model.reactions.get_by_id(i).gene_reaction_rule)
    # ID=i
        reID=str(model.reactions.get_by_id(i)).split(':')[0]
        # if ID=='R38_num':
        #     print(ID)
        #     print(reac[i])
        # reID= reID.split('_')[0]
        reactionlist.append(reID)
        fluxlist.append(reac[i])
        reactionlist2.append(model.reactions.get_by_id(i))
    df2=pd.DataFrame({'reactionid':reactionlist,'equation':reactionlist2,'flux':fluxlist})
    merged_df = pd.merge(df1, df2, on='reactionid', suffixes=('_df1', '_df2'))
    return(merged_df)

   
   
    
def main():

    # read args
    inputfile1,inputfile2,outputdir = read_config()
    df=bing_transcriptome(inputfile1,inputfile2)
    df.to_csv(outputdir,index=0,sep='\t')

if __name__ == "__main__":
    main()
