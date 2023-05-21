import pandas as pd
import os
import argparse
import json


def json_load(path):
    """
    Description: 读取json文件,转为字典
    """
    with open(path) as f:
        dict = json.load(f)

    return dict


def pathway_name_json(json_path):
    """
    Description: 根据json文件内容，提取pathway的名字
    """
    pathway_name = json_load(json_path)['pathway']

    return pathway_name


def pathway_gene_list_json(json_path, kegg_filter_path):
    """
    Description: 根据pathway名称提取相应基因list
    Arguments:
        json_path: json文件的路径
        kegg_filter_path: kegg_filter表的路径
    Returns:
        gene_list: Pathway对应的基因的list
    """
    df_kegg_filter = pd.read_csv(kegg_filter_path, sep='\t')

    # 根据json文件内容，提取pathway的名字
    pathway_name = pathway_name_json(json_path)

    # 根据KEGG_filter表的表头，找到对应的基因的list
    df = df_kegg_filter[df_kegg_filter['Pathway'] == pathway_name]
    gene_list = df['Gene'].tolist()

    return gene_list


def gene_exp_df_json(json_path, kegg_filter_path, exp_path):
    """
    Description: 根据json文件解析出的pathway的名字，找到对应的基因的list，然后根据df_exp，返回对应基因的表达量
    Arguments:
        json_path: json文件的路径
        kegg_filter_path: kegg_filter表的路径
        exp_fname_path: exp表的路径
    Returns:
        json_str: 表达谱的json字符串
    """
    df_exp = pd.read_csv(exp_path, sep='\t')

    # 根据pathway的名字，找到对应的pathway中基因的list
    gene_list = pathway_gene_list_json(json_path, kegg_filter_path)

    # 根据基因的list，找到对应的基因的表达量
    df = df_exp[df_exp['Gene id'].isin(gene_list)]

    json_str = df.to_json(orient='records')
    
    return json_str


def gene_exp_fname_df_json(json_path, kegg_filter_path, exp_fname_path):
    """
    Description: 根据json文件解析出的pathway的名字，找到对应的基因的list，然后根据df_exp_fname，返回对应基因的表达量
    Arguments:
        json_path: json文件的路径
        kegg_filter_path: kegg_filter表的路径
        exp_fname_path: exp_fname表的路径
    Returns:
        json_str: 表达谱的json字符串
    """
    df_exp_fname = pd.read_csv(exp_fname_path, sep='\t')
    
    # 根据pathway的名字，找到对应的pathway中基因的list
    gene_list = pathway_gene_list_json(json_path, kegg_filter_path)

    # 根据基因的list，找到对应的基因的表达量
    df = df_exp_fname[df_exp_fname['Gene id'].isin(gene_list)]

    json_str = df.to_json(orient='records')

    # with open(json_output_path,'w') as f:
    #     f.write(json_str)
    
    return json_str


if __name__ == '__main__':

    workdir = '/home/dongjc/mtd/transcriptome/search_exp/'
    kegg_filter_path = workdir + 'Mt_KEGG_filter.tsv'
    exp_path = workdir + 'Mt_Exp_gsm.tsv'
    exp_fname_path = workdir + 'Mt_Exp_fname.tsv'

    input_json_path = workdir + 'pathway.json'

    # 测试
    gene_list = pathway_gene_list_json(input_json_path, kegg_filter_path)
    json_exp = gene_exp_df_json(input_json_path, kegg_filter_path, exp_path)
    json_exp_fname = gene_exp_fname_df_json(input_json_path, kegg_filter_path, exp_fname_path)
    print(gene_list)
    print(json_exp)
    print(json_exp_fname)


    """ 命令行调用
    usage: search_exp_by_pathway.py pathwat_name [--gene_list] [--exp] [--exp_fname]
    1. search_exp_by_pathway.py --gene_list
    2. search_exp_by_pathway.py --exp
    3. search_exp_by_pathway.py --exp_fname
    """
    # parser = argparse.ArgumentParser()
    # parser.add_argument('--gene_list', action='store_true', help='output gene list of the pathway')
    # parser.add_argument('--exp', action='store_true', help='output gene expression change dataframe')
    # parser.add_argument('--exp_fname', action='store_true',help='output gene expression condition dataframe')
    # args = parser.parse_args()

    # if args.gene_list:
    #     gene_list = pathway_gene_list_json(input_json_path, kegg_filter_path)
    #     print(gene_list)
    # if args.exp:
    #     json_exp = gene_exp_df_json(input_json_path, kegg_filter_path, exp_path)
    #     print(json_exp)
    # if args.exp_fname:
    #     json_exp_fname = gene_exp_fname_df_json(input_json_path, kegg_filter_path, exp_fname_path)
    #     print(json_exp_fname)
