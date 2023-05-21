import pandas as pd
import numpy as np
import os
import json
import subprocess

from search_info_by_gene import get_gene_go_info_by_gene_id



def get_gene_list_by_pathway(workdir: str, species_name: str, pathway_name: str) -> list:
    """根据菌种pathway的名字，找到对应的pathway中基因的list

    Args:
    * species_name: 菌种的名字
    * pathway_name: pathway的名字
    * workdir: 工作目录
    ----------
    Returns: 
        gene_list: Pathway对应的基因的list
    """
    # 空格换成_，大写转小写
    species_name = species_name.replace(' ', '_').lower()

    # 读取kegg_filter表
    kegg_filter_path = os.path.join(workdir, species_name, 'kegg.tsv')
    df_kegg_filter = pd.read_csv(kegg_filter_path, sep='\t')

    # 根据KEGG_filter表的表头，找到对应的基因的list
    df = df_kegg_filter[df_kegg_filter['Pathway'] == pathway_name]
    gene_list = df['Gene'].tolist()
    gene_list_info = []

    for gene in gene_list:
        # 调用get_gene_go_info_by_gene_id函数，获取基因的go注释
        gene_info = get_gene_go_info_by_gene_id(workdir, species_name, gene)
        # 将每一个基因的go注释，添加到gene_list_info中
        gene_list_info.append(gene_info)

    return gene_list_info


def get_gene_list_exp_by_pathway(workdir: str, species_name: str, pathway_name: str) -> pd.DataFrame:
    """ 根据pathway的名字，依据exp_gsm表，返回对应基因的表达谱

    Args:
    * workdir: 工作目录
    * species_name: 菌种的名字
    * pathway_name: pathway的名字
    ----------
    Returns:
        df: 以GSM ID为表头的表达谱的dataframe
    """
    # 空格换成_，大写转小写
    species_name = species_name.replace(' ', '_').lower()

    # 读取exp_gsm表
    exp_path = os.path.join(workdir, species_name, 'exp_gsm.tsv')
    df_exp = pd.read_csv(exp_path, sep='\t')

    # 调用函数，找到对应的基因的list
    gene_list = get_gene_list_by_pathway(workdir, species_name, pathway_name)

    # 根据基因的list,找到对应的表达谱,标准化为log2(n+1),保留两位小数
    df = df_exp[df_exp['Gene id'].isin(gene_list)]
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: np.log2(x + 1)).round(2)

    return df


def get_gene_list_exp_fname_by_pathway(workdir: str, species_name: str, pathway_name: str) -> pd.DataFrame:
    """ 根据pathway的名字，依据exp_fname表，返回对应基因的表达谱
    """
    species_name = species_name.replace(' ', '_').lower()

    exp_fname_path = os.path.join(workdir, species_name, 'exp_fname.tsv')
    df_exp_fname = pd.read_csv(exp_fname_path, sep='\t')

    gene_list = get_gene_list_by_pathway(workdir, species_name, pathway_name)

    df = df_exp_fname[df_exp_fname['Gene id'].isin(gene_list)]
    df.iloc[:, 1:] = df.iloc[:, 1:].applymap(lambda x: np.log2(x + 1)).round(2)

    return df


def get_pathway_id_by_pathway(workdir: str, species_name: str, pathway_name: str) -> str:
    """根据pathway的名字，依据kegg_id表，返回在相应pathway的id，生成pathway的通路图，保存在kegg-viewers目录下

    Args:
    * workdir: 工作目录
    * species_name: 菌种的名字
    * pathway_name: pathway的名字
    ----------
    Returns:
        pathway_id: pathway的id

    Usage: 
    安装：pip install kegg_viewer
    :   生成pathway的通路图的命令为：kegg_viewer -p pathway_id -O output -c cache
    :   pathway_id为pathway在kegg数据库中的id，output_dir为输出的通路图的目录，cache为缓存的目录
    :   通路图的svg矢量图输出目录为：workdir/species_name/kegg-viewers
    :   缓存的目录为：workdir/species_name/cache
    """
    # 空格换成_，大写转小写
    species_name = species_name.replace(' ', '_').lower()

    # 读取kegg_id表
    kegg_path = os.path.join(workdir, species_name, 'pathway_id.tsv')
    df_kegg_id = pd.read_csv(kegg_path, sep='\t')

    # 设置表头名称，提取pathway的id
    df_kegg_id.columns = ['Id', 'Pathway'] 
    df = df_kegg_id[df_kegg_id['Pathway'] == pathway_name]
    pathway_id = df['Id'].tolist()[0]

    # try:
    #     os.system('kegg_viewer -p ' + pathway_id + 
    #               ' -O ' + os.path.join(workdir, species_name, 'kegg_viewers') + 
    #               ' -c ' + os.path.join(workdir, species_name, 'cache'))

    # except Exception as e:
    #     raise Exception("Failed to run kegg_viewer command. Error message: " + str(e))
    
    # 生成pathway的通路图
    try:
        output_dir = os.path.join(workdir, species_name, 'kegg_viewers')
        cache_dir = os.path.join(workdir, species_name, 'cache')
        
        cmd = f"kegg_viewer -p {pathway_id} -O {output_dir} -c {cache_dir}"
        subprocess.run(cmd, check=True, shell=True)

    except subprocess.CalledProcessError as e:
        raise Exception(f"Failed to run kegg_viewer command. Error message: {e}")

    return pathway_id


def get_dataset_info_by_species_name(workdir: str, species_name: str) -> pd.DataFrame:
    """根据菌种的名字，返回dataset表
    """
    species_name = species_name.replace(' ', '_').lower()
    df_dataset = pd.read_csv(os.path.join(workdir, species_name, 'dataset.tsv'), sep='\t')
    
    return df_dataset


def get_sample_info_by_species_name(workdir: str, species_name: str) -> pd.DataFrame:
    """根据菌种的名字，返回sample表
    """
    species_name = species_name.replace(' ', '_').lower()
    df_sample = pd.read_csv(os.path.join(workdir, species_name, 'sample.tsv'), sep='\t')
    
    return df_sample


if __name__ == '__main__':

    workdir = 'mtd/transcriptome/search_exp'

    # 测试
    species_name = 'Myceliophthora thermophila'
    pathway_name = 'Autophagy - yeast'
    print(get_gene_list_by_pathway(workdir, species_name, pathway_name))
    # print(get_gene_list_exp_by_pathway(workdir, species_name, pathway_name))
    # print(get_gene_list_exp_fname_by_pathway(workdir, species_name, pathway_name))
    # print(get_pathway_id_by_pathway(workdir, species_name, pathway_name))
    # print(get_dataset_info_by_species_name(workdir, species_name))
    # print(get_sample_info_by_species_name(workdir, species_name))







