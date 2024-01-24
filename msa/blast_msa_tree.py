
import os
import pandas as pd
import subprocess
import sys



def run_blastp(blast_input_path, blast_output_path, db_prot_path, evalue=1e-6):
    """
    根据输入序列信息，对蛋白库进行blastp比对，得到比对结果。
    
    Args:
        blast_input_path (str): 输入txt文件的路径。
        blast_output_path (str): 输出txt文件的路径。
        db_prot_name (str): 使用的蛋白数据库名称
        evalue (float): evalue值。
    """
    # 构建命令行参数
    cmd = [
        'blastp',
        '-query', blast_input_path,
        '-out', blast_output_path,
        '-db', db_prot_path,
        '-outfmt', '6',
        '-evalue', str(evalue)
    ]
    # 执行命令并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr



def run_mafft(workdir, blast_seq_path, mafft_result_path):
    """使用mafft进行多序列比对: MAC版
    
    Args:
        workdir: 工作目录
        blast_result_path: 输入文件路径
        mafft_result_path: 输出文件路径
    """
    mafft_path = os.path.join(workdir, 'mafft-mac/mafft.bat')
    # print(mafft_path)

    os.system(mafft_path+" --auto "+blast_seq_path+" > "+mafft_result_path)



def run_mafft_linux(blast_seq_path, mafft_result_path):
    """使用mafft进行多序列比对: linux版

    Args:
        blast_result_path: 输入文件路径
        mafft_result_path: 输出文件路径
    """

    os.system("mafft"+" --auto "+blast_seq_path+" > "+mafft_result_path)
    


def run_fasttree(workdir, mafft_result_path, fasttree_result_path):
    """
    根据多序列比对结果，使用fasttree进行进化树构建。: MAC版

    Args:
        mafft_result_path: 输入文件路径
        fasttree_result_path: 输出文件路径
    """
    fasttee_path = os.path.join(workdir, 'fasttree-mac/FastTree')
    # 构建命令行参数
    os.system(fasttee_path+' '+mafft_result_path+" > "+fasttree_result_path)


def run_fasttree_linux(mafft_result_path, fasttree_result_path):
    """
    根据多序列比对结果，使用fasttree进行进化树构建。: linux版

    Args:
        mafft_result_path: 输入文件路径
        fasttree_result_path: 输出文件路径
    """
    os.system("FastTree"+' '+mafft_result_path+" > "+fasttree_result_path)
    

def blast_msa_tree(workdir, blast_input_path, blast_result_path, db_prot_path, blast_seq_path, mafft_result_path, fasstree_result_path, evalue=1e-6):
    """运行blastp、mafft、fasttree，将输入蛋白序列信息进行分析，得到蛋白序列文件、多序列比对结果、进化树结果。

    Args:
        workdir (str): 工作目录
        blast_input_path (str): 运行blastp的输入文件路径
        blast_result_path (str): blastp的结果文件路径
        db_prot_name (str): 使用的蛋白数据库名称
        blast_seq_path (str): blastp的输入与结果中匹配到的蛋白序列文件路径
        mafft_result_path (str): 输出多序列比对结果的路径
        fasstree_result_path (str): 输出进化树结果的路径
        evalue (float): evalue值。
    """
    # 调用blastp
    run_blastp(blast_input_path, blast_result_path, db_prot_path, evalue)
    print("blastp over!")

    # 读取blastp结果,并设置列名
    df_blast = pd.read_csv(blast_result_path, sep='\t', header=None).copy()
    df_blast.columns = ['Query ID', 'Subject Id', 'Identity', 'Alignment length', 'Mismatches',
                        'Gap', 'q.start', 'q.end', 's.start', 's.end', 'evalue', 'bit score']
    
    # 保存比对结果
    df_blast.to_csv(blast_result_path, index=False)
    
    # df_blast中，subject_id为蛋白id，有重复，需要去重
    df_blast = df_blast.drop_duplicates(subset=['Subject Id'])
    
    # 读取记录了所有菌种蛋白id与序列对应关系的表
    df_all_fungi_seq = pd.read_csv(os.path.join(workdir, 'All_Fungi.tsv'), sep='\t')

    # 数据匹配
    merged_df = pd.merge(df_blast, df_all_fungi_seq, left_on='Subject Id', right_on='Protein ID', how='inner')
    fasta_sequences = ""

    # 将匹配到的序列信息与用户输入的序列信息合并
    with open(blast_input_path, 'r') as file:
        for line in file:
            fasta_sequences += line

    # 将数据写入FASTA格式
    for index, row in merged_df.iterrows():
        # fasta_sequences += ">{}\n{}\n".format(row['Protein ID'], row['Sequence'])
        fasta_sequences += ">{}-{}\n{}\n".format(row['Protein ID'], row['Species'], row['Sequence']) # 添加了菌种名

    try:
        with open(blast_seq_path, 'w') as file:
            file.write(fasta_sequences)
    except Exception as e:
        print("wrong:", str(e))
        sys.exit(1)
    print("Protein sequence information saved!")

    # 调用mafft进行多序列比对
    run_mafft(workdir, blast_seq_path, mafft_result_path)
    print("Multiple sequence alignment completed!")

    # 调用fasttree进行进化树构建
    run_fasttree(workdir, mafft_result_path, fasstree_result_path)
    print("Evolutionary tree construction completed")
    
    return None



if __name__ == '__main__':

    # 示例调用
    blast_msa_tree("/Users/dongjiacheng/Desktop/Github/msa/",
                    "/Users/dongjiacheng/Desktop/Github/msa/input_file/blast_input.txt",
                    "/Users/dongjiacheng/Desktop/Github/msa/output_file/blast_result.txt",
                    "/Users/dongjiacheng/Desktop/Github/msa/blast/db_prot/Myceliophthora_thermophila_ATCC_42464",
                    "/Users/dongjiacheng/Desktop/Github/msa/output_file/blast_seq.fasta",
                    "/Users/dongjiacheng/Desktop/Github/msa/output_file/mafft_result.fasta",
                    "/Users/dongjiacheng/Desktop/Github/msa/output_file/fasttree_result.nwk",
                    1e-6)