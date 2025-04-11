#-*- coding:utf-8 -*-

import pandas as pd
import numpy as np
import subprocess
import plotly.express as px
import os


def run_enrich(workdir, input_path, output_path, species, padjust, enrich_type):
    """
    运行GO或KEGG富集分析R脚本。

    Args:
        workdir (str): 工作目录。
        input_path (str): 输入文件的路径。
        output_path (str): 输出文件的路径。
        species (str): 菌种名称。
        padjust (float): P值阈值。
        enrich_type (str): 分析类型，"GO" 或 "KEGG"。
    """
    # 根据enrich_type确定脚本路径
    if enrich_type == 'GO':
        script_name = 'go_enrich.R'
    elif enrich_type == 'KEGG':
        script_name = 'kegg_enrich.R'
    else:
        raise ValueError("enrich_type must be 'GO' or 'KEGG'")

    script_path = os.path.join(workdir, script_name)

    cmd = [
        'Rscript', script_path,
        '--workdir', workdir,
        '--input', input_path,
        '--output', output_path,
        '--species', species,
        '--padjust', str(padjust),
    ]

    # 执行R脚本并捕获输出
    try:
        result = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        return e.stderr


def plot_kegg_chart(workdir, kegg_result, width=1280, height=720, p_adjust=0.05, font_size=15, chart_num=30, chart_size=30, color='rdbu_r',pic_type='bubble', funciton_type='All'):
    """根据输入的kegg富集结果，绘制气泡图

    Args:
        workdir (str): 工作目录
        kegg_result (pd.DataFrame): kegg富集结果
        width (int): 图表宽度. 
        height (int): 图表高度. 
        p_adjust (float): P值阈值. 
        font_size (int): 字体大小. 
        chart_num (int): 最多显示的富集通路数量.
        chart_size (int): 气泡大小. 
        pic_type (str): 图表类型.有两种选择：'bubble'和'bar'
        color (str): 颜色. 

    Returns: 输出图表、html文件的路径.
    """
    # 读取数据
    kegg_result = pd.read_csv(kegg_result, sep='\t')

    # 数据预处理
    kegg_result = kegg_result.copy()
    kegg_result = kegg_result[['ID', 'Description', 'GeneRatio', 'p.adjust', 'Count']]
    kegg_result.columns = ["ID", "Pathway", "GeneRatio","P.adjust", 'Count']
    kegg_result["GeneRatio"] = kegg_result["GeneRatio"].apply(lambda x: round(eval(x), 3))  # GeneRatio列输出处理为浮点数
    kegg_result['P.adjust'] = kegg_result['P.adjust'].apply(lambda x: round(x, 6))  # 控制P.adjust列的小数位数

    # 数据筛选
    kegg_result = kegg_result[kegg_result['P.adjust'] < p_adjust]  # 过滤P.adjust值
    kegg_result = kegg_result.sort_values(by='Count', ascending=False)  # 按照Count列降序排列
    kegg_result = kegg_result.iloc[:chart_num]  # 取前chart_num个数据

    if funciton_type == 'All':
        pass

    # 图表公共布局设置
    layout_args = {
        'title': "KEGG Enrichment Analysis",
        'yaxis_title': "Pathway",
        'yaxis': dict(autorange="reversed"),
        'font': dict(family="Arial", size=font_size),
        'template': "plotly_white",
        'width': width,
        'height': height
    }
    # 自定义颜色轴设置
    color_axis_args = {
        'colorbar_title': "P.adjust",
        'colorbar_tickformat': ".3f",
        'colorbar': dict(
            dtick=0.005,
            tickfont=dict(
                size=12  # 在这里设置你想要的字体大小
            )    
        )
    }

    # 确定颜色轴的tickvals和ticktext
    unique_p_adjust_values = sorted(kegg_result['P.adjust'].unique())
    if len(unique_p_adjust_values) == 1:
        # 只有一个唯一的P值
        color_axis_args['colorbar'].update({
            'tickvals': unique_p_adjust_values,
            'ticktext': [f'{unique_p_adjust_values[0]:.3f}']
        })
    elif len(unique_p_adjust_values) == 2:
        # 有两个唯一的P值
        color_axis_args['colorbar'].update({
            'tickvals': [unique_p_adjust_values[0], unique_p_adjust_values[-1]],
            'ticktext': [f'{unique_p_adjust_values[0]:.3f}', f'{unique_p_adjust_values[-1]:.3f}']
        })

    # 根据pic_type绘制不同类型的图表
    if pic_type == 'bubble':
        fig = px.scatter(
            kegg_result,
            x='GeneRatio',
            y='Pathway',
            size='Count',
            color='P.adjust',
            color_continuous_scale=color,
            opacity=0.85,
            hover_data=["ID",'P.adjust', 'Count'],
            size_max=chart_size
        )
    elif pic_type == 'bar':
        fig = px.bar(
            kegg_result,
            x='Count',
            y='Pathway',
            color='P.adjust',
            color_continuous_scale=color,
            opacity=0.85,
            hover_data=['ID','P.adjust', 'Count']
        )

    # 应用颜色轴设置
    fig.update_layout(**layout_args)
    if 'tickvals' in color_axis_args['colorbar']:  # 只有在tickvals存在时才更新颜色轴
        fig.update_coloraxes(**color_axis_args)

    # 设置交互式图表下载按钮导出高分辨率图片的配置
    config = {
    'toImageButtonOptions': {
        'format': 'png',  # one of png, svg, jpeg, webp
        'filename': 'Picture',
        'height': height,
        'width': width,
        'scale': 4  # 缩放参数，增加这个参数可以提高图像的分辨率
        }
    }

    # 保存为png，scale设置为4
    output_image_path = os.path.join(workdir, "kegg.png")
    fig.write_image(output_image_path, scale=4)
    # 保存为html
    output_html_path = os.path.join(workdir, "kegg.html")
    fig.write_html(output_html_path,config=config)
    
    # 测试用
    return output_image_path, output_html_path


def plot_go_chart(workdir, go_relust, width=1280, height=720, p_adjust=0.05, font_size=15, chart_num=30, chart_size=30, pic_type='bubble', color='rdbu_r', funciton_type='All'):
    """根据输入的GO富集分析结果，根据用户选择绘制气泡图或柱状图。

    Args:
        workdir (str): 工作目录。
        go_relust (DataFrame): GO富集分析结果。
        width (int): 图表宽度. 
        height (int): 图表高度. 
        p_adjust (float): P值阈值. 
        chart_num (int): 最多显示富集功能数量. 
        chart_size (int): 图表大小. 
        pic_type (str): 图表类型，可选bubble或bar. 
        color (str): 颜色. Defaults to 'Geyser'.
        funciton_type (str): GO类型，可选BP, CC, MF, All. 

    Returns: 输出图表、html文件的路径.
    """
    # 读取数据
    go_relust = pd.read_csv(go_relust, sep='\t')

    # 数据处理
    go_relust = go_relust.copy()
    go_relust = go_relust[["category", "ID", "Description", "Count", 'GeneRatio', "p.adjust"]]
    go_relust.columns = ["Class", "ID", "Description", "Count", "GeneRatio", "P.adjust"]

    # 数据列处理
    go_relust["GeneRatio"] = go_relust["GeneRatio"].apply(lambda x: round(eval(x), 3))
    go_relust['P.adjust'] = go_relust['P.adjust'].apply(lambda x: round(x, 6))
    go_relust = go_relust.sort_values(by='Count', ascending=False)
    go_relust = go_relust[go_relust["P.adjust"] < p_adjust]
    go_relust = go_relust.iloc[:chart_num]

    # 过滤GO类型
    if funciton_type in ["BP", "CC", "MF"]:
        go_relust = go_relust[go_relust["Class"].str.contains(funciton_type)]

    # 图表公共布局设置
    layout_args = {
        'title': "GO Enrichment Analysis",
        'yaxis_title': "Description",
        'yaxis': dict(autorange="reversed"),
        # 'yaxis': dict(autorange="reversed", showgrid=False), # y轴反向，不显示网格线
        # 'xaxis': dict(showgrid=False), # 不显示网格线
        'font': dict(family="Arial", size=font_size),
        'template': "plotly_white", # 主题背景
        'width': width,
        'height': height
    }

    # 自定义颜色轴设置
    color_axis_args = {
        'colorbar_title': "P.adjust",
        'colorbar_tickformat': ".3f",
        'colorbar': dict(
            dtick=0.005,
            tickfont=dict(
                size=12  # 在这里设置你想要的字体大小
            )    
        )
    }

    # 确定颜色轴的tickvals和ticktext
    unique_p_adjust_values = sorted(go_relust['P.adjust'].unique())
    if len(unique_p_adjust_values) == 1:
        # 只有一个唯一的P值
        color_axis_args['colorbar'].update({
            'tickvals': unique_p_adjust_values,
            'ticktext': [f'{unique_p_adjust_values[0]:.3f}']
        })
    elif len(unique_p_adjust_values) == 2:
        # 有两个唯一的P值
        color_axis_args['colorbar'].update({
            'tickvals': [unique_p_adjust_values[0], unique_p_adjust_values[-1]],
            'ticktext': [f'{unique_p_adjust_values[0]:.3f}', f'{unique_p_adjust_values[-1]:.3f}']
        })

    # 根据pic_type绘制不同类型的图表
    if pic_type == "bubble":
        fig = px.scatter(
            go_relust,
            x="GeneRatio",
            y="Description",
            size="Count",
            color="P.adjust",
            color_continuous_scale=color,
            opacity=0.85,
            hover_name="Class",
            hover_data=["ID", "Description", "Count", "GeneRatio", "P.adjust"],
            size_max=chart_size,
        )
        
    elif pic_type == "bar":
        fig = px.bar(
            go_relust,
            x='Count',
            y='Description',
            color='P.adjust',
            color_continuous_scale=color,
            opacity=0.85,
            hover_data=["Class", "ID", "GeneRatio", "P.adjust"],
        )

    # 应用颜色轴设置
    fig.update_layout(**layout_args)
    if 'tickvals' in color_axis_args['colorbar']:  # 只有在tickvals存在时才更新颜色轴
        fig.update_coloraxes(**color_axis_args)

    # 设置交互式图表下载按钮导出高分辨率图片的配置
    config = {
    'toImageButtonOptions': {
        'format': 'png',  # one of png, svg, jpeg, webp
        'filename': 'Picture',
        'height': height,
        'width': width,
        'scale': 4  # 缩放参数，增加这个参数可以提高图像的分辨率
        }
    }

    # 保存为png，scale设置为4
    output_image_path = os.path.join(workdir, "go.png")
    fig.write_image(output_image_path, scale=4)
    # 保存为html
    output_html_path = os.path.join(workdir, "go.html")
    fig.write_html(output_html_path,config=config)

    return output_image_path, output_html_path


if __name__ == "__main__":
    
    # 运行GO富集分析
    run_enrich(
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/input_file/gene_list.txt",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file/go.tsv",
        "Myceliophthora thermophila",
        0.05,
        "GO"
    )

    # 绘制GO富集分析气泡图
    plot_go_chart(
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file",
        "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file/go.tsv",
        width=1000,
        height=800,
        p_adjust=0.05,
        font_size=15,
        chart_num=30,
        chart_size=30,
        pic_type='bar',
        color='rdbu_r',
        funciton_type='All'
    )

    # # 运行KEGG富集分析
    # run_enrich(
    #     "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis",
    #     "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/input_file/gene_list.txt",
    #     "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file/kegg.tsv",
    #     "Myceliophthora thermophila",
    #     0.05,
    #     "KEGG"
    # )

    # # 绘制KEGG富集分析气泡图
    # plot_kegg_chart(
    #     "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file",
    #     "/Users/dongjiacheng/Desktop/Github/omic_analysis/enrichment_analysis/output_file/kegg4.tsv",
    #     width=1000,
    #     height=800,
    #     p_adjust=0.05,
    #     font_size=15,
    #     chart_num=30,
    #     chart_size=30,
    #     color='rdbu_r',
    #     pic_type='bubble',
    #     funciton_type='All'
    # )