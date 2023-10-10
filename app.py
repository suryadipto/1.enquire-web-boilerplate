from api_entrance_point import api_entrance_point
from api_entrance_point2 import api_entrance_point2
from api_entrance_point_scanpy_plots import api_entrance_point_scanpy_plots
from robust_bias_aware.data.study_bias_scores.update_study_bias_scores import update_study_bias_scores
from robust_bias_aware.data.networks.update_networks import update_networks
import flask
from flask import Flask, render_template, request, redirect, url_for
from celery import Celery
from flask_sqlalchemy import SQLAlchemy
import json
import pandas as pd
import json
import shortuuid
from flask_apscheduler import APScheduler
import os
import time
import anndata
import pickle
import numpy as np
import itertools

########################## SET MODE: ############################################
mode=1 # Manually set: 0 (runtime testing) or 1 (actually running the app)
#################################################################################

app = Flask(__name__)
app.config["CELERY_BROKER_URL"] = "redis://localhost:6379"
celery = Celery(app.name, broker=app.config["CELERY_BROKER_URL"])
celery.conf.update(app.config)
####################### Links: #######################
host_url = os.environ.get('HOST', '127.0.0.1:5000')
# robust_home_url = f'http://{host_url}'
# robust_about_url = f'http://{host_url}/robust_about'
# robust_documentation_url = f'http://{host_url}/robust_documentation'
# run_robust_url = f'http://{host_url}/run_robust'

# host_url = os.environ.get('HOST', '0.0.0.0')
# robust_home_url = f'https://{host_url}'
# robust_about_url = f'https://{host_url}/robust_about'
# robust_documentation_url = f'http://{host_url}/robust_documentation'
# run_robust_url = f'https://{host_url}/run_robust'

######################################################


# --------------------------------------------------------------
# data_uploaded=0

fileFormatIncorrect=0
not_enough_data=0

app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///db/robust.db'
app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
app.config['SECRET_KEY'] = 'secret'

@app.route('/run_scAnalyzer', methods=['POST'])
def run_scAnalyzer():
    # return str(error_Flag)
    if flask.request.method == 'POST':
        
        #####################################
        # GET = 0; POST = 1.
        # # empty_upload = 0; non_empty_upload (path exists, path does not exist) = 1.1, 1.0; non_empty_upload (path exists, nuclear_expression.cytoplasm_expression.cell_expression.extracellular_expression) = 1.1.(error_a.error_b.error_c.error_d)
        # # # # uploaded_data_positive = 1, uploaded_data_negative = 0.
        #####################################

        path_to_pe_data = request.form.get("path_to_pe_data")
        alpha = request.form.get("alpha")

        if str(path_to_pe_data)=='':
            error_run_scAnalyzer='1.0'
            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, data_uploaded=data_uploaded)
            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer)
            return render_template('empty_pe_data_path_error.html')
        else:
            if os.path.exists(str(path_to_pe_data)):
                if not str(path_to_pe_data).endswith('/'):
                    path_to_pe_data=str(path_to_pe_data)+'/'
                path_to_relativeAreaWise_nuclear_expression_count=str(path_to_pe_data)+'df___relativeAreaWise_nuclear_expression_count.csv'
                path_to_relativeAreaWise_cytoplasm_expression_count=str(path_to_pe_data)+'df___relativeAreaWise_cytoplasm_expression_count.csv'
                path_to_relativeAreaWise_cell_expression_count=str(path_to_pe_data)+'df___relativeAreaWise_cell_expression_count.csv'
                path_to_relativeAreaWise_extracellular_expression_count=str(path_to_pe_data)+'df___relativeAreaWise_extracellular_expression_count.csv'
                ##############################################################################
                if os.path.exists(str(path_to_relativeAreaWise_nuclear_expression_count)):
                    error_a=1
                else:
                    error_a=0
                ##############################################################################
                if os.path.exists(str(path_to_relativeAreaWise_cytoplasm_expression_count)):
                    error_b=1
                else:
                    error_b=0
                ##############################################################################
                if os.path.exists(str(path_to_relativeAreaWise_cell_expression_count)):
                    error_c=1
                else:
                    error_c=0
                ##############################################################################
                if os.path.exists(str(path_to_relativeAreaWise_extracellular_expression_count)):
                    error_d=1
                else:
                    error_d=0
                ##############################################################################

                errors_a_b_c_d = '('+str(error_a)+'.'+str(error_b)+'.'+str(error_c)+'.'+str(error_d)+')'
                if errors_a_b_c_d=='(0.0.0.0)':
                    error_run_scAnalyzer='1.1.(0.0.0.0)'
                    # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, data_uploaded=data_uploaded)
                    # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer)
                    # # return render_template('protein_expression_absent_error.html', error_run_scAnalyzer=error_run_scAnalyzer)
                    # ######################## Testing: ###############################
                    # error_run_scAnalyzer='fileFormatIncorrect'
                    # fileFormatIncorrect=1
                    # return render_template('incorrect_file_format_error.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect)
                    # ######################## Testing: ###############################
                    error_run_scAnalyzer='notEnoughRows'
                    not_enough_data=1
                    # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                    return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                else:
                    if error_a==1:
                        try:
                            path_to_relativeAreaWise_nuclear_expression_count=pd.read_csv('df___relativeAreaWise_nuclear_expression_count.csv')
                        except:
                            error_run_scAnalyzer='fileFormatIncorrect'
                            fileFormatIncorrect=1
                            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect, data_uploaded=data_uploaded)
                            return render_template('incorrect_file_format_error.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect)
                    if error_b==1:
                        try:
                            path_to_relativeAreaWise_cytoplasm_expression_count=pd.read_csv(f'df___relativeAreaWise_cytoplasm_expression_count.csv')
                        except:
                            error_run_scAnalyzer='fileFormatIncorrect'
                            fileFormatIncorrect=2
                            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect, data_uploaded=data_uploaded)
                            return render_template('incorrect_file_format_error.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect)
                    if error_c==1:
                        try:
                            path_to_relativeAreaWise_cell_expression_count=pd.read_csv(f'df___relativeAreaWise_cell_expression_count.csv')
                        except:
                            error_run_scAnalyzer='fileFormatIncorrect'
                            fileFormatIncorrect=3
                            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect, data_uploaded=data_uploaded)
                            return render_template('incorrect_file_format_error.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect)
                    if error_d==1:
                        try:
                            path_to_relativeAreaWise_extracellular_expression_count=pd.read_csv(f'df___relativeAreaWise_extracellular_expression_count.csv')
                        except:
                            error_run_scAnalyzer='fileFormatIncorrect'
                            fileFormatIncorrect=4
                            # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect, data_uploaded=data_uploaded)
                            return render_template('incorrect_file_format_error.html', error_run_scAnalyzer=error_run_scAnalyzer, fileFormatIncorrect=fileFormatIncorrect)

                empty_file_flag=[error_a, error_b, error_c, error_d]

                for i in range(len(empty_file_flag)):
                    if empty_file_flag[i]==0:
                        continue
                    else:
                        if i==0:
                            if path_to_relativeAreaWise_nuclear_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughRows'
                                not_enough_data=1
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                            if path_to_relativeAreaWise_nuclear_expression_count.shape[1] <= 2:
                                error_run_scAnalyzer='notEnoughColumns'
                                not_enough_data=1
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                        elif i==1:
                            if path_to_relativeAreaWise_cytoplasm_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughRows'
                                not_enough_data=2
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                            if path_to_relativeAreaWise_cytoplasm_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughColumns'
                                not_enough_data=2
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                        elif i==2:
                            if path_to_relativeAreaWise_cell_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughRows'
                                not_enough_data=3
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                            if path_to_relativeAreaWise_cell_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughColumns'
                                not_enough_data=3
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                        elif i==3:
                            if path_to_relativeAreaWise_extracellular_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughRows'
                                not_enough_data=4
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                            if path_to_relativeAreaWise_extracellular_expression_count.shape[0] <= 2:
                                error_run_scAnalyzer='notEnoughColumns'
                                not_enough_data=4
                                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data, data_uploaded=data_uploaded)
                                return render_template('insufficient_data_error.html', error_run_scAnalyzer=error_run_scAnalyzer, not_enough_data=not_enough_data)
                # relativeAreaWise_nuclear_expression_count=pd.read_csv(f'df___relativeAreaWise_nuclear_expression_count.csv')
                # relativeAreaWise_cytoplasm_expression_count=pd.read_csv(f'df___relativeAreaWise_cytoplasm_expression_count.csv')
                # relativeAreaWise_cell_expression_count=pd.read_csv(f'df___relativeAreaWise_cell_expression_count.csv')
                # relativeAreaWise_extracellular_expression_count=pd.read_csv(f'df___relativeAreaWise_extracellular_expression_count.csv')
                # centroid_list=pd.read_csv(f'centroid_list.csv')
                return render_template('run_scAnalyzer.html')
            else:
                error_run_scAnalyzer='1.1.0'
                # return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer, data_uploaded=data_uploaded)
                return render_template('index.html', error_run_scAnalyzer=error_run_scAnalyzer)
                

        # # if os.path.exists(seeds):
        # #     seeds = read_terminals(seeds)
        # # else:
        # #     raise ValueError(f'Illegal value {seeds} for parameter "seeds".\n'
        # #                      f'Must be a valid path.')
        
        
        # df___relativeAreaWise_nuclear_expression_count=pd.read_csv(f'df___relativeAreaWise_nuclear_expression_count.csv')
        # df___relativeAreaWise_cytoplasm_expression_count=pd.read_csv(f'df___relativeAreaWise_cytoplasm_expression_count.csv')
        # df___relativeAreaWise_cell_expression_count=pd.read_csv(f'df___relativeAreaWise_cell_expression_count.csv')
        # df___relativeAreaWise_extracellular_expression_count=pd.read_csv(f'df___relativeAreaWise_extracellular_expression_count.csv')



# # def _get_path_to_study_bias_scores(study_bias_scores, namespace):
# #     if isinstance(study_bias_scores, pd.DataFrame):
# #         study_bias_scores.to_csv(f'robust_bias_aware/data/study_bias_scores/{namespace}/custom_study_bias_scores.csv', index=False)
# #         study_bias_scores=f'robust_bias_aware/data/study_bias_scores/{namespace}/custom_study_bias_scores.csv'
# #         return study_bias_scores
# #     else:
# #         if os.path.exists(study_bias_scores):
# #             return study_bias_scores
# #         if study_bias_scores not in ['None','BAIT_USAGE', 'STUDY_ATTENTION']:
# #             warnings.warn(f'Illegal value {study_bias_scores} for parameter "study_bias_scores".\n'
# #                         f'==> Setting parameter "study_bias_scores" to "BAIT_USAGE"')
# #             study_bias_scores = 'BAIT_USAGE'
# #         else:
# #             if study_bias_scores=='None':
# #                 return study_bias_scores
# #     return f'robust_bias_aware/data/study_bias_scores/{namespace}/{study_bias_scores}.csv'



#         if patient_id=="20751_TCL":
#             path="./spatial_proteomics/feature_extraction/"
#             # # df___extracellular_expression_count, df___relativeAreaWise_extracellular_expression_count, df___cytoplasm_expression_count, df___relativeAreaWise_cytoplasm_expression_count, df___nuclear_expression_count, df___relativeAreaWise_nuclear_expression_count, df___cell_expression_count, df___relativeAreaWise_cell_expression_count, df___extracellular_expression_count_binary, df___relativeAreaWise_extracellular_expression_count_binary, df___cytoplasm_expression_count_binary, df___relativeAreaWise_cytoplasm_expression_count_binary, df___nuclear_expression_count_binary, df___relativeAreaWise_nuclear_expression_count_binary, df___cell_expression_count_binary, df___relativeAreaWise_cell_expression_count_binary=api_entrance_point2(path)
#             # dataframes, adata_objects=api_entrance_point2(path)
#             # # return str(1)
#             # # return os.getcwd()
#             try:
#                 df___relativeAreaWise_cell_expression_count, adata___relativeAreaWise_cell_expression_count, df___relativeAreaWise_nuclear_expression_count, adata___relativeAreaWise_nuclear_expression_count, df___relativeAreaWise_cytoplasm_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, df___relativeAreaWise_extracellular_expression_count, adata___relativeAreaWise_extracellular_expression_count=api_entrance_point2(path)
                
#                 error_flag=0
#                 return render_template('run_scAnalyzer.html', error_flag=error_flag)
#             except:
#                 error_flag=1
#                 return render_template('run_scAnalyzer.html', error_flag=error_flag)
#     return render_template('run_scAnalyzer_get_error.html', error_flag=error_flag)



@app.route('/scAnalyzer_show_gene_expression', methods = ['GET', 'POST'])
def scAnalyzer_show_gene_expression():
    if request.method=='POST':
        gene_symbols = request.form.get("textbox_gene_symbols")
        gene_symbols_list=split_seeds(gene_symbols)
        gene_symbols_list_=gene_symbols_list.copy()

        relativeAreaWise_nuclear_expression_count=pd.read_csv(f'df___relativeAreaWise_nuclear_expression_count.csv')
        relativeAreaWise_cytoplasm_expression_count=pd.read_csv(f'df___relativeAreaWise_cytoplasm_expression_count.csv')
        relativeAreaWise_cell_expression_count=pd.read_csv(f'df___relativeAreaWise_cell_expression_count.csv')
        relativeAreaWise_extracellular_expression_count=pd.read_csv(f'df___relativeAreaWise_extracellular_expression_count.csv')
        centroid_list=pd.read_csv(f'centroid_list.csv')
        centroidList_x=centroid_list['centroid_x'].values.tolist()
        centroidList_y=centroid_list['centroid_y'].values.tolist()
        cell_area_distributions=pd.read_csv('demo_cell_areas.csv')
        nucleusArea_list=cell_area_distributions['nucleus_area'].values.tolist()
        cytoplasmArea_list=cell_area_distributions['cytoplasm_area'].values.tolist()
        extracellularArea_list=cell_area_distributions['extracellular_area'].values.tolist()
        # cell_area_list=cell_area_distributions.values.tolist()
        cell_area_list=[nucleusArea_list, cytoplasmArea_list, extracellularArea_list]
        

        col_names=relativeAreaWise_nuclear_expression_count.columns.tolist()

        cell_names_column_name=col_names[0]
        relativeAreaWise_nuclear_expression_count.rename(columns = {cell_names_column_name:'Cell_names'}, inplace = True)

        cell_names_list=relativeAreaWise_nuclear_expression_count.Cell_names.values.tolist()

        # # return cell_names_list[0]
        
        # display_column_names=gene_symbols_list.append(cell_names_column_name)
        # display_column_names=gene_symbols_list.append('Cell_names')
        DISPLAY_relativeAreaWise_nuclear_expression_count=relativeAreaWise_nuclear_expression_count.loc[:, gene_symbols_list]

        # return str(max(centroidList_x))

        cellNames_coordinatesX_dict=dict(zip(cell_names_list, centroidList_x))
        cellNames_coordinatesY_dict=dict(zip(cell_names_list, centroidList_y))

        ###################

        two_element_combinations=list(itertools.combinations(gene_symbols_list,2))
        two_element_combinations_str=[]
        for i in range(len(two_element_combinations)):
            two_element_combinations_str.append(' '.join(list(two_element_combinations[i])))

        # return (two_element_combinations_str[0])

        three_element_combinations=list(itertools.combinations(gene_symbols_list,3))
        four_element_combinations=list(itertools.combinations(gene_symbols_list,4))

        ###################


        # cellAreas_dict=dict(zip(cell_names_list, cell_area_list))

        RETURN_LIST=[]

        return_dict_template=['x', 'y', 'genes_expressed', 'gene_expression', 'cell_name', 'CELL_AREA_LIST_DICT']

        cell_area_list_pie_chart_template=['name', 'y']

        CELL_AREA_LIST_DICT=[]

        for i in range(len(cell_area_list[0])):
            cell_area_list_dict=[]
            # return cell_area_list[i][0]
            cell_area_list_dict.append(dict(zip(cell_area_list_pie_chart_template, ['nuclear area', cell_area_list[0][i]])))
            cell_area_list_dict.append(dict(zip(cell_area_list_pie_chart_template, ['cytoplasmic area', cell_area_list[1][i]])))
            cell_area_list_dict.append(dict(zip(cell_area_list_pie_chart_template, ['extracellular area', cell_area_list[2][i]])))
            # return cell_area_list_dict[2]
            CELL_AREA_LIST_DICT.append(cell_area_list_dict)
        
        # return CELL_AREA_LIST_DICT[2570][2]
        
        # return str(len(CELL_AREA_LIST_DICT))
        
        # return CELL_AREA_LIST_DICT[0][2]



        count=0
        for index, row in DISPLAY_relativeAreaWise_nuclear_expression_count.iterrows():
            # return DISPLAY_relativeAreaWise_nuclear_expression_count.Cell_names.iloc[index]
            # # return_list=[]
            # # print(row['c1'], row['c2'])
            genes_expressed=[]
            gene_expression=[]
            flag=1
            for i in gene_symbols_list_:
                # return row[i]
                if row[i]>0.5:
                    genes_expressed.append(i)
                    gene_expression.append(row[i])
                    flag=0
            if flag==0:
                return_list=[cellNames_coordinatesX_dict[DISPLAY_relativeAreaWise_nuclear_expression_count.Cell_names.iloc[index]], cellNames_coordinatesY_dict[DISPLAY_relativeAreaWise_nuclear_expression_count.Cell_names.iloc[index]], genes_expressed, gene_expression, DISPLAY_relativeAreaWise_nuclear_expression_count.Cell_names.iloc[index], CELL_AREA_LIST_DICT[count]]
                return_dict=dict(zip(return_dict_template, return_list))
                RETURN_LIST.append(return_dict)
            count+=1
        
        # return RETURN_LIST[0]
        
        return render_template('scAnalyzer_show_gene_expression.html', RETURN_LIST=RETURN_LIST)
                


# def getCombinations(seq):
#     combinations = list()
#     for i in range(0,len(seq)):
#         for j in range(i+1,len(seq)):
#             combinations.append([seq[i],seq[j]])
#     return combinations


        # return str(relativeAreaWise_nuclear_expression_count.columns.tolist()[1])

        
        # # cell_names=df.Courses.values.tolist()

        # cellNames_coordinates_dict=centroid_list.tolist()



@app.route('/scAnalyzer', methods = ['GET', 'POST'])
def scAnalyzer():
    error_flag=404
    if request.method == 'POST':
        option = request.form.get('analysis-option-selected')
        if option=="scanpy":
            error_flag=1
            return render_template('scanpy.html', error_flag=error_flag)
        elif option=="squidpy":
            error_flag=2
            return render_template('squidpy.html', error_flag=error_flag)
        elif option=="scViz":
            error_flag=3
            return render_template('scViz.html', error_flag=error_flag)
        elif option=="scML":
            error_flag=4
            return render_template('scML.html', error_flag=error_flag)
        elif option=="":
            error_flag=5
            return render_template('run_scAnalyzer.html', error_flag=error_flag)
    return render_template('scAnalyzer_get_error2.html')


@app.route('/scanpy',  methods = ['GET', 'POST'])
def scanpy():
    # date = request.args.get('date', None)
    error_flag=404
    # if request.method == 'POST':
    patient_id=20751
    with open(f'Protein_Expression_{patient_id}.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        Protein_Expression_, df___extracellular_expression_count, df___relativeAreaWise_extracellular_expression_count, df___cytoplasm_expression_count, df___relativeAreaWise_cytoplasm_expression_count, df___nuclear_expression_count, df___relativeAreaWise_nuclear_expression_count, df___cell_expression_count, df___relativeAreaWise_cell_expression_count, df___extracellular_expression_count_binary, df___relativeAreaWise_extracellular_expression_count_binary, df___cytoplasm_expression_count_binary, df___relativeAreaWise_cytoplasm_expression_count_binary, df___nuclear_expression_count_binary, df___relativeAreaWise_nuclear_expression_count_binary, df___cell_expression_count_binary, df___relativeAreaWise_cell_expression_count_binary, adata___extracellular_expression_count, adata___relativeAreaWise_extracellular_expression_count, adata___cytoplasm_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___nuclear_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___cell_expression_count, adata___relativeAreaWise_cell_expression_count, adata___extracellular_expression_count_binary, adata___relativeAreaWise_extracellular_expression_count_binary, adata___cytoplasm_expression_count_binary, adata___relativeAreaWise_cytoplasm_expression_count_binary, adata___nuclear_expression_count_binary, adata___relativeAreaWise_nuclear_expression_count_binary, adata___cell_expression_count_binary, adata___relativeAreaWise_cell_expression_count_binary = pickle.load(f)
    adata___nuclear_expression_count=api_entrance_point_scanpy_plots(adata___relativeAreaWise_extracellular_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___relativeAreaWise_cell_expression_count)

    return render_template('scanpy.html', error_flag=error_flag)
    # return render_template('scanpy_get_error.html')

@app.route('/scanpy_direct',  methods = ['GET', 'POST'])
def scanpy_direct():
    error_flag = request.args.get('error_flag', None)
    patient_id=request.args.get('patient_id', None)
    with open(f'Protein_Expression_{patient_id}.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
        Protein_Expression_, df___extracellular_expression_count, df___relativeAreaWise_extracellular_expression_count, df___cytoplasm_expression_count, df___relativeAreaWise_cytoplasm_expression_count, df___nuclear_expression_count, df___relativeAreaWise_nuclear_expression_count, df___cell_expression_count, df___relativeAreaWise_cell_expression_count, df___extracellular_expression_count_binary, df___relativeAreaWise_extracellular_expression_count_binary, df___cytoplasm_expression_count_binary, df___relativeAreaWise_cytoplasm_expression_count_binary, df___nuclear_expression_count_binary, df___relativeAreaWise_nuclear_expression_count_binary, df___cell_expression_count_binary, df___relativeAreaWise_cell_expression_count_binary, adata___extracellular_expression_count, adata___relativeAreaWise_extracellular_expression_count, adata___cytoplasm_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___nuclear_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___cell_expression_count, adata___relativeAreaWise_cell_expression_count, adata___extracellular_expression_count_binary, adata___relativeAreaWise_extracellular_expression_count_binary, adata___cytoplasm_expression_count_binary, adata___relativeAreaWise_cytoplasm_expression_count_binary, adata___nuclear_expression_count_binary, adata___relativeAreaWise_nuclear_expression_count_binary, adata___cell_expression_count_binary, adata___relativeAreaWise_cell_expression_count_binary = pickle.load(f)
    adata___nuclear_expression_count=api_entrance_point_scanpy_plots(adata___relativeAreaWise_extracellular_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___relativeAreaWise_cell_expression_count)
    return render_template('scanpy_direct.html', error_flag=error_flag)
    # return render_template('scanpy_get_error.html')

# @app.route('/home_', methods = ['GET', 'POST'])
# def home_():
#     if request.method == 'POST':
#         date = request.form.get('date')
#         return redirect(url_for('booking', date=date))
#     return render_template('main/index.html')


# @app.route('/booking')
# def booking():
#     date = request.args.get('date', None)
#     return render_template('main/booking.html', date=date)

@app.route('/scanpy_direct_link',  methods = ['GET', 'POST'])
def scanpy_direct_link():
    # date = request.args.get('date', None)
    error_flag=200
    # if request.method == 'POST':
    patient_id=20751
    return redirect(url_for('scanpy_direct', error_flag=error_flag, patient_id=patient_id))
    # # with open(f'Protein_Expression_{patient_id}.pkl', 'rb') as f:  # Python 3: open(..., 'rb')
    # #     Protein_Expression_, df___extracellular_expression_count, df___relativeAreaWise_extracellular_expression_count, df___cytoplasm_expression_count, df___relativeAreaWise_cytoplasm_expression_count, df___nuclear_expression_count, df___relativeAreaWise_nuclear_expression_count, df___cell_expression_count, df___relativeAreaWise_cell_expression_count, df___extracellular_expression_count_binary, df___relativeAreaWise_extracellular_expression_count_binary, df___cytoplasm_expression_count_binary, df___relativeAreaWise_cytoplasm_expression_count_binary, df___nuclear_expression_count_binary, df___relativeAreaWise_nuclear_expression_count_binary, df___cell_expression_count_binary, df___relativeAreaWise_cell_expression_count_binary, adata___extracellular_expression_count, adata___relativeAreaWise_extracellular_expression_count, adata___cytoplasm_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___nuclear_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___cell_expression_count, adata___relativeAreaWise_cell_expression_count, adata___extracellular_expression_count_binary, adata___relativeAreaWise_extracellular_expression_count_binary, adata___cytoplasm_expression_count_binary, adata___relativeAreaWise_cytoplasm_expression_count_binary, adata___nuclear_expression_count_binary, adata___relativeAreaWise_nuclear_expression_count_binary, adata___cell_expression_count_binary, adata___relativeAreaWise_cell_expression_count_binary = pickle.load(f)
    # # adata___nuclear_expression_count=api_entrance_point_scanpy_plots(adata___relativeAreaWise_extracellular_expression_count, adata___relativeAreaWise_cytoplasm_expression_count, adata___relativeAreaWise_nuclear_expression_count, adata___relativeAreaWise_cell_expression_count)

    # return render_template('scanpy.html', error_flag=error_flag)
    # # return render_template('scanpy_get_error.html')


@app.route('/', methods=['GET', 'POST'])
def index():
    return render_template('index.html')


@app.route('/scAnalyzer_about', methods=['GET'])
def robust_about():
    return render_template('scAnalyzer_about.html')


@app.route('/scAnalyzer_documentation', methods=['GET'])
def robust_documentation():
    return render_template('scAnalyzer_documentation.html')


@app.route('/run_robust', methods=['POST', 'GET'])
def run_robust():
    return render_template('run_robust.html')


@app.route('/results', methods=['POST', 'GET'])
def results():
    if flask.request.method == 'POST':
        custom_id = _generate_custom_id()
        NETWORK, NAMESPACE, STUDY_BIAS_SCORE = _initialize_dropdown_params()
        study_bias_score_data = 'BAIT_USAGE'
        uploaded_network = ''
        namespace, alpha, beta, n, tau, study_bias_score, gamma, path_to_graph, uploaded_network, seeds = _initialize_input_params(
            NETWORK, NAMESPACE, STUDY_BIAS_SCORE)
        # return str(uploaded_network)
        is_graphml = False
        in_built_network = request.form.get("network_selection")
        ppi_network_contents_df = pd.DataFrame()
        provided_network, ppi_network_contents_df, is_graphml = _get_network_contents(is_graphml, in_built_network,
                                                                                      ppi_network_contents_df,
                                                                                      uploaded_network, NETWORK)
        custom_studybiasdata_input_df = pd.DataFrame()
        study_bias_score_data, custom_studybiasdata_input_df = _get_study_bias_data_contents(
            custom_studybiasdata_input_df, study_bias_score)
        error_statement = _empty_seeds_error(seeds)
        if not error_statement == 'None':
            return render_template('run_robust.html', error=error_statement)
        # return(str(ppi_network_contents_df))
        ppi_network_contents_df, error_statement = _network_error(in_built_network, is_graphml, ppi_network_contents_df)
        if not error_statement == 'None':
            return render_template('run_robust.html', error=error_statement)
        custom_studybiasdata_input_df, error_statement = _custom_study_bias_data_error(study_bias_score,
                                                                                       custom_studybiasdata_input_df)
        if not error_statement == 'None':
            return render_template('run_robust.html', error=error_statement)
        
        seeds_list=split_seeds(seeds)
        seeds_list_alphabetical=arrange_alphabetical(seeds_list)
        seeds_str_alphabetical=list_to_str(seeds_list_alphabetical)

        if path_to_graph in ['BioGRID', 'APID', 'STRING'] and study_bias_score in ['No', 'BAIT_USAGE', 'STUDY_ATTENTION']:
            param_str=path_to_graph+seeds_str_alphabetical+namespace+str(alpha)+str(beta)+str(n)+str(n)+str(tau)+str(study_bias_score)+str(gamma)
            found_added_task = _is_saved_results(param_str)
            # # if not found_added_task==0:
                # # # return found_added_task
                # # return redirect(f'/saved_results/{found_added_task}')

        input_dict = _make_input_dict(path_to_graph, seeds, namespace, alpha, beta, n, tau, study_bias_score,
                                      study_bias_score_data, gamma, in_built_network, provided_network, is_graphml, param_str)
        if mode==1:
            celery_running_task.delay(custom_id, input_dict)
            record_taskAdded = Task_Added(custom_id)
            db.session.add(record_taskAdded)
            db.session.commit()
            return render_template('running_celery_task.html', custom_id=custom_id)
        elif mode==0:
            robust_run_time, total_celery_running_task_time=celery_running_task(custom_id, input_dict)
            return str(robust_run_time)
    else:
        return render_template('results_get_error.html')

def list_to_str(list_):
    str_ =" ".join(str(item) for item in list_)
    return str_


def split_seeds(seeds):
    seeds_=str(seeds)
    seeds_ = seeds_.split()
    # print(seeds_)
    return seeds_

def arrange_alphabetical(seeds_list):
    seeds_list.sort()
    return seeds_list

@app.route('/saved_results/<int:saved_id>', methods=['POST', 'GET'])
def retrieve(saved_id):
    found_added_task = _is_added_task(saved_id)
    if found_added_task == 1:
        found_done_task = 0
        done_tasks = Robust.query.all()
        for done_task in done_tasks:
            if str(done_task.custom_id) == str(saved_id):
                retrievedRecord = Robust.query.get(done_task.id)
                found_done_task = 1
        if found_done_task == 1:
            custom_id, path_to_graph, seeds, namespace, alpha, beta, n, tau, study_bias_score, study_bias_score_data, gamma, in_built_network, provided_network, is_graphml, nodeData_str, edgeDataSrc_str, edgeDataDest_str, is_seed_str, parameter_str = query_Robust_database(
                retrievedRecord)
            input_network = _check_input_network(provided_network)
            if nodeData_str == "":
                return render_template('empty_output_returned.html', retrievedRecord=retrievedRecord,
                                       input_network=input_network)
            input_dict = _make_input_dict(path_to_graph, seeds, namespace, alpha, beta, n, tau, study_bias_score,
                                          study_bias_score_data, gamma, in_built_network, provided_network, is_graphml, parameter_str)
            _nodes = _convert_comma_separated_str_to_list(nodeData_str)
            src = _convert_comma_separated_str_to_list(edgeDataSrc_str)
            dest = _convert_comma_separated_str_to_list(edgeDataDest_str)
            _edges = zip(src, dest)
            _is_seed = _convert_comma_separated_str_to_list(is_seed_str)
            is_seed_int = _convert_strList_to_intList(_is_seed)
            node_data = _make_node_data(_nodes, is_seed_int)
            edge_data = _make_edge_data(_edges)
            outputData_dict = _make_dict(node_data, edge_data)
            OutputData_json = _convert_dict_to_json(outputData_dict)
            accessLink = _make_access_link(custom_id)
            input_seeds = _split_data_to_list(seeds)
            return render_template('saved_results.html', retrievedRecord=retrievedRecord, input_dict=input_dict,
                                   OutputData_json=OutputData_json, namespace=namespace, accessLink=accessLink,
                                   input_network=input_network, input_seeds=input_seeds, n=n)
        else:
            return render_template('running_celery_task.html', custom_id=saved_id)
    else:
        return render_template('no_such_task_exists.html', custom_id=saved_id)


def _initialize_input_params(NETWORK, NAMESPACE, STUDY_BIAS_SCORE):
    try:
        namespace = NAMESPACE[int(request.form.get("namespace"))]  # dropdown list
    except:
        namespace = 'GENE_SYMBOL'
    try:
        alpha = float(request.form.get('alpha'))  # number field
    except:
        alpha = 0.25
    try:
        beta = float(request.form.get('beta'))  # number field
    except:
        beta = 0.9
    try:
        n = int(request.form.get('n'))  # number field
    except:
        n = 30
    try:
        tau = float(request.form.get('tau'))  # number field
    except:
        tau = 0.1
    try:
        study_bias_score = STUDY_BIAS_SCORE[int(request.form.get("study_bias_score"))]  # dropdown list
    except:
        study_bias_score = 'BAIT_USAGE'
    try:
        gamma = float(request.form.get('gamma'))  # number field
    except:
        gamma = 1.0
    try:
        path_to_graph = NETWORK[int(request.form.get('inbuilt_network_selection'))]  # number field
    except:
        path_to_graph = 'BioGRID'
    try:
        uploaded_network = str(request.form.get("uploaded_ppi_network_filename"))
    except:
        pass
    seeds = ''
    try:
        seeds = request.form.get("textbox_seeds")
    except:
        pass
    return namespace, alpha, beta, n, tau, study_bias_score, gamma, path_to_graph, uploaded_network, seeds


def _convert_list_to_comma_separated_string(_list):
    _str = ','.join(_list)
    return _str


def query_Robust_database(retrievedRecord):
    custom_id = retrievedRecord.custom_id
    path_to_graph = retrievedRecord.path_to_graph
    seeds = retrievedRecord.seeds
    namespace = retrievedRecord.namespace
    alpha = retrievedRecord.alpha
    beta = retrievedRecord.beta
    n = retrievedRecord.n
    tau = retrievedRecord.tau
    study_bias_score = retrievedRecord.study_bias_score
    study_bias_score_data = retrievedRecord.study_bias_score_data
    gamma = retrievedRecord.gamma
    in_built_network = retrievedRecord.in_built_network
    provided_network = retrievedRecord.provided_network
    is_graphml = retrievedRecord.is_graphml
    nodeData_str = retrievedRecord.nodeData_str
    edgeDataSrc_str = retrievedRecord.edgeDataSrc_str
    edgeDataDest_str = retrievedRecord.edgeDataDest_str
    is_seed_str = retrievedRecord.is_seed_str
    parameter_str = retrievedRecord.parameter_str
    return custom_id, path_to_graph, seeds, namespace, alpha, beta, n, tau, study_bias_score, study_bias_score_data, gamma, in_built_network, provided_network, is_graphml, nodeData_str, edgeDataSrc_str, edgeDataDest_str, is_seed_str, parameter_str

def _is_saved_results(param_str):
    found_saved_task=0
    saved_tasks=Robust.query.all()
    for saved_task in saved_tasks:
        if str(saved_task.parameter_str)==str(param_str):
            found_saved_task=saved_task.custom_id
    return found_saved_task

def _is_added_task(saved_id):
    found_added_task = 0
    added_tasks = Task_Added.query.all()
    for added_task in added_tasks:
        if str(added_task.custom_id) == str(saved_id):
            # retrievedRecord_addedTask = Task_Added.query.get(added_task.id)
            found_added_task = 1
    return found_added_task


def _get_network_contents(is_graphml, in_built_network, ppi_network_contents_df, uploaded_network, NETWORK):
    if in_built_network == "Yes":
        try:
            provided_network = NETWORK[int(request.form.get("inbuilt_network_selection"))]
        except:
            provided_network = 'BioGRID'
    elif in_built_network == "No":
        if not uploaded_network.endswith('.graphml'):
            try:
                provided_network = request.form.get("network_contents")
                data2 = list(map(lambda x: x.split(' '), provided_network.split("\r\n")))
                ppi_network_contents_df = pd.DataFrame(data2[1:], columns=data2[0])
            except:
                pass
        elif uploaded_network.endswith('.graphml'):
            is_graphml = True
            try:
                provided_network = request.form.get("network_contents")
                data2 = list(map(lambda x: x.split(' '), provided_network.split("\r\n")))
                ppi_network_contents_df = pd.DataFrame(data2[1:], columns=data2[0])
            except:
                pass
    return provided_network, ppi_network_contents_df, is_graphml


def _make_input_dict(path_to_graph, seeds, namespace, alpha, beta, n, tau, study_bias_score, study_bias_score_data,
                     gamma, in_built_network, provided_network, is_graphml, param_str):
    input_dict = {
        "path_to_graph": path_to_graph,
        "seeds": seeds,
        "namespace": namespace,
        "alpha": alpha,
        "beta": beta,
        "n": n,
        "tau": tau,
        "study_bias_score": study_bias_score,
        "study_bias_score_data": study_bias_score_data,
        "gamma": gamma,
        "in_built_network": in_built_network,
        "provided_network": provided_network,
        "is_graphml": is_graphml,
        "param_str": param_str
    }
    return input_dict


def _custom_study_bias_data_error(study_bias_score, custom_studybiasdata_input_df):
    if study_bias_score == 'CUSTOM':
        if custom_studybiasdata_input_df.empty:
            error_statement = 'Custom study bias data has to be uploaded!'
        else:
            numRows_df = custom_studybiasdata_input_df.shape[0]
            if numRows_df == 0:
                error_statement = 'Custom study bias data with zero rows uploaded. Please add atleast one row excluding the column headers.'
            else:
                if custom_studybiasdata_input_df.shape[1] < 2:
                    error_statement = 'Custom study bias data with less than two columns uploaded. Please add two columns.'
                elif custom_studybiasdata_input_df.shape[1] >= 2:
                    custom_studybiasdata_input_df = custom_studybiasdata_input_df.iloc[:, :2]
                    error_statement = 'None'
    else:
        error_statement = 'None'
    return custom_studybiasdata_input_df, error_statement


def _network_error(in_built_network, is_graphml, ppi_network_contents_df):
    if in_built_network == "No":
        if not is_graphml == True:
            if ppi_network_contents_df.empty:
                error_statement = 'Custom network has to be uploaded!'
            else:
                numRows_df2 = ppi_network_contents_df.shape[0]
                if numRows_df2 == 0:
                    error_statement = 'Custom network with zero rows uploaded. Please add atleast one row excluding the column headers.'
                else:
                    if ppi_network_contents_df.shape[1] < 2:
                        error_statement = 'Custom network with less than two columns uploaded. Please add two columns.'
                    elif ppi_network_contents_df.shape[1] >= 2:
                        # Custom network with more than two columns uploaded. First two columns retained:
                        ppi_network_contents_df = ppi_network_contents_df.iloc[:, :2]
                        error_statement = 'None'
        else:
            error_statement = 'None'
    else:
        error_statement = 'None'
    return ppi_network_contents_df, error_statement


def _empty_seeds_error(seeds):
    if seeds == '' or seeds == None:
        error_statement = 'Seeds cannot be empty'
    else:
        error_statement = 'None'
    return error_statement


def _get_study_bias_data_contents(custom_studybiasdata_input_df, study_bias_score):
    if study_bias_score == 'CUSTOM':
        try:
            study_bias_score_data = request.form.get("custom_studybiasdata_contents_textbox")
            data = list(map(lambda x: x.split(' '), study_bias_score_data.split("\r\n")))
            custom_studybiasdata_input_df = pd.DataFrame(data[1:], columns=data[0])
        except:
            pass
    else:
        study_bias_score_data = study_bias_score
    return study_bias_score_data, custom_studybiasdata_input_df


def _initialize_dropdown_params():
    network = ['BioGRID', 'APID', 'STRING']
    namespace = ['GENE_SYMBOL', 'ENTREZ', 'ENSEMBL', 'UNIPROT']
    study_bias_score = ['No', 'BAIT_USAGE', 'STUDY_ATTENTION', 'CUSTOM']
    return network, namespace, study_bias_score


def _generate_custom_id():
    s = shortuuid.ShortUUID(alphabet="0123456789")
    custom_id = s.random(length=5)
    return custom_id


def _convert_strList_to_intList(str_list):
    intList = []
    for i in str_list:
        intList.append(int(i))
    return intList


def _convert_comma_separated_str_to_list(str_data):
    list_data = str_data.split(",")
    return list_data


def _make_node_data(_nodes, is_seed_int):
    node_data = []
    for i in range(len(_nodes)):
        if is_seed_int[i] == 1:
            node_dict = {"id": _nodes[i], "group": "important"}
        else:
            node_dict = {"id": _nodes[i], "group": "gene"}
        node_data.append(node_dict)
    return node_data


def _make_edge_data(_edges):
    edge_data = []
    for i, j in _edges:
        edge_dict = {"from": i, "to": j, "group": "default"}
        edge_data.append(edge_dict)
    return edge_data


def _make_dict(node_data, edge_data):
    outputData_dict = {"nodes": node_data, "edges": edge_data}
    return outputData_dict


def _convert_dict_to_json(outputData_dict):
    OutputData_json = json.dumps(outputData_dict)
    return OutputData_json


def _make_access_link(id):
    accessLink = f'{host_url}/saved_results/' + str(id)
    # accessLink='127.0.0.1:5000/saved_results/'+str(id)
    return accessLink


def _check_input_network(provided_network):
    if provided_network in ['BioGRID', 'APID', 'STRING']:
        input_network = provided_network
    else:
        input_network = 'custom'
    return input_network


def _split_data_to_list(data):
    str_data = str(data)
    list_data = str_data.split()
    return list_data


if __name__ == '__main__':
    app.run(debug=os.environ.get('DEBUG', '1') == '1', host='0.0.0.0')
