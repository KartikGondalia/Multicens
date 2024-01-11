import os
import re
import time
import datetime
import uuid
import sqlite3
import json
import shutil
import io
import csv
import logging
import itertools
import numpy as np
import pandas as pd
from pathlib import Path
from memory_profiler import profile
from util.ranking import getRanking
from werkzeug.utils import secure_filename
from datetime import datetime, timedelta, timezone
from util.util import allowed_file, create_corr_matrix
from flask import Flask, render_template, request, session, redirect, flash, make_response, send_from_directory, send_file,jsonify

from util.centrality import global_centrality, local_centrality, right_target_global_centrality_t
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.mime.application import MIMEApplication


tissue_names = {
    'adipose_visceral': 'Adipose_Visceral_Omentum',
    'adipose_subutaneous' : 'Adipose_Subcutaneous',
    'adrenal': 'Adrenal_Gland',
    'artery_aorta' : 'Artery_Aorta',
    'artery_coronary' : 'Artery_Coronary',
    'artery_tibial' : 'Artery_Tibial',
    'brain_amygdala' : 'Brain_Amygdala',
    'brain_anterior' : 'Brain_Anterior_cingulate_cortex_BA24',
    'brain_caudate' : 'Brain_Caudate_basal_ganglia',
    'brain_cerebellar' : 'Brain_Cerebellar_Hemisphere',
    'brain_cerebellum' : 'Brain_Cerebellum',
    'brain_cortex' : 'Brain_Cortex',
    'brain_frontal_cortex' : 'Brain_Frontal_Cortex_BA9',
    'brain_hippocampus' : 'Brain_Hippocampus',
    'brain_hypothalamus' : 'Brain_Hypothalamus',
    'brain_nucleus' : 'Brain_Nucleus_accumbens_basal_ganglia',
    'brain_putamen' : 'Brain_Putamen_basal_ganglia',
    'brain_spinal_cord' : 'Brain_Spinal_cord_cervical_c-1',
    'brain_substantia' : 'Brain_Substantia_nigra',
    'breast': 'Breast_Mammary_Tissue',
    'cells_cultured' : 'Cells_Cultured_fibroblasts',
    'cells_evb_transformed' : 'Cells_EBV-transformed_lymphocytes',
    'colon_sigmoid' : 'Colon_Sigmoid',
    'colon_transverse' : 'Colon_Transverse',
    'esophagus_gastroesophageal' : 'Esophagus_Gastroesophageal_Junction',
    'esophagus_mucosa' : 'Esophagus_Mucosa',
    'esophagus_muscularis' : 'Esophagus_Muscularis',
    'heart_atrial': 'Heart_Atrial_Appendage',
    'heart_left_ventricle' : 'Heart_Left_Ventricle',
    'kidney': 'Kidney_Cortex',
    'liver': 'Liver',
    'lung' : 'Lung',
    'salivary_gland' : 'Minor_Salivary_Gland',
    'muscle': 'Muscle_Skeletal',
    'nerve_tibial' : 'Nerve_Tibial',
    'ovary': 'Ovary',
    'pancreas': 'Pancreas',
    'pituitary': 'Pituitary',
    'prostate' : 'Prostate',
    'skin_not_sun' : ' Skin_Not_Sun_Exposed_Suprapubic',
    'skin_sun' : 'Skin_Sun_Exposed_Lower_leg',
    'intestine': 'Small_Intestine_Terminal_Ileum',
    'spleen' : 'Spleen',
    'stomach': 'Stomach',
    'testis' : 'Testis',
    'thyroid': 'Thyroid',
    'uterus': 'Uterus',
    'vagina' : 'Vagina',
    'Blood' : 'Whole_Blood',
    'fp_ctl':'MSBB_BM_10_Control',
    'stg_ctl':'MSBB_BM_22_Control',
    'phg_ctl':'MSBB_BM_36_Control',
    'ifg_ctl':'MSBB_BM_44_Control',
    'fp_ad':'MSBB_BM_10_Definite_AD',
    'stg_ad':'MSBB_BM_22_Definite_AD',
    'phg_ad':'MSBB_BM_36_Definite_AD',
    'ifg_ad':'MSBB_BM_44_Definite_AD'
}

gtissue_list = ['adipose_visceral','adipose_subutaneous','artery_aorta', 'artery_coronary', 'artery_tibial', 'brain_amygdala', 'brain_anterior', 'brain_caudate', 'brain_cerebellar', 'brain_cerebellum', 'brain_cortex', 'brain_frontal_cortex', 'brain_hippocampus', 'brain_hypothalamus', 'brain_nucleus', 'brain_putamen', 'brain_spinal_cord', 'brain_substantia', 'cells_cultured', 'cells_evb_transformed', 'colon_sigmoid', 'colon_transverse', 'esophagus_gastroesophageal', 'esophagus_mucosa', 'esophagus_muscularis', 'heart_atrial', 'heart_left_ventricle', 'kidney', 'liver', 'lung', 'salivary_gland', 'muscle', 'nerve_tibial', 'ovary', 'pancreas', 'pituitary', 'prostate', 'skin_not_sun', 'skin_sun', 'intestine', 'spleen', 'stomach', 'testis', 'thyroid', 'uterus', 'vagina', 'Blood']
mtissue_list = ['FP_CTL','STG_CTL','PHG_CTL','IFG_CTL','FP_AD','STG_AD','PHG_AD','IFG_AD']
mtissue_list = [tissue.upper() for tissue in mtissue_list]
hormone_list = ['adrenaline', 'aldosterone', 'angiotensin', 'cortisol', 'estradiol', 'glucagon', 'insulin', 'leptin', 'norepinephrine', 'progesterone', 'somatotrophin', 'thyroxin', 'vitamind']
measures = {
    'local': 'Local Centrality',
    'global': 'Global Centrality',
    'query': 'Query Set Centrality'
}

dataset_list = ['GTex', 'MSBB', 'Others']
gtex_path = 'dataset/GTEx_Analysis_v8_eQTL_expression_matrices'
gtex_cov_path = 'dataset/GTEx_Analysis_v8_eQTL_covariates'
msbb_path = 'dataset/MSBB'

hormone_list.sort()
# tissue_list.append('other')

app = Flask(__name__)
app.secret_key = '$#ghFkGH56-1aRGFGHGtrfJH'
app.config["SESSION_PERMANENT"] = False
# Example to set MAX_CONTENT_LENGTH in Flask
app.config['MAX_CONTENT_LENGTH'] = 512 * 1024 * 1024  # for 16MB limit
app.config['SESSION_TYPE'] = 'filesystem'
app.config['PERMANENT_SESSION_LIFETIME'] = timedelta(minutes=60)  # Session timeout of 1 hour
session_id = ''



# Session management
@app.before_request
def before_request():
    check_session_expiry()

def check_session_expiry():
    if 'last_activity' in session:
        last_activity_time = session['last_activity']
        session_timeout = app.config['PERMANENT_SESSION_LIFETIME']
        # Check if session has expired
        if datetime.now(timezone.utc) - last_activity_time > session_timeout:            
            shutil.rmtree(os.path.join('uploads', session_id))
            file_path = os.path.join('downloads', session_id)
            if os.path.exists(file_path):
                os.remove(file_path)         
            session.clear()

    # Update last activity time
    session['last_activity'] = datetime.now(timezone.utc)
    
# Home page
@app.route('/')
def index():
    cookies = request.cookies
    print("Cookies:", cookies)
    global session_id
    session_id = str(uuid.uuid4().hex)
    session['id'] = session_id
    return render_template('home.html', gtissue_list=gtissue_list,mtissue_list=mtissue_list, hormone_list=hormone_list, measures=measures, dataset_list=dataset_list)


@app.route('/tool')
def tool():
    return render_template('tool.html')

@app.route('/help')
def help():
    return render_template('help.html')


# Download the Stati Files
@app.route('/download/<filename>', methods=['GET'])
def download_file(filename):
    # Make sure the filename is safe
    filename = secure_filename(filename)
    # Define the full path to the file
    filepath = os.path.join(app.root_path, 'static/sampleData', filename)
    
    # Send the file
    return send_file(filepath, as_attachment=True)

@app.route('/download_csv/<uniqueID>/<type>', methods=['GET'])
def download_csv(uniqueID, type):
    conn = sqlite3.connect('results.db')
    c = conn.cursor()
    c.execute("SELECT common_genes, common_sample FROM results WHERE uniqueID=?", (uniqueID,))
    row = c.fetchone()
    conn.close()

    if row:
        data = row[0] if type == 'common_genes' else row[1]
        data_list = json.loads(data) if data else []

        # Convert data_list to CSV, each item in a new row
        si = io.StringIO()
        cw = csv.writer(si)
        for item in data_list:
            cw.writerow([item])  # Write each item as a new row

        output = make_response(si.getvalue())
        output.headers["Content-Disposition"] = f"attachment; filename={uniqueID}.csv"
        output.headers["Content-type"] = "text/csv"
        return output
    else:
        return "No data available", 404



@app.route('/contact')
def contact():
    return render_template('contact.html')


# Tool v1.0
@app.route('/get_ranking', methods=['POST'])
def get_ranking():
    if request.method == 'POST':
        start_time = time.time()
        if len(request.files)<=0:
            flash('No selected file')
            return redirect(request.url)

        file1 = request.files['input1']
        file2 = request.files['input2']
        source_tissue = request.form.get('source_tissue')
        target_tissue = request.form.get('target_tissue')
        hormone = request.form.get('hormone')
        
        source_path = 'static/uploads/' + uuid.uuid4().hex + Path(file1.filename).suffix
        target_path = 'static/uploads/' + uuid.uuid4().hex + Path(file2.filename).suffix
        file1.save(source_path)
        file2.save(target_path)

        # result = pd.read_csv(source)
        result = pd.DataFrame()
        
        try:
            result = getRanking(source_path, target_path, hormone, source_tissue, target_tissue)
        except:
            pass

        # delete file(s) from the server
        os.remove(source_path)
        os.remove(target_path)
        
        result['rank'] = range(1, len(result.index)+1, 1)
        result.columns = ["Gene Name", "Centrality", "Rank"] 
        
        result_csv = 'static/downloads/' + uuid.uuid4().hex + '.csv'
        
        result.to_csv(result_csv, index=False)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f'{execution_time} min(s)'
        print("Execution time:", time_taken, "\n")

        return render_template('ranking.html', columns=result.columns, rows=list(result.values.tolist()), path=result_csv, time_taken=time_taken)
    
@app.route('/run_sample', methods=['GET'])
def run_sample():
    if request.method == 'GET':
        start_time = time.time()
        source_tissue = 'pancreas'
        target_tissue = 'muscle'
        hormone = 'insulin'
        
        source_path = 'static/downloads/insulin_source_genes_full.csv'
        target_path = 'static/downloads/insulin_target_genes_full.csv'
        
        result = pd.DataFrame()
        
        try:
            result = getRanking(source_path, target_path, hormone, source_tissue, target_tissue)
        except Exception as e:
            print(str(e))
        
        result['rank'] = range(1, len(result.index)+1, 1)
        result.columns = ["Gene Name", "Centrality", "Rank"]   
        
        result_csv = 'static/downloads/' + uuid.uuid4().hex + '.csv'
        
        result.to_csv(result_csv, index=False)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f'{execution_time} min(s)'
        print("Execution time:", time_taken, "\n")

        return render_template('ranking.html', columns=result.columns, rows=list(result.values.tolist()), path=result_csv, time_taken=time_taken)


# Tool v2.0 (beta)
@app.route('/example_run/<uniqueID>', methods=['GET'])
def example_run(uniqueID):
    #print("In side")
    conn = sqlite3.connect('results.db')
    c = conn.cursor()

    # Fetch the data from the database
    c.execute("""SELECT result_data, tNames, submissionTime, expirationTime, csv_path, 
                 time_taken, measure, common_genes, common_sample, query_genes, 
                 no_common_sample, no_common_gene
                 FROM example WHERE uniqueID=?""", (uniqueID,))
    result_data = c.fetchone()
    conn.close()

    if result_data:
        columns = ['Tissue Name', 'Gene Name', 'Centrality', 'Rank']
        result_list = json.loads(result_data[0])
        tNames = result_data[1]
        submissionTime = result_data[2]
        expirationTime = result_data[3]
        path = result_data[4]
        time_taken = result_data[5]
        measure = result_data[6]
        # Capitalize the first letter
        measure = measure[0].upper() + measure[1:]

        common_genes = json.loads(result_data[7]) if result_data[7] else []
        common_sample = json.loads(result_data[8]) if result_data[8] else []
        query_genes = json.loads(result_data[9]) if result_data[9] else []
        no_common_sample = result_data[10]
        no_common_gene = result_data[11]

        # Render the template with the data
        
        return render_template('example.html', columns=columns, 
                                rows=result_list, tNames=tNames, uniqueID=uniqueID,
                                path=path,
                                measure=measure, common_genes=common_genes, 
                                common_sample=common_sample, query_genes=query_genes, 
                                no_common_sample=no_common_sample, no_common_gene=no_common_gene)
    else:
        # Handle the case where no data is found
        return jsonify({"error": "No data found for the provided uniqueID"}), 404


# Create a dictionary to store the results
results = {}

# Create a function to generate a unique ID
def generate_unique_id():
    return str(uuid.uuid4())


# Create a function to store the result for a particular ID
def store_result(id, result):
    results[id] = result

# Create a function to get the result for a particular ID
def get_result(id):
    return results.get(id)

def convert_expiration_time_to_date_time(expiration_time):
    expiration_time_datetime = datetime.datetime.fromtimestamp(expiration_time)
    return expiration_time_datetime.strftime("%Y/%m/%d %H:%M")

def send_email(to_addr, subject, body):
    # Gmail SMTP server settings
    smtp_server = 'smtp.gmail.com'
    smtp_port = 587  # TLS Port

    # Your Gmail account credentials
    gmail_user = 'birdslab.iitm@gmail.com'  # Your Gmail email address
    gmail_password = 'osecabqlhdbrkwbd '  # Your Gmail password or an app-specific password

    # Create a MIMEText object to represent the email body
    msg = MIMEMultipart()
    msg['From'] = gmail_user
    msg['To'] = to_addr
    msg['Subject'] = subject

    # Attach the body of the email
    msg.attach(MIMEText(body, 'html'))  # Assuming the body is HTML

    # Create an SMTP connection and send the email
    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.starttls()  # Upgrade the connection to TLS
        server.login(gmail_user, gmail_password)
        server.sendmail(gmail_user, to_addr, msg.as_string())
        server.close()
        print("Email sent successfully")
    except Exception as e:
        print(f"Error sending email: {str(e)}")
        

def store_initial_entry(uniqueID, tNames, submissionTime, expirationTime, measure, query_genes):
    conn = sqlite3.connect('results.db')
    c = conn.cursor()
    
    c.execute("""INSERT INTO results 
                 (uniqueID, result_data, tNames, submissionTime, expirationTime, csv_path, time_taken, measure, query_genes) 
                 VALUES (?, ?, ?, ?, ?, NULL, NULL, ?, ?)""", 
              (uniqueID, None, tNames, submissionTime, expirationTime, measure, json.dumps(query_genes)))
    conn.commit()
    conn.close()

def update_result_in_db(uniqueID, result, time_taken, result_csv):
    conn = sqlite3.connect('results.db')
    c = conn.cursor()
    
    c.execute("""UPDATE results 
                 SET result_data = ?, time_taken = ?, csv_path = ?
                 WHERE uniqueID = ?""", 
              (json.dumps(result.values.tolist()), time_taken, result_csv, uniqueID))
    conn.commit()
    conn.close()





@app.route('/get_result/<uniqueID>', methods=['GET'])
def display_result(uniqueID):
    conn = sqlite3.connect('results.db')
    c = conn.cursor()
    app.logger.info(f"Fetching result for uniqueID: {uniqueID}")

    try:
        c.execute("""SELECT result_data, tNames, submissionTime, expirationTime, csv_path, 
                     time_taken, measure, common_genes, common_sample, query_genes, no_common_sample, no_common_gene
                     FROM results WHERE uniqueID=?""", (uniqueID,))
        result_data = c.fetchone()
    except sqlite3.Error as e:
        app.logger.error(f"Database error for uniqueID: {uniqueID} - {e}")
        return "Database error", 500
    finally:
        conn.close()

    if not result_data:
        return "No result found for this uniqueID", 404

    if result_data[0] is None:
        # Ensure the format string matches your database's format
        submission_datetime_str = result_data[2]
        # Remove the time zone information for simplicity
        submission_datetime_str = re.sub(r' GMT[+-]\d{4} \(.+\)$', '', submission_datetime_str)

        try:
            submission_datetime = datetime.strptime(submission_datetime_str, '%a %b %d %Y %H:%M:%S')
            # Check if more than 4 hours have elapsed since submission
            if datetime.now() - submission_datetime > timedelta(hours=4):
                # More than 4 hours have passed, return a timeout error
                return "Your result calculation has exceeded the maximum time limit", 400
        except ValueError as e:
            # Error in parsing the submission time
            app.logger.error(f"Error parsing submissionTime for uniqueID: {uniqueID} - {e}")
            return "Error processing submission time", 400

        # Less than 4 hours have passed, continue waiting for the result
        return render_template('waiting.html')
    else:
        columns = ['Tissue Name', 'Gene Name', 'Centrality', 'Rank']
        result_list = json.loads(result_data[0])
        tNames = result_data[1]
        submissionTime = result_data[2]
        expirationTime = result_data[3]
        path = result_data[4]
        time_taken = result_data[5]
        measure = result_data[6]
        # Capitalize the first letter
        measure = measure[0].upper() + measure[1:]

        common_genes = json.loads(result_data[7]) if result_data[7] else []
        common_sample = json.loads(result_data[8]) if result_data[8] else []
        query_genes = json.loads(result_data[9]) if result_data[9] else []
        no_common_sample = result_data[10]
        no_common_gene = result_data[11]

        return render_template('ranking.html', columns=columns, 
                                rows=result_list, tNames=tNames, uniqueID=uniqueID,
                                submissionTime=submissionTime, expirationTime=expirationTime, 
                                time_taken=time_taken, path=path,
                                measure=measure, common_genes=common_genes, 
                                common_sample=common_sample, query_genes=query_genes, 
                                no_common_sample=no_common_sample, no_common_gene=no_common_gene)




# You can uncomment and use these calls after initializing the required parameters
# store_initial_entry(uniqueID, tNames, submissionTime, expirationTime)
# update_result_in_db(uniqueID, result, time_taken, result_csv)
    
def cleanup_expired_results():
    current_time = datetime.now().strftime("%Y/%m/%d %H:%M")
    
    conn = sqlite3.connect('results.db')
    c = conn.cursor()

    # Assuming the expirationTime is in the format "%Y/%m/%d %H:%M"
    c.execute("DELETE FROM results WHERE expirationTime < ?", (current_time,))

    conn.commit()
    conn.close()

# You can run this function periodically (e.g., using a scheduler or cron job).



@app.route('/get_centrality', methods=['POST', 'GET'])
@profile
def get_centrality():
    
    
    if request.method == 'POST':
        
        start_time = time.time()
        # submition_time_datetime = datetime.fromtimestamp(start_time).strftime("%Y/%m/%d %H:%M")

        measure = request.form.get('measure')

        submissionTime = request.form.get('submissionTime')
        expirationTime = request.form.get('expirationTime')
        uniqueID = request.form.get('uniqueID')
        tNames = request.form.get('tNames')
        
        store_initial_entry(uniqueID, tNames, submissionTime, expirationTime, measure, 'Null')
        email = request.form.get('email')
        

        tissues = json.loads(request.form.get('tissues'))

        #print("Emergency")
        #print(tissues)
        # tissues = sorted(tissues, key=lambda t: t['name'])
                # Store the tissue names with unique IDs and expiration times

        
        # Generate a URL for the stored result
        url = f'www.birdslab.in/get_result/{uniqueID}'

        if email:
            table = pd.DataFrame({
                'Unique ID': uniqueID,
                'URL': url,
                'Submission time': submissionTime,
                'Expiration time': expirationTime,
                'Tissue name': tNames,
                'Measure name': measure
            }, index=[0])

            table_html = table.to_html()

            # Add your custom message
            message = """
            Hi, <br>

            Your Job status is given below: <br>
            <br>
            {}<br><br>
            
            if you have any inquiries, please email us at birdslab.iitm@gmail.com <br>
            <br>
            Thanks, <br>
            Birds Lab <br>
            Department of Computer Science and Engineering <br>
            IIT Madras, Chennai <br>
            India - 600036 <br>
            """.format(table_html)

            send_email(email, 'Multicens', message)

        result = pd.DataFrame()
        #print(1)        
        try:
            folder_path = os.path.join('static/sampleData', uniqueID)
            os.makedirs(folder_path, exist_ok=True)
            files = request.files.getlist('files')
            for file in files:
                #print(2)
                if file and allowed_file(file.filename):
                    #print("2a")
                    filename = secure_filename(file.filename)
                    file_path = os.path.join(folder_path, filename)
                    file.save(file_path)
            #print("3")        
            for tissue in tissues:
                #print("3b")
                if tissue['dataset'] == 'GTex' :
                    #print("3aa")
                    tissue_name = tissue_names[tissue['name'].lower()]
                    tissue['exp'] = os.path.join(gtex_path, f'{tissue_name}.v8.normalized_expression.bed.gz')
                    tissue['cov'] = os.path.join(gtex_cov_path,f'{tissue_name}.v8.covariates.txt')
                elif tissue['dataset'] == 'MSBB' :
                    #print("3ab")
                    tissue_name = tissue_names[tissue['name'].lower()]
                    #print(tissue_name)
                    tissue['exp'] = os.path.join(msbb_path, f'{tissue_name}.bed.gz')
                    #print("GK")
                else:
                    #print("3ab")
                    tissue['exp'] = os.path.join(folder_path, tissue['file'])
            #print("You")
            #print(tissues)
            A, gene_count = create_corr_matrix(uniqueID, tissues)
            #print("Hii")
            #print(A)
            #print("got")
            n = len(A)/len(tissues)
            #print("4")
            if measure == 'query':
                #print("4a")              
                query_tissue = request.form.get('query-tissue')
                query_file = request.files.get('query-file')
                query_index = int(request.form.get('tissue-index'))                
                s = int(query_index * n)
                e = int(s + n)
                query_gene_list = pd.read_csv(query_file, header=None)[0].tolist()
                query_gene_list = [f"{gene}.{query_index+1}" for gene in query_gene_list]
                genes = A.columns[s:e]
                common_target_genes = np.intersect1d(genes, query_gene_list)
                genes_indices = [i for i, e in enumerate(A.columns) if e in common_target_genes]
                _, g = right_target_global_centrality_t(A.values, num_layers=len(tissues), target_tissue = query_index, target_gene_indices = genes_indices, start=s, end=e, p=0.9)
                result['Centrality'] = g
            elif measure == 'global':
                #print("4b")
                result['Centrality'] = global_centrality(A.values, len(tissues), p=0.9)
            else:
                #print("4c")               
                result['Centrality'] = local_centrality(A.values, len(tissues), p=0.9)
            name_list = [[tissues[i]['name'] for _ in range(gene_count[i])] for i in range(len(tissues))]
            result['Tissue Name'] =  list(itertools.chain(*name_list))
            result['Gene Name'] = [re.sub(r'\.\d+$', '', s) for s in A.columns]
            #print("5")
            if measure == 'local':
                #print("5a")
                ranking = pd.DataFrame()
                for tissue in tissues:
                    df = result.loc[result['Tissue Name']==tissue['name']]
                    df = df.sort_values(by='Centrality', ascending=False)
                    df['Rank']  = [i+1 for i in range(len(df))]
                    ranking = pd.concat([ranking, df], ignore_index=True)
                result = ranking
            else:
                #print("5b")
                result = result.sort_values(by='Centrality', ascending=False)
                result['Rank'] = [i+1 for i in range(len(result))]
            result = result[['Tissue Name', 'Gene Name', 'Centrality', 'Rank']]
            #print("5c")
        except Exception as e:
            return jsonify({"error": str(e)}), 500
        
        # print(result)
        result_csv = os.path.join('static/sampleData', uniqueID + '.csv')    
        result.to_csv(result_csv)
        result['Centrality'] = np.round(result['Centrality'], 4)
        end_time = time.time()
        execution_time = end_time - start_time
        execution_time = round(execution_time/60, 2)
        # unit = 'sec(s)' if execution_time<1 else 'min(s)'
        # time_taken = f'{execution_time} {unit}'
        time_taken = f'{execution_time} min(s)'
        print("Execution time:", time_taken, "\n")
        #print(result_csv)

        update_result_in_db(uniqueID, result, time_taken, result_csv)
        
                
        return render_template(
    'ranking.html', 
    columns=result.columns, 
    rows=list(result.values.tolist()), 
    path=result_csv, 
    time_taken=time_taken,
    uniqueID=uniqueID,  # Added this
    tNames=tNames,      # Added this
    submissionTime=submissionTime,  # Added this
    expirationTime=expirationTime   # Added this
)



if __name__ == '__main__':
    app.run(debug=True)
    