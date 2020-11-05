import os,subprocess,base64
import dash_html_components as html
import dash_table
import pandas as pd


def copy2_to_datadir(contents,filename):
    content_type, content_string = contents.split(',')
    decoded = base64.b64decode(content_string)
    try:
        c_str = decoded.decode('utf-8')
        with open('fasta_data/{}'.format(filename),'w') as o:
            o.write(c_str)
    except Exception as e:
        print(e)
        return html.Div([
            '上传过程中出现错误.'
        ])
    return html.Div('上传 {} 成功'.format(filename))

def execute_shell_command(commandline):
    ret = subprocess.Popen(commandline,shell = True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    result = ret.communicate()
    ret.wait()
    if ret.returncode == 0:
        return 'ok'
    else:
        stderr = str(result[0],encoding = 'utf8')
        return stderr

def execute_kaptive(filenames,ref,temp_id):
    command_kaptive = 'conda run -n abricate ~/bin/kaptive.py --threads 1 -a {} -k {} -o output/kaptive_{}'.format(filenames,ref,temp_id)
    res = execute_shell_command(command_kaptive)
    if res != 'ok':
        return [html.Strong('Excute Kaptive Failed!!!'),html.P(res)]
    df = pd.read_table('output/kaptive_{}_table.txt'.format(temp_id))
    tab = dash_table.DataTable(
        id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
        style_table={'overflowX': 'auto'},
        style_cell={
            'height': 'auto',
            'lineHeight': '15px',
            # all three widths are needed
            'minWidth': '50px', 'width': '150px', 'maxWidth': '250px',
            'whiteSpace': 'normal'
        },
         
    )
    return html.Div([html.Div('Known Loci 结果:'),tab])

def execute_abricate(filenames,temp_id):
    commandline_abricate = 'conda run -n abricate abricate {} >output/results_{}.tab;conda run -n abricate abricate --summary output/results_{}.tab >output/summary_{}.tab'.format(filenames,temp_id,temp_id,temp_id)
    res = execute_shell_command(commandline_abricate)
    if res != 'ok':
        return [html.Strong('Excute abricate Failed!!!'),html.P(res)]
    df = pd.read_table('output/summary_{}.tab'.format(temp_id))
    df['#FILE']=[os.path.basename(i) for i in df['#FILE']]
    tab = dash_table.DataTable(
        id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
    )
    return html.Div([html.Div('耐药基因/毒力基因结果:'),tab])

def execute_mlst(filenames,temp_id):
    command_mlst = 'conda run -n mlst mlst -q {} >output/{}.tsv'.format(filenames,temp_id)
    res = execute_shell_command(command_mlst)
    if res != 'ok':
        return [html.Strong('Excute MLST Failed!!!'),html.P(res)]
    df = pd.read_table('output/{}.tsv'.format(temp_id),names=['File','scheme','ST']+['gene_'+str(i) for i in range(1,8)])
    df['File']=[os.path.basename(i) for i in df['File']]
    tab = dash_table.DataTable(
        id='table',
        columns=[{"name": i, "id": i} for i in df.columns],
        data=df.to_dict('records'),
    )
    return html.Div([html.Div('ST分型结果:'),tab])