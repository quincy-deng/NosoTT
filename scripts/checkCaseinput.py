import re
import datetime
import numpy as np

time = lambda x :datetime.datetime.strptime(x,'%Y-%m-%d').date() 
def filter_col(df):
    itemx = re.compile(r'[0-9a-zA-Z_]+')
    style = [
        {
            'if': {
                'filter_query': '{{item}} = {}'.format(i),
                'column_id': 'item',
            },
            'backgroundColor': '#FF4136',
            'color': 'white'
        } for i in df['item'] if i is None or not itemx.match(str(i))
    ] + [
            {
                'if': {
                    'filter_query': '{{Admdate}} = {} && {{Sampdate}} = {} && {{Chgdate}} = {}'.format(i,j,k),
                    'column_id': ['Admdate','Sampdate','Chgdate'],
                },
                'backgroundColor': '#FF4136',
                # 'color': 'white'
            } for i,j,k in  zip(df['Admdate'],df['Sampdate'],df['Chgdate']) if i is None or time(i) > time(j)  or time(j) > time(k)
        ] + [
                {
                    'if': {
                        'filter_query': '{{{}}} is blank'.format(col),
                        'column_id': col
                    },
                    'backgroundColor': '#FF4136',
                    'color': 'white'
                } for col in df.columns
            ]
    return style
