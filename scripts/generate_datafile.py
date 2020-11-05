import pandas as pd
import datetime

def main():
    args = get_argument_parser().parse_args()

def generate_csv(casedata):
    dates = ['Admdate','Sampdate','Chgdate']

    # casedata = pd.read_csv(casedate,names=['item','name','Admdate','Sampdate','Chgdate'])
    baseDate = min(casedata['Sampdate'])

    for ind,row in casedata.iterrows():
        bDate = datetime.datetime.strptime(baseDate,"%Y-%m-%d")
        Adate,Sdate,Cdate = [datetime.datetime.strptime(row[i],"%Y-%m-%d") for i in dates]
        casedata['Admdate'][ind], casedata['Sampdate'][ind],casedata['Chgdate'][ind] = \
            (Adate - bDate).days,(Sdate - bDate).days, (Cdate - bDate).days

    dates = pd.DataFrame({'A':casedata['name'],'B':casedata['Sampdate']})
    hosts = pd.DataFrame({'A':casedata['name'],'B':casedata['item']})
    hostTimes = pd.DataFrame({'A':casedata['item'],'B':casedata['Admdate'],'C':casedata['Chgdate']})

    try:
        dates.to_csv('{}dates.csv'.format(basedir),index=False,header=False)
        hosts.to_csv('{}hosts.csv'.format(basedir),index=False,header=False)
        hostTimes.to_csv('{}hostTimes.csv'.format(basedir),index=False,header=False)
    except IOError as e:
        print("Save datefile failed!:%s", e)
    
    print("Generate date file success")

def get_argument_parser():
    """Specifies the command line arguments required by the script."""
    parser = argparse.ArgumentParser(description='trimmomatic',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    add_arguments_to_parser(parser)
    return parser

def add_arguments_to_parser(parser):
    parser.add_argument('-s','--sample', type=list, required=True, help='Samples list')
