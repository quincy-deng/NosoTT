import glob
import sys
import argparse
import os
import subprocess

suffixs = ["1.fastq.gz","1.fq.gz"]

def main():
    args = get_argument_parser().parse_args()

    check_for_program()

    fastqs = find_raw_file(args)
    names = find_name(fastqs)
    trimmed_fastqs = trim_fastq(args,fastqs,names)
    assemble(trimmed_fastqs,names)

def assemble(trimmed_fastqs,names):
    x = r"spades.py -o ${path}/assembled_data/$j -1 ${path}/trimmed_data/${j}_forward_paired.fq.gz -2 ${path}/trimmed_data/${j}_reverse_paired.fq.gz -k 41,49,57,65,77,85,93 --cov-cutoff auto"
    os.mkdir("assembled_fasta")
    for (f1,f2),name in zip(trimmed_fastqs,names):
        spades_command = ['spades.py','-o', "assembled_fasta"+name, '-1', f1, '-2', f2, '-k', 
                          '41,49,57,65,77,85,93','--cov-cutoff', 'auto']
        spades_process = subprocess.Popen(spades_command,stdout=subprocess.PIPE, 
                                          stderr=subprocess.PIPE)
        _, err = spades_process.communicate()
        
        if err:
            quit_with_error("trimmomatic running error:\n" + convert_bytes_to_str(err))

def trim_fastq(args,fastqs,names):
    fls = []
    for (f1,f2),name in zip(fastqs,names):
        path = os.path.abspath(args.input)
        f_pair, f_unpair, r_pair, r_unpair =path, path, path, path
        os.mkdir("trimmed_data")
        f_pair = r"trimmed_data/{}_forward_paired.fq.gz".format(name)
        f_unpair = r"trimmed_data/{}_forward_unpaired.fq.gz".format(name)
        r_pair = r"trimmed_data/{}_reverse_paired.fq.gz".format(name)
        r_unpair = r"trimmed_data/{}_reverse_unpaired.fq.gz".format(name)
        
        trimmomatic_command = ['trimmomatic', 'PE', '-phred33', f1, f2, f_pair, f_unpair, r_pair,r_unpair,'LEADING:3', 'TRAILING:3', 'SLIDINGWINDOW:4:15','MINLEN:36']
        trimmomatic_process = subprocess.Popen(trimmomatic_command, stdout=subprocess.PIPE,
                                               stderr=subprocess.PIPE)
        _, err = trimmomatic_process.communicate()

        if err:
            quit_with_error("trimmomatic running error:\n" + convert_bytes_to_str(err))
        fls.append([f_pair,r_pair])
    return fls

def find_name(fastqs):
    names = []
    for f1,_ in fastqs:
        f1 = os.path.basename(f1)
        temp =  f1.split("_")
        if len(temp) > 2:
            names.append("_".join(temp[1:-1]))
        else:
            quit_with_error("fastq file name is not nomal")
    return names  

def find_raw_file(args):
    fastq_fls = []
    if not os.path.exists(args.input):
        quit_with_error("directory not correct ")
    path = os.path.abspath(args.input)
    for suffix in suffixs:
        suffix2 = suffix.replace('1','2')
        fls = glob.glob("{}/*{}".format(path,suffix))
        if len(fls) > 0:
            for i in fls:
                i2 = i.replace(suffix,suffix2)
                if os.path.exists(i2):
                    fastq_fls.append([i,i2])
    if fastq_fls:
        return fastq_fls
    else:
        quit_with_error("not found fastq file,please check directory")

def get_argument_parser():
    """Specifies the command line arguments required by the script."""
    parser = argparse.ArgumentParser(description='trimmomatic',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    add_arguments_to_parser(parser)
    return parser

def add_arguments_to_parser(parser):
    parser.add_argument('-i','--input', type=str, required=False, default='.',
                        help='fastq directory')

def quit_with_error(message):
    """Displays the given message and ends the program's execution."""
    print('Error:', message, file=sys.stderr)
    sys.exit(1)

def check_for_program():
    """Checks to make sure the required BLAST+ tools are available."""
    if not find_program('trimmomatic'):
        quit_with_error('could not find trimmomatic')
    if not find_program('spades.py'):
        quit_with_error('could not find spades')

def find_program(name):
    """Checks to see if a program exists."""
    process = subprocess.Popen(['which', name], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    out, err = process.communicate()
    return bool(out) and not bool(err)

def convert_bytes_to_str(bytes_or_str):
    """
    This function is for both Python2 and Python3. If the input is a str, it just returns that
    same str. If not, it assumes its bytes (and we're in Python3) and it returns it as a str.
    """
    if isinstance(bytes_or_str, str):
        return bytes_or_str
    else:
        return bytes_or_str.decode()


if __name__ == '__main__':
    main()
