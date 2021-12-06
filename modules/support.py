import os 
import sys
#from termcolor import colored
def str2bool(v):
#Adapted from https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1','True'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0','False'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def handle_booleans(argument,default_argument):
    if argument == None:
        argument = default_argument
    elif argument.isnumeric() == True:
        return int(argument)
    else:
        argument = str2bool(argument)
    return argument

def which(program):
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file
    return None


def check_samtools():
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    cmd_path=which('samtools')
    if cmd_path:
        print('SAMTOOLS in path')
        return True
    else:
        sys.stdout.write('\nERROR: DRUMMER requires Samtools')
        sys.stdout.write('\n\nPlease install samtools and add it to your path following the instructions at: http://www.htslib.org/download/\n')
        return False
def check_bedtools():
#Adapted from https://github.com/pinellolab/CRISPResso2/blob/master/CRISPResso2/CRISPRessoPooledCORE.py
    cmd_path=which('bedtools')
    if cmd_path:
        print('Bedtools in path\n')
        return True
    else:
        sys.stdout.write('\nERROR: DRUMMER requires Bedtools\n')
        sys.stdout.write('\n\nPlease install Bedtools and add it to your path following the instructions at: https://bedtools.readthedocs.io/en/latest/content/installation.html\n')

        return False
        
def print_logo(mode):
    pooled_string = '''
     ________________________________
    | __   __                 __  __ |
    ||  \ |_/ |  | |\/| |\/| |_  |_/ |
    ||__/ | \ |__| |  | |  | |__ | \ |
    |________________________________|
    \n'''
    print(pooled_string)
    print('RUNNING in {}'.format(mode))






class MyException(Exception):
    pass