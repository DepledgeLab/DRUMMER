import os
import argparse
def create_output(output_dir_path:str,input_file:str,output_dir_name:str,subdir=False):
	"""Takes in an input file and returns an output file
	"""
	# Making the new output directory
	input_file = input_file.split('/')[-1]
# 	print('output_dir_path',output_dir_path)
	output_dir = ''
	if subdir != False:
	    output_dir = output_dir_path + '/' + output_dir_name + '/' + subdir
	else:
		output_dir = output_dir_path + '/' + output_dir_name
	os.makedirs(output_dir, exist_ok = True)
	# Name of the file
	split_input_file = input_file.split('.')
	split_input_file.insert(-1,output_dir_name)
	input_file = '.'.join(split_input_file)
	return output_dir + '/' + input_file
		
if __name__ == '__main__':
	ap = argparse.ArgumentParser(description = 'Takes in the output of bam-readcount \
and returns a text file containing the count of each nucleotide at the position and the \
fraction of the reference nucleotide among all reads.')
	requiredGrp = ap.add_argument_group('required arguments')
	requiredGrp.add_argument("-i","--input", required=True, help="input file location")
	requiredGrp.add_argument("-o","--output", required=True, help="input file location")
	args = vars(ap.parse_args())
	input = args['input']
	output = args['output']
	create_output(input,output)
	
	
