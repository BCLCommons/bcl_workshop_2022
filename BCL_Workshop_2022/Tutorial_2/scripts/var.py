#! /usr/bin/env python
import numpy as np
from argparse import ArgumentParser

##################################################################################################################
#
# Computes the row variance for all columns separated by a user-specified delimiter
#
##################################################################################################################


# Compute the variance with gratuitous use of numpy
def main():
	# Read in command-line arguments
	args = cmdlineparse()
	data = np.loadtxt(args.input, dtype=np.float64, delimiter=args.delimiter)
	var = np.var(np.transpose(data),axis=1)
	np.savetxt(args.output, var)
	return 0


# Parse command-line options
def cmdlineparse():
	parser = ArgumentParser(description="command line arguments")
	parser.add_argument("-input", dest="input", required=True,
						help="Input file on which the row-wise variance will be computed", metavar="<input>")
	parser.add_argument("-output", dest="output", required=True, help="Name of output file", metavar="<output file>")
	parser.add_argument("-delimiter", dest="delimiter", default=",", required=False, help="Set column delimiter",
						metavar="<delimiter>")
	args = parser.parse_args()
	return args


# Execute
if __name__ == '__main__':
	main()
