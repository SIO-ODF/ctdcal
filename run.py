import sys
import os
import argparse
import converter_scaffolding as cnv

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Main function of the script should it be run as a stand-alone utility.
# -------------------------------------------------------------------------------------
def main(argv):
    parser = argparse.ArgumentParser(description='Convert SBE CTD data into engineering units')
    parser.add_argument('hexFile', metavar='hex File', help='the .hex data file to process')
    parser.add_argument('xmlconFile', metavar='XMLCON File', help='the .XMLCON data file to process')
    parser.add_argument('-o', metavar='output file', dest='outFile', help='name and location of output file')
    parser.add_argument('-d', '--debug', action='store_true', help=' display debug messages')

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    if not os.path.isfile(args.hexFile):
        errPrint('ERROR: Input hex file:', args.hexFile, 'not found\n')
        sys.exit(1)

    if not os.path.isfile(args.xmlconFile):
        errPrint('ERROR: Input xmlcon file:', args.xmlconFile, 'not found\n')
        sys.exit(1)

    if args.outFile:
        if os.path.isdir(os.path.dirname(args.outFile)):
            outputFile = args.outFile
        else:
            errPrint('ERROR: Unable to save output to:', args.outFile, '  Destination is unreachable.\n')
            sys.exit(1)
    else:
	    outputFile = os.path.splitext(args.xmlconFile)[0] + '.converted'

    output = cnv.cnv_handler_2(args.hexFile, args.xmlconFile, DEBUG)

    if output == None:
        errPrint('No data returned from converter')
    else:
        debugPrint('Saving output to:', outputFile)	
        with open(outputFile, 'w') as f:
            f.write(output)

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])