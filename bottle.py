import sys
import os
import argparse
import converter_scaffolding as cnv
import bottle_lib as btl

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
    parser = argparse.ArgumentParser(description='Convert .converted file to bottle file')
    parser.add_argument('convertedFile', metavar='converted File', help='the .converted data file to process')
    parser.add_argument('-o', metavar='output file', dest='outFile', help='name and location of output file')
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')

    args = parser.parse_args()
    if args.debug:
        global DEBUG
        DEBUG = True
        debugPrint("Running in debug mode")

    if not os.path.isfile(args.convertedFile):
        errPrint('ERROR: Input .converted file:', args.convertedFile, ' not found\n')
        sys.exit(1)

    if args.outFile:
        if os.path.isdir(os.path.dirname(args.outFile)):
            outputFile = args.outFile
        else:
            errPrint('ERROR: Unable to save output to:', args.outFile, '  Destination is unreachable.\n')
            sys.exit(1)
    else:
        outputFile = os.path.splitext(args.convertedFile)[0] + '.bottle'

    output = btl.handler(args.convertedFile, DEBUG)

    if output == None:
        errPrint('No data returned from converter')
    else:
        debugPrint('Saving output to:', outputFile)
        with open(outputFile, 'w') as f:
            f.write(output)

        debugPrint('Success!')

DEBUG = False

def debugPrint(*args, **kwargs):
    if DEBUG:
        errPrint(*args, **kwargs)

def errPrint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

# -------------------------------------------------------------------------------------
# Required python code for running the script as a stand-alone utility
# -------------------------------------------------------------------------------------
if __name__ == '__main__':
	main(sys.argv[1:])
