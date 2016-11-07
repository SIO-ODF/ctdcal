import sys
import os
import argparse
import converter_scaffolding as cnv
import web_viewer as wv

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
    parser.add_argument('-d', '--debug', action='store_true', help='display debug messages')
    parser.add_argument('-w', metavar='dest dir', dest='wvDir', nargs='?', const='default', help='build web viewer, optionally specify the destination directory for the webviewer files')

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

    if args.wvDir:
        debugPrint('Adding webviewer')
        webView = wv.WebViewer()
        if DEBUG:
            webView._setDebug()

        debugPrint('\tSetting webviewer directory to: ', end='')
        if args.wvDir == 'default':
            wvDir = os.path.join(os.path.dirname(outputFile),'webviewer')
            debugPrint(wvDir + '... ', end='')
        else:
            wvDir = args.wvDir
            debugPrint(wvDir + '... ', end='')

        try:
            os.mkdir(wvDir, 0o755);
        except FileExistsError as exception:
            pass
        except PermissionError as exception:
            errPrint('ERROR: do not have write permission at destination directory')
            sys.exit(1)

        debugPrint('Success!')

        webView._setWebviewerFolder(wvDir)

        debugPrint('\tCreating webviewer directory structure... ', end='')

        if webView._buildScaffolding():
            debugPrint('Success!')
        else:
            debugPrint('ERROR!')
            errPrint('Unable to created webviewer directories')
            sys.exit(1)

        debugPrint("Building data file for webviewer... ", end='')
        output = webView._buildData(outputFile)
        debugPrint('Success!')

        debugPrint('Saving webviewer data to file... ', end='')
        webView._saveData(output)
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
