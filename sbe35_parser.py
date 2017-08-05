'''
WE'LL DO IT LIVE

Turn a screendump of SBE-35RT to a csv file that can be easily loaded.

Files are expected to be as clean as possible (for now). Need to fix later.

Joseph Gum
December 1, 2016
'''
import sys
import io
import csv

def parse(f, f_out):
    #second half of final output
    array = []
    reader = csv.reader(f, delimiter=' ')

    #break row down to a standardized size
    for row in reader:
        row2 = []
        for x in row:
            if x != '':
                row2.append(x)

        #if not a fixed size, assume it's a comment and to be ignored
        if len(row2) != 17:
            continue

        #start to build list to be written out
        row3 = []
        #hardcoded madness because it probably won't change anytime soon.
        row3.append(row2[0])
        #concatenate date/time, convert to ISO8601 later
        row3.append(row2[1] + ' ' + row2[2] + ' ' + row2[3] + ' ' + row2[4])
        row3.append(row2[7])
        row3.append(row2[10])
        row3.append(row2[13])
        row3.append(row2[16])
        array.append(row3)

    #create header to be plased in
    header = ['index_memory', 'datetime', 'bottle_number', 'diff', 'raw_value', 'T90']

    #fix here to output according to new name
    with open(f_out + '.btl_t', 'w') as csvf:
        csvw = csv.writer(csvf)
        csvw.writerow(header)
        for x in array:
            csvw.writerow(x)

def main(argv):
    with open(argv[0], 'r') as filename:
        parse(filename, argv[0])

if __name__ == '__main__':
    main(sys.argv[1:])
