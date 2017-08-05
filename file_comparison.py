import pandas as pd
import logging

logging.basicConfig(filename='errors/compare_summary.txt', format='%(asctime)s:%(levelname)s:%(message)s', level=logging.INFO)
#0 for diffs of columns that stand out, 1 for full summaries and diffs
DEBUG_FLAG = 1
ssscc = ['00101', '00102', '00201', '00202', '00301', '00401', '00402', '00501', '00601', '00602', '00701', '00801', '00901', '01001', '01101', '01201', '01301', '01401', '01402', '01501', '01502', '01601','01602','01701','01702','01801','01802','01901','01904','02001','02002','02101','02201','02301','02401','02402','02501','02601','02602','02801','02901','03001','03101','03201']

def compare_all(ssscc, sequence='pressure', verbosity=0):
    for x in ssscc:
        print('\nCast number: ' + x)
        #Start with time first, then pressure
        if sequence == 'time':
            df1 = pd.read_csv('data/time/' + x + '_time.csv', header = [1,2])
            df2 = pd.read_csv('data/initial/time/' + x + '_time.csv', header = [1,2])
        if sequence == 'pressure':
            df1 = pd.read_csv('data/pressure/' + x + '.csv', header = [12,13])
            df2 = pd.read_csv('data/initial/pressure/' + x + '.csv', header = [12,13])

        #Find the difference of the values.
        #TODO - try to determine in df3 where the biggest discrepancies happened
        try:
            df3 = df1.subtract(df2)
            series_sum = df3.sum()

            #show results to user
            if verbosity >= 1:
                print(series_sum)

            if verbosity >= 0:
                #if the difference of values is greater than 0.0, throw it up for checking
                unchanged_file_counter = 0
                for i, y in enumerate(series_sum):
                    if y != 0.0:
                        print(str(series_sum.index[i]) + ' difference: ' + str(y))
                        unchanged_file_counter += 1
                if unchanged_file_counter == 0:
                    logging.info('Cast {} {} is unchanged after fitting'.format(x,sequence))
        except TypeError:
            print('Not all data points are same type, skip for now')

if __name__ == "__main__":
    compare_all(ssscc, 'time', DEBUG_FLAG)
    print('Start pressure differences')
    compare_all(ssscc,'pressure', DEBUG_FLAG)
