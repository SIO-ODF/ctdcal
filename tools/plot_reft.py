#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pandas as pd
import matplotlib.pyplot as plt

from common import Cast


class Reft(Cast):
    """
    Simple class to store cast details.
    """
    def __init__(self, cast_no):
        self.cast_no = cast_no
        self.reft = pd.read_csv('../data/reft/%s_reft.csv' % cast_no)
        super().__init__()

    def parse_data(self):
        self.data['pressure'] = self.btl['CTDPRS']
        self.data['temperature 1'] = self.btl['CTDTMP1']
        self.data['temperature 2'] = self.btl['CTDTMP2']
        data = pd.merge(self.btl, self.reft, how='left', on='btl_fire_num').set_index(pd.Index([n for n in range(1, 37)]))
        self.data['bottle number'] = data['btl_fire_num'].astype(int)
        self.data['T ref'] = data['T90']

    def _flag_temp(self):
        # quick and dirty pre-fit flagging...
        self.data['DIFF1'] = self.data['T ref'] - self.data['temperature 1']
        self.data['DIFF2'] = self.data['T ref'] - self.data['temperature 2']
        self.data['Flag'] = False
        self.data.loc[(self.data['pressure'] <= 500) & (abs(self.data['DIFF1'] >= 0.02)), 'Flag'] = True
        self.data.loc[(self.data['pressure'] <= 1000) & (self.data['pressure'] > 500) & (abs(self.data['DIFF1'] >= 0.01)), 'Flag'] = True
        self.data.loc[(self.data['pressure'] <= 2000) & (self.data['pressure'] > 1000) & (abs(self.data['DIFF1'] >= 0.005)), 'Flag'] = True
        self.data.loc[(self.data['pressure'] > 2000) & (abs(self.data['pressure'] >= 0.0025)), 'Flag'] = True

    def qc_plot_temp(self):
        """
        Plot primary, secondary and reference temperatures. Label bottle stops and indicate
        any identified by flagging.
        """
        self._flag_temp()
        # data = self.data.loc[self.data['T ref'].
        fig, ax = plt.subplots(figsize=(6, 8))
        ax.plot(self.data['T ref'], self.data['pressure'], 'go', label='T Ref')
        for x, y, l, f in zip(self.data['T ref'], self.data['pressure'], self.data['bottle number'], self.data['Flag']):
            if f:
                ax.annotate(l, (x, y), xytext=(5, -10), textcoords='offset points', color='red')
            else:
                ax.annotate(l, (x, y), xytext=(5, -10), textcoords='offset points')
        ax.plot(self.data['temperature 1'], self.data['pressure'], 'k-', label='CTD 1')
        ax.plot(self.data['temperature 2'], self.data['pressure'], 'r:', label='CTD 2')
        ax.set_title('Cast %s' % self.cast_no)
        ax.set_xlabel('Temperature (T90)')
        ax.set_ylabel('Pressure (dbar)')
        ax.legend(loc='lower right')
        ax.invert_yaxis()
        ax.grid()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="""Plot bottle temperatures
                                     for quick and dirty QC.""")
    parser.add_argument('casts', metavar='SSSCC', type=str, nargs='+',
                        help='casts to plot')

    args = parser.parse_args()
    casts = [arg for arg in args.casts]

    for cast_no in casts:
        cast = Reft(cast_no)
        cast.parse_data()
        cast.qc_plot_temp()


if __name__ == '__main__':
    main()
