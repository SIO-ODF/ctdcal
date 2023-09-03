#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import matplotlib.pyplot as plt
import pandas as pd


class Cast(object):
    """
    Simple class to store cast details.
    """
    def __init__(self, cast_no):
        self.cast_no = cast_no
        self.t_raw = pd.read_csv('../data/reft/%s_reft.csv' % self.cast_no)
        self.btl = pd.read_pickle('../data/bottle/%s_btl_mean.pkl' % self.cast_no)
        cols = ['CTDPRS', 'CTDTMP1', 'CTDTMP2', 'btl_fire_num', 'T90']
        data = pd.merge(self.btl, self.t_raw, how='left', on='btl_fire_num')
        self.data = data[cols]

    def flag_temp(self):
        # quick and dirty pre-fit flagging...
        self.data['DIFF1'] = self.data['T90'] - self.data['CTDTMP1']
        self.data['DIFF2'] = self.data['T90'] - self.data['CTDTMP2']
        self.data['flag_1'] = False
        self.data.loc[(self.data['CTDPRS'] <= 500) & (abs(self.data['DIFF1'] >= 0.02)), 'flag_1'] = True
        self.data.loc[(self.data['CTDPRS'] <= 1000) & (self.data['CTDPRS'] > 500) & (abs(self.data['DIFF1'] >= 0.01)), 'flag_1'] = True
        self.data.loc[(self.data['CTDPRS'] <= 2000) & (self.data['CTDPRS'] > 1000) & (abs(self.data['DIFF1'] >= 0.005)), 'flag_1'] = True
        self.data.loc[(self.data['CTDPRS'] > 2000) & (abs(self.data['DIFF1'] >= 0.0025)), 'flag_1'] = True

    def qc_plot_temp(self):
        """
        Plot primary, secondary and reference temperatures. Label bottle stops and indicate
        any identified by flagging.
        """
        fig, ax = plt.subplots(figsize=(6, 8))
        ax.plot(self.data['T90'], self.data['CTDPRS'], 'go', label='RefT')
        for x, y, l, f in zip(self.data['T90'], self.data['CTDPRS'], self.data['btl_fire_num'], self.data['flag_1']):
            if f:
                ax.annotate(int(l), (x, y), xytext=(5, -10), textcoords='offset points', color='red')
            else:
                ax.annotate(int(l), (x, y), xytext=(5, -10), textcoords='offset points')
        ax.plot(self.data['CTDTMP1'], self.data['CTDPRS'], 'k-', label='CTD1')
        ax.plot(self.data['CTDTMP2'], self.data['CTDPRS'], 'r:', label='CTD2')
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
        cast = Cast(cast_no)
        cast.flag_temp()
        cast.qc_plot_temp()


if __name__ == '__main__':
    main()
