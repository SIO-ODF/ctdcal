#!/usr/bin/env python
# -*- coding: utf-8 -*-


import argparse
import matplotlib.pyplot as plt

from common import Cast, Coefficients

_parameter_names_rinko = [
        'salinity',
        'analog_oxygen',
        'analog_temperature',
        'oxygen temperature',
        'oxygen saturation'
    ]


class Rinko(Cast):
    """
    Simple class to store cast details.
    """
    def __init__(self, cast_no):
        self.cast_no = cast_no
        super().__init__()

    def parse_data(self):
        self.data['pressure'] = self.raw['CTDPRS']
        self.data['salinity'] = self.raw['CTDSAL']
        self.data['analog_oxygen'] = self.raw['U_DEF_poly1']
        self.data['analog_temperature'] = self.raw['U_DEF_poly2']

    def load_calibration(self, cal_file):
        cal = Coefficients(cal_file)
        cal.load_csv()
        self.coeffs = cal.coeffs

    def proc_ox(self):
        self.data['oxygen temperature'] = self._proc_t()
        self.data['oxygen saturation'] = self._proc_ox_ptcomp()

    def _proc_t(self):
        """
        Optode temperature used to correct oxygen saturation. From Rinko III manual.
        """
        v = self.data['analog_temperature']
        return self.coeffs['CC_ta'] + self.coeffs['CC_tb'] * v + self.coeffs['CC_tc'] * v ** 2 + self.coeffs['CC_td'] * v ** 3

    def _proc_ox_ptcomp(self):
        """
        Temperature and pressure compensated oxygen saturation (%), from Rinko III manual.
        """
        v = self.data['analog_oxygen']
        tval = (self.data['oxygen temperature'] - 25)
        do = ((self.coeffs['CC_oa'] / (1 + (self.coeffs['CC_od'] * tval) + (self.coeffs['CC_of'] * tval**2)))
              +
              (self.coeffs['CC_ob'] / (v * (1 + (self.coeffs['CC_od'] * tval) + (self.coeffs['CC_of'] * tval**2)) + self.coeffs['CC_oc']))
              )
        p = (self.coeffs['CC_og'] + self.coeffs['CC_oh'] * do)
        pd = p * (1 + self.coeffs['CC_oe'] * self.data['pressure'] * 0.01)
        return pd

    def qc_plot_ox(self):
        """
        Plot Rinko oxygen saturation vs depth for qc.
        """
        fig, ax = plt.subplots(figsize=(6, 8))
        ax.plot(self.data['oxygen saturation'], self.data['pressure'], 'b-', label='Oxy Sat')
        ax.set_title('Cast %s' % self.cast_no)
        ax.set_xlabel('Oxygen Saturation (%)')
        ax.set_ylabel('Pressure (dbar)')
        ax.legend(loc='lower right')
        ax.invert_yaxis()
        ax.grid()
        plt.show()


def main():
    parser = argparse.ArgumentParser(description="""Plot Rinko oxygen data
                                     for quick and dirty QC.""")
    parser.add_argument('cal_file', type=str,
                        help='calibration file')
    parser.add_argument('casts', metavar='SSSCC', type=str, nargs='+',
                        help='casts to plot')

    args = parser.parse_args()
    if not args.cal_file:
        parser.print_help()
        return
    cal_file = args.cal_file
    casts = [arg for arg in args.casts]

    for cast_no in casts:
        cast = Rinko(cast_no)
        cast.parse_data()
        cast.load_calibration(cal_file)
        cast.proc_ox()
        cast.qc_plot_ox()


if __name__ == '__main__':
    main()
