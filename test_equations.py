'''Copied and pasted from a ipython workspace, needs to turn into full unittest later'''

import libODF_sbe_equations_dict as sbe_eq
import unittest

print('Temp:')
'''
SBE 3+ temperature sensor, file 030867-120625c.pdf
SN: 0867
Calib date: 25-Jun-2012
'''
CALIB_TEMP_EX_1 = {'G': 4.84983043E-3, 'H' : 6.85991891E-4, 'I' : 3.07523322E-5, 'J' : 2.72771980E-6, 'F0' : 1000.0}
CALIB_TEMP_FREQ_1 = [6227.9391,6587.5975,6586.5119,7113.0853,7670.1725,7669.3400,8255.9574,8873.9097,8872.1426,9522.0840,10203.5835,10918.6965,11666.8669,12449.9590]
CALIB_TEMP_RES_1 = [-1.5075,1.0012,0.9937,4.4943,8.0008,7.9957,11.4971,15.0000,14.9902,18.4944,21.9937,25.4954,28.9939,32.4952]

list_temp = sbe_eq.temp_its90_dict(CALIB_TEMP_EX_1, CALIB_TEMP_FREQ_1)
#for x, y in zip(list_temp, CALIB_TEMP_RES_1):
    #print(x, y)

print('Cond:')

CALIB_COND_EX_1 = {'G' : -4.08627809e+000, 'H' : 4.63770989e-001, 'I' : -8.22760889e-004, 'J' : 6.15523357e-005, 'CPcor' : -9.5700e-008, 'CTcor' : 3.2500e-006}
CALIB_COND_FREQ_1 = [2974.44,3859.84,8327.42,8546.18,10044.17,10407.82,11468.11,11810.09]
CALIB_COND_T_1 = [0.0000,-1.0000,-1.0000,1.0000,15.0000,18.5000,29.0000,32.5000]
CALIB_COND_P_1 = [0.0000,34.6104,34.6104,34.6110,34.6125,34.6125,34.6117,34.6068]
CALIB_COND_RES_1 = [0.00000,2.78953,2.78953,2.96010,4.24921,4.59425,5.67261,6.04354]

list_cond = sbe_eq.cond_dict(CALIB_COND_EX_1,CALIB_COND_FREQ_1,CALIB_COND_T_1,CALIB_COND_P_1)
for x, y in zip(list_cond, CALIB_COND_RES_1):
    print(x, y)


'''
class SBEEquationsTest(unittest.TestCase):

    def test_temp_1(self):
        self.assertEqual()
'''
