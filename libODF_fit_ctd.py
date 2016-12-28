#!/usr/bin/env python
import numpy as np
import conversions as convert
import density_enthalpy_48 as density

def offset(offset, col, inMat):
    """offset column of data 

    Input:
        - inMat, 1d numpy array with return True np.isnans()
    Output:
        - Mat with offset applied to column of data 
    Example:
        >>> outArray = offset(offset, col, inMat) 
    """
    for i in range(0, len(inMat)):
        inMat[col][i] = inMat[col][i] + offset

    return inMat
