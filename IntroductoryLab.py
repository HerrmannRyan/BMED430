import numpy as np
import pandas as pandas

a = 3.45
b = 12.525
sigfigs1 = 3
sigfigs2 = 5

L_input = ['%.*g' % (sigfigs1, a), '%.*g' %(sigfigs2,b)]
L_isigs = ['3','5']
L_iheader = ['a','b']

# lists for math operations

L_results = []
L_rindex = []
L_rsigs = []
#
#
# case 1 = b-a
#
sigfigs = 3