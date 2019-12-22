import numpy as np
n           = 20
amin        = 1   # 0.1 micron
amax        = 20   # 10 micron
agraini     = np.linspace(amin-1,amax,n+1)
agrain      = 0.5 * ( agraini[1:] + agraini[:-1] )

print(agraini[1:])
print(agraini[:-1])
