# Jets
Jet algorithms are fun!

Comparison with fastjet
0.4ms

Test1 is my dumbest possible approach, using TLorentzVector extensively
850ms

Test2 uses TLorentzVector only temporarily, storing data as doubles
50ms

Test3 uses a multidimensional cache of all kt^2 values
12ms

Test4 uses a flat cache of the best kt^2 pair for each object
9.2ms

Test5 uses arrays rather than vectors, with masking rather than deletion
8.2ms
