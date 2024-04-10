
ContStruc_REBiStar= [ 1 1 1 1 1
                      1 1 0 0 1
                      1 1 1 0 0
                      1 0 1 1 0
                      1 0 0 1 1];

ContStruc_OPT= [ 1 1 0 0 1
                 0 1 1 0 0
                 1 1 1 1 1
                 0 0 0 1 0  %%PROBLEM: optimization outputs a matrix with no inputs on the 4th system! 
                 1 0 0 1 1];

ContStruc_cyc= [ 1 1 0 0 1
                 1 1 1 0 0
                 0 1 1 1 0
                 0 0 1 1 1   
                 1 0 0 1 1];

ContStruc_Pij= [ 1 1 0 0 0
                 1 1 1 0 1
                 0 1 1 1 0
                 0 0 1 1 1   
                 0 1 0 1 1];

ContStruc_NotPij= [ 1 0 1 1 1
                      0 1 0 1 0
                      1 0 1 0 1
                      1 1 0 1 0   
                      1 0 1 0 1]; 
