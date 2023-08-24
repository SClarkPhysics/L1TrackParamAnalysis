## Track Parameter Uncertainty Analysis Code

Code corrsponding to Internal Note: ADD IN NUMBER

Code to analyze NTuples from the L1TrackNTupleMaker and produce the Covariance Matrix Elements of the KF Track Parameters. This takes output from the 5 param Kalman Filter and reconstructs the
resolution from the 4 param Kalman Filter. The parametrization has been conducted separately and is based on the track's radius and Hit Pattern. Here, we read in the values from the covariance matrix stored in ParamFiles/radius.csv. 

Run as:
```
python parametrizer.py myInputFile.root outputFileName
```

This will run over myInputFile.root, copy the tree named `eventTree` and add two branches: `param_slope` and `param_pt`. `param_slope` is the covariance matrix element corresponding to that track,
based on it's radius and hit pattern. `param_pt` is the reconstructed pT calculated as  $p_T - slope \times d_0$, where slope is the covariance matrix element. The output is saved in OutputTrees/outputFileName.root

To see an example, run: 
```
python parametrizer.py rootFiles/SingleMu_pt2to100_100k_D49_5Params.root newTree
```

This will make OutputTrees/newTree.root which is a copy of eventTree from the input file plus the two new branches.


Steven Clark
svclark96@gmail.com

