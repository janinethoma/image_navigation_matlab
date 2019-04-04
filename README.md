 CVPR'2019 submission, matlab demo as part of supplementary material.  
 Paper ID: 4087
 Paper title: "Image-based Navigation using Visual Features and Map"
========================================================================================

This code was tested in MATLAB R2017b and with mosek 8.1.0.64 (https://www.mosek.com/downloads/).

Please, run the code in 'demo.m' to view our demo.
You may need to adjust the path to your mosek directory inside 'demo.m'.


This demo does the following using the concepts introduced in our paper:
1) Find landmarks in a subsequence of the Oxford Robotcar run from  2015-10-29 12:18:17
2) Match a short query sequence from 2014-11-18 13:20:12
Both, reference and query, had to be kept short due to supplementary
material data limit.

If you do not have mosek installed, you can have a look at the saved
figures in the results folder instead.

The produced outputs are:
- Scatter plot of original reference and query sequences
- Topology of reference sequence used for finding landmarks with network flow
- Selected landmarks
- Accuracy vs. distance plot of the final matching

