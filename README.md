# trapping-sets-enumeration

Upd: 20 January 2022. Gift for the New Year)  Dr. Chad A. Cole kindly share his Matlab source code for Trapping Sets analysis, famous research in detail described at report
Chad A. Cole, Eric. K. Hall A General Method for Finding Low Error Rates of LDPC Codes https://arxiv.org/abs/cs/0605051 and in paper Enumerating of Trapping sets using Cole Impulse tree method "C. A. Cole, S. G. Wilson, E. K. Hall and T. R. Giallorenzi, “A general method for finding low error rates of LDPC codes, ” submitted to IEEE Trans. on Inform. Theory, June 2006.". 


Thank you very much Dr. Chad A. Cole!



Impulse tree have size (3,5), for bigger tree it take too lot of time and require use multi-core CPU or GPU implementation.
To use first convert from sparse (alist) or quasi-cyclic (qc) representation to graph represetation using sparse2graph matlab script after run enumerating trapping sets. Source code shall be upload after publishing paper about probabalistical relaxation of Cole's method. Due to some limitation relaxation version and related articles shall be published after 12.10.2025. However, when I have free time, I will publish an article and source files to Trapping set search in a less effective method based on the ideas outlined in my Phd Thesis, "Topologically Driven Methods for Construction Of Multi-Edge Type (Multigraph with nodes puncturing) Quasi-Cyclic Low-density Parity-check Codes for Wireless Channel, WDM Long-Haul and Archival Holographic Memory", https://arxiv.org/abs/2011.14753. 



Example 1 (QC representation):
H=qc2sparse('tanner.qc'); 
sparse2graph(H, 'tanner.graph');


TS_enum.exe -maxit 25 -fast  -qc 31 -x  2  tanner.graph tanner.trap


Example 2 (sparse represenation):
H=alist2sparse('PEGirReg252x504'); 
sparse2graph(H, 'PEGirReg252x504.graph');


TS_enum.exe -maxit 25 -fast   PEGirReg252x504.graph PEGirReg252x504.trap


Example 3 from 


Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004.


load 9_3.mat


sparse2graph(a, '9_3.graph')


run bat file '9_3.cmd'
