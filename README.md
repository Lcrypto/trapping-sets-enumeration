# trapping-sets-enumeration
Enumerating of Trapping sets using Cole Impulse tree method "C. A. Cole, S. G. Wilson, E. K. Hall and T. R. Giallorenzi, “A general
method for finding low error rates of LDPC codes, ” submitted to IEEE Trans. on Inform. Theory, June 2006." https://arxiv.org/abs/cs/0605051
Impulse tree have size (3,5), for bigger tree it take too lot of time and require use multi-core CPU or GPU implementation.
To use first convert from sparse (alist) or quasi-cyclic (qc) representation to graph represetation using sparse2graph matlab script after run enumerating trapping sets. Source code shall be upload after publishing paper about probabalistical relaxation of Cole's method.


Example 1 (QC representation):
H=qc2sparse('tanner.qc'); 
sparse2graph(H, 'tanner.graph');


TS_enum.exe -maxit 25 -fast  -qc 31 -x  2  tanner.graph tanner.trap

Example 2 (sparse represenation):
H=alist2sparse('PEGirReg252x504'); 
sparse2graph(H, 'PEGirReg252x504.graph');


TS_enum.exe -maxit 25 -fast   PEGirReg252x504.graph PEGirReg252x504.trap
