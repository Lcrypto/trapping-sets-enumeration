# Importance Sampling Cole's  LDPC trapping sets enumeration and weighing method
Subfolder ["Article An Modified Cole's Importance Sampling Method For Low Error Floor QC-LDPC Codes Construction"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction) contained parity-check matrixes, EMD Spectrum Enumerator and Trapping set enumerator (w\ot ordering) in zip compresed files for five QC-LDPC codes constructed for article An Modified Cole's Importance Sampling Method For Low Error Floor QC-LDPC Codes Construction.

![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction/SZ_Lin_Weigh_Sim.png)

Subfolder ["Chad Cole LDPC Error Floor Archive"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Chad%20Cole%20LDPC%20Error%20Floor%20Archive)
contained  Dr. Chad A. Cole Importance Sampling for Trapping set enumerating and weighing Matlab source code.

Root folder of repo contain instruction to use binary file TS_enum.exe implementation of Cole's Trapping sets enumeration method, [1,2]:


Impulse tree have size (3,5), for bigger tree it take too lot of time and require use multi-core CPU or GPU implementation.
To use first convert from sparse (alist) or quasi-cyclic (qc) representation to graph represetation using sparse2graph matlab script after run enumerating trapping sets. 



Example 1 (QC representation):
H=qc2sparse('tanner.qc'); 
sparse2graph(H, 'tanner.graph');


TS_enum.exe -maxit 25 -fast  -qc 31 -x  2  tanner.graph tanner.trap


Example 2 (sparse represenation):
H=alist2sparse('PEGirReg252x504'); 
sparse2graph(H, 'PEGirReg252x504.graph');


TS_enum.exe -maxit 25 -fast   PEGirReg252x504.graph PEGirReg252x504.trap


Example 3 code from [3]: 


load 9_3.mat

sparse2graph(a, '9_3.graph')

run bat file '9_3.cmd'



Parameters of TS_enum.exe:



TS_enum.exe   tanner_graph.alist_like_format  output_list_TS 




  tanner_graph.alist_like_format  Define parity-check matrix file

  output_list_TS                  Define output file with list of TS.



Optional parameters:

  -x channel_factor               Channel factor. By default equal 4.0.

  -e impulse_magnitude            Set the magnitude of error impulse. The likelihoods
                                  of bits corrupted by error impulse will be defined as
                                  channel_factor * (1.0 - impulse_magnitude).
                                  By default equal 1.

  -f Start_from_VN                Start search for trapping sets from Start_from_VN variable.
                                  By default equal 1.

  -n Cardinal                     Enumerate TS only at Cardinal variable nodes.

  -qc Authomorphism               Search every (Start_from_VN + j * Authomorphism), j = 0, 1, 2, ...
                                  variable for trapping set. Decrease seach space when finding trapping
                                  sets under quasi-cyclic and cyclic codes.
                                  By default equal 1 (search every variable for trapping set)

  -gamma gamma                    Set the value of input bits not touched
                                  by error impulse. Their likelihoods will be defined as
                                  channel_factor * gamma
                                  By default equal 0.6.

  -rl likelihood                  Set reliable likelihood level for the stopping
                                  criterion. By default equal 20.

  -clip saturate                  Set clipping level for message values to saturate.
                                  By default equal 60.

  -fast                           Use quantized BP (integer SSE/AVX instead float) .

  -maxit iters                    Set maximum number of iterations to iters.
                                  By default equal 100.

  -critical                       Save trapping set to file when it was found. If program
                                  is unexpectedly terminated all found trapping sets
                                  would not be lost.

  -quiet                          Don't print a message for each found trapping set.

  -maxodd N                       Consider only trapping sets of form (a, b) with b <= N.
                                  By default equal limeted only by impulse tree size settings.


 When I have free time, I will publish an source files to Trapping set search method based on the LP outlined in my Phd Thesis, "Topologically Driven Methods for Construction Of Multi-Edge Type (Multigraph with nodes puncturing) Quasi-Cyclic Low-density Parity-check Codes for Wireless Channel, WDM Long-Haul and Archival Holographic Memory", https://arxiv.org/abs/2011.14753. 


P.S.
Upd: 20 January 2022. Gift for the New Year)  Dr. Chad A. Cole kindly share his Matlab source code (folder Chad Cole LDPC Error Floor Archive) for Trapping Sets analysis, famous research in detail described at report [1] and in paper [2]. 


Thank you very much Dr. Chad A. Cole!


References


[1].  Chad A. Cole, et al. Hall A General Method for Finding Low Error Rates of LDPC Codes, May 2006 https://arxiv.org/abs/cs/0605051


[2].  "C. A. Cole, S. G. Wilson, E. K. Hall and T. R. Giallorenzi, “A general method for finding low error rates of LDPC codes, ” submitted to IEEE Trans. on Inform. Theory, June 2006."


[3]. Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004.


