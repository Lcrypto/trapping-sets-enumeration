# Importance Sampling using modified and original Cole's LDPC trapping sets enumeration and weighing method for BER/FER error-floor estimation
Subfolder ["Article An Modified Cole's Importance Sampling Method For Low Error Floor QC-LDPC Codes Construction"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction) contained parity-check matrices, EMD Spectrum Enumerators(estimated using my tool [Extrinsic Message Degree Spectrum](https://github.com/Lcrypto/EMD-Spectrum-LDPC) )  and Trapping set enumerator (w\ot ordering) in zip compresed files for five QC-LDPC codes constructed for article An Modified Cole's Importance Sampling Method For Low Error Floor QC-LDPC Codes Construction.

![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction/Table_1.png)

![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction/Table_2.png)

![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction/Table_3_5.png)

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




In order to showcase the superiority of the proposed modified Cole method over linear programming methods, we will examine a basic search example TS(10,6) at Margulis Code. Employing CPLEX 12.8 with prior knowledge of 1,381,471 variable nodes takes approximately 5.30 seconds (1401.17 ticks) in parallel mode employing up to 32 threads. In contrast, our proposed modified Cole method uncovers 66 Trapping sets, including the target TS(10,6), in less than one second utilizing a single process.


TS variable nodes, (x_0)...(x_end) on Fig: 1 381 471 1691 1935 1988 2263 2412 2587 2617 


TS check nodes, (c_0)...(c_end) on Fig: 5 60 370 382 460 568 572 611 646 714 782 850 888 1100 1125 1137 1223 1253 



Trapping_set_submatrix =



     1     0     0     0     0     0     0     1     0     0
     0     0     0     0     0     0     0     0     1     1
     0     1     0     0     0     0     0     0     0     0
     0     0     0     1     0     0     0     1     0     0
     0     0     1     0     0     0     1     0     0     0
     0     0     0     0     0     1     0     0     0     1
     0     0     1     0     0     0     0     0     0     0
     0     0     0     1     1     0     0     0     0     0
     0     0     0     0     0     1     0     0     0     0
     0     0     1     0     1     0     0     0     0     0
     0     0     0     1     0     1     0     0     0     0
     0     0     0     0     1     0     0     0     1     0
     0     0     0     0     0     0     1     0     0     0
     1     1     0     0     0     0     0     0     0     0
     0     0     0     0     0     0     0     1     0     0
     1     0     0     0     0     0     1     0     0     0
     0     1     0     0     0     0     0     0     1     0
     0     0     0     0     0     0     0     0     0     1



![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction/Margulis_TS(10%2C6).png)



Second example to showcase the superiority of the our method over Velasquez-Subramanilinear programming method, we will examine a larges TS(a,b), where b is minimal search example  TS(62,16) at Margulis Code [4896 2448] from http://www.inference.org.uk/mackay/codes/data.html#l28 (Pascal Vontobel and Joachim Rosenthal Constructed using Ramanujan graphs and ideas from Margulis. (q=17, p=5, Allerton 2000) . Employing  with prior knowledge of 6, 41, 1132, 1140, 1508, 1561, 1988, 2506, 3154, 3182, 3444, 3575, 4141, 4710 variable nodes (involved in cycles, codes girth 12) 
and 2 step of cycles increase 
step_1  6 41 364 719 769 880 970 1097 1132 1140 1164 1274 1385 1508 1561 1625 1681 1716 1988 2194 2216 2369 2506 2738 2765 2950 3154 3182 3444 3575 3940 4012 4109 4141 4162 4355 4463 4547 4581 4710 4800 4864 
step_2  6 41 65 159 364 719 769 880 970 981 1097 1132 1140 1164 1274 1279 1385 1508 1561 1625 1681 1716 1823 1979 1988 2194 2216 2369 2460 2506 2738 2765 2795 2855 2950 2976 3024 3154 3182 3271 3444 3566 3575 3836 3940 4012 4109 4141 4162 4355 4396 4401 4463 4547 4581 4671 4710 4800 4859 4864 



We found TS(62,16) with  TS variable nodes: 6 41 65 159 364 546 719 769 880 970 981 1097 1132 1140 1164 1274 1279 1385 1508 1561 1625 1681 1716 1819 1823 1979 1988 2194 2216 2369 2460 2506 2738 2765 2795 2855 2950 2976 3024 3154 3182 3271 3444 3566 3575 3836 3940 4012 4109 4141 4162 4355 4396 4401 4463 4547 4581 4671 4710 4800 4859 4864 





and it take approximately  1.25 sec. (292.92 ticks) in parallel mode employing up to 32 threads. For example in article https://www.sciencedirect.com/science/article/pii/S030439752200127X Alvaro Velasquez, K. Subramani, Piotr Wojciechowski,
On the complexity of and solutions to the minimum stopping and trapping set problems, Theoretical Computer Science, Volume 915,
2022,Pages 26-44 was claimed as one of our main results is finding a stopping set of size 48 in the 
Margulis Code [4896 2448] (see Fig. 11 in https://www.cs.ucf.edu/~velasquez/StoppingSets/) in 700451 seconds. We need less than one minutes overall (search of variable nodes in cycles and after search of TS(a,b)) to found much more complex TS(62,16) in Margulis Code(4896 2448).




Detailed classification of harmfull Trapping sets with figure and submatrix for QC-LDPC code 1 and code 5 you can see at https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction


When I find some free time, I plan to make available the source files for a Trapping set search method that is based on the linear programming outlined in my PhD thesis titled "Topologically Driven Methods for Construction of Multi-Edge Type Quasi-Cyclic Low-Density Parity-Check Codes (Multigraph with Nodes Puncturing) for Wireless Channel, WDM Long-Haul, and Archival Holographic Memory". You can find my thesis at https://arxiv.org/abs/2011.14753.


P.S.
Update as of January 20th, 2022: As a gift for the New Year, Dr. Chad A. Cole has kindly shared his MATLAB source code for Trapping Set analysis. The source code can be found in the "Chad Cole LDPC Error Floor Archive" folder. Trapping Sets are an important research topic that is thoroughly described in Dr. Cole's report [1] and paper [2].

We would like to express our gratitude to Dr. Chad A. Cole for sharing his work with the community. You can find more information about Dr. Cole's research on his GitHub page at https://github.com/chadac8j.




References


[1].  Chad A. Cole, et al. Hall A General Method for Finding Low Error Rates of LDPC Codes, May 2006 https://arxiv.org/abs/cs/0605051


[2].  "C. A. Cole, S. G. Wilson, E. K. Hall and T. R. Giallorenzi, “A general method for finding low error rates of LDPC codes, ” submitted to IEEE Trans. on Inform. Theory, June 2006."


[3]. Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004.


