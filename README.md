# Linear Programming, Importance Sampling using modified and original Cole's LDPC trapping sets enumeration and weighing method for BER/FER error-floor estimation


The subfolder titled ["LP"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/LP) contains a Jupyter Notebook script named [`main_MLP_Trapping_Set_Enumeration_with_VNs.ipynb`](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/LP/main_MLP_Trapping_Set_Enumeration_with_VNs.ipynb), which implements a method using Mixed Integer Linear Programming (MILP) for enumerating trapping sets. This approach is detailed in the following articles:
- Usatyuk V.S., Egorov S.I. "Trapping Sets Search Using the Method of Mixed Integer Linear Programming with a Priori List of Variable Nodes." Proceedings of the Southwest State University. 2024; 27(4):79-97. (In Russ.) [DOI: 10.21869/2223-1560-2023-27-4-79-97](http://dx.doi.org/10.21869/2223-1560-2023-27-4-79-97)
- V. S. Usatyuk and S. I. Egorov, "Mixed Integer Linear Programming Method for Energy-Based Model Trapping Sets Enumerating," presented at the 2024 26th International Conference on Digital Signal Processing and its Applications (DSPA), Moscow, Russian Federation, pp. 1-6. [DOI: 10.1109/DSPA60853.2024.10510058](https://ieeexplore.ieee.org/document/10510058)




### [Colab ready to work script url](https://colab.research.google.com/gist/Lcrypto/c589fac512fc821024708dd27c91ae63/trapping_set_enumeration_with_vns.ipynb), just upload to your files folder parity-check matrix file [Mackay_408.33.864.alist](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/LP/Mackay_408.33.864.alist)   (or Mackay_96.3.967.alist, tanner.qc from ["LP"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/LP) folder), or  your parity-check matrix in alist/qc format. GPU not necessary to use.



 The method was applied to search for trapping sets in LDPC codes, utilizing the mathematical linear programming package IBM CPLEX Optimization Studio version 22.1.0.0 The computations ran on a 16-core AMD Ryzen 3950X processor with 128GB KF3200 DDR4 RAM, utilizing 32 threads. In the Margulis code (2640, 1320), the proposed method identified the trapping set TS(6,6) in just 0.29 seconds. Notably, compared to the Velasquez-Subramani method, the proposed method achieved a speedup of 15082 times. Leveraging its high speed and thorough search capabilities, the method successfully identified trapping sets TS(62,16) and TS(52,14) for the first time in the Margulis code (4896, 2474). 
Finally, employing the proposed method with CPLEX Community Edition version 22.1.1.0, we discovered trapping sets, specifically $TS(102,2)$ and $TS(108,4)$, within the Mackay LDPC code $(408,204)$. These findings necessitated a brute-force trapping set search, demanding more than $10^{101}$ operations.


![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/LP/ExtrLargeTS.png)





The subfolder titled  ["Article An Modified Cole's Importance Sampling Method For Low Error Floor QC-LDPC Codes Construction"](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction) contains multiple txt file and compressed zip file for five different QC-LDPC codes that were constructed for the article titled "Low Error Floor QC-LDPC Codes Construction Using Modified Cole's Trapping Sets Enumerating" by V. S. Usatyuk at the 2023 25th International Conference on Digital Signal Processing and its Applications (DSPA) in Moscow, Russian Federation, https://ieeexplore.ieee.org/document/10113442. These zip file contain parity-check matrices, EMD Spectrum Enumerators (estimated using my Extrinsic Message Degree Spectrum tool), and Trapping set enumerators (without ordering). In addition to the previously mentioned contents, the subfolder also includes a detailed classification of the most harmful Trapping sets for QC-LDPC codes 1 and 5. This classification is presented with [figures and submatrix](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction) to provide a clear understanding of these Trapping sets. Additionally, there is other relevant information available in this folder.



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



To demonstrate the effectiveness of [our LP method](https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/LP) (  from [my Phd thesis](https://github.com/Lcrypto/Phd-thesis-LDPC-code-construction-method-using-cycles-topology-optimization-code-pseudo-codewords) ) over the Velasquez-Subramani linear programming method, we will analyze a large Tanner graph TS(a,b) with minimal search example TS(62,16) at Margulis Code (4896, 2448) from http://www.inference.org.uk/mackay/codes/data.html#l28. This particular Tanner graph was constructed using Ramanujan graphs and ideas from Margulis (q=17, p=5, J. Rosenthal and P. O. Vontobel, "Constructions of regular and irregular LDPC codes using Ramanujan graphs and ideas from Margulis," Proceedings. 2001 IEEE International Symposium on Information Theory,  2001, pp. 4-). For search we use 14 variable nodes (involved in cycles, codes girth 12) namely 6, 41, 1132, 1140, 1508, 1561, 1988, 2506, 3154, 3182, 3444, 3575, 4141, and 4710 with prior knowledge. Additionally, we will use the two-step cycle increase method:


step_1 VN:  6 41 364 719 769 880 970 1097 1132 1140 1164 1274 1385 1508 1561 1625 1681 1716 1988 2194 2216 2369 2506 2738 2765 2950 3154 3182 3444 3575 3940 4012 4109 4141 4162 4355 4463 4547 4581 4710 4800 4864 


step_2 VN:  6 41 65 159 364 719 769 880 970 981 1097 1132 1140 1164 1274 1279 1385 1508 1561 1625 1681 1716 1823 1979 1988 2194 2216 2369 2460 2506 2738 2765 2795 2855 2950 2976 3024 3154 3182 3271 3444 3566 3575 3836 3940 4012 4109 4141 4162 4355 4396 4401 4463 4547 4581 4671 4710 4800 4859 4864 

after use LP method.

We found TS(62,16) with 

TS variable nodes, (x_0)...(x_end) on Fig: 6 41 65 159 364 546 719 769 880 970 981 1097 1132 1140 1164 1274 1279 1385 1508 1561 1625 1681 1716 1819 1823 1979 1988 2194 2216 2369 2460 2506 2738 2765 2795 2855 2950 2976 3024 3154 3182 3271 3444 3566 3575 3836 3940 4012 4109 4141 4162 4355 4396 4401 4463 4547 4581 4671 4710 4800 4859 4864 



TS check nodes, (c_0)...(c_end) on Fig: 8 15 23 53 86 113 122 123 132 133 204 225 234 261 377 460 486 514 540 541 559 587 706 723 744 751 756 774 777 814 849 874 876 880 882 934 935 963 995 1009 1073 1083 1098 1106 1119 1133 1148 1225 1233 1244 1269 1342 1348 1377 1390 1402 1465 1502 1525 1536 1562 1586 1589 1625 1646 1658 1666 1670 1708 1746 1763 1764 1766 1774 1801 1806 1807 1832 1893 1894 1896 1908 1912 1942 1959 1967 1975 1994 2003 2013 2078 2090 2190 2242 2271 2281 2320 2336 2339 2370 


![alt text](https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/Margulis%20Ramanujan%2017.5%20code/TS(62%2C16)Margulis.png)

Trapping_set_submatrix =


     1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
     0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
     0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
     0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
     0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
     0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	1	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	1	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
     0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
     0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0



In parallel mode using up to 32 threads, our method takes approximately 3.67 sec (292.92 ticks) to complete. This is significantly faster than the results reported in Alvaro Velasquez, K. Subramani, and Piotr Wojciechowski's article "On the complexity of and solutions to the minimum stopping and trapping set problems" published in Theoretical Computer Science (Volume 915, 2022, Pages 26-44), https://www.sciencedirect.com/science/article/pii/S030439752200127X. In their study, they claimed that one of their main results was finding a stopping set of size 48 in the Margulis Code (4896, 2448) in 700451 seconds. Our method, on the other hand, took less than one minute overall to find a much more complex TS(62,16) in Margulis Code (4896, 2448) after searching variable nodes in cycles followed by the search of TS(a,b). These results demonstrate the superiority of our method over Velasquez-Subramani linear programming method.


Detailed classification of harmfull Trapping sets with figure and submatrix for QC-LDPC code 1 and code 5 you can see at https://github.com/Lcrypto/trapping-sets-enumeration/tree/master/Article%20An%20Modified%20Cole's%20Importance%20Sampling%20Method%20For%20Low%20Error%20Floor%20QC-LDPC%20Codes%20Construction


When I find some free time, I plan to make available the source files for a Trapping set search method that is based on the linear programming outlined in my PhD thesis titled "Topologically Driven Methods for Construction of Multi-Edge Type Quasi-Cyclic Low-Density Parity-Check Codes (Multigraph with Nodes Puncturing) for Wireless Channel, WDM Long-Haul, and Archival Holographic Memory". You can find my thesis at https://arxiv.org/abs/2011.14753.


You can find the C++ implementation of the stopping sets enumerator using linear programming at the following link:  https://github.com/Lcrypto/Stopping-sets-enumerator 


P.S.
Update as of January 20th, 2022: As a gift for the New Year, Dr. Chad A. Cole has kindly shared his MATLAB source code for Trapping Set analysis. The source code can be found in the "Chad Cole LDPC Error Floor Archive" folder. Trapping Sets are an important research topic that is thoroughly described in Dr. Cole's report [1] and paper [2].

We would like to express our gratitude to Dr. Chad A. Cole for sharing his work with the community. You can find more information about Dr. Cole's research on his GitHub page at https://github.com/chadac8j.




References


[1].  Chad A. Cole, et al. Hall A General Method for Finding Low Error Rates of LDPC Codes, May 2006 https://arxiv.org/abs/cs/0605051


[2].  "C. A. Cole, S. G. Wilson, E. K. Hall and T. R. Giallorenzi, “A general method for finding low error rates of LDPC codes, ” submitted to IEEE Trans. on Inform. Theory, June 2006."


[3]. Tao Tian, C. R. Jones, J. D. Villasenor and R. D. Wesel, "Selective avoidance of cycles in irregular LDPC code construction," in IEEE Transactions on Communications, vol. 52, no. 8, pp. 1242-1247, Aug. 2004.


[4]. Rosnes E., Ytrehus O.  "An Algorithm to Find All Small-Size Stopping Sets of Low-Density Parity-Check Matrices," 2007 IEEE International Symposium on Information Theory, Nice, France, 2007, pp. 2936-2940


[5]. Margulis G. A.  "Explicit Constructions of Graphs Without Short Cycles and Low Density Codes", Combinatorica, vol. 2, no. 1, pp. 71-78, 1982


[6]. Rosenthal J., Vontobel P. O.  "Constructions of regular and irregular LDPC codes using Ramanujan graphs and ideas from Margulis," Proceedings. 2001 IEEE International Symposium on Information Theory, 2001, pp. 4


[7]. David J.C. MacKay Encyclopedia of Sparse Graph Codes   https://www.inference.org.uk/mackay/codes/data.html


[8]. Velasquez A., et al. Finding Minimum Stopping and Trapping Sets: An Integer Linear Programming Approach // In: Lee J., Rinaldi G., Mahjoub A. (eds) Comb. Optim. ISCO 2018. Lect. Notes in Comp. Science, V. 10856., pp. 402–415


[9]. Velasquez A., Subramani K., Wojciechowski P. On the complexity of and solutions to the minimum stopping and trapping set problems // Theor. Comput. Sci. – 2022. – Vol. 915. pp. 26-44.


[10]. Usatyuk V.S., Egorov S.I. Construction of LDPC Codes Using Cole's Modified Importance Sampling Method. Proceedings of the Southwest State University. 2023;27(1):92-110. (In Russ.) https://doi.org/10.21869/2223-1560-2023-27-1-92-110


[11]. Usatyuk V.S., Egorov S.I. Trapping Sets Search Using the Method of Mixed Integer Linear Programming with a Priori List of Variable Nodes. Proceedings of the Southwest State University. 2024; 27(4):79-97. (In Russ.) http://dx.doi.org/10.21869/2223-1560-2023-27-4-79-97


[12]. V. S. Usatyuk and S. I. Egorov, "Mixed Integer Linear Programming Method for Energy-Based Model Trapping Sets Enumerating," 2024 26th International Conference on Digital Signal Processing and its Applications (DSPA), Moscow, Russian Federation, 2024, pp. 1-6, doi: 10.1109/DSPA60853.2024.10510058. https://ieeexplore.ieee.org/document/10510058


