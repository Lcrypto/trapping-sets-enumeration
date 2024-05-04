# This repositories contained  Mixed Linear Programing method for Trappings set search

The folder contains a Jupyter Notebook script named (`main_MLP_Trapping_Set_Enumeration_with_VNs.ipynb`)[https://github.com/Lcrypto/trapping-sets-enumeration/blob/master/LP/main_MLP_Trapping_Set_Enumeration_with_VNs.ipynb], which implements a Mixed Integer Linear Programming (MILP) method for enumerating trapping sets as described in detail in articles [11-12].

The method introduced in the paper offers an efficient way to identify trapping sets in codes represented on graphs, using a mixed linear programming approach. This method ensures a comprehensive search, which is crucial for nonlinear communication channels, energy-based models, natural language processing deep neural networks, metric learning, and other complex dynamic systems. By analyzing the symmetry and asymmetry properties of dynamic systems through trapping sets (TS(a, 0) for codewords and TS(a, b) for pseudocodewords), the method provides a detailed understanding of their behavior.

Implemented in Python, this method is publicly available on GitHub and supports both the Community Edition (with limited conditions) and the Commercial Edition of CPLEX. The technique involves solving a mixed integer linear programming problem using a predefined list of variable nodes participating in the shortest cycles with small Extrinsic Message Degree values within the code graph. These cycles, similar to topological invariants, represent multidimensional voids formed by code and pseudocode words. Therefore, this method can be seen as an approach to construct topological complexes and calculate topological invariants.

The method was applied to search for trapping sets in LDPC codes using IBM CPLEX Optimization Studio version 22.1.0.0. The computations were performed on a 16-core AMD Ryzen 3950X processor with 128GB of DDR4 RAM, utilizing 32 threads. For example, in the Margulis code (2640, 1320), the method identified the trapping set TS(6,6) in just 0.29 seconds, showing a significant speed improvement compared to the Velasquez-Subramani method.

Additionally, the method successfully identified trapping sets TS(62,16) and TS(52,14) for the first time in the Margulis code (4896, 2474), showcasing its thorough search capabilities. Furthermore, in the Mackay LDPC code (408,204), the method discovered trapping sets TS(102,2) and TS(108,4) with the help of CPLEX Community Edition version 22.1.1.0. Notably, these findings required a brute-force trapping set search, which would have demanded over $10^{101}$ operations without this efficient method.


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

[12] V. S. Usatyuk and S. I. Egorov, "Mixed Integer Linear Programming Method for Energy-Based Model Trapping Sets Enumerating," 2024 26th International Conference on Digital Signal Processing and its Applications (DSPA), Moscow, Russian Federation, 2024, pp. 1-6, doi: 10.1109/DSPA60853.2024.10510058. https://ieeexplore.ieee.org/document/10510058
