%%%%%%%%%%%%%
%% This software was developed by Chad Cole as a grad student while funded by L-3
%% Communications West.
%%
%% TSearch.m - Requires input variable named H_sparse, which is the LDPC
%%   code parity-check matrix.
%%
%%  This program implements the first step of the 3-step procedure outlined
%%  in "A General Method for Finding Low Error Rates of LDPC Codes."
%%
%%  Most types of moderate block length LDPC codes are supported - regular,
%%  irregular, and high-rate for example.  The core decoder implementation
%%  is the same as in the program 'irregularSPA.m,' another piece of
%%  software in this suite, so the input variable describing the code should
%%  be named H_sparse.
%%%%%%%%%%%%

clear all
% randomly set seed based on processor clock
randn('state',sum(100*clock));
format compact
%function IS_cwfinder()

% load H2000_no_8_cycle 
% load QCH2331
% load QCH135 
% load QCH8160 
% load QR_H1008
% load QR_H2048 
% load QR_H8192 
% load H504_4_8_no4cycle
% load H504_4_8_no6_fallback4_best
% load H1000_4_8_gallager 
% load test_1200H_3_6_no6
load H500_rate08_3_15
% load H504goodcode_no82
% load H4000_wesel_test
% load H2048_8023an
% load H1008_nox2
% load H1000_rate08_4_20
% load H1200_4_8_no6cycle
% load H96_firstmackay
% load H96-964_H_sparse
% H = [1 1 1 1 0 0 0 0 0 0;
%      1 0 0 0 1 1 1 0 0 0;
%      0 1 0 0 1 0 0 1 1 0;
%      0 0 1 0 0 1 0 1 0 1;
%      0 0 0 1 0 0 1 0 1 1]
% H_sparse = sparse(H);
               
n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

H_col_sum = sum(H_sparse,1);
H_row_sum = sum(H_sparse,2);
%B will be sorted in ascending order
B = unique(H_col_sum);
num_v_deg = length(B);
R = unique(H_row_sum);
num_c_deg = length(R);
deg_c_prof = zeros(length(R), 2);
deg_v_prof = zeros(length(B), 2);

row_perm = zeros(1,n);

H_new = spalloc(n-k,n,length(find(H_sparse)));
running_index = 0;
for ii = 1:num_v_deg
    change_loc = find(H_col_sum == B(num_v_deg + 1 - ii));
    deg_v_prof(ii, 1) = length(change_loc);
    deg_v_prof(ii, 2) = B(num_v_deg + 1 - ii);
%This step orders columns from highest degree to lowest
    for jj = 1:length(change_loc)
        running_index = running_index + 1;
        H_new(:,running_index) = H_sparse(:, change_loc(jj));
    end
end
H_sparse = H_new;
clear H_new;

H_new = spalloc(n-k,n,length(find(H_sparse)));
running_index = 0;
for ii = 1:num_c_deg
    change_loc = find(H_row_sum == R(num_c_deg + 1 - ii));
    deg_c_prof(ii, 1) = length(change_loc);
    deg_c_prof(ii, 2) = R(num_c_deg + 1 - ii);
%This step orders rows from highest degree to lowest
    for jj = 1:length(change_loc)
        running_index = running_index + 1;
        H_new(running_index,:) = H_sparse(change_loc(jj),:);
    end
end
H_sparse = H_new;
clear H_new;

Vlist=zeros(n-k,max(sum(H_sparse,2))+1);  % list of V-nodes
Clist=zeros(n,max(sum(H_sparse,1))+1);    % list of C-nodes
for jj=1:n-k
    Vlist(:,1)=sum(H_sparse,2);
    icnt=0;
    for ii=1:n
        if H_sparse(jj,ii)==1
            icnt=icnt+1;
            Vlist(jj,icnt+1)=ii;
        end
    end
end

for ii=1:n
    Clist(:,1)=(sum(H_sparse,1))';
    jcnt=0;
    for jj=1:n-k
        if H_sparse(jj,ii)==1
            jcnt=jcnt+1;
            Clist(ii,jcnt+1)=jj;
        end
    end
end

num_edges_v = sum(Clist(:,1));
num_edges_c = sum(Vlist(:,1));
if (num_edges_c ~= num_edges_v)
    text = ['error - row/col mismatch in H']
end
num_edges = num_edges_c;
clear num_edges_c num_edges_v

Lr_ind = zeros(1, num_edges);
Lq_ind = zeros(1, num_edges);

tot_count = 0;
for kk = 2:length(Clist(1,:))
    for jj=1:n
        if(Clist(jj,kk) ~= 0)
             tot_count = tot_count + 1;
             Lq_ind(tot_count) = (kk-2)*n + jj;
        end
    end
end

row_index = zeros(n-k, 1);
tot_count = 0;
for kk=2:length(Clist(1,:))
    for jj=1:n
      if(Clist(jj,kk) ~= 0)
        tot_count = tot_count + 1;
        row=Clist(jj, kk); 
        row_index(row) = row_index(row) + 1;
        Lr_ind(tot_count) = row + (row_index(row)-1)*(n-k);
      end
    end
end

Lq=zeros(n,max(Clist(:,1)));  % messages up
Lr=zeros(n-k,max(Vlist(:,1)));  % messages down
LQ=zeros(1,n);
Lc=zeros(1,n);

% TS/cw_unique_count(i) count how many times the i^{th} TS or codeword was
% found with the search program.  The most dominant TS will be found many
% times whereas less dominant ones may only be found once or twice.  This
% is one guage for how dominant a specific TS may be.
TS_unique_count = zeros(1,1);
cw_unique_count = zeros(1,1);
num_cw = 0;
num_miscor = 0;
itenum=50;  % maximum # of iterations

num_USC = zeros(1, itenum);
USC_vec = zeros(1, n-k);
c_hat_hist = zeros(itenum, n);
num_TS = 0;
TS = spalloc(1000000,n,30000*30);
cw = spalloc(10000,n,3000*30);

ite_hist = int8(zeros(1,1200000));

small_num = eps;

%%%%%%%%%%%%%%%5
% The selection of error impulse parameters is more art than science, so I
% will list a few examples.  The variable names are the same as used in the
% paper "A General Method for Finding Low Error Rates of \textsc{LDPC}
% Codes".
%  INCLUDE SOME EXAMPLE H MATRICES
% The H1008_nox2 code is a good test case because the search algorithm
% works very well for this type of regular (3,6) code.  A good choice of
% parameters for this code are epsilon_1 = 4, epsilon_2 = gamma, gamma = 0.7, ebnodb = 7.
% It takes 44 minutes to run the search with these settings and 761 TS are
% found.


% H500_rate08_3_15 is a (500,400) (3,15) regular code.  Since this is a very densely-packed code,
% the impulse values can be small, and only 2 variable nodes below the root
% need to used as impulse candidates:
%  epsilon_1 = 1.5, gamma = 0.8, epsilon_2 = gamma, ebnodb = 7 dB, v_num =
%  2.

% The TS_unique_count structure keeps a tabulation of how many times each
% TS has been found with the search.  This is a very good indicator of how
% dominant a specific TS is.  For example, the TS with max(TS_unique_count)
% is probably the most dominant TS in the code.  This should be verified
% with step two of the 3-step procedure of course, but it usually is true.

% In conclusion, this method takes much time to get familiar with; it is a
% difficult concept.  I have spent many hours working with many different codes.
% But i can assure you that the method is very powerful and if it is your
% objective to find error floors of LDPC codes, it is very worthwhile to
% make the required investment to learn how to use this program/procedure
% effectively.

epsilon_1 = 1.5;
gamma = 0.8*ones(1,length(deg_v_prof(:,1)));
ebnodb = [7];
epsilon_2 = gamma(1);
% epsilon_2 = -0.7;
v_num = 2;

ebno = 10^(ebnodb/10);
esno = ebno*(rate);%*3;  %rate 1/2, 3 bits/sym
Lc_coef=4*esno;
sigma_2 = 1./(2*esno);
sigma = sqrt(sigma_2);

num_checks = 0;
num_combos = 0;
num_combos_save = zeros(1, length(deg_v_prof(:,1)));

%Start building trees from which to apply error impulses starting at the
%highest numbered variable node because the H matrix has been put into a
%form where the smallest degree variable nodes, and thus the ones most
%likely to cause the error floor, are at the highest numbered positions.
v_count = n+1;
% v_count = n-3;

v_node_hit_hist = zeros(1,n);

%Check syndrome weight of error impulse bits in case we use the pruning
%procedure to limit the number of impulse candidates that actually use the
%decoder to attempt correction.
syn_save_hist = zeros(1, 81);
syn_save_hist_hit = zeros(1, 81);

%Send all-zeros vector, which is BPSK modulated to all-ones.
a = ones(1,n);
tic
for iii = 1:length(deg_v_prof(:,1))
    %Need to get the indices of (d_v choose v_num) error impulse candidates
    %at variable tier one.
    deg_current = deg_v_prof(length(deg_v_prof(:,1))+1-iii,2);
    num_combos = 0;
    if (v_num >= deg_current)
       num_combos = 1;
       indices = [1:deg_current]';
    else    
      indices = zeros(v_num, nchoosek(deg_v_prof(length(deg_v_prof(:,1))+1-iii,2), v_num));
      for pp = 1:2^deg_current
       binary = int2vec(pp, deg_current);
       if (sum(binary) == v_num)
           num_combos = num_combos + 1;
           indices(:, num_combos) = find(binary);
       end
      end
    end
    num_combos_save(iii) = num_combos;

    for jj = 1:deg_v_prof(length(deg_v_prof(:,1))+1-iii,1)
      v_count = v_count - 1
     for kk = 1:num_combos
        
% The c_nodes structure contains the check-nodes involved in the tree rooted 
%  at the v_count^{th} variable node
    c_nodes = Clist(v_count, indices(:,kk)+1);
% The v_tier_1 structure contains all possible variable nodes that are in
% variable tier 1 of the tree rooted at the v_count^{th} variable node.
    v_tier_1 = zeros(length(c_nodes), max(Vlist(c_nodes, 1)));
    for ppp = 1:length(c_nodes)
        v_tier_1(ppp, 1:Vlist(c_nodes(ppp),1)-1) = setxor(Vlist(c_nodes(ppp),2:Vlist(c_nodes(ppp),1)+1),v_count);
    end
    % v_tier_1_ind contains indices of nonzero elements of v_tier_1
    v_tier_1_ind = find(v_tier_1);
    %store candidates for epsilon_2 impulses in this matrix which only has
    %to be calculated once for each of our n trees
    epsilon_mat = zeros(length(v_tier_1_ind), (max(Vlist(:,1))-1)*(max(Clist(v_tier_1(v_tier_1_ind),1))-1));
    for ppp = 1:length(v_tier_1_ind)
        tot = 0;
        for qqq = 1:Clist(v_tier_1(v_tier_1_ind(ppp)),1)
            if sum(Clist(v_tier_1(v_tier_1_ind(ppp)),qqq+1)*ones(1,length(c_nodes)) == c_nodes)==0
                epsilon_mat(ppp, tot+1:tot+Vlist(Clist(v_tier_1(v_tier_1_ind(ppp)),qqq+1),1)-1) ...
                    = setxor(Vlist(Clist(v_tier_1(v_tier_1_ind(ppp)),qqq+1), 2: ...
                    Vlist(Clist(v_tier_1(v_tier_1_ind(ppp)),qqq+1),1)+1), ...
                    v_tier_1(v_tier_1_ind(ppp)));
                tot = tot + Vlist(Clist(v_tier_1(v_tier_1_ind(ppp)),qqq+1),1)-1;
            end
        end
    end
                    
                    
    %tot_perms counts the number of permutations of impulse bits for each
    %  tree
    tot_perms = 1;
    for ppp = 1:length(c_nodes)
        tot_perms = tot_perms*(Vlist(c_nodes(ppp),1)-1);
    end

    %bit_combos enumerates all tot_perms length(c_nodes)-bit error impulse
    %  candidate sets
    bit_combos = zeros(length(c_nodes),tot_perms);
    count = zeros(length(c_nodes),1);
    for ppp = 1:tot_perms
        v_prod = 1;
        for zzz = 1:length(c_nodes)
            if (mod(ppp, v_prod)==0)
                count(zzz) = mod(count(zzz)+1,Vlist(c_nodes(zzz),1)-1);
            end
            v_prod = v_prod*(Vlist(c_nodes(zzz),1)-1);
        end
        for zzz = 1:length(indices(:,1))
            bit_combos(zzz, ppp) =  v_tier_1(zzz, count(zzz)+1);
        end
    end
    
    for ll = 1:tot_perms
    
%%%%The following code can be 
%%%%uncommented for finding TS in codes where searching the first layer
%%%%of variable nodes falls short, typically for larger variable degree
%%%%nodes, such as {4,8} codes with n > 2000 or so.  Don't forget to also
%%%%uncomment the end of the for loop below this block of code.

%%%% Find c_nodes connected to leftmost v_node at tier 1.  

%     v_root=intersect(Vlist(c_nodes(1),2:Vlist(c_nodes(1),1)+1), ...
%               bit_combos(:,ll));
%     v_root_cnodes = setxor(Clist(v_root,2:Clist(v_root,1)+1), c_nodes(1));
%     v_tier_2 = zeros(length(v_root_cnodes), max(Vlist(v_root_cnodes, 1)));
%     for ppp = 1:length(v_root_cnodes)
%         v_tier_2(ppp, 1:Vlist(v_root_cnodes(ppp),1)-1) = setxor(Vlist(v_root_cnodes(ppp),2:Vlist(v_root_cnodes(ppp),1)+1),v_root);
%     end
% 
%     tot_perms_2 = 1;
%     for ppp = 1:length(v_root_cnodes)
%         tot_perms_2 = tot_perms_2*(Vlist(v_root_cnodes(ppp),1)-1);
%     end
% 
%     bit_combos_2 = zeros(length(v_root_cnodes),tot_perms_2);
%     count_2 = zeros(length(v_root_cnodes),1);
%     for ppp = 1:tot_perms_2
%         v_prod = 1;
%         for zzz = 1:length(v_root_cnodes)
%             if (mod(ppp, v_prod)==0)
%                 count_2(zzz) = mod(count_2(zzz)+1,Vlist(v_root_cnodes(zzz),1)-1);
%             end
%             v_prod = v_prod*(Vlist(v_root_cnodes(zzz),1)-1);
%         end
%         for zzz = 1:length(v_root_cnodes)
%             bit_combos_2(zzz, ppp) =  v_tier_2(zzz, count_2(zzz)+1);
%         end
%     end
%           
%     for mm = 1:tot_perms_2
% 
        num_checks = num_checks + 1;
          d = zeros(1,n);
          d(v_count) = 1;
          d(bit_combos(:, ll)) = 1;
%           d(bit_combos_2(:, mm)) = 1;
          v_node_hit_hist(bit_combos(:, ll)) = v_node_hit_hist(bit_combos(:, ll)) + 1;
          v_node_hit_hist(v_count) = v_node_hit_hist(v_count) + 1;
          syn_sum = sum(mod(d*H_sparse',2));
    if (syn_sum <= 14)
          syn_save_hist(syn_sum+1) = syn_save_hist(syn_sum+1) + 1;

          y = gamma(iii)*a - d*epsilon_1;
          %Only want nonzero values of epsilon_mat to index into y,
          % this not problem for (check) regular codes.  also a problem
          % here if girth < 8
          [garbage,garbage2,useful] = intersect(bit_combos(:, ll),v_tier_1(v_tier_1_ind));
%           y(epsilon_mat(useful,:)) = epsilon_2;

%      MPA
% step 1: initialization

Lc=Lc_coef*y;

num_total = 0;
for ii = 1:num_v_deg
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = Lc(num_total+1:num_total + deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2));
    num_total = num_total + deg_v_prof(ii, 1);
end

stopsig=0;  % check whether a valid codeword has been found
%itenum=20;  % maximum # of iterations
itestep=0;  

while (stopsig==0) && (itestep<itenum)
%while (itestep<itenum)
itestep=itestep+1;

% step 2: update Lr
%     Lr(Lr_ind) = Lq(Lq_ind);
%     num_total = 0;
%     for ii = 1:num_c_deg
%         aterm = sign(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)));
%         atermprod=prod(aterm,2)*ones(1,deg_c_prof(ii,2));
%         Lr_min_mat=sort(abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))');
%         Lr_min_vals=Lr_min_mat(1,:)'*ones(1,length(Lr_min_mat(:,1)));
%         Lr_min_vals_2 = Lr_min_mat(2,:)'*ones(1,length(Lr_min_mat(:,2)));
%         Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
%             (atermprod.*aterm).*((((Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2))))-1)*-1).*Lr_min_vals + ...
%             (Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))).*Lr_min_vals_2);
% %        bterm=-log(tanh(abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))*0.5) + small_num); %try speed-up from built-in matlab tanh?
% %        btermsum=sum(bterm,2)*ones(1,deg_c_prof(ii,2));
% %         Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
% %           (atermprod.*aterm).*-log(tanh(abs(btermsum-bterm)*0.5) + small_num);
%         num_total = num_total + deg_c_prof(ii, 1);
%     end

    Lr(Lr_ind) = Lq(Lq_ind);
    num_total = 0;
    for ii = 1:num_c_deg
        aterm = sign(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)));
        atermprod=prod(aterm,2)*ones(1,deg_c_prof(ii,2));
        bterm=-log(tanh(abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))*0.5) + small_num); %try speed-up from built-in matlab tanh?
        btermsum=sum(bterm,2)*ones(1,deg_c_prof(ii,2));
        Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
         (atermprod.*aterm).*-log(tanh(abs(btermsum-bterm)*0.5) + small_num);
%              1.05*(atermprod.*aterm).*-log(tanh(abs(btermsum-bterm)*0.5) + small_num);
%          (atermprod.*aterm).*-log(tanh(abs(btermsum-bterm)*0.5));
        num_total = num_total + deg_c_prof(ii, 1);
    end


Lq(Lq_ind) = Lr(Lr_ind);
num_total = 0;
for ii = 1:num_v_deg
    Lq_sum = sum(Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)),2);
    LQ(num_total+1:num_total+deg_v_prof(ii, 1)) = Lq_sum' + Lc(num_total+1:num_total+deg_v_prof(ii, 1));
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = ...
        LQ(num_total+1:num_total+deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2)) - ...
        Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2));
    num_total = num_total + deg_v_prof(ii, 1);
end

% step 5: validity check
c_hat=(1-sign(LQ))/2;

num_USC(itestep) = sum(mod(c_hat*H_sparse',2));
if num_USC(itestep)==0
  stopsig=1;
end %if

c_hat_hist(itestep,:) = c_hat;

end  %iterations while

ite_hist(num_checks) = itestep;

if (sum(c_hat) ~= 0)
    
syn_save_hist_hit(syn_sum+1) = syn_save_hist_hit(syn_sum+1) + 1;

if(stopsig==1) %miscorrected block
      num_miscor = num_miscor + 1
      not_done_outer = 1;
      lll = 0;
      while (lll < num_cw) && not_done_outer
        lll = lll + 1;
        if (sum(cw(lll,:) == c_hat) == n)
            not_done_outer = 0;
            cw_unique_count(lll) = cw_unique_count(lll) + 1;
        end
    end
    if (not_done_outer) %add new cw
        num_cw = num_cw + 1
        cw(num_cw,:) = c_hat;
        cw_unique_count(num_cw) = 1;
    end
%      cw(num_miscor,:) = c_hat;
else %look for TS
    [val, location] = min(num_USC);
    loc = find(c_hat_hist(location,:));
    % print current TS info to screen
    poop = [length(loc), val]
    TS_temp = zeros(1,n);
    TS_temp(loc) = 1;
    not_done_outer = 1;
    lll = 0;
    %make sure we have a "good" TS
    if (length(loc) > val)
    while (lll < num_TS) && not_done_outer
        lll = lll + 1;
        if (sum(TS(lll,:) == TS_temp) == n)
            not_done_outer = 0;
            TS_unique_count(lll) = TS_unique_count(lll) + 1;
        end
    end
    if (not_done_outer) %add new TS
        num_TS = num_TS + 1
        TS(num_TS,:) = TS_temp;
        TS_unique_count(num_TS) = 1;
    end
    end
end

end

end %if we decode

% weight = 1;
% num_frame_errors = num_frame_errors + weight;
% variance = variance + weight^2;

%      end % tot_perms_2 loop - uncomment along with above stuff
end %tot_perms loop


      end % num_combos loop
    end %v_count loop
end % num deg prof while
toc

TS = TS(1:num_TS,:);
cw = cw(1:num_cw,:);
%save Hfilename_TS TS cw
