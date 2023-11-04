%%%%%%%%%
% This file contains the implementation of Step-2 of the 3-Step error floor
%   analysis procedure.  The input is an H matrix in sparse format called
%   H_sparse.  The output variable is d_e_2 which should be saved and used
%   in step 3, which is implemented in ISldpc.m.

clear all

% load H4000_wesel_test

% load H504_4_8_no6_fallback4_best
% load test_1200H_3_6_no6

% load H1000_4_8_gallager
% load H1008_nox2
% load QR_H8192 
% load Margulis11qHCV
% load H1200_4_8_no6cycle
% load H504goodcode_no82
% load H1008_peg
load H96_firstmackay

n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

%%%%%%
% These are MPA messages:
%  Lq = Variable->Check
%  Lr = Check->Variable
%  LQ = Marginal LLR used in making a hard decision.
Lq=zeros(n,max(sum(H_sparse,1)));
Lr=zeros(n-k,max(sum(H_sparse,2)));
LQ=zeros(1,n);

%%%% Organize H columns from highest degree -> lowest.  this is necessary
%%%% for the decoding algorithm to make use of Matlab vectorized commands
%%%% which significantly speed up program run-time.
H_col_sum = sum(H_sparse,1);
H_row_sum = sum(H_sparse,2);
%B will be sorted in ascending order
B = unique(H_col_sum);
num_v_deg = length(B);
R = unique(H_row_sum);
num_c_deg = length(R);
deg_c_prof = zeros(length(R), 2);
deg_v_prof = zeros(length(B), 2);

%%%%
% Set up degree profile info for V-nodes & C-nodes and make sure H is
%   ordered from highest to lowest degree in columns and rows.
%
% Example - n-k = 500, half parities have deg 6, half 8
% always list deg profiles in DESCENDING ORDER
%         : deg_c_prof = [250, 8;
%                         250, 6]
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
    
%%%%%%% The Clist Vlist Data Structures
% Clist is an n by (max(deg_v)+1) matrix with the first column having the
% variable node degree (deg_v) of each column of H.  The next deg_v entries
% for that row are the check node numbers which that variable node is
% connected to.  Any extra columns following are set to zero.
% Vlist is an (n-k) by (max(deg_c)+1) matrix with the first column having the
% check node degree (deg_c) of each row of H.  The next deg_c entries
% for that row are the variable node numbers which that check node is
% connected to.  Any extra columns following are set to zero.
 
Vlist=zeros(n-k,max(sum(H_sparse,2))+1); 
Clist=zeros(n,max(sum(H_sparse,1))+1); 

for jj=1:n-k
    Vlist(jj,1)=sum(H_sparse(jj,:));
    icnt=0;
    for ii=1:n
        if H_sparse(jj,ii)==1
            icnt=icnt+1;
            Vlist(jj,icnt+1)=ii;
        end
    end
end

for ii=1:n
    Clist(ii,1)=sum(H_sparse(:,ii));
    jcnt=0;
    for jj=1:n-k
        if H_sparse(jj,ii)==1
            jcnt=jcnt+1;
            Clist(ii,jcnt+1)=jj;
        end
    end
end

%Maximum number of MPA iterations before we declare a failure.
itenum = 50;

small_num = eps; %how expensive is eps()? 

% Simple H matrix consistency check.
num_edges_v = sum(Clist(:,1));
num_edges_c = sum(Vlist(:,1));
if (num_edges_c ~= num_edges_v)
    text = ['error - row/col mismatch in H']
end
num_edges = num_edges_c;
clear num_edges_c num_edges_v

Lr_ind = zeros(1, num_edges);
Lq_ind = zeros(1, num_edges);

%%%%%
% Set up the edge connections between Lr (C->V messages) and Lq (V->C).
%  Remember, Matlab matrices are indexed by columns, so if Lr is of
%  dimension (n-k)X max(d_c)+1, then Lr(n-k+1) is the same as Lr(1,2).
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



%load H1008_peg_del4_op6_ebno6_TS
%load H504_peg_del3_op8_ebno6_TS
%load H504_no82_del36_op8_ebno5_TS
%load H1200_4_8_no6cycles_5bitimpulse_TS
% load H1008_nox2_4delta_07op_7db_TS 
% load QR_8192_4bit_TS
% load Margulis_TS_tot

% load H504_4_8_no6_fallback4_best_TS
% load test_1200H_3_6_no6_TS
% load dinoi_TS_save
% load H204_TS
% load H1200_4_8_no6cycles_5bitimpulse_TS_extra_2
% load H1200_4_8_no6cycles_5bitimpulse_TS
% load H504_no82_TS_temp
% load H4000_wesel_test_TS_cw
% load H1008_peg_del3_op6_ebno6_TS
load H96-964_5_1_TS

num_TS = length(TS(:,1));
temp = zeros(num_TS,2);
temp(:,1) = sum(TS(1:num_TS,:),2);
temp(:,2) = sum(mod(TS(1:num_TS,:)*H_sparse',2),2)

% TS = TS(find(d_e < 30),:);
TS = TS(intersect(find(temp(:,1)==5), find(temp(:,2)==1)),:);


e_max = 3.5;
e_min = 1;
num_intervals = 10;

mu_coef = 1;

num_bits = zeros(length(ebnodb), length(TS(:,1)));

ebnodb = 6;
% ebnodb = [3 4 5 6 7];
a = ones(1,n);
e_vec = zeros(length(TS(:,1)), length(ebnodb), 2);


for outer = 1:length(TS(:,1))

    mu = TS(outer,:);

    for ebnoloop = 1:length(ebnodb)

    upper = e_max
    lower = e_min;

    ebno = 10^((ebnodb(ebnoloop))/10);
    esno = ebno*(rate);%*3;  %rate 1/2, 3 bits/sym
    sigma_2 = 1/(2*esno);
    sigma = sqrt(sigma_2);
    Lc_coef = 4*esno;
    
    variance = 0;
    num_frame_errors = 0; %don't forget to reset these to 0 for each ebno

 for e_loop = 1:num_intervals
     mu_coef = lower + (upper - lower)/2;

     y = a - mu_coef*mu;

  %      MPA
% step 1: initialization
Lc=Lc_coef*y;
num_total = 0;
for ii = 1:num_v_deg
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = Lc(num_total+1:num_total + deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2));
    num_total = num_total + deg_v_prof(ii, 1);
end


stopsig=0;  % check whether a valid codeword has been found
itestep=0;  

while (stopsig==0) && (itestep<itenum)

itestep=itestep+1;

% step 2: update Lr
%      Lr(Lr_ind) = Lq(Lq_ind);
%     num_total = 0;
%     for ii = 1:num_c_deg
%         aterm = sign(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)));
%         atermprod=prod(aterm,2)*ones(1,deg_c_prof(ii,2));
%         Lr_min_mat=sort(abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))');
%         Lr_min_vals=Lr_min_mat(1,:)'*ones(1,length(Lr_min_mat(:,1)));
%         Lr_min_vals_2 = Lr_min_mat(2,:)'*ones(1,length(Lr_min_mat(:,1)));
%         Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
%             Lc_alpha(alpha_loop)*(atermprod.*aterm).*((((Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2))))-1)*-1).*Lr_min_vals + ...
%             (Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))).*Lr_min_vals_2);
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
        num_total = num_total + deg_c_prof(ii, 1);
    end

% step 3 & 4: update Lq, LQ
Lq(Lq_ind) = Lr(Lr_ind);
num_total = 0;
for ii = 1:num_v_deg
    Lq_sum = sum(Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)),2);
    LQ(num_total+1:num_total+deg_v_prof(ii, 1)) = Lq_sum' + Lc(num_total+1:num_total+deg_v_prof(ii, 1));
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = ...
        (LQ(num_total+1:num_total+deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2)) - ...
        Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)));
    num_total = num_total + deg_v_prof(ii, 1);
end


% step 5: validity check
c_hat=(1-sign(LQ))/2;
USC_vec = mod(c_hat*H_sparse',2);
num_USC(itestep) = sum(USC_vec);
if num_USC(itestep) ==0
  stopsig=1;
end %if

% LQ_hist_save_exact{outer}(itestep,:) = c_hat;
% [val, q] = sort(LQ);
% LQ_sort_hist{outer}(itestep,:) = q;

c_hat_hist(itestep,:) = c_hat;
end  %while

if (sum(c_hat) ~= 0)

num_bits(ebnoloop, outer) = sum(c_hat);

[val, location] = min(num_USC);
loc = find(c_hat_hist(location,:));

TS_temp = zeros(1,n);
TS_temp(loc) = 1;

%%%%%If want the distance to the error boundary corresponding to the
%%%%%specific TS we're shifting to, then uncomment the next if statement.  
%%%%%When commented, it will return the distance to ANY error, which is 
%%%%%probably the most practically important.

% if (sum(TS_temp == mu)==n)
   upper = lower + (upper - lower)/2;
% else
%    lower = lower + (upper - lower)/2;
% end    
else
   lower = lower + (upper - lower)/2;
end

end %intervals loop
  e_vec(outer, ebnoloop, 1) = upper;
  e_vec(outer, ebnoloop, 2) = lower;
end %ebno loop

end %outer


d_e_2 = zeros(length(TS(:,1)), length(ebnodb));
for jj = 1:length(ebnodb)
   for ii=1:length(TS(:,1))
      temp = full(sum(TS(ii,:)));
      d_e_2(ii, jj) = temp*e_vec(ii,jj,2)^2;
   end
end

% save H4000_wesel_test_d_e_2 d_e_2