%%%%%%%%%%%%%
% This program is a general Matlab LDPC decoder simulator for an AWGN channel.   
% Both regular and irregular LDPC code performance can be simulated.  The
% full belief propagation and min-sum approximation, using full double
% precision messages, are the two decoding algorithms currently
% implemented.
%
% The full belief propagation version of this code can simulate roughly 100
% Kb/s (coded bits).  Although this is not near as much throughput as most
% hardware simulations, the dominant TS making up the error floor are 
% collected is this software implementation.
%
% The codes to be simulated should be in H_sparse form


clear all
% function irregularSPA()
randn('state',sum(100*clock));

%%%%%
% These would be some example files containing H_sparse
% load H1000_4_8_gallager 
% load H504goodcode_no82
% load QCH2331
% load QCH8160 
% load ramanujanq17p5
% load H504_4_8_no4cycle
% load H96_4_8_no4cycle
% load 128x256regular.mat H
% load H1008_nox2
% load H96_4_8_70cw
% load H504_4_8_no6_fallback4_best
% load H26_reallysmall_no4cycle
% load test_1200H_3_6_no6
% load H1200_4_8_no6cycle
load H2048_8023an
% load H4000_wesel_test
% load H128_highrate_test
% load H500_rate08_3_15 
% load H500_rate08_3_15_best
% load H1000_rate08_4_20

% H_sparse = sparse(H);

n = length(H_sparse(1,:));
n_k = length(H_sparse(:,1));
k = n - n_k;
rate = k/n;

%%%%%%%% To generate G matrix for encoding uncomment the following. 
% This is generally not necessary because sending all zeros codeword ok for
% linear codes.  Requires files rearrange_cols, inv_GF2

% H = full(H_sparse);
% H = rearrange_cols(H);
% %swap full rank cols to last part
% temp = H(:, 1:n-k);
% H(:, 1:n-k) = H(:, n-k+1:n);
% H(:, n-k+1:n) = temp;
% C2=H(:, n-k+1:n);
% %only need to do this N^3 operation once
% C2_inv = inv_GF2(C2);
% H_system = mod(C2_inv*H, 2);
% P_trans = H_system(:, 1:n-k);
% G = cat(2, eye(k), P_trans');

%Max number of noisy messages sent per SNR
num_trials = 20000000;

%%%%%%
% These are MPA messages:
%  Lq = Variable->Check
%  Lr = Check->Variable
%  LQ = Marginal LLR used in making a hard decision.
Lq=zeros(n,max(sum(H_sparse,1)));
Lr=zeros(n-k,max(sum(H_sparse,2)));
LQ=zeros(1,n);

num_frame_errors = 0;
num_bit_err = 0;

%SNR values in DeciBels
% ebnodb = [1.5:0.5:3.0];
ebnodb = [3:.5:6];

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
%Keep a history of number of iterations necessary for each decoding.
ite_hist = zeros(1,num_trials);
%Number of UnSatisfied Checks at each iteration.
num_USC = zeros(1, itenum);
%Store which checks unsatisfied.
USC_vec = zeros(1, n-k);
c_hat_hist_save = cell(1,itenum);
USC_vec_save = cell(1,itenum);
c_hat_hist = zeros(itenum, n);

%Store Trapping Set, noise, and codeword information
num_TS = 0;
TS = zeros(100, n);
noise_save = zeros(100, n);
cw = zeros(10,n);

%small_num defines the magnitude of largest messages from C->V
small_num = eps;

%Number of miscorrections (MPA converges to incorrect valid codeword)
num_miscor = 0;

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

temp1 = zeros(size(Lr));
temp2 = zeros(size(Lq));
tic
%%%%%
% Main SNR loop
for ebnoloop = 1:length(ebnodb)
    ebno = 10^((ebnodb(ebnoloop))/10)
    esno = ebno*(rate);
    sigma_2 = 1/(2*esno);
    sigma = sqrt(sigma_2);
    Lc_coef = 4*esno;
    num_frame_errors = 0;
    num_bit_errors = 0;
    t = 0;
    %Simulate till we collect 20 errors or num_trials reached w/o getting
    % 20 errors.  This can be adjusted to affect quality of estimate.
while ((t < num_trials) && (num_frame_errors < 20))
  t = t + 1;
%%% If encoding necessary
%  u = (1-sign(rand(1,k)-0.5))/2;
%   a = mod(u*G, 2);
%   a_mod = -1*(a*2 - 1);

 noise = randn(1,n);
  %for easy encoding, send all 0's message.
 a_mod = ones(1,n);

   y = a_mod + sigma*noise;  %regular MC

%      MPA
% Initialization
Lc=Lc_coef*y;
num_total = 0;
for ii = 1:num_v_deg
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = Lc(num_total+1:num_total + deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2));
    num_total = num_total + deg_v_prof(ii, 1);
end


stopsig=0;  % If a valid codeword has been found, set this to `1'
itestep=0;  % Iteration number


while (stopsig==0) && (itestep<itenum)

    itestep=itestep+1;

% Update Lr

%%%%%% This code is for Min-Sum Algorithm
%% To switch between the min-sum and BP algorithm, just uncomment the
%% following lines and comment the lines below corresponding to the BP
%% implementation.

%     Lr(Lr_ind) = Lq(Lq_ind);
%     num_total = 0;
%     for ii = 1:num_c_deg
%         aterm = sign(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)));
%         atermprod=prod(aterm,2)*ones(1,deg_c_prof(ii,2));
%         Lr_min_mat=sort(abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))');
%         Lr_min_vals=Lr_min_mat(1,:)'*ones(1,length(Lr_min_mat(:,1)));
%         Lr_min_vals_2 = Lr_min_mat(2,:)'*ones(1,length(Lr_min_mat(:,1)));
%         Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
%             (atermprod.*aterm).*((((Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2))))-1)*-1).*Lr_min_vals + ...
%             (Lr_min_vals == abs(Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))).*Lr_min_vals_2);
%         num_total = num_total + deg_c_prof(ii, 1);
%     end
            
%%%%%% This code is for full Belief Propagation Algorithm
    temp1(Lr_ind) = Lq(Lq_ind);
    num_total = 0;
    for ii = 1:num_c_deg
        aterm = sign(temp1(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)));
        atermprod=prod(aterm,2)*ones(1,deg_c_prof(ii,2));
        bterm=-log(tanh(abs(temp1(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)))*0.5) + small_num); %try speed-up from built-in matlab tanh?
        btermsum=sum(bterm,2)*ones(1,deg_c_prof(ii,2));
        Lr(num_total+1:num_total+deg_c_prof(ii, 1), 1:deg_c_prof(ii, 2)) = ...
          (atermprod.*aterm).*-log(tanh(abs(btermsum-bterm)*0.5) + small_num);
        num_total = num_total + deg_c_prof(ii, 1);
    end

% Update Lq
temp2(Lq_ind) = Lr(Lr_ind);
num_total = 0;
for ii = 1:num_v_deg
    Lq_sum = sum(temp2(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)),2);
    LQ(num_total+1:num_total+deg_v_prof(ii, 1)) = Lq_sum' + Lc(num_total+1:num_total+deg_v_prof(ii, 1));
    Lq(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2)) = ...
        LQ(num_total+1:num_total+deg_v_prof(ii, 1))'*ones(1,deg_v_prof(ii, 2)) - ...
        temp2(num_total+1:num_total+deg_v_prof(ii, 1), 1:deg_v_prof(ii, 2));
    num_total = num_total + deg_v_prof(ii, 1);
end

% Codeword check

c_hat=(1-sign(LQ))/2;
USC_vec = mod(c_hat*H_sparse',2);
num_USC(itestep) = sum(USC_vec);

if num_USC(itestep)==0
  stopsig=1;
end

%Keep a history of hard decision for the current decoding.
c_hat_hist(itestep,:) = c_hat;

end

ite_hist(t) = itestep;

if (sum(c_hat) ~= 0)
    if(stopsig==1) %miscorrected block
      num_miscor = num_miscor + 1;
      cw(num_miscor,:) = c_hat;
    else %look for TS, as defined in TCOM paper
      [val, location] = min(num_USC);
      loc = find(c_hat_hist(location,:)~=0);
      TS_temp = zeros(1,n);
      TS_temp(loc(1:length(loc))) = 1;
      num_TS = num_TS + 1
      noise_save(num_TS,:) = noise;
      TS(num_TS,:) = TS_temp;
    end

    weight = 1;
    num_frame_errors = num_frame_errors + weight;
    num_bit_errors = num_bit_errors + sum(abs(c_hat-a_mod)); 

    noise_save(num_frame_errors,:) = noise;
    c_hat_hist_save{num_frame_errors} = c_hat_hist;
    USC_vec_save{num_frame_errors}  = USC_vec;
end

end %trials while

Pf_vec(ebnoloop) = (num_frame_errors/t);
Pb_vec(ebnoloop) = num_bit_errors/(t*n);

end %ebno loop

temp = zeros(num_TS,2);
temp(:,1) = sum(TS(1:num_TS,:),2);
temp(:,2) = sum(mod(TS(1:num_TS,:)*H_sparse',2),2)

toc
%%%%%% To plot Frame Error Rate and Bit Error Rate
% semilogy(ebnodb, Pf_vec, 'k',ebnodb, Pb_vec, 'b');
% title('LDPC (504, 252) Monte Carlo (20 errors collected)');
% xlabel ('E_b/N_o (dB)');
% ylabel ('FER/BER');

