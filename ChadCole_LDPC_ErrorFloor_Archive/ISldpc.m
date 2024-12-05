%%%%%%%%%%%%%
%% ISldpc_ms_new_doc.m - Requires input H_sparse matrix.
%%
%%  This program implements the last step of the 3-step program outlined
%%  in "A General Method for Finding Low Error Rates of LDPC Codes."  The
%%  last step uses the technique of importance sampling (IS) to get a
%%  better estimate than Monte Carlo (MC) can provide of the FER/BER at 
%%  higher SNR's.
%%
%%  The core decoder implementation
%%  is the same as in the program 'irregularSPA.m,' another piece of
%%  software in this suite.
%%%%%%%%%%%%

clear all
randn('state',sum(100*clock));
format compact
%function ISldpc_ms_new_doc()

% load H1008_nox2
% load H1200_4_8_no6cycle
% load H2048_8023an
% load H504_no82
H = [1 1 1 1 0 0 0 0 0 0;
     1 0 0 0 1 1 1 0 0 0;
     0 1 0 0 1 0 0 1 1 0;
     0 0 1 0 0 1 0 1 0 1;
     0 0 0 1 0 0 1 0 1 1]
H_sparse = sparse(H);


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


small_num = eps;
%%% Set SNR's for a desired IS performance estimate.
% ebnodb = [4];
% ebnodb = [4 4.5 5 5.5 6];
ebnodb = 4:1:24;

%% The following variables collect some information which can help
%% determine how accurate the IS estimate is.
num_new_error_event = 0;
%% The error_event variable keeps track of new errors that occur throughout
%% the simulation that were not included in the initial list obtained from
%% the TSearch.m program.
error_event = spalloc(1000000,n,10000*30);
%% This counts how many occurrences of each new error events occurs.
error_event_count = zeros(1,1000);

itenum=50;  % maximum # of iterations


% load d_eH504_no82_del36_1
% load H504_no82_del36_op8_ebno5_TS
% load H2048_8023an_cw_TS
load H10_test

%% Only use the TS with d_e^2 < some threshold
% TS = TS(find(d_e_2 < 20),:);

if (exist('TS') == 0) || isempty(TS) 
    num_TS = 0;
    TS = [];
else
    num_TS = length(TS(:,1));
    temp = zeros(num_TS,2);
    temp(:,1) = sum(TS(1:num_TS,:),2);
    temp(:,2) = sum(mod(TS(1:num_TS,:)*H_sparse',2),2);
end

if exist('cw') == 0
    cw = [];
end
    

%% Another way to pick certain TS based on their (a,b) value, for example
%% in the following we just use (10,2) TS.
% TS = TS(intersect(find(temp(:,1)==10), find(temp(:,2)==2)),:);

%% mu is the list of mean-shift candidates in the f^* used in IS
mu = [TS; cw];
mu = sparse(mu);
num_mu = length(mu(:,1));

%% Keep track of any new codewords found throughout simulation.
num_miscor = 0;
if isempty(cw)
    num_cw = 0;
else
    num_cw = length(cw(:,1));
end

%% mu_coef specifies the amount we mean shift in each bit of a TS.  The
%% optimal value is not `1' unless we're shifting towards a valid codeword,
%% but for most TS, `1' will suffice.
mu_coef = 1.0;

%% We typically generate noise realizations shifted towards each of the mu
%% bias points an equal amount of time, for example 5000 in this case.
num_trials = 400*num_mu*ones(1,length(ebnodb));

%% The FER/BER estimates for each SNR.
Pf_vec = zeros(1, length(ebnodb));
Pb_vec = zeros(1, length(ebnodb));

%% The hit_rate counts the number of `hits' (errors) for each shift point at each SNR.
%% We can use this info as an indicator of how good a certain shift point
%% is.  If no hits were made in 5000 realizations of noise shifted towards
%% a certain TS, then that TS is probably not a significant contributor to
%% the error rate.  `intended hits' occur when the error is the exact same one
%% as that which we shifted towards.  These are the most common types of
%% hits and as SNR increases, the number of `intended hits' should approach
%% the number of total hits.
hit_rate = zeros(length(ebnodb), num_mu);
hit_rate_intended = zeros(length(ebnodb), num_mu);

%% ite_hist keeps a history of the number of iterations needed for each
%% decoding.
ite_hist = zeros(1,num_trials(1));
%% ite_hist_mean contains average number of iterations required at each SNR
ite_hist_mean = zeros(1, length(ebnodb));

for ebnoloop = 1:length(ebnodb)
    ebno = 10^((ebnodb(ebnoloop))/10)
    esno = ebno*(rate);%*3;  %rate 1/2, 3 bits/sym
    sigma_2 = 1/(2*esno);
    sigma = sqrt(sigma_2);
    Lc_coef = 4*esno;
    a = ones(1,n);
    variance = 0;
    num_frame_errors = 0; %don't forget to reset these to 0 for each ebno
    num_bit_errors = 0;
    t = 0;
    coef = -1/(2*sigma_2);
    %% s_f is the scale factor necessary to prevent exp(M) evaluating to
    %% zero (in Matlab) when we calculate the weight function.
    s_f = (n/2);

while (t < num_trials(ebnoloop))
 
  t = t + 1;
  %% For easy encoding, send all 0's (s0) message.
  %% The initial point in decoding space is biased towards one of the
  %% vectors in the mu matrix.
  y = a - mu_coef*mu(mod(t,num_mu)+1,:) + sigma*randn(1,n);
  %% This is different than the typical Monte Carlo simulation which would
  %% send instead:
  %% y = a + sigma*randn(1,n);

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

%%% Min-sum decoding:

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
%         num_total = num_total + deg_c_prof(ii, 1);
%     end

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

c_hat_hist(itestep,:) = c_hat;
end  %while

ite_hist(t) = itestep;

%% We have a `hit' (error) if final decoding state not the all-zeros
%% codeword.
if (sum(c_hat) ~= 0)
    %% Update hit_rate variable.
  hit_rate(ebnoloop, mod(t,num_mu)+1) = hit_rate(ebnoloop, mod(t,num_mu)+1) + 1;
  %% Determine error event which occurred.
  [val, location] = min(num_USC);
  loc = find(c_hat_hist(location,:));
  TS_temp = zeros(1,n);
  TS_temp(loc) = 1;
%% Calculate weight function (both for frame and bit error) for each `hit'
  denom = 0;
  for p = 1:num_mu
     denom = denom + exp(coef*norm(y - a + mu_coef*mu(p,:))^2 + s_f);
  end
  weight = num_mu*exp(coef*norm(y - a)^2 + s_f)/denom;
  bit_weight = weight*(sum(c_hat)/n);

  num_frame_errors = num_frame_errors + weight;
  num_bit_errors = num_bit_errors + bit_weight;
  %% This is on-line variance estimate.
  variance = variance + weight^2;

  if(stopsig==1) %miscorrected block
      num_miscor = num_miscor + 1;
      %% Build a list of new codewords discovered with IS that were not
      %% included in initial mu list.
    yes = 1;
    for jj = 1:num_cw
        if sum(cw(jj,:) == c_hat) == n
            yes = 0;
        end
    end
    if (yes)
            num_cw = num_cw + 1;
            cw(num_cw,:) = c_hat;
    end
  end
  

if (sum(TS_temp == mu(mod(t,num_mu)+1,:))==n)
     hit_rate_intended(ebnoloop, mod(t,num_mu)+1) = hit_rate_intended(ebnoloop,mod(t,num_mu)+1) + 1;

else
    %% Keep track of all new error events and their frequency of
    %% occurrence.  The following lines of code in this else clause can be
    %% commented out if this information is not required.  It does add some
    %% overhead, especially for longer codes.
    yes = 1;
    for jj = 1:num_mu
        if sum(mu(jj,:) == TS_temp) == n
            yes = 0;
        end
    end
    if (yes)
    for jj = 1:num_new_error_event 
        if sum(error_event(jj,:) == TS_temp) == n
            yes = 0;
            error_event_count(jj) = error_event_count(jj) + 1;
        end
    end
    if (yes)
            num_new_error_event = num_new_error_event + 1;
            error_event(num_new_error_event,:) = TS_temp;
            error_event_count(num_new_error_event) = 1;
    end
    end
end

end

end %trials while
 
Pf_vec(ebnoloop) = (num_frame_errors/t);
Pb_vec(ebnoloop) = num_bit_errors/t;
var_vec(ebnoloop) = variance/t;
ite_hist_mean(ebnoloop) = mean(ite_hist(1:t));

end %ebno loop

%% std is the online (sample) estimate for the standard deviation of our
%% estimate.
std = sqrt((var_vec - Pf_vec.^2)/t);

% save(['H1200_4_8_IS'],'Pf_vec','var','ebnodb','num_trials');

%% The next lines plot the FER estimate with circles above and below one standard
%% deviation (using the online calculation of variance).
% semilogy(ebnodb, Pf_vec, 'k-', ebnodb, Pf_vec + std, 'o', ebnodb, Pf_vec - std, 'o');
% title('LDPC (204, 102), Trials/SNR - 4000 for m-s; 10^8 for MC');
% legend('Mackay (204,102)');%, +/-1 std conf interval','Q(sqrt(2*esno*dminH))','Location','SouthWest');
% % axis([3 8 1e-14 1]);
% xlabel ('E_b/N_o (dB)');
% ylabel ('FER');
% grid on;

