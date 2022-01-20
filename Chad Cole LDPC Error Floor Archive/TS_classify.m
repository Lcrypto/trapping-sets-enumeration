%% This software was developed by Chad Cole as a grad student while funded by L-3
%% Communications West.
%%% Must have H_sparse for the proper TS in the workspace for this to work.
%%%  Returns in 'temp' the (a,b) characteristics of the TS in TS structure
num_TS = length(TS(:,1));
temp = zeros(num_TS,2);
temp(:,1) = sum(TS(1:num_TS,:),2);
temp(:,2) = sum(mod(TS(1:num_TS,:)*H_sparse',2),2)
