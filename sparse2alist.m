function sparse2alist(H,filename)

A = full(H);

M = size(A,1);
N = size(A,2);

ratebyrows = zeros(M,1);
ratebycols = zeros(N,1);

for i = 1:M
    rateinrow = find(A(i,:));
    ratebyrows(i) = size(rateinrow,2);
end


for i = 1:N
    rateincol = find(A(:,i));
    ratebycols(i) = size(rateincol,1);
end

maxinrows = max(ratebyrows);
maxincols = max(ratebycols);

fID = fopen(filename, 'w');
fprintf(fID,'%d\t', N, M);
fprintf(fID,'\n');
fprintf(fID,'%d\t', maxincols, maxinrows);
fprintf(fID,'\n');
fprintf(fID,'%d\t', ratebycols);
fprintf(fID,'\n');
fprintf(fID,'%d\t', ratebyrows);
fprintf(fID,'\n');

for i = 1:N
    rateincol = find(A(:,i));
    fprintf(fID,'%d\t', rateincol);
    fprintf(fID,'\n');
end

for i = 1:M
    rateinrow = find(A(i,:));
    fprintf(fID,'%d\t', rateinrow);
    fprintf(fID,'\n');
end

fclose(fID);


