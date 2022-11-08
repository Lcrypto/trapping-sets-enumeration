 function [M] = alist2sparse2(fname)


fid = fopen(fname); 
Cols = fscanf(fid,'%d',1); 
Rows = fscanf(fid,'%d',1);

A = zeros(Rows, Cols);

maxincols = fscanf(fid,'%d',1);  
maxinrows = fscanf(fid,'%d',1);

rateincols = fscanf(fid,'%d',[1 Cols]);
rateinrows = fscanf(fid,'%d',[1 Rows]);

for i = 1:Cols

    for j = 1:rateincols(i)

            k = fscanf(fid,'%d',1);
            A(k,i) = 1;

    end

end

M = sparse(A);
