function H = weigthtwoqc2sparse(qcHFileName)

fid = fopen(qcHFileName, 'r');

n = fscanf(fid, '%d', [1 1]);
m = fscanf(fid, '%d', [1 1]);
z = fscanf(fid, '%d', [1 1]);

I = sparse(eye(z));
Z = sparse(zeros(z));
H = sparse([]);

for i = 1:m
    lH = sparse([]);
    for j = 1:n
        shift2=0;
        shift = fscanf(fid, '%d', [1 1]);
        shift2 = fscanf(fid, '&%d', [1 1]);
        if shift == -1
            lH = [lH Z];
        else    
            if shift2~= 0
                lH = [lH circshift(I, [0 shift])+circshift(I, [0 shift2])];
            else
                lH = [lH circshift(I, [0 shift])];
            end
        end
    end
    H = [H; lH];
end

fclose(fid);

end
