function sparse2graph(H, filename)
    file = fopen(filename, 'w');
    
    [M, N] = size(H);
    fprintf(file, '%d %d\n', N, M);
    
    column_weights = full(sum(H,1));
    max_col = max(column_weights);
        
    row_weights = full(sum(H,2))';
    max_row = max(row_weights);
   
    for i = 1:M
        indexes = full(find(H(i, :)));
        pad_size = max_row - length(indexes);
        if  pad_size > 0
            %indexes = [indexes zeros(1, pad_size)];
            indexes = [indexes];
        end
        fprintf(file, '%d ', indexes);
        fprintf(file, '\n');
    end
    
    fclose(file);
end