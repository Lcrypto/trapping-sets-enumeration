function [ a, b] = parse_file_a_b( filename )
    a = [];
    b = [];
    file = fopen(filename);
    
    while 1
        str = fgets(file);
        if str == -1
            break;
        end
        
        start = find(str == '(');
        start = start + 1;
        end_ = find(str == ',');
        end_ = end_ - 1;
        str_value = str(start:end_);
        a = [a str2num(str_value)];
        
        start = find(str == ',');
        start = start + 1;
        end_ = find(str == ')');
        end_ = end_ - 1;
        str_value = str(start:end_);
        b = [b str2double(str_value)];
    end
end