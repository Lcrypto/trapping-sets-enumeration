[a, b] = parse_file_a_b( '20_4_128girth4emd2.MS.trap' );
TS_count_array=zeros(100,100);
for a_value=2:30
    for b_value=2:30
        I = find(a==a_value);
        J = find(b==b_value);
        Common=intersect(I,J);
        TS_count_array(a_value, b_value)=size(Common,2);
    end
end