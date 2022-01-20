   function [binvec] = int2vec(val, len)
  %set a maximum input value 2^{len}  - 1
  binvec = zeros(1,len);
  for i = len-1:-1:0
      if (val >= 2^i)
          val = val - 2^i;
          binvec(i+1) = 1;
      else
          binvec(i+1) = 0;
      end
  end
%  return binvec;
  
          
      