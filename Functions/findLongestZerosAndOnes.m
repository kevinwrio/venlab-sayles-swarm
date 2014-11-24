function [found,from,to]=findLongestZerosAndOnes(c,whatToFind)
%hack by Nasser M. Abbasi, sept 2011

%INPUT: c, arrray c of logicals
% whatToFind, true or false (to look for 1's or 0's sequence)
%OUTPUT:
% found: true if found at least one sequence
% from : sequence starts at
% to : sequence ends at
%

% error checking on arguments goes here

found = true;
from=[]; to=[];

r = arrayfun(@(i) countFrom(c,i,whatToFind),1:length(c));
if any(r)
     z = find(r==max(r));
     from = z(1);
     to = r(z(1))+z(1)-1;
else
     found = false;
end

end

function k = countFrom(c,i,good)
k=0;
if c(i)==good
     while c(i)==good
         k=k+1;
         i=i+1;
         if i>length(c)
             break
         end
     end
end
k;

end
%-------------------------- 