function SP = mat_gen(matsize,numspec)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function created to fill a matrix with a number of (plant) populations
%
% input is:
% matsize (1D size of matrix)
% numspec (number of populations)
%
% output is:
% SP ( matrix with populations)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


boundaries=linspace(0,1,(numspec+1));
SP2 = rand(matsize);
SP = zeros(matsize);
if numspec == 2
SP(SP2 <= boundaries(2)) = 1;
SP(SP2 > boundaries(2) & SP2 <= boundaries(3)) = 2;
end
if numspec == 3
SP(SP2 <= boundaries(2)) = 1;
SP(SP2 > boundaries(2) & SP2 <= boundaries(3)) = 2;
SP(SP2 > boundaries(3) & SP2 <= boundaries(4)) = 3;
end
if numspec == 4
SP(SP2 <= boundaries(2)) = 1;
SP(SP2 > boundaries(2) & SP2 <= boundaries(3)) = 2;
SP(SP2 > boundaries(3) & SP2 <= boundaries(4)) = 3;
SP(SP2 > boundaries(4) & SP2 <= boundaries(5)) = 4;
end
if numspec == 5
SP(SP2 <= boundaries(2)) = 1;
SP(SP2 > boundaries(2) & SP2 <= boundaries(3)) = 2;
SP(SP2 > boundaries(3) & SP2 <= boundaries(4)) = 3;
SP(SP2 > boundaries(4) & SP2 <= boundaries(5)) = 4;
SP(SP2 > boundaries(5) & SP2 <= boundaries(6)) = 5;
end