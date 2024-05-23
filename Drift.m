function M = Drift(L)
% transfer matrix for a drift of length L
M = [1 L 0 0 0 0; 
     0 1 0 0 0 0; 
     0 0 1 L 0 0; 
     0 0 0 1 0 0;
     0 0 0 0 1 0;
     0 0 0 0 0 1];