function M = Poleface(angle,r)
%angle = -abs(angle);

M = [     1         0         0            0  0  0 ;
     tand(angle)/r  1         0            0  0  0 ;
          0         0         1            0  0  0 ;
          0         0  -tand(angle)/r      1  0  0 ;
          0         0         0            0  1  0 ;
          0         0         0            0  0  1];