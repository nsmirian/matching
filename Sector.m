function M = Sector(angle,r)
% horizontal sector dipole magnet
arad = angle*pi/180; 
l = abs(r*arad);

M = [cosd(angle)         r*sind(angle) 0  0  0 r*(1-cosd(angle))   ;
     -1/r*sind(angle)     cosd(angle)  0  0  0   sind(angle)       ;
          0                  0         1  l  0      0              ;
          0                  0         0  1  0      0              ;
    -sind(angle)    -r*(1-cosd(angle)) 0  0  1 -(arad-sin(arad))*r ;
          0                  0         0  0  0      1              ];

