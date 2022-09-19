function [xx,yy]=utm_geo(rlon,rlat,UTM_PROJECTION_ZONE)

% convert to Cartesian

% rlat=29;
% rlon=105;
% UTM_PROJECTION_ZONE=48;


PI = 3.141592653589793d0;
degrad=PI/180.d0;
raddeg=180.d0/PI;
semimaj=6378206.4d0;
semimin=6356583.8d0;
scfa=0.9996d0;
north=0.d0;
east=500000.d0;
rx = rlon;
 ry = rlat;
 % save original parameters
  rlon_save = rlon;
  rlat_save = rlat;
  rx_save = rx;
  ry_save = ry;
  
  % define parameters of reference ellipsoid
  e2=1.0-(semimin/semimaj)^2.0;
  e4=e2*e2;
  e6=e2*e4;
  ep2=e2/(1.-e2);
 
  dlon = rlon;
  dlat = rlat;
    
  zone = double(UTM_PROJECTION_ZONE);
  cm = zone*6.0 - 183.0;
  cmr = cm*degrad;
  
    rlon = degrad*dlon;
  rlat = degrad*dlat;

  delam = dlon - cm;
  if (delam < -180.) 
      delam = delam + 360.;
  end
  if (delam > 180.) 
      delam = delam - 360.;
  end
  delam = delam*degrad;

  f1 = (1. - e2/4. - 3.*e4/64. - 5.*e6/256)*rlat;
  f2 = 3.*e2/8. + 3.*e4/32. + 45.*e6/1024.;
  f2 = f2*sin(2.*rlat);
  f3 = 15.*e4/256.*45.*e6/1024.;
  f3 = f3*sin(4.*rlat);
  f4 = 35.*e6/3072.;
  f4 = f4*sin(6.*rlat);
  rm = semimaj*(f1 - f2 + f3 - f4);
  if (dlat == 90. | dlat == -90.) then
    xx = 0.;
    yy = scfa*rm;
  else
    rn = semimaj/sqrt(1. - e2*sin(rlat)^2);
    t = tan(rlat)^2;
    c = ep2*cos(rlat)^2;
    a = cos(rlat)*delam;

    f1 = (1. - t + c)*a^3/6.;
    f2 = 5. - 18.*t + t^2 + 72.*c - 58.*ep2;
    f2 = f2*a^5/120.;
    xx = scfa*rn*(a + f1 + f2);
    f1 = a^2/2.;
    f2 = 5. - t + 9.*c + 4.*c^2;
    f2 = f2*a^4/24.;
    f3 = 61. - 58.*t + t^2 + 600.*c - 330.*ep2;
    f3 = f3*a^6/720.;
    yy = scfa*(rm + rn*tan(rlat)*(f1 + f2 + f3));
  end
  xx = xx + east;
  yy = yy + north;
