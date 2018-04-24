function [ rs,ts ] = f_coefficients(ep1,ep2,kz1,kz2)
%   caculate frensnel coefficients of each layer
rs=(kz1-kz2)/(kz1+kz2);
ts=2*kz1/(kz1+kz2);
rp=(ep2*kz1-ep1*kz2)/(ep1*kz2+ep2*kz1);
tp=2*ep1*ep2*kz1/(ep1*kz2+ep2*kz1);


end

