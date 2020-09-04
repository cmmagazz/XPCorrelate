function [x,y] = transform_frm_img(fwd,x,y,T,transtype)
T = maketform(transtype,T);

if fwd == 1
    [x,y] = tformfwd(T,x,y);
    x = x - x(1);
    y = y - y(1);
elseif fwd==0
    [x,y] = tforminv(T,x,y);
    x = x - x(1);
    y = y - y(1);
elseif fwd==2
    [x,y] = tformfwd(T,x,y);
end

end
