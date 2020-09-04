
function mat = gridify_vector(vec,lenx,leny)
%     mat = zeros(lenx,leny);
    for y = 1:leny
        mat(y,:) = vec(lenx*(y-1) + 1:lenx*y); %%%%%%%%%%%%%%%%%
    end
end