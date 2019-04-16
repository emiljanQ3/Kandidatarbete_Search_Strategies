function [pos_res, indice] = splitPositionDataPartitioned(pos_a, n, myCircle)
    pos_res = zeros([1,2,n]);
    indice  = zeros([1,2]);

    k = 0;
    while k + n <= length(pos_a)
       [dPos, dIndice] = splitPositionData(pos_a(:,:,k + 1: k + n), myCircle);
       if(size(dIndice,1)>1)
            pos_res = [pos_res ; dPos];
       end
       indice  = [indice ; k + dIndice(1:end-1,:)];
       k = k + n;
    end
   
    %Last interation
    [dPos, dIndice] = splitPositionData(pos_a(:,:,k + 1:end), myCircle);
    if(size(dIndice,1)>1)
        dPos = cat(3, dPos , zeros(size(dPos,1), size(dPos,2), n - size(dPos,3)));
        pos_res = [pos_res ; dPos];
    end 
    indice  = [indice ; k + dIndice];
    
    %Remove first row
    pos_res = pos_res(2:end,:,:);
    indice  = indice(2:end,:);

    close all
end