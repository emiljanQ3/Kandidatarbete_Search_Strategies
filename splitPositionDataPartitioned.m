function [pos_res, indice] = splitPositionDataPartitioned(pos_a, n, myCircle)
   pos_res = zeros([1,2,n]);
   indice  = zeros([1,2]);
   
   k = 1;
   while k + n - 1 <= length(pos_a)
       [dPos, dIndice] = splitPositionData(pos_a(:,:,k: k + n - 1), myCircle);
       pos_res = [pos_res ; dPos];
       indice  = [indice ; dIndice(1:end-1,:)];
       
       k = k + n;
   end
   
   %Last interation
   [dPos, dIndice] = splitPositionData(pos_a(:,:,k:end), myCircle);
   dPos = cat(3, dPos , zeros(size(dPos,1), size(dPos,2), n - size(dPos,3)));
   pos_res = [pos_res ; dPos];
   indice  = [indice ; dIndice];
   
   %Remove first row
   pos_res = pos_res(2:end,:,:);
   indice  = indice(2:end,:);
   
end