function [M,globalM] = magnitude_of_vels(U,V,W)

M = zeros(size(U));
globalM = zeros(size(U,4), 1);
for t=1:size(U,4)
    M(:,:,:,t) = sqrt(U(:,:,:,t).^2 + V(:,:,:,t).^2 + W(:,:,:,t).^2); 
    globalM(t) = sum( M(:,:,:,t) , 'all') / sum( (M(:,:,:,t) > 0), 'all'); 
end

end