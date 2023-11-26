% 基于密度排序
function rho_mean = kernel (temp_x,temp_pop)

    h = 1.1;

    % 阈值
    rho_min = exp(-h^2 / (2*h^2));
 
    if isempty(temp_pop)
        rho_mean = zeros(size(temp_x,1),1);
    else
        dist = pdist2(temp_x,temp_pop);

        rho = exp(-dist.^2 ./ (2 * h^2));
        rho(rho < rho_min) = 0;
        rho_mean = mean(rho, 2);
    end
end