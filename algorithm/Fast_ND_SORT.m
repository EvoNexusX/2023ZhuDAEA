function [index,F] = Fast_ND_SORT(obj)
    num = size(obj,2);
    obj(1,:) = -obj(1,:);
    obj1 = obj(:,obj(2,:)==0);
    
    obj2 = obj(:,obj(2,:)~=0);
    F = ones(1,num);
    [F2,~] = NDSort(obj2',num);
    idx2 = obj(2,:)~=0;
    F(idx2) = F(idx2) + F2;

    obj_new = zeros(num,2);
    
    obj_new(:,1) = F';
    obj_new(:,2) = obj(1,:)';

    [~,index] = sortrows(obj_new);

end