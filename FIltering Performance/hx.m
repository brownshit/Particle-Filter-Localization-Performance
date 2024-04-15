function [pt_observ] = hx(pt_pred)
    %measurement equation.
    pt_observ = zeros(3,1);
    ankx = 630; anky = 540;
    ANK = ...
        [0,ankx,0;
        0,0,anky];
    for ind1 = 1:1:3
        pt_observ(ind1,1) = sqrt((ANK(1,ind1)-pt_pred(1,1))^2+(ANK(2,ind1)-pt_pred(2,1))^2);
    end
end