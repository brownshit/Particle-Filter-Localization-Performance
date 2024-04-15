function H = H_linearized(pos)
    % return h vector.
    % Anc | 2x3 : state x Anc Num
    % pos | 2x1 : state
    %       For EKF
    H = zeros(3,2);
    ankx = 630; anky = 540;
    Anc = ...
        [0,ankx,0;
        0,0,anky];
    for ind1 = 1:1:3
        H(ind1,1) = (pos(1,1)-Anc(1,ind1))...
            /sqrt((Anc(1,ind1)-pos(1,1))^2+(Anc(2,ind1)-pos(2,1))^2);
        H(ind1,2) = (pos(2,1)-Anc(2,ind1))...
            /sqrt((Anc(1,ind1)-pos(1,1))^2+(Anc(2,ind1)-pos(2,1))^2);
    end
end