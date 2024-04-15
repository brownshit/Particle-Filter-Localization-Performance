function [pt_pred] = fx(pt_prev,dt,prev_pt_vel,E_w)
    % size      : 4xNeff | p_x; v_x; p_y; v_y
    % linear system dynamic equation.
    
    F = eye(2);
    pt_pred = F*pt_prev + dt*prev_pt_vel + E_w;
end