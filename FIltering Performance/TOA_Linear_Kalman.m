function [esti_pos, esti_cov,Kalman_Gain, pred_cov] =...
    TOA_Linear_Kalman...
    (esti_pos_prev, esti_cov_prev,...
    E_w,Q,noise_var, measu_dist,esti_vel_prev)
    %this function is for Kalman Filter

    %codes implementation.
    %init_pos = [0;0];
    dt = 0.01;
    A = [1,0;0,1];  %for modified state ttransient matrix.
    ankx = 630; anky = 540;
    Anchor_1 = [0;0];
    Anchor_2 = [ankx;0];
    Anchor_3 = [0;anky];
    
    %measurement step should be implemented by TOA algorithm base.

    z = zeros(2,1);
    H = zeros(2,2);
    R = zeros(2,2);
    %should be renewal by each step.

    %%
    %set H, z, R first
    z(1,1) = (measu_dist(2,1)^2-measu_dist(1,1)^2)+...
        (Anchor_1(1,1)^2+Anchor_1(2,1)^2)-...
        (Anchor_2(1,1)^2+Anchor_2(2,1)^2);

    z(2,1) = (measu_dist(3,1)^2-measu_dist(1,1)^2)+...
        (Anchor_1(1,1)^2+Anchor_1(2,1)^2)-...
        (Anchor_3(1,1)^2+Anchor_3(2,1)^2);

    H(1,:) = 2*[(Anchor_1(1,1)-Anchor_2(1,1)),...
        (Anchor_1(2,1)-Anchor_2(2,1))];

    H(2,:) = 2*[(Anchor_1(1,1)-Anchor_3(1,1)),...
        (Anchor_1(2,1)-Anchor_3(2,1))];
                

    %{
    
    % involve error for calculate Kalman Gain

    z(3,1) = (measu_dist(3,1)^2-measu_dist(2,1)^2)+...
        (Anchor_2(1,1)^2+Anchor_2(2,1)^2)-...
        (Anchor_3(1,1)^2+Anchor_3(2,1)^2);
    H(3,:) = 2*[(Anchor_2(1,1)-Anchor_3(1,1)),...
        (Anchor_2(2,1)-Anchor_3(2,1))];
    %}

    
    %%
    %prediction step
    %use modified term.
    pred_pos = A*esti_pos_prev + esti_vel_prev.*dt+E_w;     %should get Exppectation of omega; process noise. by varience
    %pred_vel = esti_vel_prev;           %should calc vel after Filter algorithm ends.(After escaping function)
    pred_cov = A*esti_cov_prev*A' + Q;  %Q : modified Q(zero-mean)

    %%
    R(1,1) = measu_dist(1,1)^2+measu_dist(2,1)^2;
    R(1,2) = measu_dist(1,1)^2;     R(2,1) = R(1,2);
    %R(1,3) = -1*measu_dist(2,1)^2;     R(3,1) = R(1,3);

    R(2,2) = measu_dist(1,1)^2+measu_dist(3,1)^2;
    %R(2,3) = measu_dist(3,1)^2;     R(3,2) = R(2,3);

    %R(3,3) = measu_dist(2,1)^2+measu_dist(3,1)^2;

    R = R.*4.*noise_var;

    %%
    %calc Kalman Gain
    Kalman_Gain = pred_cov*H'*inv(H*pred_cov*H'+R);

    %%
    %estiamtion step
    esti_pos = pred_pos+Kalman_Gain*(z-H*pred_pos);
    esti_cov = pred_cov-Kalman_Gain*H*pred_cov;

    %% After function ends,
    % should calc Q and R and vel of estimation numerically by steps.(outside of the function.)
    % by each step, should approxi. vel by estimation pos. is it optimal...?
    % -> should compare with vel_LPF.

    % idea #1, constant Q
    % with const Q, should implement 
    % diagonal Q 
    % decrease of Q by param by gamma

    %idea #2
    % with time-varient Q.. should compare these.

    % idea #3
    % n square as Random Variable
    %{
    
        %results store by txt files.
        fileID = fopen('TOA_Kalman_pred_movPnt.txt', 'w');
        fprintf(fileID, '%.7f\n', error_memory);
        fclose(fileID);
    %}
end