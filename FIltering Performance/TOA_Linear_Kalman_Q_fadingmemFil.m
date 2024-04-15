function [esti_pos, esti_cov,Kalman_Gain, pred_cov, Q_pres, pred_pos, sub_opt_pos,Alpha] =...
    TOA_Linear_Kalman_Q_fadingmemFil...
    (esti_pos_prev, esti_cov_prev,...
    Q_prev,noise_var, measu_dist,esti_vel_prev,alpha,a,E_w,...
    Q_max, Q_mean)
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
       
    %%
    R(1,1) = measu_dist(1,1)^2+measu_dist(2,1)^2;
    R(1,2) = measu_dist(1,1)^2;     R(2,1) = R(1,2);
    %R(1,3) = -1*measu_dist(2,1)^2;     R(3,1) = R(1,3);

    R(2,2) = measu_dist(1,1)^2+measu_dist(3,1)^2;
    %R(2,3) = measu_dist(3,1)^2;     R(3,2) = R(2,3);

    %R(3,3) = measu_dist(2,1)^2+measu_dist(3,1)^2;

    R = R.*4.*noise_var;
   
    %%
    %prediction step
    %use modified term.
    pred_pos = A*esti_pos_prev + esti_vel_prev.*dt;     %should get Exppectation of omega; process noise. by varience
    %pred_vel = esti_vel_prev;           %should calc vel after Filter algorithm ends.(After escaping function)

    %%
    %before we step into pred. cov, we need to calculate Q matrix.
    %instantaneous process noise

    %=============================================================
    % Revised | 2024.02.11

    % For testing
    TOA_pos = TOA(measu_dist);

    %Q matrix
    %define; sub optimal state x
    
    %operator = inv(pinv(H)*R*(pinv(H)')) + Q_prev;
    
    %lower equation could take more comp. cost of our alg.
    % gotten by Maximum likelihood Estimation. could try
    % Relative Maximum likelihood
            
    %state_sub_op = (pred_pos' + TOA_pos'*operator)*(eye(2)+operator);

%%
    %Proposed
    %instant_w = state_sub_op - pred_pos;

    %Substitude state_sub_op as TOA_pos
    %esti_epsilon_{k-1}
    
    esti_epsilon_prev = zeros(size(pred_pos,1),1);
    for ind = 1:1:size(pred_pos,1)
        %draw samples from prev. information.
        esti_epsilon_prev(ind,1) = esti_pos_prev(ind,1)+sqrt(esti_cov_prev(ind,ind))*randn;
    end
%{

    esti_epsilon_prev = zeros(size(pred_pos,1),2^size(pred_pos,1)+1);
    esti_epsilon_prev(1,2) = sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(1,3) = sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(1,4) = -1*sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(1,5) = -1*sqrt(esti_cov_prev(ind,ind));
    
    esti_epsilon_prev(2,2) = sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(2,3) = -1*sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(2,4) = sqrt(esti_cov_prev(ind,ind));
    esti_epsilon_prev(2,5) = -1*sqrt(esti_cov_prev(ind,ind));
    
    instant_w_sub = zeros(2,2^size(pred_pos,1)+1);
    Alpha = zeros(1,2^size(pred_pos,1)+1);
    Q_pres_sub = zeros(2,2,2^size(pred_pos,1)+1);

    for vari = 1:1:2^size(pred_pos,1)+1 %5
        instant_w_sub(:,vari) = TOA_pos - pred_pos - A*esti_epsilon_prev(:,vari);       %revision/ 24.02.19
                                       %A*esti_epsilon_prev; 's appropriate processing is needed     
        %this TOA_pos value could be substitude into sub optimal
        %position value.
        %=============================================================
        if alpha == 100
            Alpha(1,vari) = 1/(1+exp(-a*(abs(instant_w_sub(:,vari))-E_w_max/2)));      %adaptible Q_varient
            % abs(instant_w) means, we need to imply the w's amplitude on
        else
            Alpha(1,vari) = alpha;
        end
        Q_pres_sub(:,:,vari) = Alpha(1,vari).*instant_w_sub(:,vari)*instant_w_sub(:,vari)' + (1-Alpha(1,vari)).*Q_prev;
        
    end
%}        
    instant_w_sub = TOA_pos - pred_pos;% - A*esti_epsilon_prev;       %revision/ 24.02.19
                                       %A*esti_epsilon_prev; 's appropriate processing is needed     
    %this TOA_pos value could be substitude into sub optimal
    %position value.
    %=============================================================
    
    if alpha == 101
        Alpha = 1/(1+exp(-a*(abs(instant_w_sub)-E_w/2)));      %adaptible Q_varient
        Q_pres_sub = Alpha.*Q_max + (1-Alpha).*Q_mean;
        % Q should be determined by...
        %   high curvature : big Q
        %   low  curvature : small Q
    else
        if alpha == 100
            Alpha = 1/(1+exp(-a*(abs(instant_w_sub)-E_w/2)));      %adaptible Q_varient
            % abs(instant_w) means, we need to imply the w's amplitude on
        else
            Alpha = alpha;
        end
        Q_pres_sub = Alpha.*instant_w_sub*instant_w_sub' + (1-Alpha).*Q_prev;
        
    end



    %if Q is ruined.
    if Q_pres_sub(1,1)<0
        disp('(1,1) is less than 0 : Err')
        if Q_pres_sub(2,2)<0
            disp('(2,2) is less than 0 : Err')
        end
    end

    %new method.
    if alpha == 102
        inv_R_revis = H'*pinv(R)*H;
    
        %disp(Q_pres_sub)
    
        sub_opt_pos = (Q_pres_sub + inv(inv_R_revis))*...
            (pinv(Q_pres_sub)*pred_pos...
            + inv_R_revis*TOA_pos);
    
        instant_w = sub_opt_pos - pred_pos - A*esti_epsilon_prev;
    
        if alpha == 100
            Alpha = 1/(1+exp(-a*(abs(instant_w)-E_w/2)));      %adaptible Q_varient
        end
        Q_pres = Alpha.*instant_w*instant_w' + (1-Alpha).*Q_prev;
    else
        Q_pres = Q_pres_sub;
        sub_opt_pos = NaN;
    end

%%
    
%{
    %Revision/ 24.02.16
    inv_R_revis = H'*pinv(R)*H;
    
    %disp(Q_pres_sub)

    sub_opt_pos = (Q_pres_sub + inv(inv_R_revis))*...
        (inv(Q_pres_sub)*pred_pos...
        + inv_R_revis*TOA_pos);

    instant_w = sub_opt_pos - pred_pos - A*esti_epsilon_prev;

    if alpha == 100
        Alpha = 1/(1+exp(-a*(abs(instant_w)-E_w_max/2)));      %adaptible Q_varient
    end
    Q_pres = Alpha.*instant_w*instant_w' + (1-Alpha).*Q_prev;


%}
    %%
    %pred error cov mat
    pred_cov = A*esti_cov_prev*A' + Q_pres;  %Q : modified Q(zero-mean)

    %%
    %calc Kalman Gain
    Kalman_Gain = pred_cov*H'*pinv(H*pred_cov*H'+R);

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