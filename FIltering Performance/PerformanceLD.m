function PerformanceLD()    
    close all;
    %we'll use 'cm' scale to prevent floating point error
    format long e;
        
    %to generate noise by different 
    sigma_square = [0.25;1;4;16;64];
    alpha_1 = 0.2;
    alpha_2 = 1.0;        %   hyper parameter for Fading Q algorithm.
    alp_trigger = 100;

    %sigma should be 7.5[cm]

    %our given conditions of error should be lower than 30cm.

%%
    %basic setting for Kalman filter

    iter_Num = 1e4;     %iteration num.     to get mean value
    point_num = 100;
    dt = 0.01;

%   should be get      a 01 / 14

    %Get Q and E_w from txt files.
    Q_mean = zeros(2,2,5);
    fileID = fopen('Q_mean.txt', 'r');
    Q_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            for i22=1:1:2
               Q_mean(i2,i22,i1) = Q_raw((i1-1)*4+(i2-1)*2+i22,1);
            end
        end
    end
    
    E_w_mean = zeros(2,5);
    fileID = fopen('E_w_mean.txt', 'r');
    E_w_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            E_w_mean(i2,i1) = E_w_raw((i1-1)*2+i2,1);
        end
    end

%%
    Q_max = zeros(2,2,5);
    fileID = fopen('Q_max.txt', 'r');
    Q_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            for i22=1:1:2
               Q_max(i2,i22,i1) = Q_raw((i1-1)*4+(i2-1)*2+i22,1);
            end
        end
    end
    
    E_w_max = zeros(2,5);
    fileID = fopen('E_w_mean.txt', 'r');
    E_w_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            E_w_max(i2,i1) = E_w_raw((i1-1)*2+i2,1);
        end
    end

%%    
    ankx = 630; anky = 540;
    Anchor_1 = [0;0];
    Anchor_2 = [ankx;0];
    Anchor_3 = [0;anky];
    
    %measurement step should be implemented by TOA algorithm base.

    z = zeros(2,1);
    H = zeros(2,2);
    
    TOA_error = zeros(1,5);
    Kalman_F_error = zeros(1,5);
    
    TOA_esti = zeros(2,point_num,size(sigma_square,1));

    %%
    %%
    %%
    %%
    % Naive Estimation Storage
    esti_pos_LK = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_LK = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_LK = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_LK = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_LK = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_LK = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_LK = zeros(1,point_num,size(sigma_square,1));

    %%
    %%
    % Fading Q Estimation Storage
    esti_pos_FDQ = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_FDQ = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_FDQ = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_FDQ = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_FDQ = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_FDQ = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_FDQ = zeros(1,point_num,size(sigma_square,1));

    Fd_Q = zeros(2,2,point_num,size(sigma_square,1));
    %%
    % Fading Q Estimation Storage
    esti_pos_FDQ_2 = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_FDQ_2 = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_FDQ_2 = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_FDQ_2 = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_FDQ_2 = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_FDQ_2 = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_FDQ_2 = zeros(1,point_num,size(sigma_square,1));

    Fd_Q_2 = zeros(2,2,point_num,size(sigma_square,1));
    %%
    % Fading Q Estimation Storage       |   Adap
    esti_pos_FDQ_Adap = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_FDQ_Adap = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_FDQ_Adap = zeros(2,2,point_num,size(sigma_square,1));
    sub_opt_pos = zeros(2,2,point_num,size(sigma_square,1));

    Tr_esti_cov_FDQ_Adap = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_FDQ_Adap = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_FDQ_Adap = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_FDQ_Adap = zeros(1,point_num,size(sigma_square,1));

    Fd_Q_Adap = zeros(2,2,point_num,size(sigma_square,1));
    %%
    %%
    % (Naive)EKF with 2-iter Estimation Storage
    esti_pos_EKF_2iter = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_EKF_2iter = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_EKF_2iter = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_EKF_2iter = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_EKF_2iter = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_EKF_2iter = zeros(2,3,point_num,size(sigma_square,1));
    pos_KG_EKF_2iter = zeros(1,point_num,size(sigma_square,1));
%%
%%
%%
%%
    % Naive Estimation Storage
    esti_pos_LK_2 = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_LK_2 = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_LK_2 = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_LK_2 = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_LK_2 = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_LK_2 = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_LK_2 = zeros(1,point_num,size(sigma_square,1));
    
    %%
    % (Naive)EKF with 2-iter Estimation Storage
    esti_pos_EKF_2iter_2 = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_EKF_2iter_2 = zeros(2,2,point_num,size(sigma_square,1));
    pred_cov_EKF_2iter_2 = zeros(2,2,point_num,size(sigma_square,1));
    
    Tr_esti_cov_EKF_2iter_2 = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_EKF_2iter_2 = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_EKF_2iter_2 = zeros(2,3,point_num,size(sigma_square,1));
    pos_KG_EKF_2iter_2 = zeros(1,point_num,size(sigma_square,1));

    %%
    %%
    %%
    %%
    for var = 1:1:size(sigma_square,1)
        %varience step.
        disp(sigma_square(var,1))
        %%
        
        % should get error cov.matrix of step 3; [2,2]
        %   2023/11/26
        % for computational effort, I just decrease size of iter_num for cov
        % calculator.
        P3 = init_err_cov_calc(sigma_square(var,1),iter_Num/10);

        %for iter = 1:1:iter_Num
        %%
        %for Q

        % 01 / 14 : revision from here

        esti_vel_prev_LK = zeros(2,1);
        esti_vel_prev_LK_2 = zeros(2,1);
        esti_vel_prev_FDQ = zeros(2,1);
        esti_vel_prev_FDQ_2 = zeros(2,1);
        esti_vel_prev_FDQ_Adap = zeros(2,1);
        esti_vel_prev_EKF_2iter = zeros(2,1);
        esti_vel_prev_EKF_2iter_2 = zeros(2,1);
        %%
        %Moving Sequence

        % outputs are 100*3 matrices...
        [MD, TP] = ArbitraryPoint3D(sigma_square(var,1));

        measu_dist = MD';
        true_pos = TP';
        
        for mov_pnt = 1:1:point_num
            
            % measu_dist(mov_pnt,1)으로 대체해야 한다. 2024.01.12
            z(1,1) = (measu_dist(2,mov_pnt)^2-measu_dist(1,mov_pnt)^2)+...
                (Anchor_1(1,1)^2+Anchor_1(2,1)^2)-...
                (Anchor_2(1,1)^2+Anchor_2(2,1)^2);

            z(2,1) = (measu_dist(3,mov_pnt)^2-measu_dist(1,mov_pnt)^2)+...
                (Anchor_1(1,1)^2+Anchor_1(2,1)^2)-...
                (Anchor_3(1,1)^2+Anchor_3(2,1)^2);

            H(1,:) = 2*[(Anchor_1(1,1)-Anchor_2(1,1)),...
                (Anchor_1(2,1)-Anchor_2(2,1))];

            H(2,:) = 2*[(Anchor_1(1,1)-Anchor_3(1,1)),...
                (Anchor_1(2,1)-Anchor_3(2,1))];
                        
%{
            z(3,1) = (measu_dist(3,mov_pnt)^2-measu_dist(2,mov_pnt)^2)+...
                (Anchor_2(1,1)^2+Anchor_2(2,1)^2)-...
                (Anchor_3(1,1)^2+Anchor_3(2,1)^2);
        
            H(3,:) = 2*[(Anchor_2(1,1)-Anchor_3(1,1)),...
                (Anchor_2(2,1)-Anchor_3(2,1))];
%}
            

            TOA_esti(:,mov_pnt,var) = pinv(H)*z;      %Left pseudo inv
%%
%%            
            % Naive LKF /   Q_mean
            if mov_pnt >=4
                [esti_pos_LK(:,1,mov_pnt,var), esti_cov_LK(:,:,mov_pnt,var),...
                    Kalman_G_LK(:,:,mov_pnt,var),pred_cov_LK(:,:,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman(...
                    esti_pos_LK(:,1,mov_pnt-1,var), esti_cov_LK(:,:,mov_pnt-1,var),...
                    E_w_mean(:,var), Q_mean(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_LK(:,1));
            else
                % Before KF alg. starts
                esti_pos_LK(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_LK(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_LK(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_LK(1,mov_pnt,var) = trace(esti_cov_LK(:,:,mov_pnt,var));
            Tr_pred_cov_LK(1,mov_pnt,var) = trace(pred_cov_LK(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_LK(1,mov_pnt,var) = ...
                Kalman_G_LK(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_LK(:,1) = (esti_pos_LK(:,1,mov_pnt,var)-esti_pos_LK(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_LK(:,:,mov_pnt) = P3;     %should get by prev. simul.
                end
            end
%%            
            % Naive LKF /   Q_max
            if mov_pnt >=4
                [esti_pos_LK_2(:,1,mov_pnt,var), esti_cov_LK_2(:,:,mov_pnt,var),...
                    Kalman_G_LK_2(:,:,mov_pnt,var),pred_cov_LK_2(:,:,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman(...
                    esti_pos_LK_2(:,1,mov_pnt-1,var), esti_cov_LK_2(:,:,mov_pnt-1,var),...
                    E_w_max(:,var), Q_max(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_LK_2(:,1));
            else
                % Before KF alg. starts
                esti_pos_LK_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa_2 = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_LK_2(:,:,mov_pnt,var) = err_toa_2*err_toa_2';
                pred_cov_LK_2(:,:,mov_pnt,var) = err_toa_2*err_toa_2';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_LK_2(1,mov_pnt,var) = trace(esti_cov_LK_2(:,:,mov_pnt,var));
            Tr_pred_cov_LK_2(1,mov_pnt,var) = trace(pred_cov_LK_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_LK_2(1,mov_pnt,var) = ...
                Kalman_G_LK_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_LK_2(:,1) = (esti_pos_LK_2(:,1,mov_pnt,var)-esti_pos_LK_2(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_LK_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                end
            end
%%
%%          
            % LKF with time varient Q with diagonal | alpha : 0.2
            if mov_pnt >=4
                [esti_pos_FDQ(:,1,mov_pnt,var), esti_cov_FDQ(:,:,mov_pnt,var),...
                    Kalman_G_FDQ(:,:,mov_pnt,var),pred_cov_FDQ(:,:,mov_pnt,var),...
                    Fd_Q(:,:,mov_pnt,var),...
                    pred_pos_FDQ(:,1,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman_Q_fadingmemFil(...
                    esti_pos_FDQ(:,1,mov_pnt-1,var), esti_cov_FDQ(:,:,mov_pnt-1,var),...
                    Fd_Q(:,:,mov_pnt-1,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_FDQ(:,1),alpha_1,1,1);
            else
                % Before KF alg. starts
                esti_pos_FDQ(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_FDQ(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_FDQ(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_FDQ(1,mov_pnt,var) = trace(esti_cov_FDQ(:,:,mov_pnt,var));
            Tr_pred_cov_FDQ(1,mov_pnt,var) = trace(pred_cov_FDQ(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_FDQ(1,mov_pnt,var) = ...
                Kalman_G_FDQ(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_FDQ(:,1) = (esti_pos_FDQ(:,1,mov_pnt,var)-esti_pos_FDQ(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_FDQ(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    Fd_Q(:,:,mov_pnt,var) = Q_max(:,:,var); %init Q matrix (using Naive KF's Q)
                end
            end
%%          
            % LKF with time varient Q with diagonal | alpha : 0.9
            if mov_pnt >=4
                [esti_pos_FDQ_2(:,1,mov_pnt,var), esti_cov_FDQ_2(:,:,mov_pnt,var),...
                    Kalman_G_FDQ_2(:,:,mov_pnt,var),pred_cov_FDQ_2(:,:,mov_pnt,var),...
                    Fd_Q_2(:,:,mov_pnt,var),...
                    pred_pos_FDQ(:,1,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman_Q_fadingmemFil(...
                    esti_pos_FDQ_2(:,1,mov_pnt-1,var), esti_cov_FDQ_2(:,:,mov_pnt-1,var),...
                    Fd_Q_2(:,:,mov_pnt-1,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_FDQ_2(:,1),alpha_2,1,1);
            else
                % Before KF alg. starts
                esti_pos_FDQ_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_FDQ_2(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_FDQ_2(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_FDQ_2(1,mov_pnt,var) = trace(esti_cov_FDQ_2(:,:,mov_pnt,var));
            Tr_pred_cov_FDQ_2(1,mov_pnt,var) = trace(pred_cov_FDQ_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_FDQ_2(1,mov_pnt,var) = ...
                Kalman_G_FDQ_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_FDQ_2(:,1) = (esti_pos_FDQ_2(:,1,mov_pnt,var)-esti_pos_FDQ_2(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_FDQ_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    Fd_Q_2(:,:,mov_pnt,var) = Q_max(:,:,var); %init Q matrix (using Naive KF's Q)
                end
            end
            
%%          
            % LKF with time varient Q with diagonal | adaptible alpha 
                                                    % _with huristic algorithm
            if mov_pnt >=4
                %[temp_alp] determination======
                a = 2*log(9)*pinv(E_w_max(:,var));
                %==============================
                [esti_pos_FDQ_Adap(:,1,mov_pnt,var), esti_cov_FDQ_Adap(:,:,mov_pnt,var),...
                    Kalman_G_FDQ_Adap(:,:,mov_pnt,var),pred_cov_FDQ_Adap(:,:,mov_pnt,var),...
                    Fd_Q_Adap(:,:,mov_pnt,var),...
                    pred_pos_FDQ(:,1,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman_Q_fadingmemFil(...
                    esti_pos_FDQ_Adap(:,1,mov_pnt-1,var), esti_cov_FDQ_Adap(:,:,mov_pnt-1,var),...
                    Fd_Q_Adap(:,:,mov_pnt-1,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_FDQ_Adap(:,1),alp_trigger,...
                    a,E_w_max(:,var));
            else
                % Before KF alg. starts
                esti_pos_FDQ_Adap(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_FDQ_Adap(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_FDQ_Adap(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_FDQ_Adap(1,mov_pnt,var) = trace(esti_cov_FDQ_Adap(:,:,mov_pnt,var));
            Tr_pred_cov_FDQ_Adap(1,mov_pnt,var) = trace(pred_cov_FDQ_Adap(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_FDQ_Adap(1,mov_pnt,var) = ...
                Kalman_G_FDQ_Adap(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_FDQ_Adap(:,1) = (esti_pos_FDQ_Adap(:,1,mov_pnt,var)-esti_pos_FDQ_Adap(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_FDQ_Adap(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    Fd_Q_Adap(:,:,mov_pnt,var) = Q_max(:,:,var); %init Q matrix (using Naive KF's Q)
                end
            end
            
%%
%%
            % EKF 2-iterative.  /   Q_mean
            % TOA_Extended_Kalman

            % EKF Algorithm is ruined / 2024.02.05 / divegence issues.

            if mov_pnt >=4
                [esti_pos_EKF_2iter(:,1,mov_pnt,var), esti_cov_EKF_2iter(:,:,mov_pnt,var),...
                    Kalman_G_EKF_2iter(:,:,mov_pnt,var),pred_cov_EKF_2iter(:,:,mov_pnt,var)]...
                    =...
                    TOA_Extended_Kalman(...
                    esti_pos_EKF_2iter(:,1,mov_pnt-1,var), esti_cov_EKF_2iter(:,:,mov_pnt-1,var),...
                    E_w_mean(:,var), Q_mean(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_EKF_2iter(:,1));
            else
                % Before KF alg. starts
                esti_pos_EKF_2iter(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_EKF_2iter(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_EKF_2iter(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_EKF_2iter(1,mov_pnt,var) = trace(esti_cov_EKF_2iter(:,:,mov_pnt,var));
            Tr_pred_cov_EKF_2iter(1,mov_pnt,var) = trace(pred_cov_EKF_2iter(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_EKF_2iter(1,mov_pnt,var) = ...
                Kalman_G_EKF_2iter(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_EKF_2iter(:,1) = (esti_pos_EKF_2iter(:,1,mov_pnt,var)-esti_pos_EKF_2iter(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_EKF_2iter(:,:,mov_pnt) = P3;     %should get by prev. simul.
                end
            end

%%
            % EKF 2-iterative.  /   Q_max
            % TOA_Extended_Kalman

            if mov_pnt >=4
                [esti_pos_EKF_2iter_2(:,1,mov_pnt,var), esti_cov_EKF_2iter_2(:,:,mov_pnt,var),...
                    Kalman_G_EKF_2iter_2(:,:,mov_pnt,var),pred_cov_EKF_2iter_2(:,:,mov_pnt,var)]...
                    =...
                    TOA_Extended_Kalman(...
                    esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var), esti_cov_EKF_2iter_2(:,:,mov_pnt-1,var),...
                    E_w_max(:,var), Q_max(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt), esti_vel_prev_EKF_2iter_2(:,1));
            else
                % Before KF alg. starts
                esti_pos_EKF_2iter_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                esti_cov_EKF_2iter_2(:,:,mov_pnt,var) = err_toa*err_toa';
                pred_cov_EKF_2iter_2(:,:,mov_pnt,var) = err_toa*err_toa';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_EKF_2iter_2(1,mov_pnt,var) = trace(esti_cov_EKF_2iter_2(:,:,mov_pnt,var));
            Tr_pred_cov_EKF_2iter_2(1,mov_pnt,var) = trace(pred_cov_EKF_2iter_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_EKF_2iter_2(1,mov_pnt,var) = ...
                Kalman_G_EKF_2iter_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_EKF_2iter_2(:,1) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
                if mov_pnt == 3
                    esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                end
            end


%%
            % Error calc
            % TOA alg error calc.
            TOA_error(1,var) = TOA_error(1,var) +...
                sqrt((TOA_esti(1,mov_pnt,var)-true_pos(1,mov_pnt))^2+...
                (TOA_esti(2,mov_pnt,var)-true_pos(2,mov_pnt))^2);    %error
            
            %{
            %% Kalamn FIlter

            if mov_pnt>=4
                % Kalman Filtering
                [esti_pos(:,:,mov_pnt),esti_cov(:,:,mov_pnt),R] = ...
                    TOA_Kalman(esti_pos(:,:,mov_pnt-1), esti_cov(:,:,mov_pnt-1),...
                        E_w(:,var),Q(:,:,var),sigma_square(var,1), measu_dist,esti_vel_prev);
                %calc vel of prev.
                esti_vel_prev = (esti_pos(:,:,mov_pnt)-esti_pos(:,:,mov_pnt-1))./dt;
                

                %tracking position.
                %esti_pos_global(:,1,mov_pnt,iter) = esti_pos(:,1,mov_pnt);
    
                %cov_register_bystep
                cov_register_bystep(:,:,mov_pnt) = cov_register_bystep(:,:,mov_pnt) + esti_cov(:,:,mov_pnt);
                trace_cov_register_bystep(1,mov_pnt) = trace_cov_register_bystep(1,mov_pnt) + trace(esti_cov(:,:,mov_pnt));%esti_cov(1,1,mov_pnt)*esti_cov(2,2,mov_pnt)
                
                %add error by esti_pos.
                error = esti_pos(:,:,mov_pnt)-exact_Pos';
                Kalman_F_error(1,var) = Kalman_F_error(1,var)+sqrt(error(1,1)^2 + error(2,1)^2);
    
            elseif mov_pnt ~= 1
                if mov_pnt == 3
                    esti_cov(:,:,3) = P2;

                    cov_register_bystep(:,:,mov_pnt) = cov_register_bystep(:,:,mov_pnt) + esti_cov(:,:,mov_pnt);
                    trace_cov_register_bystep(1,mov_pnt) = trace_cov_register_bystep(1,mov_pnt) + trace(esti_cov(:,:,mov_pnt));%esti_cov(1,1,mov_pnt)*esti_cov(2,2,mov_pnt)
    
                end
                %set by TOA values for point [1,1],[2,2]
                esti_pos(:,1,mov_pnt) = TOA(measu_dist);
                %calc vel of prev.
                esti_vel_prev = (esti_pos(:,:,mov_pnt)-esti_pos(:,:,mov_pnt-1))./dt;
    
                %for diag Q
                esti_pos_diag(:,1,mov_pnt) = TOA(measu_dist);
                %calc vel of prev.
                esti_vel_prev_diag = (esti_pos_diag(:,:,mov_pnt)-esti_pos_diag(:,:,mov_pnt-1))./dt;
            end
            %}


        end
        %end

        %%
        %TOA_error(1,var) =... 
        % TOA_error(1,var)./(iter_Num*point_num); %10 for MovPnt.
            %%
        %trace_cov_register_bystep(1,:) = trace_cov_register_bystep(1,:)./iter_Num;
        %cov_register_bystep(:,:,:) = cov_register_bystep(:,:,:)./iter_Num;
    
        
        %calc just for Kalman Filter point 3,3 to 9,9   
        %disp("[cov_register_bystep]");
        %disp(cov_register_bystep);
    
        %{
        det_R = zeros(1,point_num);
        for ind =1:1:point_num
            det_R(1,ind) = abs(det(R_average(:,:,ind)));
        end
        figure;
        plot(plt_point, det_R,"-*");
        title("step R")
        grid on
        %}
    
        %plot each cov by point
        %{
        figure;
        plot(plt_point, trace_cov_register_bystep,"--o");
        title_name = "Trace of State covarience matirx [ Var : "+sigma_square(var,1)+" ]";
        title(title_name);
        xlabel("[Point]");
        ylabel("Trace(C_x)");
        legend("Q");
        grid on
        %}
        

        TOA_point = TOA_esti(:,:,var)'; %3 x 100 -> 100 x 3'
        for ind = 1:1:100
            Naive_KAL_FIL_pos(:,ind) = esti_pos_LK(:,1,ind,var);
            Naive_KAL_FIL_pos_2(:,ind) = esti_pos_LK_2(:,1,ind,var);

            FDQ_KAL_FIL_pos(:,ind) = esti_pos_FDQ(:,1,ind,var);
            FDQ_KAL_FIL_pos_2(:,ind) = esti_pos_FDQ_2(:,1,ind,var);
            FDQ_KAL_FIL_pos_Adap(:,ind) = esti_pos_FDQ_Adap(:,1,ind,var);
            %Gotten by maximum likelihood estimation
            %sub_OPT_pos(:,ind) = sub_opt_pos(:,1,ind,var);
            pred_posi(:,ind) = pred_pos_FDQ(:,1,ind,var);

            EKF_2iter_KAL_FIL_pos(:,ind) = esti_pos_EKF_2iter(:,1,ind,var);
            EKF_2iter_KAL_FIL_pos_2(:,ind) = esti_pos_EKF_2iter_2(:,1,ind,var);
        end


        %Figure #1 : postion traking comparisions by each algorithms
        figure;
        plot(true_pos(1,:), true_pos(2,:), 'co-');
        hold on
        plot(TOA_point(:,1), TOA_point(:,2), 'ko-');

        plot(Naive_KAL_FIL_pos(1,:),Naive_KAL_FIL_pos(2,:), '^-');
        plot(Naive_KAL_FIL_pos_2(1,:),Naive_KAL_FIL_pos_2(2,:), '^-');

        plot(FDQ_KAL_FIL_pos(1,:),FDQ_KAL_FIL_pos(2,:), '^-');
        plot(FDQ_KAL_FIL_pos_2(1,:),FDQ_KAL_FIL_pos_2(2,:), '^-');
        plot(FDQ_KAL_FIL_pos_Adap(1,:),FDQ_KAL_FIL_pos_Adap(2,:), '^-');
        %plot(sub_OPT_pos(1,:),sub_OPT_pos(2,:), '^-');

        
        plot(pred_posi(1,:),pred_posi(2,:), '^-');
        
        plot(EKF_2iter_KAL_FIL_pos(1,:),EKF_2iter_KAL_FIL_pos(2,:), '^-');
        plot(EKF_2iter_KAL_FIL_pos_2(1,:),EKF_2iter_KAL_FIL_pos_2(2,:), 'm^-');

        hold off
        legend('True Position',...
            'TOA measurement',...
            'Kalman Filter_{Naive, mean Q}',...
            'Kalman Filter_{Naive, max Q}',...
            'Kalman Filter_{time varient Q, \alpha = 0.2}',...
            'Kalman Filter_{time varient Q, \alpha = 1.0}',...
            'Kalman Filter_{time varient Q, Adaptable \alpha}',...
            ...%'MLEs based Position_{sub}',...
            'pred position(for comparision)',...
            'Extended Kalman Filter_{with 2-iter, Naive, mean Q}',...
            'Extended Kalman Filter_{with 2-iter, Naive, max Q}')
        xlabel('X');
        ylabel('Y');
        titles = "[Variences:"+sigma_square(var,1)+"] Points";
        title(titles);
        grid on;
        axis equal;
%{

        ax = linspace(1,point_num,point_num);
        %Figure #2 : covarience traking
        figure;
        plot(ax,Tr_pred_cov_LK(1,:,var), '^--');
        hold on
        plot(ax,Tr_esti_cov_LK(1,:,var), 'o--');
        plot(ax,Tr_pred_cov_FDQ(1,:,var), '^--');
        plot(ax,Tr_esti_cov_FDQ(1,:,var), 'o--');
        plot(ax,Tr_pred_cov_EKF_2iter(1,:,var), '^--');
        
        plot(ax,Tr_esti_cov_EKF_2iter(1,:,var), 'o--');
        hold off
        legend('Pred Cov.(Q mean)','Esti Cov.(Q mean)','Pred Cov.(Q varient)','Esti Cov.(Q varient)','Pred Cov.(EKF, Q mean)','Esti Cov.(EKF, Q mean)')
        xlabel('k-steps[time step]');
        ylabel('Tr(P_k)');
        titles = "[Variences:"+sigma_square(var,1)+"] Tr[Error Cov.] traking";
        title(titles);
        grid on;

        %Figure #3 : Kalman Gain traking
        figure;
        plot(ax,pos_KG_LK(1,:,var), 'o--');
        hold on
        plot(ax,pos_KG_FDQ(1,:,var), 'o--');
        
        plot(ax,pos_KG_EKF_2iter(1,:,var), 'o--');
        hold off
        legend('Gain (KF)','Gain (Q varient KF)','Gain (EKF)')
        xlabel('k-steps[time step]');
        ylabel('Kalman Gain_c_o_r_r_e_s_p_o_n_d_i_n_g _P_O_S');
        titles = "[Variences:"+sigma_square(var,1)+"] Kalman Gain_c_o_r_r_e_s_p_o_n_d_i_n_g _P_O_S traking";
        title(titles);
        grid on;

%}
    end



%{

    Kalman_F_error = Kalman_F_error./(iter_Num*7);   
    
    disp("TOA error:")
    disp(TOA_error')
    disp("TOA based K.F error:")
    disp(Kalman_F_error')

    figure;
    plot(sigma_square,TOA_error,":o");
    hold on
    %for KF Q
    plot(sigma_square, Kalman_F_error,"--o");
    hold off

    xlabel("Varience");
    ylabel("Error Rate");
    legend("TOA algorithm",...
        "TOA based Kalman Filter _w_i_t_h _Q",...
        "location","northwest")
    title("TOA location determination");
    grid on
%}
    
    
end