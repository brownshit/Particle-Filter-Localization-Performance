function PerformanceLDPF()
    close all;
    %we'll use 'cm' scale to prevent floating point error
    format long e;
        
    pt_num = 1e3;
    %to generate noise by different 
    sigma_square = [0.25;1;4;16;64];
    alp_trigger_2 = 101;

    %sigma should be 7.5[cm]

    %our given conditions of error should be lower than 30cm.

%%
    %basic setting for Kalman filter

    iter_Num = 1e4;     %iteration num.     to get mean value
    point_num = 100;
    dt = 0.01;
    %%
    measu_dist = zeros(3,point_num, size(sigma_square,1));
    true_pos = zeros(2,point_num, size(sigma_square,1));

    Naive_KAL_FIL_pos_2 = zeros(2,point_num, size(sigma_square,1));
    FDQ_KAL_FIL_pos_Adap_2 = zeros(2,point_num, size(sigma_square,1));
    EKF_2iter_KAL_FIL_pos_2 = zeros(2,point_num, size(sigma_square,1));
    PF_pos = zeros(2,point_num, size(sigma_square,1));
    ASIR_PF_pos = zeros(2,point_num, size(sigma_square,1));
    RPF_pos = zeros(2,point_num, size(sigma_square,1));
    RBPF_pos = zeros(2,point_num, size(sigma_square,1));
    %%
    
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

    abs_E_w_mean = zeros(2,5);
    fileID = fopen('abs_E_w_mean.txt', 'r');
    abs_E_w_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            abs_E_w_mean(i2,i1) = abs_E_w_raw((i1-1)*2+i2,1);
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

    abs_E_w_max = zeros(2,5);
    fileID = fopen('abs_E_w_max.txt', 'r');
    abs_E_w_raw = fscanf(fileID, '%f');
    fclose(fileID);
    for i1=1:1:size(sigma_square,1)
        for i2=1:1:2
            abs_E_w_max(i2,i1) = abs_E_w_raw((i1-1)*2+i2,1);
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
    % Fading Q Estimation Storage       |   Adap mode 101
    esti_pos_FDQ_Adap_2 = zeros(2,1,point_num,size(sigma_square,1));
    esti_cov_FDQ_Adap_2 = zeros(2,2,point_num,size(sigma_square,1));
    pred_pos_FDQ_Adap_2 = zeros(2,1,point_num,size(sigma_square,1));
    pred_cov_FDQ_Adap_2 = zeros(2,2,point_num,size(sigma_square,1));
    sub_opt_pos_2 = zeros(2,1,point_num,size(sigma_square,1));

    Tr_esti_cov_FDQ_Adap_2 = zeros(1,point_num,size(sigma_square,1));
    Tr_pred_cov_FDQ_Adap_2 = zeros(1,point_num,size(sigma_square,1));

    Kalman_G_FDQ_Adap_2 = zeros(2,2,point_num,size(sigma_square,1));
    pos_KG_FDQ_Adap_2 = zeros(1,point_num,size(sigma_square,1));

    Fd_Q_Adap_2 = zeros(2,2,point_num,size(sigma_square,1));
    Alpha_2 = zeros(1,point_num,size(sigma_square,1));

%%
    %particle filter storage
    pf_x_esti = zeros(2,point_num,size(sigma_square,1));

%%
    %ASIR PF filter storage
    pf_x_esti_ASIR = zeros(2,point_num,size(sigma_square,1));

%%
    %ASIR PF filter storage
    pf_x_esti_RPF = zeros(2,point_num,size(sigma_square,1));

%%
    %RBPF LPFQ PF filter storage
    pf_x_esti_RBPF = zeros(2,point_num,size(sigma_square,1));
    pf_v_esti_RBPF = zeros(2,point_num,size(sigma_square,1));

    Q_max_RBPF = zeros(2,2,5);
    Q_max_RBPF = zeros(2,2,5);
%%
    for var = 1:1:size(sigma_square,1)
        %varience step.
        fprintf('[Noise Varience] : %3d\n',sigma_square(var,1))
        
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

        esti_vel_prev_LK_2 = zeros(2,1);
        esti_vel_prev_FDQ_Adap_2 = zeros(2,1);
        esti_vel_prev_EKF_2iter_2 = zeros(2,1);
  %%      
        % Each pt vel.
        prev_pt_vel = zeros(2,pt_num);
        abs_process_n = zeros(2,1);

        pt = zeros(2,pt_num);
        w = zeros(1,pt_num);
%%
        prev_pt_vel_ASIR = zeros(2,pt_num);
        abs_process_n_ASIR = zeros(2,1);

        pt_ASIR = zeros(2,pt_num);
        w_ASIR = zeros(1,pt_num);
%%
        prev_pt_vel_RPF = zeros(2,pt_num);
        abs_process_n_RPF = zeros(2,1);

        pt_RPF = zeros(2,pt_num);
        w_RPF = zeros(1,pt_num);
%%      RBPF
        B_f_pf_inv = eye(2).*(2/dt^2);

        pt_pf_RBPF = zeros(2,pt_num);
        pt_kf_pred_RBPF = zeros(2,pt_num);
        w_RBPF = zeros(1,pt_num);

        kf_cov_pred = zeros(2,2);

        Q_max_RBPF(:,:,var) = Q_max(:,:,var)*B_f_pf_inv;
        Q_mean_RBPF(:,:,var) = Q_mean(:,:,var)*B_f_pf_inv;
        Q_LPF_prev = zeros(2,2);
                %%
        %Moving Sequence

        % outputs are 100*3 matrices...
        [MD, TP] = ArbitraryPoint3D(sigma_square(var,1));

        measu_dist(:,:,var) = MD';
        true_pos(:,:,var) = TP';


%%
        fprintf('=> [Progress] : %3d%%\n',0)
        
        for mov_pnt = 1:1:point_num
            
            % measu_dist(mov_pnt,1)으로 대체해야 한다. 2024.01.12
            z(1,1) = (measu_dist(2,mov_pnt,var)^2-measu_dist(1,mov_pnt,var)^2)+...
                (Anchor_1(1,1)^2+Anchor_1(2,1)^2)-...
                (Anchor_2(1,1)^2+Anchor_2(2,1)^2);

            z(2,1) = (measu_dist(3,mov_pnt,var)^2-measu_dist(1,mov_pnt,var)^2)+...
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
            

            TOA_esti(:,mov_pnt,var) = inv(H'*H)*H'*z;      %Left pseudo inv
 %%            
            % Naive LKF /   Q_max
            if mov_pnt >=4
                [esti_pos_LK_2(:,1,mov_pnt,var), esti_cov_LK_2(:,:,mov_pnt,var),...
                    Kalman_G_LK_2(:,:,mov_pnt,var),pred_cov_LK_2(:,:,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman(...
                    esti_pos_LK_2(:,1,mov_pnt-1,var), esti_cov_LK_2(:,:,mov_pnt-1,var),...
                    E_w_max(:,var), Q_max(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt,var), esti_vel_prev_LK_2(:,1));
            else
                % Before KF alg. starts
                esti_pos_LK_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);

                %cov
                %err_toa_2 = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                if mov_pnt == 3
                    esti_cov_LK_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                else
                    %esti_cov_LK_2(:,:,mov_pnt,var) = err_toa_2*err_toa_2';
                end
                %pred_cov_LK_2(:,:,mov_pnt,var) = err_toa_2*err_toa_2';
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_LK_2(1,mov_pnt,var) = trace(esti_cov_LK_2(:,:,mov_pnt,var));
            Tr_pred_cov_LK_2(1,mov_pnt,var) = trace(pred_cov_LK_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_LK_2(1,mov_pnt,var) = ...
                Kalman_G_LK_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_LK_2(:,1) = (esti_pos_LK_2(:,1,mov_pnt,var)-esti_pos_LK_2(:,1,mov_pnt-1,var))/dt;
            end

%%
            % LKF with time varient Q with diagonal | adaptible alpha 
                                                    % _with huristic algorithm
                                                    % MODE 101
            if mov_pnt >=4
                %[temp_alp] determination======
                a = 2*log(9)*pinv(sigma_square(var,1)*abs_E_w_mean(:,var));
                %==============================
                [esti_pos_FDQ_Adap_2(:,1,mov_pnt,var), esti_cov_FDQ_Adap_2(:,:,mov_pnt,var),...
                    Kalman_G_FDQ_Adap_2(:,:,mov_pnt,var),pred_cov_FDQ_Adap_2(:,:,mov_pnt,var),...
                    Fd_Q_Adap_2(:,:,mov_pnt,var),...
                    pred_pos_FDQ_Adap_2(:,1,mov_pnt,var),...
                    sub_opt_pos_2(:,1,mov_pnt,var),...
                    Alpha_2(1,mov_pnt,var)]...
                    =...
                    TOA_Linear_Kalman_Q_fadingmemFil(...
                    esti_pos_FDQ_Adap_2(:,1,mov_pnt-1,var), esti_cov_FDQ_Adap_2(:,:,mov_pnt-1,var),...
                    Fd_Q_Adap_2(:,:,mov_pnt-1,var), sigma_square(var,1), measu_dist(:,mov_pnt,var), esti_vel_prev_FDQ_Adap_2(:,1),alp_trigger_2,...
                    a,sigma_square(var,1)*abs_E_w_mean(:,var),...
                    Q_max(:,:,var), Q_mean(:,:,var));
            else
                % Before KF alg. starts
                esti_pos_FDQ_Adap_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                %err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                %esti_cov_FDQ_Adap(:,:,mov_pnt,var) = err_toa*err_toa';
                %pred_cov_FDQ_Adap(:,:,mov_pnt,var) = err_toa*err_toa';
                if mov_pnt == 3
                    esti_cov_FDQ_Adap_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    Fd_Q_Adap_2(:,:,mov_pnt,var) = Q_mean(:,:,var); %init Q matrix (using Naive KF's Q)
                end
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_FDQ_Adap_2(1,mov_pnt,var) = trace(esti_cov_FDQ_Adap_2(:,:,mov_pnt,var));
            Tr_pred_cov_FDQ_Adap_2(1,mov_pnt,var) = trace(pred_cov_FDQ_Adap_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_FDQ_Adap_2(1,mov_pnt,var) = ...
                Kalman_G_FDQ_Adap_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_FDQ_Adap_2(:,1) = (esti_pos_FDQ_Adap_2(:,1,mov_pnt,var)-esti_pos_FDQ_Adap_2(:,1,mov_pnt-1,var))/dt;
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
                    E_w_max(:,var), Q_max(:,:,var), sigma_square(var,1), measu_dist(:,mov_pnt,var), esti_vel_prev_EKF_2iter_2(:,1));
            else
                % Before KF alg. starts
                esti_pos_EKF_2iter_2(:,1,mov_pnt,var) = TOA_esti(:,mov_pnt,var);
                %cov
                %err_toa = TOA_esti(:,mov_pnt,var)-true_pos(:,mov_pnt);
                %esti_cov_EKF_2iter_2(:,:,mov_pnt,var) = err_toa*err_toa';
                %pred_cov_EKF_2iter_2(:,:,mov_pnt,var) = err_toa*err_toa';
                if mov_pnt == 3
                    esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                end
            end
            
            % To visualization, get trace of error cov matrix.
            Tr_esti_cov_EKF_2iter_2(1,mov_pnt,var) = trace(esti_cov_EKF_2iter_2(:,:,mov_pnt,var));
            Tr_pred_cov_EKF_2iter_2(1,mov_pnt,var) = trace(pred_cov_EKF_2iter_2(:,:,mov_pnt,var));
            
            %       Additionally, get Frobenious norm of Kalman Gain.
            pos_KG_EKF_2iter_2(1,mov_pnt,var) = ...
                Kalman_G_EKF_2iter_2(1,1,mov_pnt,var);

            if mov_pnt ~= 1
                esti_vel_prev_EKF_2iter_2(:,1) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
            end
%%
%%
%%
            % Implementation of PF
            % pt_num
            %{
                - should determine;
                    abs_process_n : from max_Q
                    init pos
            %}
            
            if mov_pnt >=4
                abs_process_n(1,1) = sqrt(Q_max(1,1,var));
                abs_process_n(2,1) = sqrt(Q_max(2,2,var));
                pt_prev = pt;

                [pf_x_esti(:,mov_pnt,var),pt,w] = ...
                    generic_PF(pt,w,...
                    measu_dist(:,mov_pnt,var),dt,...
                    sigma_square(var,1),abs_process_n,E_w_max(:,var),...
                    pt_num,prev_pt_vel);

                % (pt_pos - pt_pos)/dt
                prev_pt_vel = (pt-pt_prev)./dt;

            else %step for {1,2,3}
                % Before PF alg. starts
                pf_x_esti(:,mov_pnt,var) = TOA_esti(:,mov_pnt,var);      %position estimation
                if mov_pnt == 3
                    %init particle generation, and set weights either.
                    pt(1,:) = pf_x_esti(1,mov_pnt,var) + sqrt(P3(1,1))*randn(1,pt_num);
                    pt(2,:) = pf_x_esti(2,mov_pnt,var) + sqrt(P3(2,2))*randn(1,pt_num);

                    w = ones(1,pt_num)*1/pt_num;

                    %special case of calc vel.
                    % (pt_pos - pos by TOA)/dt
                    prev_pt_vel = (pt-pf_x_esti(:,2,var))./dt;    %step 3th postion.
                end
            end
            %{
                %eliminate under codes.. processing in upper codes.
                if mov_pnt ~= 1
                    %vel_update    | size : 
                    prev_pt_vel(:,:) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
                    if mov_pnt == 3
                        esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    end
                end
            %}
%%
            % Implementation of ASIR PF
            %{
            prev_pt_vel_ASIR = zeros(2,pt_num);
            abs_process_n_ASIR = zeros(2,1);
    
            pt_ASIR = zeros(2,pt_num);
            w_ASIR = zeros(1,pt_num);

                %       여기부터. 변수명 바꾸면 된다.
                        현재 ASIR까지 코드만 짬

            %}
            if mov_pnt >=4
                abs_process_n_ASIR(1,1) = sqrt(Q_max(1,1,var));
                abs_process_n_ASIR(2,1) = sqrt(Q_max(2,2,var));
                pt_prev_ASIR = pt_ASIR;

                [pf_x_esti_ASIR(:,mov_pnt,var),pt_ASIR,w_ASIR] = ...
                    ASIR_PF(pt_ASIR,w_ASIR,...
                    measu_dist(:,mov_pnt,var),dt,...
                    sigma_square(var,1),abs_process_n_ASIR,E_w_max(:,var),...
                    pt_num,prev_pt_vel_ASIR);

                % (pt_pos - pt_pos)/dt
                prev_pt_vel_ASIR = (pt_ASIR-pt_prev_ASIR)./dt;

            else %step for {1,2,3}
                % Before PF alg. starts
                pf_x_esti_ASIR(:,mov_pnt,var) = TOA_esti(:,mov_pnt,var);      %position estimation
                if mov_pnt == 3
                    %init particle generation, and set weights either.
                    pt_ASIR(1,:) = pf_x_esti_ASIR(1,mov_pnt,var) + sqrt(P3(1,1))*randn(1,pt_num);
                    pt_ASIR(2,:) = pf_x_esti_ASIR(2,mov_pnt,var) + sqrt(P3(2,2))*randn(1,pt_num);

                    w_ASIR = ones(1,pt_num)*1/pt_num;

                    %special case of calc vel.
                    % (pt_pos - pos by TOA)/dt
                    prev_pt_vel_ASIR = (pt_ASIR-pf_x_esti_ASIR(:,2,var))./dt;    %step 3th postion.
                end
            end
            %{
                %eliminate under codes.. processing in upper codes.
                if mov_pnt ~= 1
                    %vel_update    | size : 
                    prev_pt_vel(:,:) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
                    if mov_pnt == 3
                        esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    end
                end
            %}

%%
            % Implementation of RPF
            %{
            prev_pt_vel_RPF = zeros(2,pt_num);
            abs_process_n_RPF = zeros(2,1);
    
            pt_RPF = zeros(2,pt_num);
            w_RPF = zeros(1,pt_num);

                %       여기부터. 변수명 바꾸면 된다.

            %}
            if mov_pnt >=4
                abs_process_n_RPF(1,1) = sqrt(Q_max(1,1,var));
                abs_process_n_RPF(2,1) = sqrt(Q_max(2,2,var));
                pt_prev_RPF = pt_RPF;
                
                %Regularized_PF

                [pf_x_esti_RPF(:,mov_pnt,var),pt_RPF,w_RPF] = ...
                    Regularized_PF(pt_RPF,w_RPF,...
                    measu_dist(:,mov_pnt,var),dt,...
                    sigma_square(var,1),abs_process_n_RPF,E_w_max(:,var),...
                    pt_num,prev_pt_vel_RPF);

                % (pt_pos - pt_pos)/dt
                prev_pt_vel_RPF = (pt_RPF-pt_prev_RPF)./dt;

            else %step for {1,2,3}
                % Before PF alg. starts
                pf_x_esti_RPF(:,mov_pnt,var) = TOA_esti(:,mov_pnt,var);      %position estimation
                if mov_pnt == 3
                    %init particle generation, and set weights either.
                    pt_RPF(1,:) = pf_x_esti_RPF(1,mov_pnt,var) + sqrt(P3(1,1))*randn(1,pt_num);
                    pt_RPF(2,:) = pf_x_esti_RPF(2,mov_pnt,var) + sqrt(P3(2,2))*randn(1,pt_num);

                    w_RPF = ones(1,pt_num)*1/pt_num;

                    %special case of calc vel.
                    % (pt_pos - pos by TOA)/dt
                    prev_pt_vel_RPF = (pt_RPF-pf_x_esti_RPF(:,2,var))./dt;    %step 3th postion.
                end
            end
            %{
                %eliminate under codes.. processing in upper codes.
                if mov_pnt ~= 1
                    %vel_update    | size : 
                    prev_pt_vel(:,:) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
                    if mov_pnt == 3
                        esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    end
                end
            %}

            %%
            % Implementation of RBPF
            % pt_num
            %{
                - should determine;
                    abs_process_n : from max_Q
                    init pos

            [x_esti,v_esti,...
                pt_pf_pres,pt_kf_pred_pres, w_pres,...
                kf_cov_pred_pres,...
                Q_LPF_return] = RBPF...
                (pt_pf_prev, pt_kf_pred_prev,w_prev,...
                kf_cov_pred_prev,...
                noise_var,Npt,dt,z,...
                Q_mean, Q_max,Q_LPF_prev, x_esti_prev,...
                a,E_w)

            also, process Q_mean and Q_max for function input.
                
            %}
            %{
            
            if mov_pnt >=4
                a = 2*log(9)*pinv(sigma_square(var,1)*abs_E_w_mean(:,var));

                [pf_x_esti_RBPF(:,mov_pnt,var),pf_v_esti_RBPF(:,mov_pnt,var),...
                pt_pf_RBPF,pt_kf_pred_RBPF, w_RBPF,...
                kf_cov_pred,...
                Q_LPF_prev] = ...
                    RBPF(pt_pf_RBPF, pt_kf_pred_RBPF,w_RBPF,...
                    kf_cov_pred,...
                    sigma_square(var,1),pt_num,dt,measu_dist(:,mov_pnt,var),...
                    Q_mean(:,:,var), Q_max(:,:,var),Q_LPF_prev,...
                    pf_x_esti_RBPF(:,mov_pnt-1,var),pf_v_esti_RBPF(:,mov_pnt-1,var),...
                    a,sigma_square(var,1)*abs_E_w_mean(:,var));
            else %step for {1,2,3}
                % Before PF alg. starts
                pf_x_esti_RBPF(:,mov_pnt,var) = TOA_esti(:,mov_pnt,var);      %position estimation
                if mov_pnt == 3
                    %init particle generation, and set weights either.
                    

                    % particle 초기 생성자체에 문제가 있는 것 같다.
                    % 실제 상황에 적절한 분산을 가질 수 있게 해야한다.
                    pt_pf_RBPF(1,:) = pf_x_esti_RBPF(1,mov_pnt,var) + sqrt(P3(1,1))*randn(1,pt_num);
                    pt_pf_RBPF(2,:) = pf_x_esti_RBPF(2,mov_pnt,var) + sqrt(P3(2,2))*randn(1,pt_num);
                    
                    pt_kf_pred_RBPF(1,:) = ...
                        (pt_pf_RBPF(1,:)-pf_x_esti_RBPF(1,mov_pnt-1,var))./dt;
                    pt_kf_pred_RBPF(2,:) = ...
                        (pt_pf_RBPF(2,:)-pf_x_esti_RBPF(2,mov_pnt-1,var))./dt;

                    %empirical cov. of kf pt.
                    mean_kf_pt = mean(pt_kf_pred_RBPF,2);
                    accum_cov = zeros(2,2);
                    for ind = 1:pt_num
                        err = pt_kf_pred_RBPF(:,ind) - mean_kf_pt;
                        accum_cov = accum_cov + (err*err');
                    end
                    kf_cov_pred = accum_cov./pt_num;
    
                    %kf_cov_pred = pt_init_vel_err_cov_calc(sigma_square(var,1),iter_Num,pt_num);

                    w_RBPF = ones(1,pt_num)*1/pt_num;

                    Q_LPF_prev = Q_mean(:,:,var);      %init process noise's err. cov
                    disp(Q_LPF_prev)
                    %disp(Q_max_RBPF(:,:,var))
                    %disp(Q_mean_RBPF(:,:,var))
                end
            end
            %}

            
            %{
                %eliminate under codes.. processing in upper codes.
                if mov_pnt ~= 1
                    %vel_update    | size : 
                    prev_pt_vel(:,:) = (esti_pos_EKF_2iter_2(:,1,mov_pnt,var)-esti_pos_EKF_2iter_2(:,1,mov_pnt-1,var))/dt;
                    if mov_pnt == 3
                        esti_cov_EKF_2iter_2(:,:,mov_pnt) = P3;     %should get by prev. simul.
                    end
                end
            %}
%%
%%
%%

            TOA_point = TOA_esti(:,:,var)'; %3 x 100 -> 100 x 3'
            for ind = 1:1:100
                Naive_KAL_FIL_pos_2(:,ind,var) = esti_pos_LK_2(:,1,ind,var);
                FDQ_KAL_FIL_pos_Adap_2(:,ind,var) = esti_pos_FDQ_Adap_2(:,1,ind,var);
                EKF_2iter_KAL_FIL_pos_2(:,ind,var) = esti_pos_EKF_2iter_2(:,1,ind,var);
                PF_pos(:,ind,var) = pf_x_esti(:,ind,var);
                ASIR_PF_pos(:,ind,var) = pf_x_esti_ASIR(:,ind,var);
                RPF_pos(:,ind,var) = pf_x_esti_RPF(:,ind,var);
                %RBPF_pos(:,ind,var) = pf_x_esti_RBPF(:,ind,var);
            end
%%    
        fprintf('\b\b\b\b\b')
        fprintf('%3d%%\n',round(mov_pnt/point_num*100))
            
        end
        fprintf('\n')

        %Figure #1 : postion traking comparisions by each algorithms
        figure;
        plot(true_pos(1,:,var), true_pos(2,:,var), 'co-');
        hold on
        plot(TOA_esti(1,:,var), TOA_esti(2,:,var), 'ko-');
        plot(Naive_KAL_FIL_pos_2(1,:,var),Naive_KAL_FIL_pos_2(2,:,var), '^-');
        plot(FDQ_KAL_FIL_pos_Adap_2(1,:,var),FDQ_KAL_FIL_pos_Adap_2(2,:,var), '^-');
        plot(EKF_2iter_KAL_FIL_pos_2(1,:,var),EKF_2iter_KAL_FIL_pos_2(2,:,var), 'm^-');
        plot(PF_pos(1,:,var),PF_pos(2,:,var), 'r^-');
        plot(ASIR_PF_pos(1,:,var),ASIR_PF_pos(2,:,var), 'b^-');
        plot(RPF_pos(1,:,var),RPF_pos(2,:,var), '^-');
        hold off
        legend('True Position',...
            'TOA measurement',...
            'Kalman Filter_{Naive, max Q}',...
            'Kalman Filter_{time varient Q, Adaptable \alpha}',...
            'Extended Kalman Filter_{with 2-iter, Naive, max Q}',...
            'Generic Particle Filter_{max Q}',...
            'ASIR Particle Filter_{max Q}',...
            'Regularized Particle Filter_{max Q}')
        xlabel('X');
        ylabel('Y');
        titles = "Trajectory Tracking [Variences:"+sigma_square(var,1)+"]";
        title(titles);
        grid on;
        axis equal;


        figure;
        ax = linspace(1,100,100);
        plot(ax, Alpha_2(1,:,var),'ko--');
        xlabel('point num')
        ylabel('Alpha')
        legend('\alpha LKF')
        titles = "\alpha for Q_{LPF} [Variences:"+sigma_square(var,1)+"]";
        title(titles)
        grid on

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
end