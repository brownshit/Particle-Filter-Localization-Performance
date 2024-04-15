function [x_esti,v_esti,...
    pt_pf_pres,pt_kf_pred_pres, w_pres,...
    kf_cov_pred_pres,...
    Q_LPF_return] = RBPF...
    (pt_pf_prev, pt_kf_pred_prev,w_prev,...
    kf_cov_pred_prev,...
    noise_var,Npt,dt,z,...
    Q_mean, Q_max, Q_LPF_prev, x_esti_prev,v_esti_prev,...
    a,E_w)
    %RBPF for each step.
    %{
        - set importance density as prior
            prior needs information of pt_kf_pred_prev 
            (which is pt of t-1|t-2)

        - 현재의 LPF Q를 다음 스텝에 넘겨줘야한다.

        - 기존의 구한 Q를 여기 안으로 들여올 때, 처리를 해줘야한다. 
            (기존Q)*eye(2)./dt를 function 내의 Q로 넣어줘야 한다.

        - Note that) kf state is velocity.

        input : 
            pt_pf_prev       : t-1's        particle of pf
            pt_kf_pred_prev  : t-1|t-2's    particle of kf      
                이걸 받아오는게 맞나? 왜냐면 resampling과정도 안거치긴 했는데. 이건 확인을 좀 해봐야할 것
                - 딱히 상관은 없을 듯. one step kalman Filter처럼 동작하기 때문.

            measu_var        : sigma square
            Npt              : particle number
            dt               : time stride
            Q_mean, Q_max    : matrices to determine present LPF Q
            Q_LPF_prev       : prev Q defined with prev LPF step.
            z                : measurement of distance target to each Anchors
            
        output : 
            x_esti           : pf_esti
            v_esti           : kf_esti
            pt_pf_pres       : pf_pt
            pt_kf_pred_pres  : kf_pt_present        
            w_pres           : (resampled)weight
            kf_cov_pred_pres : cov for next step.
            Q_LPF_return     : Q_LPF for next teim step
    
    %}

    pt_pf_pres = zeros(2,Npt);
    pt_kf_pred_pres = zeros(2,Npt);
    kf_cov_pred_pres = zeros(2,2);  %kf_pred of t|t-1
    w_pres = zeros(1,Npt);

    std = sqrt(noise_var);
    %%
    % linear transform matrices
    A_pf = eye(2).*dt;
    A_kf = eye(2);
    %B_f_pf = eye(2).*(dt^2/2);
    B_f_pf = eye(2);
    %B_f_kf = eye(2).*dt;
    B_f_kf = eye(2)./dt;
    %%
    % Q_mean and Q_max are B_f_kf*f_t-1 
    % || process outside of this function
    
    %% pred : draw sample from prior density. not using fx
    for ind_prev = 1:Npt
        mu = pt_pf_prev(:,ind_prev) + A_pf*pt_kf_pred_prev(:,ind_prev);
        sig_sq = A_pf*kf_cov_pred_prev*(A_pf') + B_f_pf*Q_LPF_prev*(B_f_pf');
        %disp(mu)
        %disp(sig_sq)
        try
            R = chol(sig_sq);
        catch ME
            disp(Q_LPF_prev)
            disp("                      ")
        end


        pt_pf_pres(:,ind_prev) = mu+chol(sig_sq)*randn(2,1);        %particle generation.
    end
    %% calc pred of /pt_kf_pres(mean)/ and /kf_cov_pred_pres(cov)/.
    A_bar_kf = A_kf-B_f_kf*pinv(B_f_pf)*A_pf;
    
    Kalman_Gain_prev = ...
        kf_cov_pred_prev*(A_pf')*...
        inv(A_pf*kf_cov_pred_prev*(A_pf')+B_f_pf*Q_LPF_prev*(B_f_pf'));

    for ind_prev = 1:Npt
        z_t_1(:,ind_prev) = ...
            pt_pf_pres(:,ind_prev)-pt_pf_prev(:,ind_prev);

        pt_kf_pred_pres(:,ind_prev) = ...
            A_bar_kf*(pt_kf_pred_prev(:,ind_prev)+...
            Kalman_Gain_prev*(z_t_1(:,ind_prev)-A_pf*pt_kf_pred_prev(:,ind_prev)))+...
            B_f_kf*pinv(B_f_pf)*z_t_1(:,ind_prev);
    end
    kf_cov_pred_pres = ...
        A_bar_kf*(kf_cov_pred_prev - Kalman_Gain_prev*A_pf*kf_cov_pred_prev)*(A_bar_kf');
    
    %figure;
    %plot(pt_pf_pres(1,:),pt_pf_pres(2,:),'o')
    %{
    nn = sum(pt_pf_pres(1,:))
    temp_sum = 0;
    for ki = 1:Npt
        temp_sum = temp_sum + (pt_pf_pres(1,ki)-nn)^2;
    end
    disp(temp_sum)
    %}

    %% weight
    % Calc. weight from likelihood of obervation, w_k = w_k-1*p(z_k|pt_k)
    
    figure;
    hold on
    
    for ind_prev = 1:Npt
        %여기가 문제.
        %여기가 문제.
        %여기가 문제.
        %여기가 문제.
        %여기가 문제.
        err = z-hx(pt_pf_pres(:,ind_prev));
        toa = TOA(z);
        %{
            
        %}
        plot(toa(1,1),toa(2,1),'ro');
        plot(pt_pf_pres(1,ind_prev),pt_pf_pres(2,ind_prev),'ko')
        
        theta = 0:0.01:2*pi; % 각도 범위
        xx = toa(1,1) + std*cos(theta); % 원의 x 좌표
        yy = toa(2,1) + std*sin(theta); % 원의 y 좌표

        plot(xx, yy,'r'); % 원 plot
        %disp(err)
        likelihood = exp((-1/(2*noise_var))*(err'*err))...
            ./((sqrt(2*pi))^size(z,1)*sqrt(noise_var^size(z,1)));
        %disp(likelihood)        %모두 0나온다.
        w_pres(ind_prev) = w_prev(ind_prev)*likelihood;    %weight revision
        %normpdf could revise into specific noise, when our noise are not
        %Gaussian.
    end
    
    hold off
    axis equal
    
    %weight normalization
    disp(sum(w_pres))       %애가 0나온다.
    w_pres = w_pres/sum(w_pres);
    disp("///////////////////////////////")
    %disp(w_pres')
    %disp(pt_pf_pres')
    disp("///////////////////////////////")
    %% estimation
    x_esti = pt_pf_pres*w_pres';        %pos_estimation by pf_particles
    v_esti = pt_kf_pred_pres*w_pres';   %vel_estimation by kf_particles
    
    %%
    % Q_LPF for step t.
    %   previous works : instant_w_sub = TOA_pos - pred_pos;
    %instant_w_sub = x_esti - (x_esti_prev+v_esti_prev.*dt);
    instant_w_sub = TOA(z) - (x_esti_prev+v_esti_prev.*dt);

    Alpha = 1/(1+exp(-a*(abs(instant_w_sub)-E_w/2)));      %adaptible Q_varient
    Q_LPF_return = Alpha.*Q_max + (1-Alpha).*Q_mean;
    %disp(a)
    disp(x_esti)
    disp(x_esti_prev)
    disp(instant_w_sub)     %이게 NaN
    disp(exp(-a*(abs(instant_w_sub)-E_w/2)))
    disp((1+exp(-a*(abs(instant_w_sub)-E_w/2))))
    disp(Alpha)
    disp("=================================")
    %%
    % Resampling step could have operation condition... N_eff < N_th = 2N/3
    var_accum = 0;
    for ind_prev = 1:Npt
        var_accum = var_accum + w_pres(ind_prev)^2; 
    end
    if 1/var_accum < 2*Npt/3
        %============================== SIR procedure. // Systematic resampling
        wtc = cumsum(w_pres);
        rpt = rand(Npt,1);
        [~,ind1] = sort([rpt; wtc']);
        ind_pres = find(ind1<= Npt) - [0:Npt-1]';
        pt_pf_pres = pt_pf_pres(:,ind_pres);
        w_pres = ones(1,Npt)*(1/Npt);
        %==============================
    end
end