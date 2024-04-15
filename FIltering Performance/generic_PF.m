function [x_esti,pt_pres,w_pres] = generic_PF(pt_prev,w_prev,z,dt,measu_var,abs_process_n,E_w,Npt,prev_pt_vel)
    %{
        Output      : pt_pres,w_pres
        - size      : 2xNeff            | p_x; p_y
                    : x_esti : estimation value

        Output      : pt_prev,w_prev
                    : z             : measurement
                    : dt            : time step
                    : measu_var     : R
                    : abs_process_n : process noise

        For this PF example, we assume that Q, R are Gaussian noise.
    %}
    
    pt_pres = zeros(2,Npt);
    w_pres = zeros(1,Npt);

    std = sqrt(measu_var);

    %% 
    % Draw sample from importance density, especailly, we assume id as prior
    %       Prediction
    for k = 1:Npt
        pt_pres(:,k) = fx(pt_prev(:,k),dt,prev_pt_vel(:,k),E_w)+...
            [abs_process_n(1,1)*randn;...
            abs_process_n(2,1)*randn];
        %Term abs_process_n.*randn(4,1); could revise into specific noise.
    end

    %%
    % Calc. weight from likelihood of obervation, w_k = w_k-1*p(z_k|pt_k)
    for k = 1:Npt
        err = z-hx(pt_pres(:,k));
        likelihood = exp((-1/(2*measu_var))*(err'*err))./((sqrt(2*pi))^size(z,1)*std);
        w_pres(k) = w_prev(k)*likelihood;    %weight revision
        %normpdf could revise into specific noise, when our noise are not
        %Gaussian.
    end
    %weight normalization
    w_pres = w_pres/sum(w_pres);

    %%
    %estimation
    x_esti = pt_pres*w_pres';

    %%
    % Resampling step could have operation condition... N_eff < N_th = 2N/3
    var_accum = 0;
    for ind = 1:Npt var_accum = var_accum + w_pres(ind)^2; end
    if 1/var_accum < 2*Npt/3
        %============================== SIR procedure. // Systematic resampling
        wtc = cumsum(w_pres);
        rpt = rand(Npt,1);
        [~,ind1] = sort([rpt; wtc']);
        ind = find(ind1<= Npt) - [0:Npt-1]';
        pt_pres = pt_pres(:,ind);
        w_pres = ones(1,Npt)*(1/Npt);
        %==============================
    end
    %disp(ind)
end