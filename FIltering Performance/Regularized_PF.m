function [x_esti,pt_pres,w_pres] = Regularized_PF(pt_prev,w_prev,z,dt,measu_var,abs_process_n,E_w,Npt,prev_pt_vel)
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
    S_k = zeros(2,2);       %empirical cov. mat
    D_k = zeros(2,2);       %empirical cov. mat

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
    for ind = 1:Npt 
        var_accum = var_accum + w_pres(ind)^2; 
    end
    if 1/var_accum < 2*Npt/3
        %============================== SIR procedure.
        for ind1 = 1:Npt
            error = pt_pres(:,k)-x_esti;
            S_k = S_k + w_pres(1,k)*(error*error'); %Emprical cov.
        end
        %D_k = chol(S_k)';       %Choelesky decomp.

%%
        n_X = 2;
        c_nx = pi;%unit vol. of hypersphere in this case; pi
        h_opt = nthroot((Npt*192),6);       % ~~7.6
        
        kernel_input_norm = zeros(1,Npt);
        dist_accum = zeros(1,Npt);

        for ind = 1:Npt
            %Calc Epanechnikov kernel.
            kernel_input_norm(1,ind) = norm(x_esti - pt_pres(:,ind),2)/h_opt;
            K_h = K_opt(c_nx,n_X,kernel_input_norm(1,ind))/(h_opt^n_X);
            %disp(K_h)
            dist_accum(1,ind) = w_pres(1,ind)*K_h;
        end
        if sum(dist_accum) == 0
            disp("[CAUTION : prob's sum is equals to 0]")
        end

        [kernel_input_norm, index] = sort(kernel_input_norm);
        sorted_dist_accum = zeros(1,Npt);

        for ind = 1:Npt
            sorted_dist_accum(1,ind) = dist_accum(1,index(1,ind));
        end
%%
        %============================== // Systematic resampling
        wtc = cumsum(w_pres);
        rpt = rand(Npt,1);
        [~,ind1] = sort([rpt; wtc']);
        ind2 = find(ind1<= Npt) - [0:Npt-1]';
        pt_pres = pt_pres(:,ind2);
        w_pres = ones(1,Npt)*(1/Npt);
        %==============================

%{
        fileID = fopen('dist_accum.txt', 'w');
        fprintf(fileID, '%.7f\n', dist_accum);
        fclose(fileID);

%}
        %fitting dist with values.
        %disp(sum(sorted_dist_accum));
        
        %sorted_dist_accum = sorted_dist_accum./sum(sorted_dist_accum);
        for ind = 1:Npt
            %disp(dist_accum(1,ind))

            epsilon = randsample(kernel_input_norm,2,true,sorted_dist_accum)';
            %for state dim is 2;
            
            %{
            epsilon(1,1) = epsilon(1,1)*sqrt(S_k(1,1));
            epsilon(2,1) = epsilon(2,1)*sqrt(S_k(2,2));
            %}
            epsilon = D_k*epsilon;

            %disp(norm(epsilon,2))
            pt_pres(:,ind) = pt_pres(:,ind) + (epsilon).*h_opt;
        end
        %==============================
    end
    %disp(ind)
end