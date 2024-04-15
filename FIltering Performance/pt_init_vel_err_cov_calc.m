function init_cov = pt_init_vel_err_cov_calc(var,iter_Num,Npt)
    %
    %
    temp_P = zeros(2,2);

    %setting Anchor
    dt = 0.01;
    ankx = 630; anky = 540;
    ANK_1 = [0,0];
    ANK_2 = [ankx,0];
    ANK_3 = [0,anky];

    exact_vel = [1;1];

    n_pt = zeros(2,Npt);

    for ind2 = 1:1:iter_Num
        %%
        for ind=1:1:4
            n(ind,1) = sqrt(var)*randn; 
        end
        %Moving Points.
        exact_Pos_prev = [2,2];
        
        Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_prev(1,1))^2+(ANK_1(1,2)-exact_Pos_prev(1,2))^2);
        Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_prev(1,1))^2+(ANK_2(1,2)-exact_Pos_prev(1,2))^2);
        Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_prev(1,1))^2+(ANK_3(1,2)-exact_Pos_prev(1,2))^2);
        
        %truely get measurement considering noise
        Noisy_Observation_1 = Observation_1+n(1,1);
        Noisy_Observation_2 = Observation_2+n(2,1);
        Noisy_Observation_3 = Observation_3+n(3,1);
        
        measu_dist = ...
            [Noisy_Observation_1;...
            Noisy_Observation_2;...
            Noisy_Observation_3]; 
        
        pos_prev = TOA(measu_dist);      %TOA of [2,2]
        %%
        for ind=1:1:4
            n(ind,1) = sqrt(var)*randn; 
        end
        %Moving Points.
        exact_Pos_pres = [3,3];
        
        Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_pres(1,1))^2+(ANK_1(1,2)-exact_Pos_pres(1,2))^2);
        Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_pres(1,1))^2+(ANK_2(1,2)-exact_Pos_pres(1,2))^2);
        Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_pres(1,1))^2+(ANK_3(1,2)-exact_Pos_pres(1,2))^2);
        
        %truely get measurement considering noise
        Noisy_Observation_1 = Observation_1+n(1,1);
        Noisy_Observation_2 = Observation_2+n(2,1);
        Noisy_Observation_3 = Observation_3+n(3,1);
        
        measu_dist = ...
            [Noisy_Observation_1;...
            Noisy_Observation_2;...
            Noisy_Observation_3]; 
        
        pos_pres = TOA(measu_dist);      %TOA of [3,3]
        %%
        for ind=1:1:2
            n_pt(ind,:) = sqrt(var)*randn(1,Npt); 
        end
        %Moving Points.
        pseudo_pt_pres = pos_pres+n_pt;
        pseudo_pt_vel_err = (pseudo_pt_pres-pos_prev)./dt - exact_vel;
        
        temp_P = temp_P + (pseudo_pt_vel_err*pseudo_pt_vel_err')./Npt;
    end
    temp_P = temp_P./iter_Num;
    init_cov = temp_P;
%{

    for ind1 = 1:1:iter_Num
        aver_error = zeros(2,1);
        for ind2 = 1:1:iter_Num
            %%
            for ind=1:1:4
                n(ind,1) = sqrt(var)*randn; 
            end
            %Moving Points.
            exact_Pos_prev = [2,2];
            
            Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_prev(1,1))^2+(ANK_1(1,2)-exact_Pos_prev(1,2))^2);
            Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_prev(1,1))^2+(ANK_2(1,2)-exact_Pos_prev(1,2))^2);
            Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_prev(1,1))^2+(ANK_3(1,2)-exact_Pos_prev(1,2))^2);
            
            %truely get measurement considering noise
            Noisy_Observation_1 = Observation_1+n(1,1);
            Noisy_Observation_2 = Observation_2+n(2,1);
            Noisy_Observation_3 = Observation_3+n(3,1);
            
            measu_dist = ...
                [Noisy_Observation_1;...
                Noisy_Observation_2;...
                Noisy_Observation_3]; 
            
            pos_prev = TOA(measu_dist);      %TOA of [2,2]
            %%
            for ind=1:1:4
                n(ind,1) = sqrt(var)*randn; 
            end
            %Moving Points.
            exact_Pos_pres = [3,3];
            
            Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_pres(1,1))^2+(ANK_1(1,2)-exact_Pos_pres(1,2))^2);
            Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_pres(1,1))^2+(ANK_2(1,2)-exact_Pos_pres(1,2))^2);
            Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_pres(1,1))^2+(ANK_3(1,2)-exact_Pos_pres(1,2))^2);
            
            %truely get measurement considering noise
            Noisy_Observation_1 = Observation_1+n(1,1);
            Noisy_Observation_2 = Observation_2+n(2,1);
            Noisy_Observation_3 = Observation_3+n(3,1);
            
            measu_dist = ...
                [Noisy_Observation_1;...
                Noisy_Observation_2;...
                Noisy_Observation_3]; 
            
            pos_pres = TOA(measu_dist);      %TOA of [3,3]
            %%
            for ind=1:1:2
                n_pt(ind,:) = sqrt(var)*randn(1,Npt); 
            end
            %Moving Points.
            pseudo_pt_pres = pos_pres+n_pt;
            pseudo_pt_vel_err = (pseudo_pt_pres-pos_prev)./dt - exact_vel;

            aver_error = aver_error + pseudo_pt_vel_err;
        end
        aver_error = aver_error./iter_Num;

        %%
        for ind=1:1:4
            n(ind,1) = sqrt(var)*randn; 
        end
        %Moving Points.
        exact_Pos_prev = [2,2];
        
        Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_prev(1,1))^2+(ANK_1(1,2)-exact_Pos_prev(1,2))^2);
        Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_prev(1,1))^2+(ANK_2(1,2)-exact_Pos_prev(1,2))^2);
        Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_prev(1,1))^2+(ANK_3(1,2)-exact_Pos_prev(1,2))^2);
        
        %truely get measurement considering noise
        Noisy_Observation_1 = Observation_1+n(1,1);
        Noisy_Observation_2 = Observation_2+n(2,1);
        Noisy_Observation_3 = Observation_3+n(3,1);
        
        measu_dist = ...
            [Noisy_Observation_1;...
            Noisy_Observation_2;...
            Noisy_Observation_3]; 
        
        pos_prev = TOA(measu_dist);      %TOA of [2,2]
        %%
        for ind=1:1:4
            n(ind,1) = sqrt(var)*randn; 
        end
        %Moving Points.
        exact_Pos_pres = [3,3];
        
        Observation_1 = sqrt((ANK_1(1,1)-exact_Pos_pres(1,1))^2+(ANK_1(1,2)-exact_Pos_pres(1,2))^2);
        Observation_2 = sqrt((ANK_2(1,1)-exact_Pos_pres(1,1))^2+(ANK_2(1,2)-exact_Pos_pres(1,2))^2);
        Observation_3 = sqrt((ANK_3(1,1)-exact_Pos_pres(1,1))^2+(ANK_3(1,2)-exact_Pos_pres(1,2))^2);
        
        %truely get measurement considering noise
        Noisy_Observation_1 = Observation_1+n(1,1);
        Noisy_Observation_2 = Observation_2+n(2,1);
        Noisy_Observation_3 = Observation_3+n(3,1);
        
        measu_dist = ...
            [Noisy_Observation_1;...
            Noisy_Observation_2;...
            Noisy_Observation_3]; 
        
        pos_pres = TOA(measu_dist);      %TOA of [3,3]
        %%
        for ind=1:1:2
            n_pt(ind,:) = sqrt(var)*randn(1,Npt); 
        end
        %Moving Points.
        pseudo_pt_pres = pos_pres+n_pt;
        pseudo_pt_vel_err = (pseudo_pt_pres-pos_prev)./dt - exact_vel;

        error = pseudo_pt_vel_err-aver_error;
        temp_P = temp_P + error*error';
    end
    init_cov = temp_P./iter_Num;
end
%}