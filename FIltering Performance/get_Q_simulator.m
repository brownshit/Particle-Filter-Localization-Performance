function [Q_mat ,Expect_of_process_noise,ABS_Expect_of_process_noise]...
    = get_Q_simulator(var, E_w,mod)
    format long e;
    %   Q simulation
    %   sigma_square = [1e-3;1e-2;1e-1;1;1e1];
    %
    %   Expect_of_process_noise is 2x1 vec; 
    %   (1,1) : E(wx) / (2,1) : E(wy)
    %var = 1e-2;
    iteration = 1e4;
    end_pnt = 100;

    dt = 0.01;
    A = [1,0;
        0,1];
    Q_mat = zeros(2,2);

    ankx = 630; anky = 540;
    Anchor_1 = [0;0];
    Anchor_2 = [ankx;0];
    Anchor_3 = [0;anky];

    %=================================
    PROCESS_NOISE = zeros(2,1);
    process_noise = zeros(2,end_pnt); %[wx;wy];

    process_noise_accumulator = [0;0];
    process_noise_accumulator_ABS = [0;0];
    process_noise_memory = zeros(2,iteration);
    %=================================

    for iter = 1:1:iteration
        [MD,TP] = ArbitraryPoint3D(var);

        meas_dist = MD';
        true_pos = TP'; %$size : 2,100
        
        TOA_memory = zeros(2,end_pnt);
        %error_memory = zeros(2,end_pnt);
        
        approxi_vel = zeros(2,end_pnt);

        for MovPnt = 2:1:end_pnt      %skip first step if we exactly know about init value.
            %{
            for ind=1:1:4
                %n(ind,1) = sqrt(sigma_square(var,1))*randn; 
                n(ind,1) = sqrt(var)*randn; 
            end
            %Moving Points.
            %exact_Pos = [d1*rand, d2*rand];
            exact_Pos = [MovPnt-1,MovPnt-1];
            true_pos(:,MovPnt) = exact_Pos';
            
            Observation_1 = sqrt((Anchor_1_Pos(1,1)-exact_Pos(1,1))^2+(Anchor_1_Pos(1,2)-exact_Pos(1,2))^2);
            Observation_2 = sqrt((Anchor_2_Pos(1,1)-exact_Pos(1,1))^2+(Anchor_2_Pos(1,2)-exact_Pos(1,2))^2);
            Observation_3 = sqrt((Anchor_3_Pos(1,1)-exact_Pos(1,1))^2+(Anchor_3_Pos(1,2)-exact_Pos(1,2))^2);
            Observation_4 = sqrt((Anchor_4_Pos(1,1)-exact_Pos(1,1))^2+(Anchor_4_Pos(1,2)-exact_Pos(1,2))^2);
            
            %truely get measurement considering noise
            Noisy_Observation_1 = Observation_1+n(1,1);
            Noisy_Observation_2 = Observation_2+n(2,1);
            Noisy_Observation_3 = Observation_3+n(3,1);
            Noisy_Observation_4 = Observation_4+n(4,1);
            
            %location detection algorithm -> by z = Hx
            %basic settings.
            
            z = [Noisy_Observation_1^2-Noisy_Observation_2^2-100;
                Noisy_Observation_1^2-Noisy_Observation_3^2;
                Noisy_Observation_1^2-Noisy_Observation_4^2+100;
                Noisy_Observation_2^2-Noisy_Observation_3^2+100;
                Noisy_Observation_2^2-Noisy_Observation_4^2+200;
                Noisy_Observation_3^2-Noisy_Observation_4^2+100];
            
            H = [0,-20;
                20,-20;
                20,0;
                20,0;
                20,20;
                0,20];
            
            measurement_pos_TOA = (inv(H'*H))*H'*z;      %Left pseudo inv
            %}
            measurement_pos_TOA = TOA(meas_dist(:,MovPnt));
            TOA_memory(:,MovPnt) = measurement_pos_TOA;

            %approxi vel.
            approxi_vel(:,MovPnt) = (TOA_memory(:,MovPnt)-TOA_memory(:,MovPnt-1))/dt;
            
            % could skip first step.
            %if we don't have confidence at our init position, should use...
            %LPF_est_temp = pos_sol_1;       %new temp
    
            %{
            error_storage_sol_1(1,var) = error_storage_sol_1(1,var) +...
                sqrt((pos_sol_1(1,1)-exact_Pos(1,1))^2 + (pos_sol_1(2,1)-exact_Pos(1,2))^2);    %error
    
            %}
            if MovPnt >= 3      
                %이때부터 process noise를 구할 수 있다.  
                % process noise의 movpnt 1,2번째 요소는 비어있다.
                process_noise(:,MovPnt) = ...
                    true_pos(:,MovPnt)...
                    -(A*TOA_memory(:,MovPnt-1)+(approxi_vel(:,MovPnt-1))*dt)...
                    -E_w(:,1);
            end
        end
        
        % I took the max value of each process noise in one procedure.
        % which could make our filter more robust.
        

        %option #1
        if mod == 1
            PROCESS_NOISE(:,1) = mean(process_noise,2);
        end

        %option #2
        if mod == 2
            PROCESS_NOISE(:,1) = max(process_noise,[],2);
        end

        %after calc 0~2, we can get process noise -> to compose Q matrix.
        %for [2;2] point, we could get process noise calc step.
        process_noise_accumulator = process_noise_accumulator + PROCESS_NOISE;
        process_noise_accumulator_ABS = process_noise_accumulator_ABS + abs(PROCESS_NOISE);

        %for figure process noise
        process_noise_memory(:,iter) = PROCESS_NOISE;

        %Q_mat accumulation.
        Q_mat = Q_mat + (PROCESS_NOISE)*(PROCESS_NOISE)';
    end
    %calc expectation of process noise.
    Expect_of_process_noise = process_noise_accumulator./iteration;
    ABS_Expect_of_process_noise = process_noise_accumulator_ABS./iteration;
    Q_mat = Q_mat./iteration;

    %=================================
    %{
    figure;
    bivariate_hist(process_noise_memory(1,:),process_noise_memory(2,:), 0.01);
    xlabel('[Wx]');
    ylabel('[Wy]');
    title('[process noise]');

    %Q_matrix Expectation
    
    disp("[E(wx)] : "+Expect_of_process_noise(1,1)+" [E(wy)] : "+Expect_of_process_noise(2,1));
    disp("we got Q matrix");
    disp(Q_mat);


    %}
end