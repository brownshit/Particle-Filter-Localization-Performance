function init_cov = init_err_cov_calc(var,iter_Num)
    %
    %
    temp_P = zeros(2,2);

    %this method is more efficient.
    
    %setting Anchor
    ankx = 630; anky = 540;
    ANK_1 = [0,0];
    ANK_2 = [ankx,0];
    ANK_3 = [0,anky];
    for ind1 = 1:1:iter_Num
        aver_error = zeros(2,1);
        for ind2 = 1:1:iter_Num
            for ind=1:1:4
                n(ind,1) = sqrt(var)*randn; 
            end
            %Moving Points.
            exact_Pos = [2,2];
            
            Observation_1 = sqrt((ANK_1(1,1)-exact_Pos(1,1))^2+(ANK_1(1,2)-exact_Pos(1,2))^2);
            Observation_2 = sqrt((ANK_2(1,1)-exact_Pos(1,1))^2+(ANK_2(1,2)-exact_Pos(1,2))^2);
            Observation_3 = sqrt((ANK_3(1,1)-exact_Pos(1,1))^2+(ANK_3(1,2)-exact_Pos(1,2))^2);
            
            %truely get measurement considering noise
            Noisy_Observation_1 = Observation_1+n(1,1);
            Noisy_Observation_2 = Observation_2+n(2,1);
            Noisy_Observation_3 = Observation_3+n(3,1);
            
            measu_dist = ...
                [Noisy_Observation_1;...
                Noisy_Observation_2;...
                Noisy_Observation_3]; 

            %accumulate error
            aver_error(:,1) = aver_error(:,1) + TOA(measu_dist) - exact_Pos';
        end
        aver_error(:,1) = aver_error(:,1)./iter_Num;

        for ind=1:1:4
            n(ind,1) = sqrt(var)*randn; 
        end
        %Moving Points.
        exact_Pos = [2,2];
        
        Observation_1 = sqrt((ANK_1(1,1)-exact_Pos(1,1))^2+(ANK_1(1,2)-exact_Pos(1,2))^2);
        Observation_2 = sqrt((ANK_2(1,1)-exact_Pos(1,1))^2+(ANK_2(1,2)-exact_Pos(1,2))^2);
        Observation_3 = sqrt((ANK_3(1,1)-exact_Pos(1,1))^2+(ANK_3(1,2)-exact_Pos(1,2))^2);
        
        %truely get measurement considering noise
        Noisy_Observation_1 = Observation_1+n(1,1);
        Noisy_Observation_2 = Observation_2+n(2,1);
        Noisy_Observation_3 = Observation_3+n(3,1);
        
        measu_dist = ...
            [Noisy_Observation_1;...
            Noisy_Observation_2;...
            Noisy_Observation_3];

        pos = TOA(measu_dist);
        error = (pos-exact_Pos')-aver_error;
        temp_P = temp_P + error*error';
    end
    init_cov = temp_P./iter_Num;
end

%{
for ind1 = 1:1:iter_Num
        aver_error = zeros(2,1);
        for ind2 = 1:1:iter_Num
            %{
            
            for ind=1:1:4
                n(ind,1) = sqrt(var)*randn; 
            end
            %Moving Points.
            exact_Pos = [2,2];
            
            Observation_1 = sqrt((ANK_1(1,1)-exact_Pos(1,1))^2+(ANK_1(1,2)-exact_Pos(1,2))^2);
            Observation_2 = sqrt((ANK_2(1,1)-exact_Pos(1,1))^2+(ANK_2(1,2)-exact_Pos(1,2))^2);
            Observation_3 = sqrt((ANK_3(1,1)-exact_Pos(1,1))^2+(ANK_3(1,2)-exact_Pos(1,2))^2);
            
            %truely get measurement considering noise
            Noisy_Observation_1 = Observation_1+n(1,1);
            Noisy_Observation_2 = Observation_2+n(2,1);
            Noisy_Observation_3 = Observation_3+n(3,1);
            
            measu_dist = ...
                [Noisy_Observation_1;...
                Noisy_Observation_2;...
                Noisy_Observation_3];
            %} 
   
            [Measu_dist, Exact_Pos] = ArbitraryPoint3D(var);
            measu_dist = Measu_dist(3,:)';
            exact_Pos = Exact_Pos(3,:)';
            %accumulate error
            aver_error(:,1) = aver_error(:,1) + ...
                TOA(measu_dist) - ...
                exact_Pos;
        end
        aver_error(:,1) = aver_error(:,1)./iter_Num;

        [Measu_dist, Exact_Pos] = ArbitraryPoint3D(var);
        measu_dist = Measu_dist(3,:)';
        exact_Pos = Exact_Pos(3,:)';

        pos = TOA(measu_dist);
        error = (pos-exact_Pos)-aver_error;
        temp_P = temp_P + error*error';
    end
    init_cov = temp_P./iter_Num;
end
%}