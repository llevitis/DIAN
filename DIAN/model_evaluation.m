function [model_solutions,model_parameters,model_times,model_CORRs,model_RMSEs] = model_evaluation(Nnodes,Nsubjects,time_range,ref_pattern,h_integration,Conn_Matrix,Long_Matrix,num_repetitions,Nparameters,BETAS0,DELTAS0,SIGMAS0,MUS0,seed_regions_1,SPEEDS)
% Yasser Iturria Medina,
% 20/06/2013, Montreal.

% Defining outputs for our model
model_solutions  = zeros(size(ref_pattern,1),Nsubjects);
model_parameters = zeros(Nparameters,Nsubjects);
model_times      = zeros(Nsubjects,1);
model_CORRs      = -10*ones(Nsubjects,1);
model_RMSEs      = 100*ones(Nsubjects,1); % root mean squared error

NodesS0    = zeros(Nnodes,1);
t0         = time_range(1);
time_final = time_range(end);
mu = 0;
I = ceil(max(time_final)/h_integration); % number of sample points (Runge-Kutta)
for repetition = 1:num_repetitions
    if exist('seed_regions_1') && ~isempty(seed_regions_1),
        disp(['Repetition ' num2str(repetition) ', seed regions 1 -> ' num2str(seed_regions_1)])
    end
    if exist('num_repetitions') && repetition > 1,
        NodesS0(seed_regions_1) = rand(length(seed_regions_1),1);
    else
        NodesS0(seed_regions_1) = 0.1;
    end
    
    % Computing future probability of Amyloidal deposition
    for mu = MUS0
        disp(['Noise (mu): ' num2str(mu)])
        for sigma = SIGMAS0
            disp(['Noise (sigma): ' num2str(sigma)])
            tic
            % profile on
            if mu ~= 0 || sigma ~= 0
                % Noise = normrnd(mu,sigma,[Nnodes I 4]);
                Noise = single(normrnd(mu,sigma,[Nnodes I 4]));
                add_noise = 1;
            else
                add_noise = 0;
            end
            for beta = BETAS0
                % profile on
                for delta = DELTAS0
                    disp(['Computing progression, beta: ' num2str(beta) ', delta: ' num2str(delta)])
                    if add_noise,
                        [prob_tf,temporal_progression,times] = RungeKutta4_version6(h_integration,NodesS0,Conn_Matrix,t0,time_final-t0,beta,delta,Noise);
                        % [prob_tf,temporal_progression,times] = RungeKutta4_version6_C(h_integration,NodesS0,Conn_Matrix,t0,time_final-t0,beta,delta,Noise);
                    else
                        [prob_tf,temporal_progression,times] = RungeKutta4_version6(h_integration,NodesS0,Conn_Matrix,t0,time_final-t0,beta,delta);
                        % [prob_tf,temporal_progression,times] = RungeKutta4_version6_C(h_integration,NodesS0,Conn_Matrix,t0,time_final-t0,beta,delta);
                    end
                    temp_corr = 1./slmetric_pw(ref_pattern,temporal_progression,'eucdist');
                    
                    % Selecting the best predictions
                    [temp_corr,temp_times] = max(temp_corr,[],2);
                    
                    ind = find(temp_corr > model_CORRs);
                    if ~isempty(ind),
                        model_CORRs(ind) = temp_corr(ind);
                        model_solutions(:,ind) = temporal_progression(:,temp_times(ind));
                        model_parameters(1,ind) = beta;
                        model_parameters(2,ind) = delta;
                        model_parameters(3,ind) = sigma;
                        if Nparameters == 4, model_parameters(4,ind) = mu; end
                        model_times(ind) = temp_times(ind)*h_integration;
                        model_RMSEs(ind) = (mean((ref_pattern(:,ind) - temporal_progression(:,temp_times(ind))).^2,1)).^0.5;
                    end
                end
                % profile viewer
            end
            % profile viewer
            toc
        end
    end
end