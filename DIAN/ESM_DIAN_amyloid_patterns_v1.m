addpath(genpath('PairWiseDistanceMetrics/')); 
addpath('results/'); 
addpath('testing_tmp'); 

% Epidemic Spreading Model (ESM; see Iturria-Medina et al., 2014, Plos Comp.
% Biol.), adapted to tau propagation/deposition. Code started by Yasser
% Iturria-Medina and Jacob Vogel, Montreal, Nov./Dec., 2017.
% Add also 'PairWiseDistanceMetrics' toolbox <------------------------- !!!

% Experiment: Predict amyloid beta deposition values for carriers in DIAN dataset 
esm_input_file = '~/projects/def-aevans/llevitis/DIAN/DIAN/testing_tmp/esm_dian_2comp_brainstem_ages_noncarriers.mat'; 
[filepath,name,ext] = fileparts(esm_input_file);
roidata = load(esm_input_file); 
output_dir = ['~/projects/def-aevans/llevitis/DIAN/DIAN/results']; % output directory <-------------------------------- !!!
%output_dir = ['/Users/liza/data/GreyMatterESM/results/']; 
%% Loading data  ---------------------------------------------------------%
N_tau_possibilities = 1; % number of different Tau measurements <------ !!!
% loading tau presence Probabilities, using a global reference
%roidata = load('/home/users/llevitis/code/DIAN/testing_tmp/esm_dian_2comp_brainstem_ages.mat');
ref_pattern(:,:,1) = roidata.test_data; % in the interval [0 1]. <----------------- !!!
% loading tau presence Probabilities, using a composite reference
%ref_pattern(:,:,2) = [rand(78,500)]; % in the interval [0 1]. <----------------- !!!
% loading tau presence Probabilities, using corresponding regional references (for "accounting" for spatial heterogeneity)
%ref_pattern(:,:,3) = [rand(78,500)]; % in the interval [0 1]. <----------------- !!!

AGEs = roidata.ages; % ages or eyo of all subjects, for calculating onset ages <----- !!!   

% Loading region-region connectivity information for defined atlas (e.g. DKT, Klein and Tourville Atlas [Frontiers in Neuroscience, 2012])
% load('/home/users/llevitis/code/GreyMatterESM/testing_tmp/ACP.mat');
Conn_Matrix = roidata.ACP; % Connectivity matrix (Nr X Nr, a common template for all subjects, in this case)  <----- !!!  
%load('/home/users/llevitis/code/GreyMatterESM/connectivity_CMU60DSI/Matrix_LONG.mat');
Long_Matrix = roidata.LONG; % distance matrix (in mm)  <----- !!!  
Conn_Matrix = Conn_Matrix + eye(size(Conn_Matrix)); % adding self connections
Nnodes      = size(Conn_Matrix,1);

%% Default parameters:
speed_Matrix = (168e-6)/(1/((1e+3)*60*60*24*365)); % For amyloid-beta. Converting 168 (nm/ms) to mm/years % Effective diff. coefficient: (1.8*1e-6/(1.6)^2)/(1/(60*60*24*365)); % mm2/years (based on Jack Waters, PLoS One, 2010)
% speed_Matrix_min and speed_Matrix_max % changing tortuosity according to S. Hrab?etov�a, C. Nicholson / Neurochemistry International 45 (2004).
% See also Nicholson et al., 2000, Progress in Brain Research, Vol 125.
Nparameters      = 3; % # parameters to estimate, i.e. production, clearance, noise level
time_range       = [0 50]; % in years
h_integration    = 7/365; % integration step in years (e.g. 1 day)
max_regions      = 2; % maximum number of initial seed regions
num_noise_levels = 10;
num_repetitions  = 1;

%[BETAS0,DELTAS0] = get_parameters_space(1e-1); % for quick exploration (delete this line and uncomment next one!!!)
[BETAS0,DELTAS0] = get_parameters_space(1e-2); % for 100 values (change this if further resolution is desired)
MUS0             = 0; % mean of the additive noise
max_sigma_value  = 0.1; % maximum std of the additive noise
SIGMAS0          = 0; % for quick exploration (delete this line and uncomment next one!!!)
%SIGMAS0          = 0:max_sigma_value/num_noise_levels:max_sigma_value; % std of the additive noise
Nsubjects        = size(ref_pattern,2); %

%% Selecting/evaluating automatically the "best" predictor model
for measure = 1:N_tau_possibilities
    disp(['Analysis for amyloid beta measurement -> ' num2str(measure)]);
    
    Nmodels     = 0; % or define also any apriori model based in literature
    counter     = 1;
    temp_model  = [];
    ant_mean_corrs = -1;
    model = Nmodels + 1;
    while counter <= max_regions,
        model_solutions  = zeros(Nnodes-length(temp_model),Nnodes,Nsubjects,'single');
        model_parameters = zeros(Nnodes-length(temp_model),Nparameters,Nsubjects,'single');
        model_times      = zeros(Nnodes-length(temp_model),Nsubjects,1,'single');
        model_CORRs      = -ones(Nnodes-length(temp_model),Nsubjects,1,'single');
        model_RMSEs      = 100*ones(Nnodes-length(temp_model),Nsubjects,1,'single');
        poss_nodes       = setdiff(1:Nnodes,temp_model);
        % matlabpool(2); % change to the desired number of cores
        % parfor node = 1:Nnodes
        for node = 1:Nnodes
            if ismember(node,poss_nodes)
                seed_regions_1 = [temp_model node];
                disp(['evaluating seed model/hyphotesis  -> ' num2str(Nmodels+1) ', seed regions: ' num2str(seed_regions_1)])
                [model_solutions(node,:,:),model_parameters(node,:,:),model_times(node,:,:),model_CORRs(node,:,:),model_RMSEs(node,:,:)] = ...
                    model_evaluation(Nnodes,Nsubjects,time_range,ref_pattern(:,:,measure),h_integration,Conn_Matrix,...
                    Long_Matrix,num_repetitions,Nparameters,BETAS0,DELTAS0,SIGMAS0,MUS0,seed_regions_1);
            end
        end
        % matlabpool('close');
        
        [i_node,j_node] = max(mean(model_CORRs,2),[],1); % across all subjects (can be optimized subject by subject)
        if mean(model_CORRs(j_node,:)) > ant_mean_corrs
            temp_model(counter) = j_node;
            ant_mean_corrs = mean(model_CORRs(j_node,:));
            Final_solutions(:,:,Nmodels + 1)  = model_solutions(j_node,:,:);
            Final_parameters(:,:,Nmodels + 1) = model_parameters(j_node,:,:);
            Final_times(:,Nmodels + 1) = model_times(j_node,:);
            Final_CORRs(:,Nmodels + 1) = model_CORRs(j_node,:);
            Final_RMSEs(:,Nmodels + 1) = model_RMSEs(j_node,:);
        else,
            break;
        end
        counter = counter + 1;                       
    end
    models(Nmodels + 1).regions = temp_model;
    
    %% Marginalizing production and clearance rates
    for model = 1:Nmodels+1
        disp(['Analyzing model ' num2str(model)])
        BETAS0_est  = Final_parameters(1,:,model)';
        DELTAS0_est = Final_parameters(2,:,model)';
        SIGMAS_est  = Final_parameters(3,:,model)';
        AGEs = reshape(AGEs, [size(roidata.ages,2),model]); %based on n
        ONSETS_est(:,model) = AGEs - Final_times(:,model);
        BETAS_est(:,model)  = (exp(-BETAS0_est)+BETAS0_est-1)./BETAS0_est;
        DELTAS_est(:,model) = (1-exp(-DELTAS0_est))./DELTAS0_est;
    end
    
    %% Saving all results
    save([output_dir filesep name '_' date '.mat'],'models','Final_solutions','Final_parameters','Final_times',...
        'Final_CORRs','Final_RMSEs','Nmodels','Nnodes','Nsubjects', 'ONSETS_est', 'time_range','ref_pattern','h_integration','Conn_Matrix',...
        'Long_Matrix','num_repetitions','Nparameters','BETAS0','DELTAS0','SIGMAS0','MUS0',...
        'AGEs','BETAS_est','DELTAS_est');
end






