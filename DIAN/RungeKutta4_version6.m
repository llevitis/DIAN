function [S,U,T] = RungeKutta4_version6(h,S0,Conn_Matrix,Time_0,Time_Points,beta0,delta0,Noise)
% Four order Runge Kutta method (without delays) to solve disease progresion model´s
% differential equations
% INPUTS
% h:  integration step (in years)
% S0: initial probability values (column vector, for Nnodes)
% Conn_Matrix: matrix (Nnodes x Nnodes).
% Time_0: initial time (years)
% Time_Points: time difference (in years) between final estimation times and initial
% time Time_0.
% beta0: parameter for birth rate.
% delta0: parameter for cure rate.
% noise_mu (optional): mean of distribution noise (gaussian).
% noise_sigma (optional): std of distribution noise (gaussian).
% OUTPUTS
% S: final solutions for Time_Points.
% U: temporal progresion from Time_0 to Time_Points.
% T: time points (from step to step) from Time_0 to Time_Points.
%-------------------------------------------------------------------------%
% Yasser Iturria Medina,
% 20/06/2013, Montreal.

x1 = Time_0;
y1 = S0;
tf = Time_Points; % difference between final estimation times and Time_0.
Nnodes = length(Conn_Matrix);

I = ceil(max(tf)/h); %number of points
U          = zeros(Nnodes,I+1);
U(:,1)     = S0;
T          = zeros(1,ceil(I)+1);
T(1)       = Time_0;
Ulag      = zeros(Nnodes,1);

if exist('Noise') && ~isempty(Noise)
    add_noise = 1;
else
    add_noise = 0;
end
Conn_Diag   = diag(Conn_Matrix); % Self connections
Conn_Matrix = Conn_Matrix - Conn_Matrix.*eye(Nnodes,Nnodes); % Eliminating self connections
for i = 1:I
    %% evaluating k1:                                 k1=h*fh(x1+(i-1)*h,u);
    Ulag = U(:,i);
    
    % Computing network dispersion factor
    VALUES_ord  = sort(Ulag');
    Disp_index  = (2/Nnodes)*sum(VALUES_ord.*(1:Nnodes))/sum(VALUES_ord) - (Nnodes + 1)/Nnodes; % Gini index
    
    % beta_t = beta0*Ulag;
    % delta_t = -delta0*U(:,i) + 1;
    beta_t   = (1-exp(-beta0*Ulag));
    delta_t  = exp(-delta0*U(:,i));
    if add_noise
        k1 = h*(-delta_t.*U(:,i) + (1-U(:,i)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*U(:,i)) + Noise(:,i,1));
    else
        k1 = h*(-delta_t.*U(:,i) + (1-U(:,i)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*U(:,i)));
    end
    
    %% evaluating k2:                    k2=h*fh(x1+(i-1)*h+(h/2),u+(k1/2));
    Ulag = U(:,i) + k1/2;
    
    % Computing network dispersion factor
    VALUES_ord  = sort(Ulag');
    Disp_index  = (2/Nnodes)*sum(VALUES_ord.*(1:Nnodes))/sum(VALUES_ord) - (Nnodes + 1)/Nnodes; % Gini index
    
    % beta_t  = beta0*Ulag;
    % delta_t = -delta0*(U(:,i)+(k1/2)) + 1;
    beta_t   = (1-exp(-beta0*Ulag));
    delta_t = exp(-delta0*(U(:,i)+(k1/2)));
    if add_noise
        k2 = h*(-delta_t.*(U(:,i)+(k1/2)) + (1-(U(:,i)+k1/2)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k1/2)) + Noise(:,i,2));
    else
        k2 = h*(-delta_t.*(U(:,i)+(k1/2)) + (1-(U(:,i)+k1/2)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k1/2)));
    end
    
    %% evaluating k3:                    k3=h*fh(x1+(i-1)*h+(h/2),u+(k2/2));
    Ulag = U(:,i) + k2/2;
    
    % Computing network dispersion factor
    VALUES_ord  = sort(Ulag');
    Disp_index  = (2/Nnodes)*sum(VALUES_ord.*(1:Nnodes))/sum(VALUES_ord) - (Nnodes + 1)/Nnodes; % Gini index
    
    % beta_t  = beta0*Ulag;
    % delta_t = -delta0*(U(:,i)+(k2/2)) + 1;
    beta_t   = (1-exp(-beta0*Ulag));
    delta_t = exp(-delta0*(U(:,i)+(k2/2)));
    if add_noise
        k3 = h*(-delta_t.*(U(:,i)+(k2/2)) + (1-(U(:,i)+k2/2)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k2/2)) + Noise(:,i,3));
    else
        k3 = h*(-delta_t.*(U(:,i)+(k2/2)) + (1-(U(:,i)+(k2/2))).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k2/2)));
    end
    
    %% evaluating k4:                        k4=h*fh(x1+(i-1)*h + h,u + k3);
    Ulag = U(:,i) + k3;
    
    % Computing network dispersion factor
    VALUES_ord  = sort(Ulag');
    Disp_index  = (2/Nnodes)*sum(VALUES_ord.*(1:Nnodes))/sum(VALUES_ord) - (Nnodes + 1)/Nnodes; % Gini index
    
    % beta_t  = beta0*Ulag;
    % delta_t = -delta0*(U(:,i)+k3) + 1;
    beta_t   = (1-exp(-beta0*Ulag));
    delta_t = exp(-delta0*(U(:,i)+k3));
    if add_noise
        k4 = h*(-delta_t.*(U(:,i)+k3) + (1-(U(:,i)+k3)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k3)) + Noise(:,i,4));
    else
        k4 = h*(-delta_t.*(U(:,i)+k3) + (1-(U(:,i)+k3)).*(Conn_Matrix*(Disp_index*beta_t.*Ulag) + (1-Disp_index)*Conn_Diag.*(U(:,i)+k3)));
    end
    
    % Final results on time point i
    T(i+1) = x1 + i*h;
    U(:,i+1) = U(:,i) + (1/6*(k1 + 2*(k2 + k3) + k4));
end
% Interpolating if needed
if length(nonzeros(Time_Points)) == 1
    S = U(:,end);
else
    for node = 1:Nnodes
        S(node,:) = interp1(T,U(node,:),Time_0 + nonzeros(Time_Points));
    end
end