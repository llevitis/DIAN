function [BETAS0,DELTAS0] = get_parameters_space(dI_beta)
% Yasser Iturria Medina,
% 20/06/2013, Montreal.

disp(['Generating parameters space...'])
tic
beta_step = 1e-4;
betas0 = 0:beta_step:100;
I = (exp(-betas0)+betas0-1)./betas0;

if nargin < 1,
    dI_beta = 1e-2;
end
Ib = eps:dI_beta:1-eps;
for i = 1:length(Ib)
    [k,ind_b] = min(abs(I - Ib(i)));
    BETAS0(i) = betas0(ind_b);
end
% figure; plot(BETAS0,Ib,'.'); grid on
BETAS0 = unique(BETAS0);
% figure; hist(BETAS0,50)
DELTAS0 = BETAS0;
toc
return
end