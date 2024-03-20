%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit Skew Exponential Power distribution parameters to RF residuals
% Created By: James Knighton (03/20/2024)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
clear all;

params = {'gsmax','P12','P50','P88','rdmax','WUE','height','SLA','LeafN'};

for m = 1:9

    disp(m);

    %Read in RF test dataset outputs
    infile = strcat('C:\PhyloTraitEst\RF_Model\RF_Test_Predict_',params{m},'.csv');
    RF_Res = readtable(infile);
    Res = table2array(RF_Res(:,2:3));
    Res(:,3) = Res(:,1) - Res(:,2);
    Res(:,1) = abs(Res(:,2));
    Res(:,2) = Res(:,3);
    Res(:,3) = [];

	% Monte Carol sampling of SEP distribution parameters
    for i = 1:1e6

        epsilon = unifrnd(0.75,1.25);     	% skew
        beta = unifrnd(0,2);       			% kurtosis
        sigma_0 = unifrnd(0,max(Res(:,1))); % heteroscedasiticty
        sigma_1 = unifrnd(0,1.25);     		% heteroscedasiticty

        L(i,1) = epsilon;
        L(i,2) = beta;
        L(i,3) = sigma_0;
        L(i,4) = sigma_1;

        Res(:,3) = Res(:,2);

        %Remove heteroscedasticity
        sigma_t = sigma_0 + sigma_1.*Res(:,1);
        Res(:,4) = Res(:,3)./sigma_t;

        %Compute SEP Log Likelihood
        n = length(Res(:,1));
        M1 = gamma(1 + beta)/(sqrt(gamma(3*(1+beta)/2))*sqrt(gamma((1 + beta)/2)));
        M2 = 1;
        mew_e = M1*(epsilon - epsilon^-1);
        sigma_e = sqrt((M2 - M1^2)*(epsilon^2 + epsilon^-2) + 2*(M1^2) - M2);
        wb = sqrt(gamma(3*(1 + beta)/2))/((1+beta)*gamma((1 + beta)/2)^(3/2));
        cb = (gamma(3*(1 + beta)/2)/gamma((1 + beta)/2))^(1/(1 + beta));
        alpha_wt = (epsilon.^sign(mew_e + sigma_e.*Res(:,4))).*(mew_e + sigma_e.*Res(:,4));
        L(i,5) = n*log((2*sigma_e*wb)/(epsilon + epsilon^-1)) - sum(log(sigma_t)) - cb*sum(abs(alpha_wt).^(2/(1 + beta)));
        L(i,6) = exp(L(i,5));

    end

	%Retain the MLE parameter estimate
    L = sortrows(L,-5);
    ModelResParams(m,1) = m;
    ModelResParams(m,2) = L(1,1); % skew
    ModelResParams(m,3) = L(1,2); % kurtosis
    ModelResParams(m,4) = L(1,3); % sigma 0
    ModelResParams(m,5) = L(1,4); % sigma x
    ModelResParams(m,6) = L(1,5); % log(L)
    ModelResParams(m,7) = L(1,6); % L

end
