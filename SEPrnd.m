function output = SEPrnd(epsilon,beta)
    
    %Generates SEP random numbers    
    M1 = gamma(1 + beta)/(sqrt(gamma(3*(1+beta)/2))*sqrt(gamma((1 + beta)/2)));
    M2 = 1;
    
    mew_e = M1*(epsilon - epsilon^-1);
    sigma_e = sqrt((M2 - M1^2)*(epsilon^2 + epsilon^-2) + 2*(M1^2) - M2);
    wb = sqrt(gamma(3*(1 + beta)/2))/((1+beta)*gamma((1 + beta)/2)^(3/2));
    cb = (gamma(3*(1 + beta)/2)/gamma((1 + beta)/2))^(1/(1 + beta));

    %Compute SEP PDF
        alpha = -10:0.01:10;
        alpha_et = (epsilon.^sign(mew_e + sigma_e.*alpha)).*(mew_e + sigma_e.*alpha);
        SEP_Dist(:,1) = alpha;
        SEP_Dist(:,2) = (2*sigma_e)*wb*exp(-1*cb*abs(alpha_et).^(2/(1 + beta)))/(epsilon + epsilon^-1);
    
    %Compute SEP CDF
    SEP_Dist(:,3) = 0.01*cumsum(SEP_Dist(:,2));
    
    %Draw random value from SEP CDF
    sample = unifrnd(0,1);
    tmp = abs(SEP_Dist(:,3) - sample);
    [idx, idx] = min(tmp);
    output = SEP_Dist(idx,1);
    
end