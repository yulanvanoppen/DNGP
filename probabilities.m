function [E, MAP] = probabilities(q_r, model, modelA, modelB, F)
    if nargin < 5
        F = eye(length(model.t));
    end
    
    [log_p, combinations] = posterior_Zn(q_r(1), q_r(2), model, modelA, modelB, F);

    [~, MAP_idx] = max(log_p);
    MAP = combinations(:, MAP_idx);

    E = exp(log_p) * combinations' / sum(exp(log_p));
end


% % % REFERENCE: post_q2()
% function log_p = posterior_q(logit_q, precomputed, Z, a, b)
%     q = 1 ./ (1 + exp(-logit_q));
% 
%     components = zeros(size(q));
%     components(q ~= 1) = log(1-q(q ~= 1));
%     components(q == 1) = log(q(q == 1)) - logit_q(q == 1);
% 
%     terms = sum(Z' .* logit_q + components, 2);
%     weights = exp(terms);
% 
%     p_data = log(sum(weights .* precomputed'));
% 
%     p_prior = sum((a-1) * log(q) + (b-1)*log(1-q));
% 
%     log_p = p_prior + p_data;
% end


% REFERENCE: Gibbs_Zn2()
function [log_p, combinations] = posterior_Zn(q, r, model, modelA, modelB, F)
    nseries = 200;
    subset_series = 1:min(nseries, min(size(modelA.y, 2), size(modelB.y, 2)));

    YA = modelA.y(:, subset_series);
    YB = modelB.y(:, subset_series);

    Si = model.flexible_covariance([model.th -Inf]);
    la = model.flexible_covariance(model.ze);

    SiA = modelA.flexible_covariance([modelA.th -Inf]);
    laA = modelA.flexible_covariance(modelA.ze);

    SiB = modelB.flexible_covariance([modelB.th -Inf]);
    laB = modelB.flexible_covariance(modelB.ze);

    if true %filter
        t_filtered = zeros(1, size(F, 1));
        t_filtered([1 end]) = model.t([1 end]);
        for idx = 2:length(t_filtered)-1
            nonzero = find(F(idx, :) ~= 0);
            t_filtered(idx) = nonzero((length(nonzero)+1)/2);
        end

        YA = F*YA;
        YB = F*YB;
        Si = F*Si*F';
        SiA = F*SiA*F';
        SiB = F*SiB*F';
        la = F*la*F';
        laA = F*laA*F';
        laB = F*laB*F';
    end

    N = length(t_filtered);
    I = eye(N);
    RA = size(YA, 2);
    RB = size(YB, 2);

    combinations = dec2bin(0:2^N-1)' - '0';
    ncomb = length(combinations);
    Z = combinations;

    [log_p0, log_p1a, log_p1b, log_p2, log_p3, log_p4, log_p5] = deal(zeros(1, ncomb));

    Sinv = svdinv(Si);
    SAinv = svdinv(SiA);
    SBinv = svdinv(SiB);

    for comb = 1:ncomb
        if mod(comb, 10000) == 0, fprintf('comb \n'), end
        dZ = diag(Z(:, comb));
        I_dZ = I - dZ;

        iLaZA = svdinv(dZ*laA*dZ + I_dZ*la*I_dZ);
        iLaZB = svdinv(dZ*laB*dZ + I_dZ*la*I_dZ);

        Sf = svdinv(Sinv + I_dZ * (RA*iLaZA + RB*iLaZB) * I_dZ);
        SfA = svdinv(SAinv + dZ * RA*iLaZA * dZ);
        SfB = svdinv(SBinv + dZ * RB*iLaZB * dZ);

        mf = Sf * I_dZ * (mean(RA*iLaZA*YA, 2) + mean(RB*iLaZB*YB, 2));
        mfA = SfA * RA * dZ * iLaZA * mean(YA, 2);
        mfB = SfB * RB * dZ * iLaZB * mean(YB, 2);

        % switchpoints = sum(abs(diff(Z(:, comb))));
        % log_p0(comb) = switchpoints * log(max(eps, q)) + (N-1 - switchpoints) * log(max(eps, 1-q));
        % log_p0(comb) = sum(Z(:, comb)) * log(r) + (N - sum(Z(:, comb))) * log(1-r);

        % switch_on = sum(max(0, diff(Z(:, comb))));
        % switch_off = sum(max(0, -diff(Z(:, comb))));
        % stay_on = sum(~diff(Z(:, comb)) & Z(2:end, comb));
        % stay_off = sum(~diff(Z(:, comb)) & ~Z(2:end, comb));
        % log_p0(comb) = switch_on * log(max(eps, q(1))) + stay_off * log(max(eps, 1-q(1))) ...
        %              + switch_off * log(max(eps, q(2))) + stay_on * log(max(eps, 1-q(2)));

        switch_on = sum(max(0, diff(Z(:, comb))));
        switch_off = sum(max(0, -diff(Z(:, comb))));
        stay_on = sum(~diff(Z(:, comb)) & Z(2:end, comb));
        stay_off = sum(~diff(Z(:, comb)) & ~Z(2:end, comb));
        log_p0(comb) = switch_on * log(max(eps, q)) + stay_off * log(max(eps, 1-q)) ...
                     + switch_off * log(max(eps, r)) + stay_on * log(max(eps, 1-r));

        log_p1a(comb) = - 1/2 * RA * svd_logdet(dZ*laA*dZ + I_dZ*la*I_dZ);
        log_p1b(comb) = - 1/2 * RB * svd_logdet(dZ*laB*dZ + I_dZ*la*I_dZ);
        log_p2(comb) = 1/2 * (svd_logdet(Sf) + svd_logdet(SfA) + svd_logdet(SfB));
        log_p3(comb) = 1/2 * (mf' * svdinv(Sf) * mf + mfA' * svdinv(SfA) * mfA ...
                                                    + mfB' * svdinv(SfB) * mfB);
        for cell = 1:RA
            log_p4(comb) = log_p4(comb) - 1/2 * YA(:, cell)' * iLaZA * YA(:, cell);
        end
        for cell = 1:RB
            log_p5(comb) = log_p5(comb) - 1/2 * YB(:, cell)' * iLaZB * YB(:, cell);
        end
    end

    log_p = log_p0 + log_p1a + log_p1b + log_p2 ...
                   + log_p3 + log_p4 + log_p5;

    [~, order] = sort(log_p);
    combinations(:, order);

    log_p = log_p - max(log_p);
    % p = exp(log_p) / sum(exp(log_p));
    % idx = mnrnd(1, p);
    % z = combinations(:, find(idx));
end



function logdet = svd_logdet(x)
    % S = svd(x);
    % values = S(S > max(S)*1e-3);
    logdet = sum(log(svd(x)));
end