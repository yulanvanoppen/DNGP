classdef PoEmodel < handle
    properties %(SetAccess = private)
        t
        y
        R
        
        bspline
        nbasis
        options
        
        prior_init
        prior_mean
        prior_resid
        
        zeta_opt
        theta_opt

        ze
        th
        
        muf
        Sf
        sf
    end
    
    methods
        function obj = PoEmodel(t, y)
            obj.t = t;
            obj.y = y;
            obj.R = size(y, 2);
            
            obj.bspline = BSpline(4, [.25 .5 .75], (t - t(1))/range(t));
            obj.nbasis = obj.bspline.card;
            
            obj.options = optimoptions('fmincon', 'Display', 'off', 'OptimalityTolerance', 1e-6);
            
                                                        % coefficients ~ N(0, 10^2)
            obj.prior_init = @(par) sum(log(normpdf(par(1:obj.nbasis), 0, 10))) ...     
                ...                                     % length_scales L(t1, t2) ~ Gamma(T/5, 1)
                + sum(log(gampdf(exp(par(2*obj.nbasis+1:end-1)), range(t)/5, 1))) ...
                ...                                     % amplitudes and noise ~ G(1, 1)
                + sum(log(gampdf(exp(par([obj.nbasis+1:2*obj.nbasis, end])), 1, 1)));
                
            % obj.prior_mean = @(par) sum(log(gampdf(exp(par(obj.nbasis+1:end)), range(t)/10, 1)))...
            %     + sum(log(gampdf(exp(par(1:obj.nbasis)), 1, 1)));
            obj.prior_mean = @(par) log(gampdf(exp(par(2)), range(t)/10, 1))...
                +                   + log(gampdf(exp(par(1)), 1, 1));

            obj.prior_resid = @(par) sum(log(gampdf(exp(par(obj.nbasis+1:end-1)), range(t)/10, 1)))...
                + sum(log(gampdf(exp(par([1:obj.nbasis end])), 1, 1)));
        end
        
        
        function fit(obj, max_iter, tol)
            if nargin < 3, tol = 1e-3; end
            if nargin < 2, max_iter = 50; end

            zeta_opt_full  = obj.fit_initial();
            obj.zeta_opt = zeta_opt_full(obj.nbasis+1:end);
            obj.theta_opt  = obj.fit_mean(obj.zeta_opt);
            
            obj.zeta_opt = [obj.zeta_opt; zeros(max_iter, length(obj.zeta_opt))];
            obj.theta_opt = [obj.theta_opt; zeros(max_iter, length(obj.theta_opt))];
            
            for iter = 1:max_iter
                Sigma = obj.flexible_covariance([obj.theta_opt(iter, :) -Inf]);
                lambda = obj.flexible_covariance(obj.zeta_opt(iter, :));
                Sf_iter = pinv(pinv(Sigma) + obj.R*(pinv(lambda)));
                muf_iter = obj.R * Sf_iter * pinv(lambda) * mean(obj.y, 2);

                try
                obj.zeta_opt(iter+1, :) = fmincon(@(par) -obj.logpost_cond(par, muf_iter), ...
                                              obj.zeta_opt(iter, :), [], [], [], [], ...
                                              -10*ones(1, size(obj.zeta_opt, 2)), ...
                                              10*ones(1, size(obj.zeta_opt, 2)), [], obj.options);

                obj.theta_opt(iter+1, :) = fmincon(@(par) -obj.logpost_full(par, obj.zeta_opt(iter+1, :)), ...
                                           obj.theta_opt(iter, :), [], [], [], [], ...
                                           -10*ones(1, size(obj.theta_opt, 2)), ...
                                           10*ones(1, size(obj.theta_opt, 2)), [], obj.options);
                catch ME
                    disp(ME)
                end

                if eucl_rel(obj.zeta_opt(iter, :), obj.zeta_opt(iter+1, :)) < tol ...
                   && eucl_rel(obj.theta_opt(iter, :), obj.theta_opt(iter+1, :)) < tol
                    break
                end
            end

            % obj.theta_opt(1, :) = obj.theta_opt(iter+1, :);
            % obj.zeta_opt(1, :) = obj.zeta_opt(iter+1, :);
            % obj.theta_opt(2:end, :) = 0;
            % obj.zeta_opt(2:end, :) = 0;
            % 
            % for iter = 1:max_iter
            %     obj.zeta_opt(iter+1, :) = fmincon(@(par) -obj.logpost_full(obj.theta_opt(iter, :), par), ...
            %                                   obj.zeta_opt(iter, :), [], [], [], [], ...
            %                                   -10*ones(1, size(obj.zeta_opt, 2)), ...
            %                                   10*ones(1, size(obj.zeta_opt, 2)), [], obj.options);
            % 
            %     obj.theta_opt(iter+1, :) = fmincon(@(par) -obj.logpost_full(par, obj.zeta_opt(iter+1, :)), ...
            %                                obj.theta_opt(iter, :), [], [], [], [], ...
            %                                -10*ones(1, size(obj.theta_opt, 2)), ...
            %                                10*ones(1, size(obj.theta_opt, 2)), [], obj.options);
            % 
            %     if eucl_rel(obj.zeta_opt(iter, :), obj.zeta_opt(iter+1, :)) < tol ...
            %        && eucl_rel(obj.theta_opt(iter, :), obj.theta_opt(iter+1, :)) < tol
            %         break
            %     end
            % end

            obj.theta_opt = obj.theta_opt(1:iter+1, :);
            obj.zeta_opt = obj.zeta_opt(1:iter+1, :);

            obj.th = obj.theta_opt(end, :);                                 % abbreviate for convenience
            obj.ze = obj.zeta_opt(end, :);
            
            Sigma = obj.flexible_covariance([obj.theta_opt(end, :) -Inf]);
            lambda = obj.flexible_covariance(obj.zeta_opt(end, :));
            obj.Sf = pinv(pinv(Sigma) + obj.R*(pinv(lambda)));
            obj.muf = obj.R * obj.Sf * pinv(lambda) * mean(obj.y, 2);
            obj.sf  = sqrt(diag(obj.Sf)); 
        end
        
        
        function par = fit_initial(obj)
            opt = Inf;
            par = zeros(1, 3*obj.nbasis+1);
            
            for n = 1:20
                fprintf('%3d', n)
                try
                    [MLE_par, MLE] = fmincon(@(par) -obj.loglik(par), ...
                                             randn(1, length(par)), [], [], [], [], ...
                                             [-5*ones(1, length(par)-1) -5], ...
                                             5*ones(1, length(par)), [], obj.options);
                catch
                    continue;
                end
                if MLE<opt
                    par = MLE_par;
                    opt = MLE;
                end
            end
            fprintf('\n')
        end
        
        
        function par = fit_mean(obj, zeta_opt)
            opt = Inf;
            % par = zeros(1, 2*obj.nbasis);
            par = zeros(1, 2);
            
            for n = 1:20
                fprintf('%3d', n)
                try
                    [MLE_par, MLE] = fmincon(@(par) -obj.logpost_full(par, zeta_opt), ...
                                             randn(1, length(par)), [], [], [], [], ...
                                             -5*ones(1, length(par)), ...
                                             5*ones(1, length(par)), [], obj.options);
                catch
                    continue;
                end
                if MLE<opt
                    par = MLE_par;
                    opt = MLE;
                end
            end
            fprintf('\n')
        end
        
        
        function p = loglik(obj, par)
            warning('off','MATLAB:singularMatrix')      % surpress warnings for singular matrices
            warning('off','MATLAB:eigs:IllConditionedA')

            N = numel(obj.y) / length(obj.t);
            mu_par = par(1:obj.nbasis);
            Sigma_par = par(obj.nbasis+1:end);

            mu = obj.bspline.cubic_basis_base * mu_par';% mean
            Sigma = obj.flexible_covariance(Sigma_par); % covariance

            constant = -N/2 * length(obj.t) * log(2*pi);% denominator terms
            det_term = -N/2 * log(det(Sigma));          % exponentiated term
            exp_term = -1/2 * sum(pagemtimes(reshape(obj.y - mu, [], 1, N), 'transpose', ...
                                             reshape(Sigma \ (obj.y - mu), [], 1, N), 'none'));

                                                        % log likelihood
            p = constant + det_term + exp_term + obj.prior_init([mu_par Sigma_par]);

            warning('on')
        end
        
        
        function p = logpost_full(obj, theta, zeta)
            warning('off','MATLAB:singularMatrix')
            warning('off','MATLAB:eigs:IllConditionedA')

            N = numel(obj.y) / length(obj.t);

            Sigma = obj.flexible_covariance([theta -Inf]);
            Sigmainv = pinv(Sigma);
            lambda = obj.flexible_covariance(zeta);
            lambdainv = pinv(lambda);
            covinv = Sigmainv + N*lambdainv;
            mu = N * pinv(covinv) * lambdainv * mean(obj.y, 2);

            logdet_covinv = log(det(covinv));
            if isinf(logdet_covinv)
                logdet_covinv = length(obj.t) * log(N) + log(det(lambdainv));
            end

            det_term = -N/2 * logdet_covinv - 1/2 * log(det(Sigma));
            exp_term = 1/2 * mu' * covinv * mu;

            p = det_term + exp_term + obj.prior_mean(theta);

            warning('on')
        end
        
        
        function p = logpost_cond(obj, zeta, mu)
            warning('off','MATLAB:singularMatrix')
            warning('off','MATLAB:eigs:IllConditionedA')

            N = numel(obj.y) / length(obj.t);

            Sigma = obj.flexible_covariance(zeta);

            constant = -N/2 * length(obj.t) * log(2*pi);
            det_term = -N/2 * log(det(Sigma));
            exp_term = -1/2 * sum(pagemtimes(reshape(obj.y - mu, [], 1, N), 'transpose', ...
                                             reshape(Sigma \ (obj.y - mu), [], 1, N), 'none'));
                                         
            p = constant + det_term + exp_term + obj.prior_resid(zeta);

            warning('on')
        end
        
        
        function K = flexible_covariance(obj, par, tp)
            if size(par, 1) < size(par, 2), par = par'; end
            if nargin < 3
                tp = obj.t;
                smooth_function = @(coef) exp(obj.bspline.cubic_basis_base * coef);
            else
                smooth_function = @(coef) exp(obj.bspline.spline(coef, (tp - tp(1)) / range(tp)));
            end

            n_coef = (length(par)-1) / 2;
            if n_coef > 1
                amplitudes = smooth_function(par(1:n_coef));    % variable GP amplitude and length scale
                length_scales = smooth_function(par(n_coef+1:end-1));
            else
                amplitudes = exp(par(1));
                length_scales = exp(par(2));
            end
            noise = exp(par(end));

            nu = 5/2;   
            denom = sqrt((length_scales.^2 + length_scales'.^2)/2);
            inner = sqrt(2*nu) * abs(tp - tp') ./ denom;
                                                                % variable length scale normalization
            normalizer = sqrt(length_scales .* length_scales') ./ denom;

                                                                % matern correlation function
            K = normalizer .* 2^(1-nu) / gamma(nu) .* inner.^nu .* besselk(nu, inner);
            K(logical(eye(size(K)))) = 1;

            K(K < 1e-12) = 0;                                   % nullify negligible elements
            K = amplitudes * amplitudes' .* K + noise^2 * eye(size(K));   % scale according to variance
        end
    end
end