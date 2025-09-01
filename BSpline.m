classdef BSpline < handle
    properties (SetAccess = immutable)
        ord                                                                 % spline order
        knots                                                               % knot locations
        card                                                                % basis cardinality
    end
    
    properties (SetAccess = private)
        const_basis
        lin_basis
        quad_basis
        cubic_basis

        t_base
        cubic_basis_base
    end
    
    methods                                                             % Constructor
        function obj = BSpline(ord, knots, t_base)
            obj.ord = ord;
            if numel(knots) == 1
                obj.knots = 1/knots : 1/knots : 1-1/knots;
                obj.card = ord + knots - 1;
            else
                obj.knots = knots(knots ~= 0 & knots ~= 1);
                obj.card = ord + sum(~ismember(knots, [0 1]));
            end
            
            obj.t_base = t_base;
            obj.make_basis();
        end
        
        
        function make_basis(obj)
            K = length(obj.knots) + 1;
                                
            obj.const_basis = @(t) [0 * t; ...                              % auxiliary zeros
                                    t >= [0; obj.knots'] & t < [obj.knots'; 1] ...  % intervals
                                    + [t .* zeros(K-1, 1); t >= 1];         % upper limit correction
                                    0 * t];                                 % auxiliary zeros
            
            comb = eye(2*(K+2));
            comb = [zeros(1, 2*(K+2)); comb(1:K+1, :) + comb(K+4:end, :); zeros(1, 2*(K+2))];
            obj.lin_basis = @(t) comb * (repmat(obj.const_basis(t), 2, 1) .* ...
                                 [[(t - [0; 0; obj.knots']) ./ max(eps, [0; obj.knots'; 1] - [0; 0; obj.knots']); 0 * t];
                                 [0 * t; ([obj.knots'; 1; 1] - t) ./ max(eps, [obj.knots'; 1; 1] - [0; obj.knots'; 1])]]);
                          
            comb = eye(2*(K+3));
            comb = [zeros(1, 2*(K+3)); comb(1:K+2, :) + comb(K+5:end, :); zeros(1, 2*(K+3))];
            obj.quad_basis = @(t) comb * (repmat(obj.lin_basis(t), 2, 1) .* ...
                                  [[(t - [0; 0; 0; obj.knots']) ./ max(eps, [0; obj.knots'; 1; 1] - [0; 0; 0; obj.knots']); 0 * t];
                                  [0 * t; ([obj.knots'; 1; 1; 1] - t) ./ max(eps, [obj.knots'; 1; 1; 1] - [0; 0; obj.knots'; 1])]]);
                              
            comb = eye(2*(K+4));
            comb = comb(1:K+3, :) + comb(K+6:end, :);
            obj.cubic_basis = @(t) (comb * (repmat(obj.quad_basis(t), 2, 1) .* ...
                                  [[(t - [0; 0; 0; 0; obj.knots']) ./ max(eps, [0; obj.knots'; 1; 1; 1] - [0; 0; 0; 0; obj.knots']); 0 * t];
                                  [0 * t; ([obj.knots'; 1; 1; 1; 1] - t) ./ max(eps, [obj.knots'; 1; 1; 1; 1] - [0; 0; 0; obj.knots'; 1])]]))';

            obj.cubic_basis_base = obj.cubic_basis(obj.t_base);
        end
        
        
        function y = spline_function(obj, beta)                         % Evaulate a spline with basis coefficients beta at x
            dim = numel(beta);
            y = @(t) reshape(obj.cubic_basis(t) * reshape(beta, obj.card, []), [length(t) dim(2:end)]);
            
        end
        
                                            
        function y = spline(obj, beta, t)                               % Evaulate a spline with basis coefficients beta at x
            if nargin == 3
                y = reshape(obj.cubic_basis(t) * reshape(beta, obj.card, []), length(t), 1);   % convert back to matrix
            else
                y = reshape(obj.cubic_basis_base * reshape(beta, obj.card, []), length(obj.t_base), 1);
            end

        end
        
                                                                        % Evaulate the derivative of BSpline.spline()
        function y = spline_diff(obj, x, beta)
            J = size(beta, 1);                                              % extract the basis cardinality and reshape beta
            beta_reshaped = reshape(beta, J, 1, []);

                                                                            % compute the linear basis combination
            y = sum(beta_reshaped .* obj.basis_diff(x, 1:J));
            y = permute(y, [3 2 1]);                                        % convert back to matrix
        end
        
        
         function y = spline_ddiff(obj, x, beta)
            J = size(beta, 1);                                              % extract the basis cardinality and reshape beta
            beta_reshaped = reshape(beta, J, 1, []);

                                                                            % compute the linear basis combination
            y = sum(beta_reshaped .* obj.basis_ddiff(x, 1:J));
            y = permute(y, [3 2 1]);                                        % convert back to matrix
         end       
        
                                                                        % B-spline basis function as in Tan and Ghosal (2019)
        function y = basis(obj, x, s)                                   % and Appendix E.2 in Ghosal and van der Vaart (2017)
            if nargin < 3, s = 1:obj.card; end
            y = obj.recursion(x, s-obj.ord, obj.ord);
        end
        
                                          
        function y = basis_diff(obj, x, s)                              % Corresponding derivatives
            if nargin < 3, s = 1:obj.card; end
            y = obj.recursion_diff(x, s-obj.ord, obj.ord);
        end

        
        function y = basis_ddiff(obj, x, s)                             % Corresponding 2nd-order derivatives
            if nargin < 3, s = 1:obj.card; end
            y = obj.recursion_ddiff(x, s-obj.ord, obj.ord);
        end
        

        function y = recursion(obj, x, j, q)                            % Recursive helper function for basis() 
            K = length(obj.knots) + 1;                                      % number of segments
                                                                            % index to segment point
            if K == 1
                t = @(k) 1 * (k > K-1);
            else
                t = @(k) (k > 0 & k <= K-1) .* obj.knots(max(1, min(k, K-1))) + (k > K-1);
            end

            if q == 1
                y = t(j)' <= x & x < t(j+1)' | x == 1 & t(j+1)' == 1;       % execute base case or recursion
            else
                y = (x - t(j)')   ./ max(eps, (t(j+q-1) - t(j))') .* obj.recursion(x, j,   q-1) ...
                  + (t(j+q)' - x) ./ max(eps, (t(j+q) - t(j+1))') .* obj.recursion(x, j+1, q-1);
            end
        end
        

        function y = recursion_diff(obj, x, j, q)                       % Recursive helper function for basis_diff()
            K = length(obj.knots) + 1;                                      % number of segments
                                                                            % index to segment point
            if K == 1
                t = @(k) 1 * (k > K-1);
            else
                t = @(k) (k > 0 & k <= K-1) .* obj.knots(max(1, min(k, K-1))) + (k > K-1);
            end

            switch q                                                        % execute base case or recursion
                case 1                          
                    y = zeros(length(j), length(x));
                case 2
                    y = (t(j)' <= x & x < t(j+1)' | x == 1 & j' == K-1) ./ max(eps, (t(j+1) - t(j))') ...
                        - (t(j+1)' <= x & x < t(j+2)' | x == 1 & j' == K-2) ./ max(eps, (t(j+2) - t(j+1))');
                otherwise
                    y = ((x - t(j)')   .* obj.recursion_diff(x, j,   q-1) + obj.recursion(x, j,   q-1)) ...
                        ./ max(eps, (t(j+q-1) - t(j))')  ...
                      + ((t(j+q)' - x) .* obj.recursion_diff(x, j+1, q-1) - obj.recursion(x, j+1, q-1)) ... 
                        ./ max(eps, (t(j+q) - t(j+1))');
            end
        end
        
        
        function y = recursion_ddiff(obj, x, j, q)                      % Recursive helper function for basis_ddiff()
            K = length(obj.knots) + 1;                                      % number of segments
                                                                            % index to segment point
            if K == 1
                t = @(k) 1 * (k > K-1);
            else
                t = @(k) (k > 0 & k <= K-1) .* obj.knots(max(1, min(k, K-1))) + (k > K-1);
            end

            switch q                                                        % execute base case or recursion
                case {1, 2}                          
                    y = zeros(length(j), length(x));
                case 3
                    y =   (j >=  0)' .* 2 .* (t(j)'   <= x & x < t(j+1)' | x == 1 & j' == K-3)       ./ max(eps, (t(j+2) - t(j)  )' .* (t(j+1) - t(j)  )') ...
                        - (j >= -1)' .* 2 .* (t(j+1)' <= x & x < t(j+2)' | x == 1 & j' == K-2) .* (1 ./ max(eps, (t(j+2) - t(j)  )' .* (t(j+2) - t(j+1))') + 1 ./ max(eps, (t(j+3) - t(j+1))' .* (t(j+2) - t(j+1))')) ...
                        + (j >= -2)' .* 2 .* (t(j+2)' <= x & x < t(j+3)')                            ./ max(eps, (t(j+3) - t(j+1))' .* (t(j+3) - t(j+2))') ...
                        +               2 .* (x == 1 & j' == K-1)                                    ./ max(eps, (t(j+1) - t(j)  )' .* (t(j+1) - t(j)  )');
                otherwise
                    y = ((x - t(j)')   .* obj.recursion_ddiff(x, j,   q-1) + 2 * obj.recursion_diff(x, j,   q-1)) ...
                        ./ max(eps, (t(j+q-1) - t(j))')  ...
                      + ((t(j+q)' - x) .* obj.recursion_ddiff(x, j+1, q-1) - 2 * obj.recursion_diff(x, j+1, q-1)) ... 
                        ./ max(eps, (t(j+q) - t(j+1))');
            end
        end
    end
end