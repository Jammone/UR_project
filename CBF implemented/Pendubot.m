classdef Pendubot < CtrlAffineSys
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        m1 = 0;
        m2 = 0;
        L1 = 0;
        Lc1 = 0;
        Lc2 = 0;
        L2 = 0;
        grav = 9.81;
        I1 = 0;
        I2 = 0;
    end
    
    methods
       function cbf = defineCbf(obj, params, symbolic_state)
            x = symbolic_state;
            amax = 0.2;
           cbf = 0.5*(amax^2 - (x(1)+x(2)-pi/2)-0.1*x(4)^2);
       end
       function [x, f, g] = defineSystem(obj, params)
            syms q1 q2 dq1 dq2
            x = [q1; q2; dq1; dq2];
            obj.m1 = params.m1;
            obj.m2 = params.m2;
            obj.L1 = params.L1;
            obj.L2 = params.L2;
            obj.grav = params.grav;
            obj.I1 = params.I1;
            obj.I2 = params.I2;
            obj.Lc1 = params.Lc1;
            obj.Lc2 = params.Lc2;
            d1= obj.L1-obj.Lc1;
            d2 =obj.L2-obj.Lc2;
            
M = [(obj.L1^2 + 2*cos(x(2))*obj.L1*obj.Lc2 + obj.Lc2^2)*obj.m2 + obj.m1*obj.Lc1^2 + obj.I1 + obj.I2, obj.Lc2*(obj.Lc2 + obj.L1*cos(x(2)))*obj.m2 + obj.I2;
     (obj.Lc2^2 + obj.L1*cos(x(2))*obj.Lc2)*obj.m2 + obj.I2,     obj.m2*obj.Lc2^2 + obj.I2];
            
       c = [x(4)*obj.L1*obj.m2*sin(x(2))*(x(3)+ x(4))*(d2 - obj.L2) + x(3)*x(4)*obj.L1*obj.m2*sin(x(2))*(d2 - obj.L2);
        -x(3)^2*obj.L1*obj.m2*sin(x(2))*(d2 - obj.L2)];
           
           e = [ obj.grav*(obj.L1*obj.m2*cos(x(1)) + obj.Lc1*obj.m1*cos(x(1))+ obj.Lc2*obj.m2*cos(x(1) + x(2)));
                    obj.grav*obj.Lc2*obj.m2*cos(x(1) + x(2))];
            
            n = c+e+[0.5*x(3);0.5*x(4)];
            
            IM = inv(M);
            
            % 4 x 1
            f = [x(3);
                x(4);
                -IM*n];

            % 4 x 2 
            g =[0,0;
                0,0;
                IM];
       end
       
      function initSys(obj, symbolic_x, symbolic_f, symbolic_g, symbolic_cbf, symbolic_clf)
            if isempty(symbolic_x) || isempty(symbolic_f) || isempty(symbolic_g)
                error('x, f, g is empty. Create a class function defineSystem and define your dynamics with symbolic expression.');
            end

            if ~isa(symbolic_f, 'sym')
                f_ = sym(symbolic_f);
            else
                f_ = symbolic_f;
            end
            if ~isa(symbolic_g, 'sym')
                g_ = sym(symbolic_g);
            else
                g_ = symbolic_g;
            end

            x = symbolic_x;
            % Setting state and input dimension.
            obj.xdim = size(x, 1);
            obj.udim = size(g_, 2);
            % Setting f and g (dynamics)
            obj.f = matlabFunction(f_, 'vars', {x});
            obj.g = matlabFunction(g_, 'vars', {x});            

            % Obtaining Lie derivatives of CBF.
            if ~isempty(symbolic_cbf)
                dcbf = simplify(jacobian(symbolic_cbf, symbolic_x));
                lf_cbf_ = dcbf * f_;
                lg_cbf_ = dcbf * g_;        
                obj.cbf = matlabFunction(symbolic_cbf, 'vars', {x});
                obj.lf_cbf = matlabFunction(lf_cbf_, 'vars', {x});
                % TODO: add sanity check of relative degree.
                obj.lg_cbf = matlabFunction(lg_cbf_, 'vars', {x});
            end

            % Obtaining Lie derivatives of CLF.    
            if ~isempty(symbolic_clf)
                dclf = simplify(jacobian(symbolic_clf, symbolic_x));
                lf_clf_ = dclf * f_;
                lg_clf_ = dclf * g_;
                obj.clf = matlabFunction(symbolic_clf, 'vars', {x});                       
                obj.lf_clf = matlabFunction(lf_clf_, 'vars', {x});
                % TODO: add sanity check of relative degree.
                obj.lg_clf = matlabFunction(lg_clf_, 'vars', {x});        
            end
      end
      function [u, B, feas, comp_time] = ctrlCbfQp(obj, x, u_ref, verbose)
            %% Implementation of vanilla CBF-QP
            % Inputs:   x: state
            %           u_ref: reference control input
            %           verbose: flag for logging (1: print log, 0: run silently)
            % Outputs:  u: control input as a solution of the CBF-CLF-QP
            %           B: Value of the CBF at current state.
            %           feas: 1 if QP is feasible, 0 if infeasible. (Note: even
            %           when qp is infeasible, u is determined from quadprog.)
            %           compt_time: computation time to run the solver.
            if isempty(obj.cbf)
                error('CBF is not defined so ctrlCbfQp cannot be used. Create a class function [defineCbf] and set up cbf with symbolic expression.');
            end

            if nargin < 3
                u_ref = zeros(obj.udim, 1);
            end
            if nargin < 4
                % Run QP without log in default condition.
                verbose = 0;
            end

            if size(u_ref, 1) ~= obj.udim
                error("Wrong size of u_ref, it should be (udim, 1) array.");
            end                

            tstart = tic;
            B = obj.cbf(x);
            LfB = obj.lf_cbf(x);
            LgB = obj.lg_cbf(x);

            %% Constraints : A * u <= b
            % CBF constraint.
            A = [-LgB];
            b = [LfB + obj.params.cbf.rate * B];                
            % Add input constraints if u_max or u_min exists.
            if isfield(obj.params, 'u_max')
                A = [A; ones(obj.udim)];
                if size(obj.params.u_max, 1) == 1
                    b = [b; obj.params.u_max * ones(obj.udim, 1)];
                elseif size(obj.params.u_max, 1) == obj.udim
                    b = [b; obj.params.u_max];
                else
                    error("params.u_max should be either a scalar value or an (udim, 1) array.")
                end
            end
            if isfield(obj.params, 'u_min')
                A = [A; -ones(obj.udim)];
                if size(obj.params.u_min, 1) == 1
                    b = [b; -obj.params.u_min * ones(obj.udim, 1)];
                elseif size(obj.params.u_min, 1) == obj.udim
                    b = [b; -obj.params.u_min];
                else
                    error("params.u_min should be either a scalar value or an (udim, 1) array")
                end
            end


            %% Cost
            if isfield(obj.params.weight, 'input')
                if size(obj.params.weight.input, 1) == 1 
                    weight_input = obj.params.weight.input * eye(obj.udim);
                elseif all(size(obj.params.weight.input) == obj.udim)
                    weight_input = obj.params.weight.input;
                else
                    error("params.weight.input should be either a scalar value or an (udim, udim) array.")
                end
            else
                weight_input = eye(obj.udim);
            end

            if verbose
                options =  optimset('Display','notify');
            else
                options =  optimset('Display','off');
            end

            % cost = 0.5 u' H u + f u
            H = weight_input;
            f_ = -weight_input * u_ref;
            [u, ~, exitflag, ~] = quadprog(H, f_, A, b, [], [], [], [], [], options);
            if exitflag == -2
                feas = 0;
                disp("Infeasible QP. CBF constraint is conflicting with input constraints.");
            else
                feas = 1;
            end
            comp_time = toc(tstart);
        end  

      
    end
end

