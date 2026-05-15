%Numerical Jacobian calculation for ODE systems
%
% Computes the Jacobian matrix J of the ODE system dydt = f(t,y) at the 
% initial condition (t0, y0) using numerical finite differences.
%
% The Jacobian J is an (neq x neq) matrix where J(i,j) = df_i/dy_j
%
% USAGE:
%   [jac, jac_sparsity] = CalcJacobian(ode, tspan, y0, params)
%
% INPUTS:
%   ode      - ODE function handle: dydt = ode(t, y, params)
%   tspan    - Time span [t_initial, t_final, ...]
%   y0       - Initial conditions (column vector)
%   params   - Parameters structure passed to ODE function
%
% OUTPUTS:
%   jac      - Jacobian matrix (neq x neq)
%   jac_sparsity - Sparsity pattern (logical matrix)

function [jac, jac_sparsity] = CalcJacobian(ode, tspan, y0, params)

    % Input validation
    if nargin < 4
        error('CalcJacobian requires at least 4 inputs: ode, tspan, y0, params');
    end
    
    % Ensure y0 is a column vector
    y0 = y0(:);
    neq = length(y0);
    
    % Initial time
    t0 = tspan(1);
    
    % Finite difference parameters
    % Use a scaled step size based on solution magnitude
    delta_rel = 1e-7;  % Relative perturbation
    delta_abs = 1e-10; % Absolute perturbation (for y_i near zero)
    
    % Initialize Jacobian matrix
    jac = zeros(neq, neq);
    
    % Compute baseline ODE evaluation f(t0, y0)
    f0 = ode(t0, y0);
    f0 = f0(:);  % Ensure column vector
    
    if length(f0) ~= neq
        error('ODE function output dimension (%d) does not match y0 dimension (%d)', ...
              length(f0), neq);
    end
    
    % Compute Jacobian using forward finite differences
    % J(:,j) = [f(t0, y0 + delta*e_j) - f(t0, y0)] / delta
    
    parfor j = 1:neq  % Use parfor for parallel computation
        % Determine perturbation size for j-th component
        delta = max(delta_abs, delta_rel * abs(y0(j)));
        
        % Create perturbed state vector
        y_perturbed = y0;
        y_perturbed(j) = y0(j) + delta;
        
        % Evaluate ODE at perturbed state
        f_pert = ode(t0, y_perturbed);
        f_pert = f_pert(:);  % Ensure column vector
        
        % Compute j-th column of Jacobian using central differences
        % For better accuracy, use: (f(y+delta) - f(y-delta)) / (2*delta)
        y_perturbed_minus = y0;
        y_perturbed_minus(j) = y0(j) - delta;
        f_pert_minus = ode(t0, y_perturbed_minus);
        f_pert_minus = f_pert_minus(:);
        
        % Central difference approximation (more accurate)
        jac(:, j) = (f_pert - f_pert_minus) / (2 * delta);
    end
    
    % Compute sparsity pattern (identify negligible Jacobian entries)
    % Entries with magnitude < 1e-12 * max(|J|) are considered zero
    jac_max = max(abs(jac(:)));
    if jac_max > 0
        tol_sparsity = 1e-12 * jac_max;
    else
        tol_sparsity = 1e-12;
    end
    
    jac_sparsity = abs(jac) > tol_sparsity;
    
    % Optionally convert to sparse if sparsity is high
    sparsity_fraction = 1 - nnz(jac_sparsity) / (neq * neq);
    if sparsity_fraction > 0.9  % If > 90% sparse
        % Could return sparse format here if needed
        % jac = sparse(jac);
    end
    
end