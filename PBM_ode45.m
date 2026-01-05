function varargout = PBM_ode45(ode,tspan,y0,options,varargin)


%% Setup

    %Solver info
    solver_name = 'PBM_ode45';

    % Stats
    nsteps  = 0;
    nfailed = 0;
    nfevals = 0;

    % Check inputs
    if nargin < 4
        options = [];
        if nargin < 3
            y0 = [];
            if nargin < 2
                tspan = [];
            end
        end
    end

    %Create ode function handle
    [ode, odeIsFuncHandle, odeTreatAsMFile] = packageAsFuncHandle(ode);

% Output
output_sol = odeIsFuncHandle && nargout == 1; % sol = odeXX(...)
output_ty  = ~output_sol && nargout > 0;      % [t,y,...] = odeXX(...)
% There might be no output requested...

sol = []; f3d = [];
if output_sol
    sol.solver = solver_name;
    sol.extdata.odefun = ode;
    sol.extdata.options = options;
    sol.extdata.varargin = varargin;
end

% Handle solver arguments
[neq, tspan, ntspan, next, t0, tfinal, tdir, y0, f0, odeArgs, ...
    options, threshold, rtol, normcontrol, normy, userhmin, hmax, htry, htspan, dataType] = ...
    odearguments(odeIsFuncHandle,odeTreatAsMFile, solver_name, ode, tspan, y0, options, varargin);
nfevals = nfevals + 1;

% Handle the output
if nargout > 0
    outputFcn = odeget(options,'OutputFcn',[]);
else
    outputFcn = odeget(options,'OutputFcn',@odeplot);
end
outputArgs = {};
if isempty(outputFcn)
    haveOutputFcn = false;
else
    haveOutputFcn = true;
    outputs = odeget(options,'OutputSel',1:neq);
    if isa(outputFcn,'function_handle')
        % With MATLAB 6 syntax pass additional input arguments to outputFcn.
        outputArgs = varargin;
    end
end
refine = max(1,odeget(options,'Refine',4));
if ntspan > 2
    outputAt = 1;          % output only at tspan points
elseif refine <= 1
    outputAt = 2;          % computed points, no refinement
else
    outputAt = 3;          % computed points, with refinement
    S = (1:refine-1) / refine;
end
printstats = strcmp(odeget(options,'Stats','off'),'on');

% Handle the event function
[haveEventFcn,eventFcn,eventArgs,valt,teout,yeout,ieout] = ...
    odeevents(odeIsFuncHandle,ode,t0,y0,options,varargin);

% Handle the mass matrix
[Mtype, M, Mfun] =  odemass(odeIsFuncHandle,ode,t0,y0,options,varargin);
if Mtype > 0  % non-trivial mass matrix
    Msingular = odeget(options,'MassSingular','no');
    if strcmp(Msingular,'maybe')
        warning(message('MATLAB:ode45:MassSingularAssumedNo'));
    elseif strcmp(Msingular,'yes')
        error(message('MATLAB:ode45:MassSingularYes'));
    end
    % Incorporate the mass matrix into ode and odeArgs.
    [ode,odeArgs] = odemassexplicit(odeIsFuncHandle,Mtype,ode,odeArgs,Mfun,M);
    f0 = ode(t0,y0,odeArgs{:});
    nfevals = nfevals + 1;
end

% Non-negative solution components
idxNonNegative = odeget(options,'NonNegative',[]);
nonNegative = ~isempty(idxNonNegative);
if nonNegative  % modify the derivative function
    [ode,thresholdNonNegative] = odenonnegative(ode,y0,threshold,idxNonNegative);
    f0 = ode(t0,y0,odeArgs{:});
    nfevals = nfevals + 1;
end

% Pull minimum step size
try
    hmin = options.MinStepSize;
catch
    userhmin = 1E-3;
end

%Initialize time and solution
t = t0;
y = y0;

% Allocate memory if we're generating output.
prototype = zeros(dataType);
nout = 0;
tout = []; yout = [];
if nargout > 0
    if output_sol
        chunk = min(max(100,50*refine), refine+floor((2^11)/neq));
        tout = zeros(1, chunk, 'like', prototype);
        yout = zeros(neq, chunk, 'like', prototype);
        f3d  = zeros(neq, 7, chunk, 'like', prototype);
    else
        if ntspan > 2                         % output only at tspan points
            tout = zeros(1, ntspan, 'like', prototype);
            yout = zeros(neq, ntspan, 'like', prototype);
        else                                  % alloc in chunks
            chunk = min(max(100,50*refine), refine+floor((2^13)/neq));
            tout = zeros(1, chunk, 'like', prototype);
            yout = zeros(neq, chunk, 'like', prototype);
        end
    end
    nout = 1;
    tout(nout) = t;
    yout(:,nout) = y;
end

% Initialize method parameters.
pow = 1/5;
A = [1/5, 3/10, 4/5, 8/9, 1, 1]; % Still used by restarting criteria
% B = [
%     1/5         3/40    44/45   19372/6561      9017/3168       35/384
%     0           9/40    -56/15  -25360/2187     -355/33         0
%     0           0       32/9    64448/6561      46732/5247      500/1113
%     0           0       0       -212/729        49/176          125/192
%     0           0       0       0               -5103/18656     -2187/6784
%     0           0       0       0               0               11/84
%     0           0       0       0               0               0
%     ];
% E = [71/57600; 0; -71/16695; 71/1920; -17253/339200; 22/525; -1/40];

% Same values as above extracted as scalars (1 and 0 values ommitted)
a2=cast(1/5, 'like', prototype);
a3=cast(3/10, 'like', prototype);
a4=cast(4/5, 'like', prototype);
a5=cast(8/9, 'like', prototype);

b11=cast(1/5, 'like', prototype);
b21=cast(3/40, 'like', prototype);
b31=cast(44/45, 'like', prototype);
b41=cast(19372/6561, 'like', prototype);
b51=cast(9017/3168, 'like', prototype);
b61=cast(35/384, 'like', prototype);
b22=cast(9/40, 'like', prototype);
b32=cast(-56/15, 'like', prototype);
b42=cast(-25360/2187, 'like', prototype);
b52=cast(-355/33, 'like', prototype);
b33=cast(32/9, 'like', prototype);
b43=cast(64448/6561, 'like', prototype);
b53=cast(46732/5247, 'like', prototype);
b63=cast(500/1113, 'like', prototype);
b44=cast(-212/729, 'like', prototype);
b54=cast(49/176, 'like', prototype);
b64=cast(125/192, 'like', prototype);
b55=cast(-5103/18656, 'like', prototype);
b65=cast(-2187/6784, 'like', prototype);
b66=cast(11/84, 'like', prototype);

e1=cast(71/57600, 'like', prototype);
e3=cast(-71/16695, 'like', prototype);
e4=cast(71/1920, 'like', prototype);
e5=cast(-17253/339200, 'like', prototype);
e6=cast(22/525, 'like', prototype);
e7=cast(-1/40, 'like', prototype);

hmin = max(16*eps(t),userhmin);
if isempty(htry)
    % Compute an initial step size h using y'(t).
    absh = min(hmax, htspan);
    if normcontrol
        rh = (norm(f0) / max(normy,threshold)) / (0.8 * rtol^pow);
    else
        rh = norm(f0 ./ max(abs(y),threshold),inf) / (0.8 * rtol^pow);
    end
    if absh * rh > 1
        absh = 1 / rh;
    end
    absh = max(absh, hmin);
else
    absh = min(hmax, max(hmin, htry));
end
f1 = f0;

% Initialize the output function.
if haveOutputFcn
    feval(outputFcn,[t tfinal],y(outputs),'init',outputArgs{:});
end

% THE MAIN LOOP

done = false;
while ~done

    % By default, hmin is a small number such that t+hmin is only slightly
    % different than t.  It might be 0 if t is 0.
    %hmin = max(16*eps(t),userhmin);
    absh = min(hmax, max(hmin, absh));    % couldn't limit absh until new hmin
    h = tdir * absh;

    %Limit minimum timestep (if option specified)
    if h < hmin
        h = hmin;
        
    end

    % Stretch the step if within 10% of tfinal-t.
    if 1.1*absh >= abs(tfinal - t)
        h = tfinal - t;
        absh = abs(h);
        done = true;
    end

    % LOOP FOR ADVANCING ONE STEP.
    nofailed = true;                      % no failed attempts
    while true
        y2 = y + h .* (b11.*f1 );
        t2 = t + h .* a2;
        f2 = ode(t2, y2, odeArgs{:});

        y3 = y + h .* (b21.*f1 + b22.*f2 );
        t3 = t + h .* a3;
        f3 = ode(t3, y3, odeArgs{:});

        y4 = y + h .* (b31.*f1 + b32.*f2 + b33.*f3 );
        t4 = t + h .* a4;
        f4 = ode(t4, y4, odeArgs{:});

        y5 = y + h .* (b41.*f1 + b42.*f2 + b43.*f3 + b44.*f4 );
        t5 = t + h .* a5;
        f5 = ode(t5, y5, odeArgs{:});

        y6 = y + h .* (b51.*f1 + b52.*f2 + b53.*f3 + b54.*f4 + b55.*f5 );
        t6 = t + h;
        f6 = ode(t6, y6, odeArgs{:});

        tnew = t + h;
        if done
            tnew = tfinal;   % Hit end point exactly.
        end
        h = tnew - t;      % Purify h.

        ynew = y + h.* ( b61.*f1 + b63.*f3 + b64.*f4 + b65.*f5 + b66.*f6 );
        f7 = ode(tnew, ynew, odeArgs{:});

        nfevals = nfevals + 6;

        % Estimate the error.
        NNrejectStep = false;
        fE = f1*e1 + f3*e3 + f4*e4 + f5*e5 + f6*e6 + f7*e7;
        if normcontrol
            normynew = norm(ynew);
            errwt = max(max(normy,normynew),threshold);
            err = absh * (norm(fE) / errwt);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ) / errwt ;
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        else
            err = absh * norm((fE) ./ max(max(abs(y),abs(ynew)),threshold),inf);
            if nonNegative && (err <= rtol) && any(ynew(idxNonNegative)<0)
                errNN = norm( max(0,-ynew(idxNonNegative)) ./ thresholdNonNegative, inf);
                if errNN > rtol
                    err = errNN;
                    NNrejectStep = true;
                end
            end
        end

        % Accept the solution only if the weighted error is no more than the
        % tolerance rtol.  Estimate an h that will yield an error of rtol on
        % the next step or the next try at taking this step, as the case may be,
        % and use 0.8 of this value to avoid failures.
        if err > rtol % Failed step

            %Enforce minimum step size, regardless of error
            if h < hmin + hmin*1E-3
                break
            end

            nfailed = nfailed + 1;
            % if absh <= hmin
            %     warning(message('MATLAB:ode45:IntegrationTolNotMet', sprintf( '%e', t ), sprintf( '%e', hmin )));
            %     solver_output = odefinalize(solver_name, sol,...
            %         outputFcn, outputArgs,...
            %         printstats, [nsteps, nfailed, nfevals],...
            %         nout, tout, yout,...
            %         haveEventFcn, teout, yeout, ieout,...
            %         {f3d,idxNonNegative}, ...
            %         t);
            %     if nargout > 0
            %         varargout = solver_output;
            %     end
            %     return;
            % end

            if nofailed
                nofailed = false;
                if NNrejectStep
                    absh = max(hmin, 0.5*absh);
                else
                    absh = max(hmin, absh * max(0.1, 0.8*(rtol/err)^pow));
                end
            else
                absh = max(hmin, 0.5 * absh);
            end
            h = tdir * absh;
            done = false;

        else                                % Successful step

            NNreset_f7 = false;
            if nonNegative && any(ynew(idxNonNegative)<0)
                ynew(idxNonNegative) = max(ynew(idxNonNegative),0);
                if normcontrol
                    normynew = norm(ynew);
                end
                NNreset_f7 = true;
            end

            break;

        end
    end
    nsteps = nsteps + 1;

    if haveEventFcn
        f = [f1 f2 f3 f4 f5 f6 f7];
        [te,ye,ie,valt,stop] = ...
            odezero(@ntrp45,eventFcn,eventArgs,valt,t,y,tnew,ynew,t0,h,f,idxNonNegative);
        if ~isempty(te)
            if output_sol || (nargout > 2)
                teout = [teout, te]; %#ok<AGROW>
                yeout = [yeout, ye]; %#ok<AGROW>
                ieout = [ieout, ie]; %#ok<AGROW>
            end
            if stop                 % Stop on a terminal event.
                % Adjust the interpolation data to [t te(end)].

                % Update the derivatives using the interpolating polynomial.
                taux = t + (te(end) - t)*A;
                [~,f(:,2:7)] = ntrp45(taux,t,y,[],[],h,f,idxNonNegative);
                f2 = f(:,2); f3 = f(:,3); f4 = f(:,4); f5 = f(:,5); f6 = f(:,6); f7 = f(:,7);

                tnew = te(end);
                ynew = ye(:,end);
                h = tnew - t;
                done = true;
            end
        end
    end

    if output_sol
        nout = nout + 1;
        if nout > length(tout)
            tout = [tout, zeros(1, chunk, 'like', prototype)];  %#ok<AGROW> % requires chunk >= refine
            yout = [yout, zeros(neq, chunk, 'like', prototype)]; %#ok<AGROW>
            f3d  = cat(3,f3d,zeros(neq, 7, chunk, 'like', prototype));
        end
        tout(nout) = tnew; %#ok<AGROW>
        yout(:,nout) = ynew; %#ok<AGROW>
        f3d(:,:,nout) = [f1 f2 f3 f4 f5 f6 f7]; %#ok<AGROW>
    end

    if output_ty || haveOutputFcn
        switch outputAt
        case 2      % computed points, no refinement
            nout_new = 1;
            tout_new = tnew;
            yout_new = ynew;
        case 3      % computed points, with refinement
            tref = t + (tnew-t)*S;
            nout_new = refine;
            tout_new = [tref, tnew];
            yntrp45 = ntrp45split(tref,t,y,h,f1,f3,f4,f5,f6,f7,idxNonNegative);
            yout_new = [yntrp45, ynew];
        case 1      % output only at tspan points
            nout_new =  0;
            tout_new = [];
            yout_new = [];
            while next <= ntspan
                if tdir * (tnew - tspan(next)) < 0
                    if haveEventFcn && stop     % output tstop,ystop
                        nout_new = nout_new + 1;
                        tout_new = [tout_new, tnew]; %#ok<AGROW>
                        yout_new = [yout_new, ynew]; %#ok<AGROW>
                    end
                    break;
                end
                nout_new = nout_new + 1;
                tout_new = [tout_new, tspan(next)]; %#ok<AGROW>
                if tspan(next) == tnew
                    yout_new = [yout_new, ynew]; %#ok<AGROW>
                else
                    yntrp45 = ntrp45split(tspan(next),t,y,h,f1,f3,f4,f5,f6,f7,idxNonNegative);
                    yout_new = [yout_new, yntrp45]; %#ok<AGROW>
                end
                next = next + 1;
            end
        end

        if nout_new > 0
            if output_ty
                oldnout = nout;
                nout = nout + nout_new;
                if nout > length(tout)
                    tout = [tout, zeros(1, chunk, 'like', tout)]; %#ok<AGROW> requires chunk >= refine
                    yout = [yout, zeros(neq, chunk, 'like', yout)]; %#ok<AGROW>
                end
                idx = nout;
                if nout_new ~= 1
                    idx = oldnout+1:nout;
                end
                tout(idx) = tout_new; %#ok<AGROW>
                yout(:,idx) = yout_new; %#ok<AGROW>
            end
            if haveOutputFcn
                stop = feval(outputFcn,tout_new,yout_new(outputs,:),'',outputArgs{:});
                if stop
                    done = true;
                end
            end
        end
    end

    if done
        break
    end

    % If there were no failures compute a new h.
    if nofailed
        % Note that absh may shrink by 0.8, and that err may be 0.
        temp = 1.25*(err/rtol)^pow;
        if temp > 0.2
            absh = absh / temp;
        else
            absh = 5.0*absh;
        end
    end

    % Advance the integration one step.
    t = tnew;
    y = ynew;
    if normcontrol
        normy = normynew;
    end
    if NNreset_f7
        % Used f7 for unperturbed solution to interpolate.
        % Now reset f7 to move along constraint.
        f7 = ode(tnew, ynew, odeArgs{:});
        nfevals = nfevals + 1;
    end
    f1 = f7;  % Already have f(tnew,ynew)

end % end of while loop

interp_data = {f3d,idxNonNegative};
solver_output = odefinalize(solver_name, sol,...
    outputFcn, outputArgs,...
    printstats, [nsteps, nfailed, nfevals],...
    nout, tout, yout,...
    haveEventFcn, teout, yeout, ieout, interp_data);

if nargout > 0
    varargout = solver_output;
end



end