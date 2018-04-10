function [bd,f,retcode,a,b] = calc_bound_rect_disj(dat,mod,al,au,bl,bu,fl,fu,sense,time_limit)

try

    %% Save general tolerance for convenience

    tol = dat.tol_gen;

    %% Update the model RHS's based on input to this subroutine

    mod.rhs(mod.a_lb) = -al; % Note negative!
    mod.rhs(mod.a_ub) =  au;
    mod.rhs(mod.b_lb) = -bl; % Note negative!
    mod.rhs(mod.b_ub) =  bu;
    mod.rhs(mod.f_lb) = -fl; % Note negative!
    mod.rhs(mod.f_ub) =  fu;

    %% Set Gurobi's options

    settings = set_solver_options;

    %% If the time_limit argument is passed in, update the solver time limit

    if nargin >= 10
        settings.gurobi.TimeLimit = time_limit;
    end

    %% If the upper bound case, then...

    if strcmp(sense,'upper')

        %% Set the model to be a maximization

        mod.modelsense = 'max';

        %% Set Gurobi's MIP focus to 3; see
        %% http://www.gurobi.com/documentation/5.6/refman/mipfocus.

        settings.gurobi.MIPFocus = 3;

        %% Solve the problem, save the diagnostics

        diags = gurobi(mod,settings.gurobi);

        %% Save the objbound if it exists...

        if isfield(diags, 'objbound')

            bd = diags.objbound;

        %% Otherwise set upper bound to +Inf

        else

            bd = +Inf;

        end

        %% If the problem came back as infeasible or unbounded or if
        %% the time limit was hit, then set the values f and retcode
        %% accordingly

        if strcmp(diags.status,'INF_OR_UNBD') || strcmp(diags.status,'TIME_LIMIT')

            f = diags.x(mod.f); % Changed by Veronica; previous line was f = [];
            retcode = -1;

        %% Otherwise, the problem should have been solved to optimality.
        %% Set f and retcode accordingly

        elseif isfield(diags,'objval')

            f = diags.x(mod.f);
            retcode = 0;

        %% Otherwise, we hit some kind of bad error. Just quit

        else
            error('Problem with upper bound');
        end

    %% Else, if the lower bound case, then...

    else

        %% Set the model to be a maximization

        mod.modelsense = 'min';

        %% Set Gurobi's MIP focus to 3; see
        %% http://www.gurobi.com/documentation/5.6/refman/mipfocus.

        settings.gurobi.MIPFocus = 3;

        %% Solve the problem, save the diagnostics

        diags = gurobi(mod,settings.gurobi);

        %% We believe that we should *not* get 'INF_OR_UNBD' as a
        %% return code (why?), but if it does occur, it is probably due
        %% to numerical 'tightness' of the lower and upper bounds on
        %% the f variable. We therefore keep resolving with loosened
        %% constraints until we do not get 'INF_OR_UNBD'

        while strcmp(diags.status,'INF_OR_UNBD')

            mod.rhs(mod.f_lb) = -(fl-tol);
            mod.rhs(mod.f_ub) =   fu+tol;
            diags = gurobi(mod,settings.gurobi);
            tol = 10*tol;

        end

        %% Save the objbound if it exists...

        if isfield(diags, 'objbound')

            bd = diags.objbound;

        %% Otherwise set upper bound to -Inf

        else

            bd = -Inf;

        end

        %% If the problem comes back as infeasible or unbounded (it
        %% shouldn't!) or if the time limit is hit, then set the return
        %% values accordingly

        if strcmp(diags.status,'INF_OR_UNBD') || strcmp(diags.status,'TIME_LIMIT')

            f = [];
            a = [];
            b = [];
            retcode = -1;

        %% Otherwise, the problem should have been solved to optimality.
        %% Set the return values accordingly

        elseif isfield(diags,'objval')

            f = diags.x(mod.f);
            a = diags.x(mod.a);
            b = diags.x(mod.b);
            retcode = 0;

        %% Otherwise, we hit some kind of bad error. Just quit

        else

            error('Problem with lower bound');

        end

    end

    %% Print the diagnostics to the screen if the problem was not
    %% solved optimally. This is for visual debugging and probably
    %% could be removed

    if ~strcmp(diags.status,'OPTIMAL')
        diags
    end

catch me

    save calc_bound_rect_disj_error
    error('Examine the file calc_bound_rect_disj_error.mat');

end

return
