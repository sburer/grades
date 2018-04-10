function [gub, a_opt, b_opt, gap_opt, feedback] = main(dat)

% Dictionary:
%
%   UB     = upper bound
%   LB     = lower bound
%   PBL    = pseudoborderline (associated with fixed a, b, f)
%   MaxPBL = maximum/worst-case PBL (associated with fixed a, b and varying f)
%
% Commenting conventions:
%
%   Comment starting with double '%%' is an action
%   Comment starting with single '%' is some other type of comment

%%%%%%%%%%%%%%%%%%%%%
%% Start CPU clock %%
%%%%%%%%%%%%%%%%%%%%%

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize B&B structure and variables %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Scalar for global UB on MaxPBL (as we are minimizing MaxPBL). An UB
% corresponds to the evaluation of MaxPBL at a fixed a, b

gub = +Inf;
gub_verified = +Inf;

% Scalar for global LB on gap (since the gap is our secondary objective,
% which we are maximizing)

gap_opt = -Inf;

% Vector of LB's for MaxPBL, one for each node

LBs = -Inf;

% Vector for UB's for gap, one for each node

GAP_UBs = +Inf;

% List of rectangles/feasible sets, one for each node. A node is defined
% by its rectangle, but we try to use the term "node" instead of
% "rectangle" to make clear the role in the B&B algorithm

INTs = [dat.l_a, dat.u_a, dat.l_b, dat.u_b];

% Initial number of nodes in the B&B tree

nnodes = size(INTs,1);

% Initialization for number of times we have encountered a Gurobi
% problem while calculating the MaxPBL UB

times_upper_prob = 0;

% Initialization for number of times we encountered a Gurobi problem
% while calculating the MaxPBL LB

times_lower_prob = 0;

% Initialization for number of times we hit the limit for tries to find
% a verified upper bound in a node/rectangle

times_hit_limit = 0;

% Structure for saving good feasible solutions, f, in the uncertainty
% set encountered during the B&B algorithm

saved_f    = []; % A matrix with each column a feasible f
saved_f_ct = []; % A vector with how many times corresponding column f has
                 % proven valuable in improving the LB during an iteration

% Structure for saving optimal solutions

a_opt = [];
b_opt = [];

% Structure for LB's of nodes that are fathomed due to their rectangle
% being too small

LBs_too_small = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform other initializations %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Save data from dat structure locally in order to reduce typing and
%% make code easier to read

tol_bor = dat.tol_bor; % Tolerance for borderlines
tol_int = dat.tol_int; % Interval tolerance
tol_gap = dat.tol_gap; % Maxgap tolerance

%% Build Gurobi model from YALMIP model. Quite complicated, but it
%% should save overhead because we can call Gurobi repeatedly faster
%% than we can call YALMIP

fprintf('%s: Gurobi model...', mfilename);

mod = build_model_rect_disj(dat);
mod = extract_gurobi_model(dat, mod);

fprintf('done\n');

%% Initialize iteration counter to 1

iter = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Perform the main B&B algorithm loop %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

while size(INTs,1) > 0

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Select next node in tree to examine %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Re-order remaining nodes first by ascending LB, then by
    %% descending diagonal width of rectangle associated with node. In
    %% particular, we do *not* sort by the gap in any way

    clear diag
    for index = 1:length(LBs)
        al = INTs(index,1);
        au = INTs(index,2);
        bl = INTs(index,3);
        bu = INTs(index,4);
        diag(index) = sqrt((au - al)^2 + (bu - bl)^2)/2;
    end
    diag = diag';
    tmp = [LBs,diag];
    [~,perm] = sortrows(tmp,[1,-2]);
    INTs = INTs(perm,:);
    LBs  = LBs(perm);
    GAP_UBs = GAP_UBs(perm,:);

    %% Pull node off tree. Update the remaining tree

    int = INTs(1,:);
    lb = LBs(1);

    INTs = INTs(2:end,:);
    LBs = LBs(2:end);
    GAP_UBs = GAP_UBs(2:end);

    %% Save info about the current node, including associated limits and
    %% midpoints of rectangle

    al = int(1); au = int(2);
    bl = int(3); bu = int(4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate UB for minimum MaxPBL on the node %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('%s (iter=%3d): Upp bd...', mfilename, iter);

    % The algorithm's primary focus is the minimization of MaxPBL.
    % Here, we need an UB for the MaxPBL on the entire rectangle
    % associated with the current node (which also serves as a global
    % UB). Specifically, we calculate MaxPBL for the midpoint of the
    % rectangle, which serves as the UB

    % S: Have we described this correctly? Maybe the comment should say:
    %
    % The algorithm's primary focus is the minimization of MaxPBL. Here,
    % we want to compute an UB for the minimum MaxPBL using the current
    % rectangle/node. Specifically, we calculate MaxPBL for the midpoint
    % of the rectangle, which serves as the UB

    local_count = 1; retcode = -1;

    while local_count <= 25 && retcode == -1

        if local_count == 1
            am = 0.5*(al + au);
            bm = 0.5*(bl + bu);
        else
            am = al + rand*(au - al);
            bm = bl + rand*(bu - bl);
        end

        [tmp_obj_val, ~, retcode] = calc_bound_rect_disj(dat, mod, am, ...
            am, bm, bm, dat.fL, dat.fU, 'upper');

        %% Increment counter if a problem with the MIP solver was detected

        if retcode == -1
            times_upper_prob = times_upper_prob + 1;
        end

        local_count = local_count + 1;

    end

    if local_count > 25 && retcode == -1
        times_hit_limit = times_hit_limit + 1;
    end

    %% Update gub and list of optimal solutions based on the UB just
    %% calculated, but only if retcode ~= -1 (that is, only if we have a
    %% verified upper bound)

    if retcode ~= -1

        [a_opt, b_opt, gub, gub_verified, gap_opt] = ...
            update_optimal(tmp_obj_val, gub, gub_verified, retcode, tol_bor, ...
            a_opt, b_opt, am, bm, dat, gap_opt, tol_gap);

    end

    fprintf('done\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate gap UB for the node/rectangle %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % S: Veronica says, this gap is 100% accurate for the given
    % node/rectangle, i.e., it's not an UB. And remember that the gap is
    % always wrt the nominal f (i.e., we don't consider uncertainty for
    % the gap)

    [gap_ub, gap_a, gap_b] = get_maxgap_ubV(al, au, bl, bu, dat.fL);

    % If we are lucky, the pair (gap_a, gap_b) just calculated will also
    % have a low value of MaxPBL, allowing us to update the gub and
    % optimal solution. So we calculate MaxPBL and try to update things

    % S: The f produced here ("another_f") is apparently the only
    % f that gets saved to the structure saved_f. Is this right?
    %
    % V: My recall was that we added more f but maybe I am
    % wrong...
    %
    % S&V: We don't remember. Maybe we *could* save the other f from
    % calling the UB a few steps above. However, either way this should
    % not affect the correctness of the code

    [tmp_obj_val, another_f, retcode] = calc_bound_rect_disj(dat, mod, ...
      gap_a, gap_a, gap_b, gap_b, dat.fL, dat.fU, 'upper');

    %% Increment counter if a problem with the MIP solver was detected

    if retcode == -1
        times_upper_prob = times_upper_prob + 1
    end

    % V: this is an interesting pair [a,b] of this rectangle, so it
    % is worthwhile to be investigated. If ret_code is -1, we could
    % check the obj_bound, and if it is interesting with respect to the
    % upper bound on the rectangle, we could maybe increase the time
    % limit? If it is higher than the upper bound computed above we can
    % just discard the pair and avoid to call update_optimal?

    % S: I don't fully understand. Can we talk about this?

    if retcode ~= -1 

        [a_opt, b_opt, gub, gub_verified, gap_opt] = ...
            update_optimal(tmp_obj_val, gub, gub_verified, retcode, tol_bor, ...
            a_opt, b_opt, gap_a, gap_b, dat, gap_opt, tol_gap);

    end

    %% Save the f just found. We only save the f (not the G) because G
    %% does not impact the number of PBL's directly

    % S: I do not fully understand the previous statement anymore

    saved_f = [saved_f, another_f];
    saved_f_ct = [saved_f_ct, 1];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate the MaxPBL LB %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate the LB of MaxPBL (saved in value lb). We plug in
    %% different values of f (sorted in descending order wrt to the
    %% number of times they have been useful) and track whether a
    %% particular f is helpful for improving the LB

    % Value of lb has been initialized above

    fprintf('%s (iter=%3d): Low bd...', mfilename, iter);

    %% Set indicator for fathoming

    fathom_indicator = 0;

    %% Loop over the saved f in reverse order

    for k = size(saved_f, 2) : -1 : 1

        %% Save the current f as tmp_f, and calculate LB based on tmp_f

        tmp_f = saved_f(:,k);
        [next_f_val, ~, retcode] = calc_bound_rect_disj(dat, mod, al, au, ...
            bl, bu, tmp_f, tmp_f, 'lower');

        %% Increment counter if a problem with the MIP solver was detected

        if retcode == -1
            times_lower_prob = times_lower_prob + 1
        end

        %% Calculate the relative error compared to the current value of lb

        relerr = calc_rel_error(next_f_val, lb);

        %% If relerr is small, then we have found the same LB and we...

        if relerr < tol_bor

            %% Average/smooth out the lb

            lb = 0.5*(lb + next_f_val); % V: Use min here instead?

        %% Otherwise, if the new bound is better than lb, then we...

        elseif next_f_val > lb

            %% Update lb, and save tmp_f as being helpful

            lb = next_f_val;
            if k < size(saved_f, 2) % Not necessary if k == size(saved_f, 2)
                saved_f_ct(k) = saved_f_ct(k) + 1;
            end

            % If the lb > gub, then break the for loop because we can
            % fathom this node immediately, i.e., we are done with this
            % node
           
            % S: Are we certain we don't later accidentally branch on
            % this node?
            %
            % V: I am not sure I understood the question
            %
            % S: I guess my question is, does the 'break' statement take
            % us to the next iteration in the B&B algorithm, or does it
            % take us to the branching below?

            relerr = calc_rel_error(lb, gub); % V: Possibly use gub_verified instead
            if relerr >= tol_bor && lb > gub
                fathom_indicator = 1;
                break
            end

        %% Otherwise, new LB was not equal or better. So we...

        else

            %% Update tmp_f as not being helpful (as long as it wasn't
            %% *just* created during this iteration) and loop

            if k < size(saved_f, 2)
                saved_f_ct(k) = saved_f_ct(k) - 1;
            end
            
        end

    end

    %% For housekeeping, reorder saved_f based on helpfulness, and then
    %% throw away those that have not proved helpful

    [saved_f_ct,perm] = sort(saved_f_ct, 'ascend');
    saved_f = saved_f(:,perm);

    ind = find(saved_f_ct > -10);
    saved_f = saved_f(:,ind);
    saved_f_ct = saved_f_ct(ind);

    fprintf('done\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Branch (as appropriate) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf('%s (iter=%3d): Branching...  ', mfilename, iter);

    %% Calculate the a_width and b_width of the current rectangle

    a_width = (au - al)/2;
    b_width = (bu - bl)/2;

    % Based on the code above, if lb > gub for this node, we would
    % have already fathomed the current node and continued to the next
    % iteration, i.e., these next lines should not be executed. Should
    % we do some kind of sanity check here?

    %% If lb <= GUB and at least one width is at least tol_int, then we...
    
    if fathom_indicator == 0
        
        if a_width >= tol_int || b_width >= tol_int

            %% Branch on the larger width 

            if a_width >= b_width
                mid = (au + al)/2;
                INTs = [ INTs ; al,mid,bl,bu ; mid,au,bl,bu ];
            else
                mid = (bu + bl)/2;
                INTs = [ INTs ; al,au,bl,mid ; al,au,mid,bu ];
            end
            LBs = [ LBs ; lb ; lb ];
            GAP_UBs = [ GAP_UBs ; gap_ub ; gap_ub ];
            nnodes = nnodes + 2;

        else

            LBs_too_small = [LBs_too_small; lb];

        end

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Fathom nodes (as appropriate) %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Initialize list of nodes to keep

    keep = [];

    %% Loop through remaining nodes...

    for i = 1:length(LBs)

        %% Calculate relerr between gub and LBs(i)

        relerr = calc_rel_error(gub, LBs(i)); % V: Possibly use gub_verified

        %% If LBs(i) is significantly lower than gub, then we...

        if relerr >= tol_bor && LBs(i) <= gub % V: Possibly use gub_verified

            %% Keep this node because its rectangle might have a
            %% better/smaller number of MaxPBL

            keep = [keep; i];

        %% Else if LBs(i) is sufficiently close to gub, then we...

        elseif relerr < tol_bor

            %% Examine the secondary objective (i.e., the gap) and keep
            %% the node only if it could have a better/larger gap

            relerr_gap = calc_rel_error(gap_opt, GAP_UBs(i));
            if GAP_UBs(i) >= gap_opt & relerr_gap >= tol_gap
                keep = [keep; i];
            end

        end

    end

    % The result is that we keep only nodes that could have: (i) a
    % better MaxPBL; or (ii) MaxPBL = gub and a better gap

    INTs = INTs(keep,:);
    LBs = LBs(keep,:);
    GAP_UBs = GAP_UBs(keep,:);

    fprintf('done\n');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% End iteration and prepare for next %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Calculate global LB

    if length(LBs_too_small) == 0
        glob_lb = min(LBs);
    else
        glob_lb = min(min(LBs), min(LBs_too_small));
    end

    %% Print diagnostic information

    if size(INTs,1) > 0
      fprintf('%s (iter=%3d): Nodes left = %d, (glob_lb <= mean_lb <= glob_ub) = (%.2f <= %.2f <= %.8f) gap = %.4f gub_verified = %.8f\n', ...
        mfilename, iter, size(INTs,1), glob_lb, mean(LBs), gub, gap_opt, gub_verified);
    end

    %% Increment iteration counter

    iter = iter + 1;

end

%% For the first optimal pair in the list (a_opt, b_opt), recompute
%% MaxPBL just to verify before returning. This allows a total time
%% limit of one hour

[tmp_obj_val, ~, retcode] = calc_bound_rect_disj(dat, mod, a_opt(1), ...
    a_opt(1), b_opt(1), b_opt(1), dat.fL, dat.fU, 'upper', 3600);

%% Process the results and print to screen

if retcode == -1

    relerr_final = +Inf;

    fprintf('Final MaxPBL NOT verified (retcode = -1)\n');

else

    relerr_final = calc_rel_error(tmp_obj_val, gub);

    if relerr_final < tol_bor
        fprintf('Final MaxPBL verified and correct (relerr with gub = %.4f)\n', relerr_final);
    else
        fprintf('Final MaxPBL verified but NOT correct (relerr with gub = %.4f)\n', relerr_final);
    end

end

%% For the first optimal pair in the list (a_opt, b_opt), recompute the
%% gap just to verify before returning. This code has been adopted from
%% verify_grades.m where a small tolerance was subtracted, but cannot
%% remember the exact reason (probably just some numerical-accuracy
%% thing)

a_local = a_opt(1) - dat.tol_gen/10;
b_local = b_opt(1) - dat.tol_gen/10;

[gap_local, ~, ~] = get_maxgap_ubV(a_local, a_local, b_local, b_local, ...
    dat.f_nom);

%% Process the results and print to screen

abserr = abs(gap_opt - gap_local);
if abserr < tol_bor
    fprintf('Final gap correct (abserr with gap_opt = %.4f)\n', abserr);
else
    fprintf('Final gap incorrect (abserr with gap_opt = %.4f)\n', abserr);
end

%%%%%%%%%%%%%%%%%%%%%%%
%% End the CPU clock %%
%%%%%%%%%%%%%%%%%%%%%%%

time = toc;

%% Create the feedback structure, which keeps statistics and errors.
%% In particular, it contains the time, the number of iterations, the
%% total number of nodes created, and an error code if CPLEX stalls or
%% similar

feedback = struct('time', 0, 'iter', iter, 'nodes', nnodes, 'error', 0, ...
    'times_upper_prob', times_upper_prob,...
    'times_lower_prob', times_lower_prob,...
    'times_hit_limit', times_hit_limit,...
    'gub', +Inf, ...
    'gub_verified', +Inf); 

%% Fill-in feedback to return

feedback.time = time;
feedback.iter = iter;
feedback.nodes = nnodes;
feedback.times_upper_prob = times_upper_prob;
feedback.times_lower_prob = times_lower_prob;
feedback.times_hit_limit = times_hit_limit;
feedback.gub = gub;
feedback.gub_verified = gub_verified;

feedback.error = 0; % Should be set to 1 if something goes wrong. Not
                    % implemented yet as of 2015-04-30 (but see next
                    % 2016-08-22)

%% Set the error message equal to 1 (if final MaxPBL value was not
%% verified) or -1 (if final MaxPBL value was verified but not correct)

if relerr_final >= +Inf
    feedback.error = 1;
elseif relerr_final > tol_bor
    feedback.error = -1;
end

%% Set the error_gap message equal to -1 (if final gap value was
%% incorrect)

feedback.error_gap = 0;
if abserr >= tol_bor
    feedback.error_gap = -1;
end

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  [a_opt, b_opt, gub, gub_verified, gap_opt] = ...
    update_optimal(ub_cur, gub, gub_verified, retcode, tol_bor, ...
    a_opt, b_opt, a_cur, b_cur, dat, gap_opt, tol_gap)

% Update gub_verified

% S: Is this how we want to update it?

if retcode ~= -1
    gub_verified = min(gub_verified, ub_cur);
end

%% Calculate relative error between GUB and current UB

relerr = calc_rel_error(ub_cur, gub);

%% If the current UB improves the GUB, i.e., it is less than the GUB by
%% more than tol_bor, then...

if relerr >= tol_bor && ub_cur < gub

    %% Update information about the set of optimal solutions because we
    %% have found a better number of PBL's

    %% Specifically, update the optimal solution (a singleton at this
    %% point in time)

    a_opt = a_cur; 
    b_opt = b_cur;

    %% And update the GUB

    gub = ub_cur;

    %% And calculate an UB on the maxgap at the optimal solution

    gap_opt = get_maxgap_ubV(a_opt, a_opt, b_opt, b_opt, dat.fL);

% Otherwise, if the current UB equals the GUB, i.e., it has relative
% error within tol_bor, then...

elseif relerr < tol_bor

    %% Investigate this solution further because we have found a
    %% solution with the same number of PBL's. Normally, we would just
    %% append this to our (current) list of optimal solutions. However,
    %% because we are also interested in the maxgap, we first check if
    %% we also found a better gap

    %% Calculate the current gap

    gap_cur = get_maxgap_ubV(a_cur, a_cur, b_cur, b_cur, dat.fL);

    %% Calculate the relative error with the current gap

    relerr_gap = calc_rel_error(gap_cur, gap_opt);

    %% If we have found a better gap, then...

    if gap_cur > gap_opt & relerr_gap >= tol_gap

        %% Update optimal solution set
      
        a_opt = a_cur;
        b_opt = b_cur;
        gub = ub_cur;
        gap_opt = gap_cur;

    %% Otherwise, we have found the same gap...

    elseif relerr_gap < tol_gap

        %% Append to (on-going) list of optimal solutions
      
        a_opt = [a_opt; a_cur];
        b_opt = [b_opt; b_cur];
        gub = max(gub, ub_cur); % To be conservative, use the worse of the two
        gap_opt = 0.5*(gap_opt + gap_cur); % Why the average exactly??? V: Use min instead?

    end

end

return
