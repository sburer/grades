function [gap_int,a_gap,b_gap] = get_maxgap_ubV(a_l,a_u,b_l,b_u,f_nom)

  % This function computes an upper bound on the maxgap (wrt to f_nom)
  % over a rectangle defined by the intervals [a_l,a_u] and [b_l,b_u].

  % First compute those entries of f_nom that overlap [a_l,a_u] and,
  % separately, [b_l,b_u].

  a_f_nom = f_nom( find( a_l <= f_nom & f_nom <= a_u ) );
  b_f_nom = f_nom( find( b_l <= f_nom & f_nom <= b_u ) );

  % Also add the closest f_nom entries (from below the lower limits) to
  % the lists. After 'union', the lists should be sorted. I believe the
  % sort order is ascending.

  % Can this fail if there are no entries of f_nom that are below?

  a_f_nom = [a_f_nom; max(f_nom(find(f_nom < a_l))) ];
  b_f_nom = [b_f_nom; max(f_nom(find(f_nom < b_l))) ];

  % Also add the upper entries a_u and b_u to the lists. After 'union',
  % the lists should be sorted. I belive the sort order is ascending.

  a_f_nom = union(a_f_nom,a_u);
  b_f_nom = union(b_f_nom,b_u);

  % Now compute the widths/gaps of the intervals given by the lists.

  a_gaps = a_f_nom(2:end) - a_f_nom(1:end-1);
  b_gaps = b_f_nom(2:end) - b_f_nom(1:end-1);

  % Now we calculate the gaps wrt a and b separately.

  [a_gap_ub, ind_a] = max(a_gaps);
  a_gap = a_f_nom(ind_a+1);
  [b_gap_ub, ind_b ] = max(b_gaps);
  % b_gap = b_f_nom(ind_a+1); % Sam: old
  b_gap = b_f_nom(ind_b+1); % Sam: new

  % Our final upper bound is the max (min?) of the two upper bounds just
  % calculated.

  %gap_ub = max(a_gap_ub,b_gap_ub);
  %now we have the min since it is the exact gap 
   gap_int = min(a_gap_ub,b_gap_ub);

return
