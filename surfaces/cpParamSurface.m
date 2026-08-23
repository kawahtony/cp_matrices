function [cpx,cpy,cpz, dist, bdy, uu, vv] = cpParamSurface(xx,yy,zz, paramf, paramf2nd, paramfEdge, optfin, surfmesh, LB, UB, paramAdjust, How, DEBUG)
%CPPARAMSURFACE   CP representation of parameterised surfaces via optimization
%   Currently for open surfaces with a single edge (e.g., mobius strip)
%   Probably you don't want this directly, see e.g., cpMobiusStrip
%
%   [cpx,cpy,cpz dist, varargout] = cpParamSurface(p,q,r,
%                         paramf, paramf2nd, edgeparamf.
%                         optf, surfmesh, LB, UB, paramAdjust,
%                         How, Debug)
%    p,q,r: points to find closest points of.
%
%    paramf: a function of the form [x,xu,xv] = paramf(u,v), which
%            evaluates the surface x(u,v) and xu and xv are the
%            partial w.r.t. u and v.
%
%    paramf2nd: similar to above but provides [xuu,xvv,xuv], the
%               2nd partials.
%
%    paramfEdge: a function which parameterizes an edge of an open
%                surface.  TODO: what should we do for a closed
%                surface?  []?
%                Currently there can only be one of these (the surface
%                can have a single edge, e.g., hemisphere or mobius
%                strip).
%
%    optf: override the detault optimization function.  pass [] to
%          use the default.
%
%    surfmesh: a mesh of the surface used to get a good initial guess,
%              stored in a cell array as {xp,yp,zp,up,vp}.  Should
%              have a good coverage of the surface, may 500 points?
%
%    LB, UB: lower/upper bounds for parameters.  If you have a
%            periodic parameter, see the next option.  TODO
%
%    Optional inputs:
%    paramAdjust: a function that is called after each optimization
%                 step (e.g., to wrap a periodic parameter back into
%                 its fundamental domain).  CAREFUL WITH THIS, DON'T
%                 DO ARBITRARY THINGS WITH IT.
%
%    How: what technique to use, various implementations, see source.
%         How=2 (Newton) is the default workhorse; see helper_newton
%         below.
%
%    DEBUG: verbosity for How=2 (Newton).  0 (default): silent except
%           for a single summary warning if any points failed to
%           converge to "tol".  1: also prints one line per point that
%           needed the fallback (didn't converge, or converged to a
%           worse point than the initial mesh guess).  >=2: verbose
%           per-iteration trace (see helper_newton).  >=10: also plots
%           (into figure 1) and pauses after every iteration -- only
%           useful interactively on a handful of points.
%
% TODO: this code is in a state of flux...  BEWARE
%
% 2026 update: helper_newton was rewritten to use Levenberg-Marquardt
% damping (instead of an ad hoc "gradient descent with a fixed 0.05
% step" fallback whenever the plain Newton step wasn't a descent
% direction) and to fall back to the initial mesh guess if it somehow
% still ends up worse.  This makes it converge far more often --
% important because for a whole Cartesian block of query points (not
% just a narrow band near the surface), many points end up near the
% surface's medial axis / cut locus (e.g., the interior "hole" of a
% torus-like shape) where the closest point is genuinely
% ill-conditioned (near-zero gradient, indefinite Hessian, sometimes
% more than one equally-close point) and plain Newton could fail to
% converge or even wander to a worse point.  Also: non-convergence no
% longer calls `keyboard` (which would otherwise hang a batch script
% waiting at a debug prompt for every failing point) -- it's tracked
% and reported once, as a single summary warning, after all points are
% processed.  The old debug-only comparison block at the end of
% helper_newton referenced several undefined functions (xs, ys, g, gp,
% d2) and would error out with "undefined function" the moment a
% point's Newton result was worse than its initial guess -- likely the
% actual crash/hang behind "Newton's failure to converge" -- that
% block has been replaced by the (correct, harmless) fallback above.

  if (nargin < 9) || isempty(LB)
    LB = [-inf -inf];
  end
  if (nargin < 10) || isempty(UB)
    UB = [inf inf];
  end
  if (nargin < 11)
    paramAdjust = [];
  end
  if (nargin < 12) || isempty(How)
    How = 0;
  end
  if (nargin < 13) || isempty(DEBUG)
    DEBUG = 0;
  end

  if isempty(optfin)
    if How == 0
      optf = @my_fmincon_optf;
    elseif How == 1
      optf = @my_lsq_optf;
    else
      % unneeded?
    end
  else
    optf = optfin;
  end

  % seems like less accurate in the x direction is normal
  opt = optimset('tolfun', 1e-14, 'tolx', 1e-12, 'tolcon', 1e-13, ...
                 'maxfunevals', 100000, 'maxiter', 100000);
  opt = optimset(opt, 'Display', 'off');
  opt = optimset(opt, 'largescale', 'off');
  opt = optimset(opt, 'algorithm', 'active-set');
  opt = optimset(opt, 'GradObj','on');
  fmincon_opt = opt;

  % seems like less accurate in the x direction is normal
  opt = optimset('tolfun', 1e-14, 'tolx', 1e-12, 'tolcon', 1e-13, ...
                 'maxfunevals', 100000, 'maxiter', 100000);
  opt = optimset(opt, 'Display', 'off');
  opt = optimset(opt, 'Jacobian', 'on');
  lsq_opt = opt;



  %% loop over the points and find the closest point for each
  x1d = xx(:); y1d = yy(:); z1d = zz(:);
  nx = length(x1d);     % number of points

  % allocate space
  cpx = zeros(nx,1);
  cpy = zeros(nx,1);
  cpz = zeros(nx,1);
  dist = zeros(nx,1);
  bdy = zeros(nx,1);
  uu = zeros(nx,1);
  vv = zeros(nx,1);

  [xp, yp, zp, up, vp] = surfmesh{:};

  nFailed = 0;   % How==2 only: count of points that needed the fallback

  fprintf('cpParamSurface: starting to process %d points\n', nx);
  for i = 1:nx
    xpt = [x1d(i); y1d(i); z1d(i)];
    p = xpt(1);
    q = xpt(2);
    r = xpt(3);

    %% initial guess
    time_guess = cputime();
    % vectorized:
    dds = (xp - xpt(1)).^2 + (yp - xpt(2)).^2 + (zp - xpt(3)).^2;
    [mindd_guess,I] = min(dds(:));
    s_initial_guess = [up(I); vp(I)];
    %[s_initial_guess, mindd_guess] = helper_initialguess(xpt, surfmesh);
    time_guess = cputime() - time_guess;

    out = [];

    if (How == 0)
      opt_time = cputime();
      [cp, dist1, bdy1, s] = helper_fmincon(xpt, optf, s_initial_guess, paramf, fmincon_opt, LB, UB, paramAdjust);
      u1 = s(1);  v1 = s(2);
      out = [out  cputime() - opt_time];
    end


    if (How == 1)
      % in my tests this was much slower
      opt_time = cputime();
      [cp, dist1, bdy1, s] = helper_lsq(xpt, optf, s_initial_guess, paramf, lsq_opt, LB, UB, paramAdjust);
      u1 = s(1);  v1 = s(2);
      out = [out  cputime() - opt_time];
    end


    if (How == 2)
      % fastest by an order of magnitude, maybe less reliable, lots
      % of parameters to tune (!)
      opt_time = cputime();
      [cp, dist1, bdy1, s, converged1] = helper_newton(xpt, mindd_guess, s_initial_guess, paramf, paramf2nd, LB, UB, paramAdjust, DEBUG);
      if (~converged1)
        nFailed = nFailed + 1;
      end
      u1 = s(1);  v1 = s(2);
      if (bdy1 == 1)
        [cpx2,cpy2,cpz2,dist2,s2] = cpParam3DCurveClosed(xpt(1), xpt(2), xpt(3), paramfEdge, [0 4*pi]);
        cp = [cpx2; cpy2; cpz2];
        dist1 = dist2;
        u1 = s2;
        v1 = 1;   % TODO: hardcoded for mobius
      end
      out = [out  cputime() - opt_time];
    end

    %disp([time_guess out])

    cpx(i) = cp(1);
    cpy(i) = cp(2);
    cpz(i) = cp(3);
    dist(i) = dist1;
    bdy(i) = bdy1;
    uu(i) = u1;
    vv(i) = v1;
  end

  if (How == 2) && (nFailed > 0)
    warning('cpParamSurface:newtonNotConverged', ...
            ['%d of %d point(s) did not converge to the Newton tolerance ' ...
             '(the best point found was used instead). This is common ' ...
             'for points far from the surface, e.g. near a medial ' ...
             'axis/interior hole, and is usually harmless if such points ' ...
             'end up outside your band. Pass DEBUG>=1 to cpParamSurface ' ...
             'for a per-point report.'], nFailed, nx);
  end

  cpx = reshape(cpx, size(xx));
  cpy = reshape(cpy, size(xx));
  cpz = reshape(cpz, size(xx));
  dist = reshape(dist, size(xx));
  bdy = reshape(bdy, size(xx));
  uu = reshape(uu, size(xx));
  vv = reshape(vv, size(xx));

  %% helper functions placed inside: access local variables
  % careful: variables inside this aren't local either (!)
  function [d2, grad] = my_fmincon_optf(uv, p)
    % todo, nargout here and below, faster when no grad? (probably not)
    [tx, txu, txv] = paramf(uv(1), uv(2));

    d2 = sum( (tx - p).^2 );
    %d2 = (x(1) - p(1))^2 + ...
    %     (x(2) - p(2))^2 + ...
    %     (x(3) - p(3))^2;

    grad = [ 2*(tx - p)' * (txu);
             2*(tx - p)' * (txv) ];

    % TODO:
    %grad = [ 2*(x(1)-p(1)) * xu(1) + ...
    %         2*(x(2)-p(2)) * xu(2) + ...
    %         2*(x(3)-p(3)) * xu(3); ...
    %         2*(x(1)-p(1)) * xv(1) + ...
    %         2*(x(2)-p(2)) * xv(2) + ...
    %         2*(x(3)-p(3)) * xv(3) ];
  end

  function [F, J] = my_lsq_optf(uv, p)
    [tx, txu, txv] = paramf(uv(1), uv(2));
    F = tx - p;
    %if nargout > 1   % Two output arguments
    J = [txu(1)  txv(1); txu(2)  txv(2); txu(3)  txv(3)];
    %end
  end

end  % end main function



%% helper functions


function [cp, dist, bdy, s, converged] = helper_newton(xpt, mindd_guess, s_guess, paramf, paramf2nd, LB, UB, paramAdjust, DEBUG)
% Newton's method, with Levenberg-Marquardt damping for robustness.
%
% This doesn't deal with boundaries---just puts a penalty to try to
% stop it from converging too far outside.  bdy will be set to 1 if
% it finishes on or outside the boundary.  You could then do a
% search on the boundary curve.
%
% Robustness notes (see also the header of cpParamSurface.m):
%  - A plain Newton step J\(-f) is only a descent direction of the
%    squared distance d^2 when J (the Hessian of d^2) is positive
%    definite.  Off the "reach" of the surface (e.g. beyond its
%    medial axis/cut locus, which any point far enough from the
%    surface can be) J need not be positive definite, so instead of a
%    single fixed-size gradient-descent fallback, we damp J by adding
%    lambda*I (increasing lambda until the resulting step actually
%    decreases d^2).  This is the standard Levenberg-Marquardt trick
%    and is much more reliably a descent method than plain Newton.
%  - Convergence is judged from the size of the (undamped-by-later-
%    wrapping) step itself, not from comparing the pre- and
%    post-paramAdjust parameter values -- otherwise a periodic
%    paramAdjust (e.g. wrapping u into [0,2*pi)) could make a
%    converged step look huge just because it crossed the periodic
%    seam.
%  - If, despite all this, the final point is somehow worse than the
%    cheap initial mesh-based guess, we just return that initial
%    guess instead: the caller is never worse off than a naive nearest
%    -mesh-vertex search.
%  - Failing to converge in maxn iterations no longer throws up a
%    `keyboard` prompt (which would hang unattended/batch runs); it's
%    reported back via the `converged` output instead, and
%    cpParamSurface.m prints one summary warning covering all points.

  if (nargin < 9) || isempty(DEBUG)
    DEBUG = 0;
  end

  tol = 1e-10;
  maxn = 80;

  s = s_guess(:);

  [x0] = paramf(s(1), s(2));
  d2cur = sum((x0-xpt).^2);

  lambda = 0;
  lambda_growth = 10;

  converged = false;
  n = 0;
  while (n < maxn)
    n = n + 1;

    [x, xu, xv] = paramf(s(1), s(2));
    [xuu, xvv, xuv] = paramf2nd(s(1), s(2));

    f = [ 2*(x-xpt)' * xu; ...
          2*(x-xpt)' * xv ];

    f1u = 2*(x-xpt)' * xuu + 2*xu' * xu;
    f2v = 2*(x-xpt)' * xvv + 2*xv' * xv;
    f1v = 2*(x-xpt)' * xuv + 2*xu' * xv;

    % add a penalty to distance for outside the boundaries
    % TODO: lots of choice here and somewhat hardcoded for mobius
    % strip (e.g., boundary is only in v).
    gamma = 100;
    pow = 3;
    % (x(u,v) - xpt)^2 + gamma*(v-1)^3
    if (s(2) > UB(2))
      f(2) = f(2) + gamma*(s(2)-UB(2))^(pow-1);
      f2v = f2v + (pow-1)*gamma*(s(2)-UB(2))^(pow-2);
    elseif (s(2) < LB(2))
      f(2) = f(2) - gamma*(s(2)-LB(2))^(pow-1);
      f2v = f2v - (pow-1)*gamma*(s(2)-LB(2))^(pow-2);
    end

    J = [f1u f1v; f1v f2v];
    Jscale = max(1, norm(J, 'fro'));

    % Levenberg-Marquardt: find the smallest damping (starting from
    % last iteration's lambda, so well-behaved regions stay cheap)
    % that gives a genuine decrease in d^2.
    change = [];
    for lmtry = 1:40
      Jd = J + lambda*Jscale*eye(2);
      dJ = Jd(1,1)*Jd(2,2) - Jd(1,2)*Jd(2,1);
      if (abs(dJ) > eps*Jscale^2*1e4)
        trial = Jd \ (-f);
      else
        trial = [];
      end
      if ~isempty(trial) && all(isfinite(trial))
        if norm(trial) > 1.0
          trial = trial * (1.0/norm(trial));  % trust-region-ish cap
        end
        strial = s + trial;
        if (~isempty(paramAdjust))
          strial_eval = paramAdjust(strial);
        else
          strial_eval = strial;
        end
        xtrial = paramf(strial_eval(1), strial_eval(2));
        d2trial = sum((xtrial-xpt).^2);
        if (d2trial <= d2cur + 1e-13*max(1,d2cur))
          change = trial;
          d2cur = d2trial;
          break;
        end
      end
      if (lambda == 0)
        lambda = 1e-8;
      else
        lambda = lambda * lambda_growth;
      end
    end

    if (DEBUG >= 2)
      fprintf('iter: n=%d, s=[%g,%g], lambda=%g, x=[%f,%f,%f]\n', n, s(1), s(2), lambda, xpt(1), xpt(2), xpt(3));
    end

    if isempty(change)
      % Levenberg-Marquardt couldn't find any improving step (should
      % be extremely rare): stop here rather than looping uselessly.
      break;
    end

    % relax the damping a bit for next time, since this one worked
    lambda = lambda / lambda_growth^2;
    if (lambda < 1e-12)
      lambda = 0;
    end

    snew = s + change;
    if (~isempty(paramAdjust))
      snew = paramAdjust(snew);
    end

    if (DEBUG >= 10)
      cp_dbg = paramf(snew(1), snew(2));
      fprintf('iter: n=%d, snew=(%f,%f), s0=(%f,%f), x=(%f,%f,%f)\n', n, snew(1), snew(2), s_guess(1), s_guess(2), xpt(1), xpt(2), xpt(3));
      set(0, 'CurrentFigure', 1);
      plot3(cp_dbg(1), cp_dbg(2), cp_dbg(3), 'bo');
      axis equal
      drawnow();
      pause
    end

    % convergence test uses the actual step taken (norm(change)), not
    % abs(s-snew): after a periodic paramAdjust wraps snew, s and snew
    % can look far apart even for a tiny true step.
    if (norm(change) < tol)
      s = snew;
      converged = true;
      break;
    end

    s = snew;
  end

  if (~converged && DEBUG >= 1)
    fprintf('cpParamSurface: Newton did not converge in %d iters, s=(%f,%f), x=(%f,%f,%f), d2=%g\n', ...
            maxn, s(1), s(2), xpt(1), xpt(2), xpt(3), d2cur);
  end

  cp = paramf(s(1), s(2));
  dist = norm(xpt - cp, 2);

  % Safety net: never do worse than the cheap initial mesh guess (this
  % also covers the case where the loop above bailed out immediately,
  % e.g. lmtry exhausted on iteration 1).
  if (dist^2 > mindd_guess + tol)
    if (DEBUG >= 1)
      fprintf('cpParamSurface: Newton result worse than initial guess (%g > %g), falling back, x=(%f,%f,%f)\n', ...
              dist^2, mindd_guess, xpt(1), xpt(2), xpt(3));
    end
    s = s_guess(:);
    cp = paramf(s(1), s(2));
    dist = sqrt(mindd_guess);
    converged = false;
  end

  % somewhat hardcoded for mobius (no s1 here)
  if ( (s(2) < LB(2)) || (s(2) > UB(2)) )
    bdy = true;
    % need to search the boundary (we haven't found it properly yet)
  else
    bdy = false;
  end

end  % newton function




function [cp, dist, bdy, res] = helper_lsq(xpt, f, initial_guess, paramf, opt, LB, UB, paramAdjust)

  t0 = initial_guess;

  %[res, fval, fvals, flg, output,lambda]...
  %    = lsqnonlin(d2, t0, [-inf; -1], [inf; 1], opt);
  [res, fval, fvals, flg, output,lambda]...
      = lsqnonlin(@(s) f(s,xpt), t0, LB, UB, opt);

  bdy = any([lambda.lower; lambda.upper] ~= 0);
  % TODO: lambda more reliable?
  %if abs(abs(res(2)) - 1) < 1e-12
  %  bdy = true;
  %else
  %  bdy = false;
  %end
  %[bdy lambda.lower' lambda.upper']

  if (flg < 1)
    warning('cpParamSurface:lsqNotConverged', 'lsq search: possibly nonconverged CP search (x=[%g,%g,%g])', xpt(1), xpt(2), xpt(3));
  end

  if (~isempty(paramAdjust))
    newres = paramAdjust(res);
  else
    newres = res;
  end

  cp = paramf(newres(1),newres(2));
  dd = fval;
  % for lsqnonlin, need squared dist
  if (abs(dd - sum((cp - xpt).^2)) > 1e-15)
    warning('cpParamSurface:lsqDistMismatch', 'actual distance doesn''t match opt result (x=[%g,%g,%g])', xpt(1), xpt(2), xpt(3));
  end

  res = newres;
  dist = sqrt(dd);
end  % helper_lsq function




function [cp, dist, bdy, res] = helper_fmincon(xpt, f, initial_guess, paramf, opt, LB, UB, paramAdjust)
  MakePlots = 0;

  if (MakePlots)
    fn = 99;
    %figure(fn); clf;
    set(0, 'CurrentFigure', fn); clf;
    [x,y,z] = paramMobiusStrip(32,Rad,Thick);
    surf(x,y,z);
    xlabel('x'); ylabel('y'); zlabel('z');
    hold on;
    axis equal
  end

  if (MakePlots)
    cp_guess = [xp(I) yp(I) zp(I)];
    plot3([xpt(1) cp_guess(1)], [xpt(2) cp_guess(2)], ...
          [xpt(3) cp_guess(3)], 'ro--')
  end

  t0 = initial_guess;

  %[res, fval, flg, output,lambda,jacobian]...
  %    = fmincon(f, t0, ...
  %              [],[], [],[], [-inf; -1], [inf; 1], [],  opt);
  [res, fval, flg, output,lambda,jacobian]...
      = fmincon(@(s) f(s,xpt), t0, ...
                [],[], [],[], LB, UB, [],  opt);

  bdy = any([lambda.lower; lambda.upper] ~= 0);

  if (flg < 1)
    warning('cpParamSurface:fminconNotConverged', 'possibly nonconverged CP search (x=[%g,%g,%g])', xpt(1), xpt(2), xpt(3));
  end

  if (~isempty(paramAdjust))
    newres = paramAdjust(res);
  else
    newres = res;
  end
  cp = paramf(newres(1),newres(2));
  dd = fval;
  if (abs(dd - sum((cp - xpt).^2)) > 2e-15)
    warning('cpParamSurface:fminconDistMismatch', 'actual distance doesn''t match opt result (x=[%g,%g,%g])', xpt(1), xpt(2), xpt(3));
  end

  dist = sqrt(dd);
  res = newres;

  if (MakePlots)
    plot3([xpt(1) cp(1)], [xpt(2) cp(2)], [xpt(3) cp(3)], 'bo-')
    drawnow()
  end
end  % helper_fmincon function
