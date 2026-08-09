function TVSA_make_LHS_13D_dataset(opts)
% TVSA_make_LHS_13D_dataset.m
%
% What it does (one file):
%   (A) Build 120k input design X using:
%       - 11D LHS (N11 = 30k, maximin)
%       - attach K=4 feasible starting states (qCO2_start, qH2O_start)
%         sampled inside boundary(alpha=0.85) from periodic results.
%   (B) (Optional) Run expensive part:
%       - one-cycle-forward evaluation for each X row
%       - outputs Y = [qCO2_postAds, qCO2_postDes, qH2O_postAds, qH2O_postDes, E_cycle_MWh, t_cycle_s]
%   (C) Save .mat and .csv
%
% Usage (minimal):
%   opts = struct();
%   TVSA_make_LHS_13D_dataset(opts);
%
% Units (consistent with your dataset definition):
%   Inputs:
%     T_ambC (C), RH_in (-,0-1), u_feed (m/s), p_vac (bar), T_reg (C),
%     t_heat,t_des,t_ads (s), width,height (m), N_cell (-),
%     qCO2_start,qH2O_start (mol/kg_sorb)
%   Outputs:
%     qCO2_postAds,qCO2_postDes,qH2O_postAds,qH2O_postDes (mol/kg_sorb),
%     E_cycle (MWh/cycle/collector), t_cycle (s/cycle)

% ---------------- defaults ----------------
if nargin<1 || isempty(opts), opts = struct(); end

opts = def(opts, 'out_dir', fullfile(pwd,'tvsa_lhs13d_out'));
opts = def(opts, 'core_name', 'TVSA_CorePhysics_v2');   % file must be on path
opts = def(opts, 'periodic_mat', fullfile(fileparts(mfilename('fullpath')), 'tvsa_endpoints_ext2048_rand5000.mat'));

% Design size controls
opts = def(opts, 'N11', 30000);      % 11D LHS points
opts = def(opts, 'K_state', 4);      % attach K feasible states => total N = N11*K
opts = def(opts, 'lhs_iters', 25);   % lhsdesign maximin iterations
opts = def(opts, 'seed', 11);

% Feasible domain (boundary in 2D)
opts = def(opts, 'boundary_alpha', 0.85);
opts = def(opts, 'state_reject_maxTries', 200000); % for rejection sampling

% Runner toggles
opts = def(opts, 'do_run', true);       % set false if you only want X
opts = def(opts, 'use_parfor', true);
opts = def(opts, 'chunk_size', 5000);   % chunk saving for safety
opts = def(opts, 'force_rerun', false);

% Core options for "one-cycle forward via periodic(max_iters=1)"
opts = def(opts, 'dt_s', 3.0);
opts = def(opts, 'tol_rel', 1e-12);     % irrelevant when max_iters=1, but keep small
opts = def(opts, 'max_iters', 1);       % critical: one iteration only
opts = def(opts, 'omega', 1.0);         % critical: no mixing with old state

% Parameter ranges (11D)
rng(opts.seed);

ranges = get_ranges_11D(); % struct with min/max

if ~exist(opts.out_dir,'dir'), mkdir(opts.out_dir); end

% File names
tag = sprintf('N11_%dk_K_%d_N_%dk', round(opts.N11/1000), opts.K_state, round(opts.N11*opts.K_state/1000));
mat_X = fullfile(opts.out_dir, ['X_13D_',tag,'.mat']);
csv_X = fullfile(opts.out_dir, ['X_13D_',tag,'.csv']);
mat_XY= fullfile(opts.out_dir, ['XY_13D_',tag,'.mat']);
csv_XY= fullfile(opts.out_dir, ['XY_13D_',tag,'.csv']);

% ---------------- (A) Make feasible state sampler from periodic results ----------------
[state_poly, state_bounds, meta_periodic] = build_feasible_state_domain(opts.periodic_mat, opts.boundary_alpha);

% ---------------- (B) Build X (120k) ----------------
if exist(mat_X,'file')==2 && ~opts.force_rerun
    fprintf('[X] Exists, loading: %s\n', mat_X);
    S = load(mat_X,'X','X_names','opts_saved');
    X = S.X; X_names = S.X_names;
else
    fprintf('[X] Generating 11D LHS (%d) and attaching K=%d states...\n', opts.N11, opts.K_state);

    X11 = lhs_11D(opts.N11, ranges, opts.lhs_iters);

    % sample feasible states: total needed = N11*K
    Ntot = opts.N11 * opts.K_state;
    states = sample_states_in_polygon(Ntot, state_poly, state_bounds, opts.state_reject_maxTries);

    % replicate 11D rows K times and append states
    X = nan(Ntot, 13);

    % expand 11D
    for k = 1:opts.K_state
        idx = (k-1)*opts.N11 + (1:opts.N11);
        X(idx,1:11) = X11;
    end

    % append states
    X(:,12) = states(:,1); % qCO2_start
    X(:,13) = states(:,2); % qH2O_start

    X_names = { ...
        'T_ambC_C','RH_in_frac','u_feed_m_s','p_vac_bar','T_reg_C', ...
        't_heat_s','t_des_s','t_ads_s','width_m','height_m','N_cell', ...
        'qCO2_start_molkg','qH2O_start_molkg'};

    opts_saved = opts; %#ok<NASGU>
    save(mat_X,'X','X_names','opts_saved','ranges','state_poly','state_bounds','meta_periodic','-v7.3');
    writetable(array2table(X,'VariableNames', matlab.lang.makeValidName(X_names)), csv_X);
    fprintf('[X] Saved: %s and %s\n', mat_X, csv_X);
end

% If you only want inputs, stop here
if ~opts.do_run
    fprintf('[DONE] Inputs generated only (opts.do_run=false).\n');
    return;
end

% ---------------- (C) Expensive runner: X -> Y ----------------
if exist(mat_XY,'file')==2 && ~opts.force_rerun
    fprintf('[XY] Exists, skipping run: %s\n', mat_XY);
    return;
end

% core handle
if exist(opts.core_name,'file')~=2
    error('Core not found on path: %s.m', opts.core_name);
end
TVSA = @(varargin) feval(opts.core_name, varargin{:});

% templates (ctrl/design/opt) compatible with core
[amb0, ctrl0, design0, opt0] = templates_core(opts);

% finalize base design once (we will override width/height/N_cell each case)
design0 = TVSA('finalize_design', design0);
phys0   = TVSA('make_phys', design0);

% quick smoke test on 1 row
fprintf('[RUN] Smoke test...\n');
[y_smoke, ok_smoke, msg_smoke] = run_one_case(TVSA, X(1,:), amb0, ctrl0, design0, phys0, opt0);
if ~ok_smoke
    disp(msg_smoke);
    error('Smoke test failed. Fix this before running 120k.');
end
fprintf('[RUN] Smoke test OK.\n');

% parallel pool
if opts.use_parfor
    p = gcp('nocreate');
    if isempty(p), parpool; end
end

N = size(X,1);
Y = nan(N,6);
ok = false(N,1);
tsec = nan(N,1);
err = strings(N,1);

Y_names = { ...
    'qCO2_postAds_molkg','qCO2_postDes_molkg', ...
    'qH2O_postAds_molkg','qH2O_postDes_molkg', ...
    'E_cycle_MWh_cycle_collector','t_cycle_s'};

% chunk loop (safer than one giant run)
nChunk = ceil(N/opts.chunk_size);
fprintf('[RUN] Total N=%d in %d chunks (chunk_size=%d)\n', N, nChunk, opts.chunk_size);

for c = 1:nChunk
    i1 = (c-1)*opts.chunk_size + 1;
    i2 = min(c*opts.chunk_size, N);
    idx = i1:i2;

    fprintf('[RUN] Chunk %d/%d | rows %d..%d\n', c, nChunk, i1, i2);

    if opts.use_parfor
        err_cell = cell(numel(idx),1);
        Yc = nan(numel(idx),6);
        okc= false(numel(idx),1);
        tsc= nan(numel(idx),1);

        parfor j = 1:numel(idx)
            ii = idx(j);
            t0 = tic;
            [y_i, ok_i, msg_i] = run_one_case(TVSA, X(ii,:), amb0, ctrl0, design0, phys0, opt0);
            tsc(j) = toc(t0);

            Yc(j,:) = y_i;
            okc(j)  = ok_i;
            err_cell{j} = msg_i;
        end

        Y(idx,:) = Yc;
        ok(idx)  = okc;
        tsec(idx)= tsc;
        for j=1:numel(idx)
            if isempty(err_cell{j}), err(idx(j))=""; else, err(idx(j))=string(err_cell{j}); end
        end

    else
        for ii = idx
            t0 = tic;
            [y_i, ok_i, msg_i] = run_one_case(TVSA, X(ii,:), amb0, ctrl0, design0, phys0, opt0);
            tsec(ii)= toc(t0);
            Y(ii,:) = y_i;
            ok(ii)  = ok_i;
            if isempty(msg_i), err(ii)=""; else, err(ii)=string(msg_i); end
        end
    end

    % incremental save
    meta = struct();
    meta.X_names = X_names;
    meta.Y_names = Y_names;
    meta.opts    = opts;
    meta.ranges  = ranges;
    meta.state_poly = state_poly;
    meta.state_bounds = state_bounds;
    meta.meta_periodic = meta_periodic;

    save(mat_XY,'X','Y','ok','tsec','err','meta','-v7.3');
end

% Write CSV (final)
T = array2table([X Y], 'VariableNames', matlab.lang.makeValidName([X_names Y_names]));
writetable(T, csv_XY);

fprintf('[DONE] Saved:\n  %s\n  %s\n', mat_XY, csv_XY);

end

%% ============================= core runner per case =============================
function [y, okflag, emsg] = run_one_case(TVSA, Xi, amb0, ctrl0, design0, phys0, opt0)
okflag = false; emsg = '';
y = nan(1,6);

try
    % decode inputs
    amb = amb0; ctrl = ctrl0; design = design0; phys = phys0; opt = opt0;

    amb.T_ambC = Xi(1);
    amb.RH_in  = Xi(2);

    ctrl.u_feed_m_s    = Xi(3);
    ctrl.p_vac_bar     = Xi(4);
    ctrl.T_regenC      = Xi(5);
    ctrl.time_heat_s   = Xi(6);
    ctrl.time_des_s    = Xi(7);
    ctrl.time_ads_s    = Xi(8);

    design.width_cell_m  = Xi(9);
    design.height_cell_m = Xi(10);
    design.N_cell        = round(Xi(11));

    qCO2_tot0 = Xi(12);
    qH2O_tot0 = Xi(13);

    % finalize design/phys for this geometry
    design = TVSA('finalize_design', design);
    phys   = TVSA('make_phys', design);

    % driver
    drv = TVSA('get_driver', amb, ctrl, opt, opt.driver_baseline);

    % build x_init_override (6x1): [qCO2f;qCO2s;qH2Of;qH2Ob;yCO2g;TKbed]
    stage0 = uint8(drv.stage_id(1));

    x0 = build_x0_from_totals(qCO2_tot0, qH2O_tot0, amb, drv, opt, stage0);

    % ONE-CYCLE FORWARD via periodic with max_iters=1 (core has override for x_init)
    sim = TVSA('run_point_periodic', drv, amb, ctrl, design, phys, opt, x0);

    % outputs: CO2 endpoints
    if isfield(sim,'cyc')
        y(1) = sim.cyc.qCO2_ads_end;
        y(2) = sim.cyc.qCO2_des_end;
        ncap = sim.cyc.nCO2_captured_eff; % mol/cycle/collector
    else
        ncap = NaN;
    end

    % outputs: H2O endpoints using stage id indices (Ads=6 end, Des=3 end)
    if isfield(sim,'drv') && isfield(sim.drv,'stage_id') && isfield(sim,'prof') && isfield(sim.prof,'qH2O_total')
        sid = sim.drv.stage_id;
        idx_ads_end = find(sid==6,1,'last'); if isempty(idx_ads_end), idx_ads_end = numel(sid); end
        idx_des_end = find(sid==3,1,'last'); if isempty(idx_des_end), idx_des_end = 1; end
        qH2O = sim.prof.qH2O_total;
        y(3) = qH2O(idx_ads_end);
        y(4) = qH2O(idx_des_end);
    end

    % t_cycle (s/cycle)
    if isfield(sim,'drv') && isfield(sim.drv,'t_s')
        tcycle = sim.drv.t_s(end);
    else
        % fallback from ctrl
        tcycle = ctrl.time_vac_s + ctrl.time_heat_s + ctrl.time_des_s + ctrl.time_cool_s + ctrl.time_press_s + ctrl.time_ads_s;
    end
    y(6) = tcycle;

    % E_cycle (MWh/cycle/collector)
    % core provides intensity: ppi.Elec_total_MWh_tCO2 (MWh/tCO2)
    if isfield(sim,'ppi') && isfield(sim.ppi,'Elec_total_MWh_tCO2') && isfinite(ncap)
        EI = sim.ppi.Elec_total_MWh_tCO2; % MWh/tCO2
        tCO2_cycle = ncap * 44.01e-6;     % mol * (kg/mol) -> kg; /1000 => t. 44.01e-6 t/mol
        y(5) = EI * tCO2_cycle;           % MWh/cycle/collector
    end

    okflag = all(isfinite(y(1:4))); % endpoints must exist

catch ME
    okflag = false;
    emsg = getReport(ME,'basic');
end
end

%% ============================= build x0 from q totals =============================
function x0 = build_x0_from_totals(qCO2_tot0, qH2O_tot0, amb, drv, opt, stage0)
% x0 = [qCO2f;qCO2s;qH2Of;qH2Ob;yCO2g;TKbed]
qCO2_tot0 = max(qCO2_tot0, 0);
qH2O_tot0 = max(qH2O_tot0, 0);

% stage-based pool fractions (k_aw typically 0 in your defaults, so f=f0)
if opt.h2o_2pool.enable
    if stage0==6
        fb0 = opt.h2o_2pool.fbound.f0_ads;
    else
        fb0 = opt.h2o_2pool.fbound.f0_regen;
    end
    fb0 = clamp(fb0, opt.h2o_2pool.fbound.f_min, opt.h2o_2pool.fbound.f_max);
    qH2Ob0 = fb0 * qH2O_tot0;
    qH2Of0 = max(qH2O_tot0 - qH2Ob0, 0);
else
    qH2Of0 = qH2O_tot0;
    qH2Ob0 = 0;
end

if opt.co2_2pool.enable
    if stage0==6
        fs0 = opt.co2_2pool.fslow.f0_ads;
    else
        fs0 = opt.co2_2pool.fslow.f0_regen;
    end
    fs0 = clamp(fs0, opt.co2_2pool.fslow.f_min, opt.co2_2pool.fslow.f_max);
    qCO2s0 = fs0 * qCO2_tot0;
    qCO2f0 = max(qCO2_tot0 - qCO2s0, 0);
else
    qCO2f0 = qCO2_tot0;
    qCO2s0 = 0;
end

yCO2g0 = clamp(amb.yCO2_air, opt.yco2.y_floor, opt.yco2.y_max);
TKbed0 = drv.TK_cmd(1);

x0 = [qCO2f0; qCO2s0; qH2Of0; qH2Ob0; yCO2g0; TKbed0];
end

%% ============================= templates =============================
function [amb0, ctrl0, design0, opt0] = templates_core(opts)
amb0 = struct('altitude_m',0,'T_ambC',20,'RH_in',0.50,'yCO2_air',400e-6);

ctrl0 = struct();
ctrl0.u_feed_m_s   = 0.028;
ctrl0.p_vac_bar    = 0.05;
ctrl0.T_regenC     = 100;
ctrl0.time_vac_s   = 60;
ctrl0.time_heat_s  = 2400;
ctrl0.time_des_s   = 20000;
ctrl0.time_cool_s  = 400;
ctrl0.time_press_s = 60;
ctrl0.time_ads_s   = 8000;

design0 = struct();
design0.width_cell_m   = 1.43;
design0.height_cell_m  = 0.10;
design0.N_cell         = 13*88;

design0.msorb_kg       = 1002.8;
design0.eps            = 0.56;
design0.dp_m           = 7.50e-4;

design0.L_flow_fixed_m = 0.0172;
design0.A_mult         = 2.10;
design0.Ng_ref         = 35;
design0.dP_misc_Pa     = 0;
design0.N_series_pass  = 1;
design0.flow_mode      = 'Aflow';

% opt defaults (subset needed by core)
opt0 = struct();
opt0.mode      = 'S2';
opt0.dt_s      = opts.dt_s;
opt0.max_iters = opts.max_iters;
opt0.tol_rel   = opts.tol_rel;
opt0.omega     = opts.omega;

% ptraj/aw blocks (required by core safety)
opt0.ptraj = struct('enable',false,'use_locked',false,'locked_file','', 'dt_s', opts.dt_s);
opt0.aw    = struct('model','simple','eps',1e-9,'max',0.999999);

% yCO2 guards
opt0.yco2 = struct('y_floor',1e-12,'y_max',0.50);

% 2-pool configs (k_aw set to 0 so pool fractions constant stage-wise)
opt0.co2_2pool = struct();
opt0.co2_2pool.enable = true;
opt0.co2_2pool.fslow  = struct('f0_ads',0.25,'f0_des',0.55,'f0_regen',0.55,'k_aw_ads',0.0,'k_aw_reg',0.0,'f_min',0.0,'f_max',1.0);
opt0.co2_2pool.kslow  = struct('ratio',0.10,'kmin',2e-6,'kmax',1e-3,'aw_pow',0.0,'enable',true);
opt0.co2_2pool.resid  = struct('enable',true,'mode','ads_scaled','apply_in_stage_ids',uint8([2 3]),'f_res',0.08,'q_abs_min',0.05,'tol_frac',1e-4);

opt0.h2o_2pool = struct();
opt0.h2o_2pool.enable = true;
opt0.h2o_2pool.fbound = struct('f0',0.25,'f0_ads',0.20,'f0_des',0.35,'f0_regen',0.35,'k_aw_ads',0.0,'k_aw_reg',0.0,'f_min',0.0,'f_max',1.0);
opt0.h2o_2pool.kslow  = struct('ratio',0.15,'kmin',5e-7,'kmax',5e-3,'aw_pow',0.0,'enable',true);
opt0.h2o_2pool.resid  = struct('enable',true,'apply_in_stage_ids',uint8([1 2 3 4 5]),'f0',0.10,'f_max',0.55,'k_T',0.10,'k_time',0.10,'Tref_K',373.15,'tref_s',20000,'q_abs_min',0.05,'tol_frac',1e-4);

% driver baseline struct (forced-only if locked file missing)
locked_name = 'locked_drv_baseline.mat';
locked_file = '';
cands = {fullfile(pwd,locked_name), fullfile(fileparts(mfilename('fullpath')),locked_name)};
for i=1:numel(cands)
    if exist(cands{i},'file')==2, locked_file=cands{i}; break; end
end
if ~isempty(locked_file)
    opt0.driver_baseline = struct('mode','auto','locked_file',locked_file);
else
    opt0.driver_baseline = struct('mode','forced_only','locked_file','');
end

end

%% ============================= 11D ranges =============================
function ranges = get_ranges_11D()
ranges = struct();

ranges.T_ambC  = [5, 40];
ranges.RH_in   = [0.10, 0.90];
ranges.u_feed  = [0.02, 0.04];
ranges.p_vac   = [0.02, 0.10];
ranges.T_reg   = [80, 120];

ranges.t_heat  = [1200, 3000];
ranges.t_des   = [10000, 20000];
ranges.t_ads   = [4000, 10000];

ranges.width   = [1.001, 1.859];
ranges.height  = [0.07, 0.13];
ranges.N_cell  = [801, 1487];

end

%% ============================= LHS 11D =============================
function X11 = lhs_11D(N, ranges, lhs_iters)
% 11D: [T,RH,u,pvac,Treg,theat,tdes,tads,width,height,Ncell]
D = 11;
U = simple_lhs(N, D, 11);  % seed can be opts.seed if you pass it in

X11 = nan(N,D);

% continuous uniform
X11(:,1) = scale_lin(U(:,1), ranges.T_ambC);
X11(:,2) = scale_lin(U(:,2), ranges.RH_in);
X11(:,3) = scale_lin(U(:,3), ranges.u_feed);
X11(:,4) = scale_lin(U(:,4), ranges.p_vac);
X11(:,5) = scale_lin(U(:,5), ranges.T_reg);

% time log-uniform
X11(:,6) = scale_log(U(:,6), ranges.t_heat);
X11(:,7) = scale_log(U(:,7), ranges.t_des);
X11(:,8) = scale_log(U(:,8), ranges.t_ads);

% geometry uniform (continuous)
X11(:,9)  = scale_lin(U(:,9),  ranges.width);
X11(:,10) = scale_lin(U(:,10), ranges.height);

% N_cell integer
X11(:,11) = round(scale_lin(U(:,11), ranges.N_cell));
X11(:,11) = clamp(X11(:,11), ranges.N_cell(1), ranges.N_cell(2));
end

function x = scale_lin(u, mm)
x = mm(1) + (mm(2)-mm(1)).*u;
end

function x = scale_log(u, mm)
lo = mm(1); hi = mm(2);
x = 10.^(log10(lo) + (log10(hi)-log10(lo)).*u);
end

%% ============================= feasible domain from periodic =============================
function [poly, bounds, meta_periodic] = build_feasible_state_domain(periodic_mat, alpha)
if exist(periodic_mat,'file')~=2
    error('Periodic MAT not found: %s', periodic_mat);
end

S = load(periodic_mat);
meta_periodic = struct();
if isfield(S,'meta'), meta_periodic = S.meta; end

% Locate columns for qCO2_postDes and qH2O_postDes
% Your periodic dataset used Y columns:
% [qCO2_postAds, qCO2_postDes, qH2O_postAds, qH2O_postDes, ...]
Y = S.Y;
ok = true(size(Y,1),1);
if isfield(S,'ok'), ok = S.ok; end

qCO2 = Y(:,2);
qH2O = Y(:,4);

good = ok & isfinite(qCO2) & isfinite(qH2O) & (qCO2>=0) & (qH2O>=0);
P = [qCO2(good) qH2O(good)];
if size(P,1) < 50
    error('Not enough good periodic points to build boundary. good=%d', size(P,1));
end

% boundary polygon in (qCO2, qH2O) plane
k = boundary(P(:,1), P(:,2), alpha);
poly = P(k,:);

bounds = struct();
bounds.qCO2 = [min(P(:,1)) max(P(:,1))];
bounds.qH2O = [min(P(:,2)) max(P(:,2))];

fprintf('[STATE] boundary(alpha=%.2f) | nGood=%d | qCO2=[%.3g %.3g] | qH2O=[%.3g %.3g]\n', ...
    alpha, size(P,1), bounds.qCO2(1), bounds.qCO2(2), bounds.qH2O(1), bounds.qH2O(2));
end

function S = sample_states_in_polygon(N, poly, bounds, maxTries)
% Rejection sample uniform in bounding box then accept in polygon
xv = poly(:,1); yv = poly(:,2);

S = nan(N,2);
n = 0; tries = 0;

while n < N && tries < maxTries
    tries = tries + 1;

    % generate batch
    nb = min(5000, N-n);
    xr = bounds.qCO2(1) + (bounds.qCO2(2)-bounds.qCO2(1))*rand(nb,1);
    yr = bounds.qH2O(1) + (bounds.qH2O(2)-bounds.qH2O(1))*rand(nb,1);

    inside = inpolygon(xr, yr, xv, yv);
    keep = find(inside);
    nk = numel(keep);

    if nk > 0
        nk2 = min(nk, N-n);
        S(n+(1:nk2),:) = [xr(keep(1:nk2)) yr(keep(1:nk2))];
        n = n + nk2;
    end
end

if n < N
    error('State sampling failed: got %d/%d within polygon after %d tries. Increase maxTries.', n, N, tries);
end
end

%% ============================= utils =============================
function s = def(s, field, val)
if ~isfield(s,field) || isempty(s.(field))
    s.(field) = val;
end
end

function y = clamp(x, lo, hi)
y = min(max(x, lo), hi);
end

function U = simple_lhs(N, D, seed)
% Toolbox-free Latin Hypercube in [0,1]
% For each dimension: stratify into N bins, sample one point per bin, permute.

if nargin < 3, seed = 1; end
rng(seed);

U = zeros(N,D);
edges = (0:N-1)'/N;          % bin starts

for j = 1:D
    % sample uniformly within each bin
    u = edges + rand(N,1)/N; % N stratified samples
    % random permutation per dimension
    U(:,j) = u(randperm(N));
end
end
