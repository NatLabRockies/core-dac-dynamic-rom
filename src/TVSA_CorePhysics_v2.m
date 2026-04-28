% ========================= TVSA_CorePhysics_v14.m =======================
function varargout = TVSA_CorePhysics_v2(action, varargin)
% =========================================================================
% TVSA_CorePhysics_v14.m  (CORE)
% - All physics, drivers, periodic solver, and utilities live here.
% - Runner calls this file via action-dispatch.
%
% v10_9 changes vs v10_8 (MINIMAL-INVASIVE patch set):
%   (A) opt.corr.ppi_RH:
%       - Single RH-dependent multiplier applied inside PPI energy calc
%         (default: applied to heating-side only; configurable).
%   (B) opt.corr.h2o_eff:
%       - "Effective water" handling for H2O-only energy pathway:
%         affects ONLY the H2O term used in the regeneration heat (rxn power)
%         via dqH2O/dt (does NOT change qH2O state or CO2 coupling).
%   (C) Vacuum/pump RH double-count audit:
%       - Adds internal audit diagnostics and an optional comparison against
%         a "total-gas" vacuum work surrogate (OFF by default).
%
% NOTE:
% - Interface and main physics remain identical to v10_7 to be drop-in.
% - Corrections are guarded by opt.corr.* enable flags (safe when absent).
% =========================================================================
if nargin < 1 || isempty(action), error('TVSA_CorePhysics_final requires an action string.'); end
act = lower(strtrim(string(action)));

switch act
    case "finalize_design"
        varargout{1} = finalize_design(varargin{1});

    case "make_phys"
        varargout{1} = make_phys_params_SendiNoteS2(varargin{1});

    case "get_driver"
        varargout{1} = get_driver_explicit(varargin{1}, varargin{2}, varargin{3}, varargin{4});

    case "run_point_periodic"
        varargout{1} = run_point_periodic(varargin{:});

    case "run_s2_suite"
        varargout{1} = run_s2_suite(varargin{:});

    case "ref_12pt"
        varargout{1} = ref_sendi_12pt();

    otherwise
        error('Unknown action: %s', act);
end
end


%% =========================================================================
% MODE APPLIER (S0/S1/S2)
% =========================================================================
function opt = apply_mode(opt, mode)
m = upper(strtrim(string(mode)));
switch m
    case "S0"
        opt.co2_2pool.enable = false;
        opt.h2o_2pool.enable = false;
        opt.co2_2pool.resid.enable = false;
        opt.h2o_2pool.resid.enable = false;
    case "S1"
        opt.co2_2pool.enable = true;
        opt.h2o_2pool.enable = true;
        opt.co2_2pool.resid.enable = false;
        opt.h2o_2pool.resid.enable = false;
    otherwise % S2
        opt.co2_2pool.enable = true;
        opt.h2o_2pool.enable = true;
        opt.co2_2pool.resid.enable = true;
        opt.h2o_2pool.resid.enable = true;
end
end

%% =========================================================================
% Step-2 suite runner (baseline + optional 12-pt validation)
% Called by wrappers as: core('run_s2_suite', amb, ctrl, design, opt)
% =========================================================================
function out = run_s2_suite(amb, ctrl, design, opt)

% ---- safety: ensure required opt sub-structs exist (wrapper may omit them) ----
opt = ensure_opt_defaults(opt);

% ---- safety: ensure opt.ptraj exists (wrapper may omit it) ----
opt = ensure_ptraj_defaults(opt);
opt = ensure_aw_defaults(opt);

% ---- mode ----
if ~isfield(opt,'mode') || strlength(string(opt.mode))==0
    opt.mode = "S2";
end
opt = apply_mode(opt, opt.mode);

% ---- phys ----
phys = make_phys_params_SendiNoteS2(design);

% Apply fixed S_cap if present (paper-defensible single scalar)
if isfield(design,'S_cap_fixed') && isfinite(design.S_cap_fixed) && design.S_cap_fixed>0
    phys.iso.mech.S_cap = design.S_cap_fixed;
end

% Optional: CO2 kinetics scale knob
if isfield(opt,'kin') && isfield(opt.kin,'kCO2_scale') && isfinite(opt.kin.kCO2_scale)
    phys.kin.kCO2_ref = phys.kin.kCO2_ref * opt.kin.kCO2_scale;
end

% ---- driver (baseline recipe) ----
driver_cfg = struct('mode',"forced_only");
if isfield(opt,'driver_baseline') && strlength(string(opt.driver_baseline))>0
    driver_cfg.mode = string(opt.driver_baseline);
end
drv_base = get_driver_explicit(amb, ctrl, opt, driver_cfg);

% ---- baseline periodic ----
sim_base = run_point_periodic(drv_base, amb, ctrl, design, phys, opt, []);

out = struct();
out.baseline = struct( ...
    'Prod_tpy', sim_base.ppi.Prod_tCO2_collector_yr, ...
    'Elec_MWh_t', sim_base.ppi.Elec_total_MWh_tCO2, ...
    'periodic', sim_base.periodic);

% Optional time-series for plotting/debug (only when requested)
if isfield(opt,'do_plots') && opt.do_plots
    out.baseline.sim = sim_base;
end

% ---- 12-pt validation (Sendi) ----
do12 = false;
if isfield(opt,'run_12pt'), do12 = logical(opt.run_12pt); end

if do12
    ref = ref_sendi_12pt();
    npt = height(ref);

    Prod_pred = nan(npt,1);
    Elec_pred = nan(npt,1);

    sims = cell(npt,1);

    % Warm-start across points for stability/speed
    x_ws = sim_base.x_end;

    for k = 1:npt
        amb_k = amb;
        ctrl_k = ctrl;

        amb_k.T_ambC = ref.T_C(k);
        amb_k.RH_in  = ref.RH(k);

        % Keep Sendi's point-specific desorption time (dominant lever)
        if isfield(ctrl_k,'time_des_s')
            ctrl_k.time_des_s = ref.t_des_s(k);
        end

        drv_k = get_driver_explicit(amb_k, ctrl_k, opt, driver_cfg);
        sim_k = run_point_periodic(drv_k, amb_k, ctrl_k, design, phys, opt, x_ws);

        x_ws = sim_k.x_end;

        Prod_pred(k) = sim_k.ppi.Prod_tCO2_collector_yr;
        Elec_pred(k) = sim_k.ppi.Elec_total_MWh_tCO2;

        if isfield(opt,'do_plots') && opt.do_plots
            sims{k} = sim_k;
        end
    end

    Prod_ref = ref.Prod_ref;
    Elec_ref = ref.Elec_ref;

    mapeP = mean(abs((Prod_pred-Prod_ref)./Prod_ref))*100;
    mapeE = mean(abs((Elec_pred-Elec_ref)./Elec_ref))*100;

    % Trend score: rank-correlation vs reference (robust to scaling bias)
    rhoP = spearman_corr_safe(Prod_pred, Prod_ref);
    rhoE = spearman_corr_safe(Elec_pred, Elec_ref);
    trend = 5.0*(rhoP + rhoE + 2.0);   % maps [-1,1]x2 -> [0,20]

    out.report_12pt = struct( ...
        'ref', ref, ...
        'Prod_pred', Prod_pred, ...
        'Elec_pred', Elec_pred, ...
        'MAPE_Prod_pct', mapeP, ...
        'MAPE_Elec_pct', mapeE, ...
        'rho_Prod_spearman', rhoP, ...
        'rho_Elec_spearman', rhoE, ...
        'trend_score', trend);

    if isfield(opt,'do_plots') && opt.do_plots
        out.report_12pt.sims = sims;
    end
end

end

function rho = spearman_corr_safe(a,b)
a = a(:); b = b(:);
ok = isfinite(a) & isfinite(b);
a = a(ok); b = b(ok);
if numel(a) < 3
    rho = 0; return;
end
ra = tiedrank(a);
rb = tiedrank(b);
ra = ra - mean(ra);
rb = rb - mean(rb);
rho = (ra.'*rb) / max(1e-12, sqrt((ra.'*ra)*(rb.'*rb)));
rho = max(min(rho,1),-1);
end



%% =========================================================================
% Geometry finalizer (thin-bed conventions fixed)
% =========================================================================
function design = finalize_design(design)
% Backward/Wrapper compatibility: accept alternative field names.
if ~isfield(design,'L_flow_fixed_m') || isempty(design.L_flow_fixed_m)
    if isfield(design,'bed_thickness_m') && ~isempty(design.bed_thickness_m)
        design.L_flow_fixed_m = design.bed_thickness_m;
    elseif isfield(design,'L_flow_m') && ~isempty(design.L_flow_m)
        design.L_flow_fixed_m = design.L_flow_m;
    else
        error('finalize_design:missing_L_flow', 'Missing design.L_flow_fixed_m (or bed_thickness_m / L_flow_m).');
    end
end

if ~isfield(design,'rb_fixed') || isempty(design.rb_fixed)
    if isfield(design,'rb_fixed_kg_m3') && ~isempty(design.rb_fixed_kg_m3)
        design.rb_fixed = design.rb_fixed_kg_m3;
    elseif isfield(design,'rb_kg_m3') && ~isempty(design.rb_kg_m3)
        design.rb_fixed = design.rb_kg_m3;
    else        % Fallback to v10_8 baseline if wrapper did not provide packing density
        design.rb_fixed = 356.38; % kg/m3
        
    end
end



% ---------------- DEFAULTS / compatibility (v10_9) ----------------
% Wrapper codes may pass a minimal `design` struct. To avoid brittle crashes,
% we back-fill common geometry/packing defaults consistent with the v10_8/Sendi
% baseline used in this project. Any missing field is filled ONCE here.

if ~isfield(design,'A_mult') || isempty(design.A_mult)
    if isfield(design,'A_mult_fixed') && ~isempty(design.A_mult_fixed)
        design.A_mult = design.A_mult_fixed;
    else
        design.A_mult = 1.0;
    end
end

if ~isfield(design,'width_cell_m')  || isempty(design.width_cell_m)
    design.width_cell_m = 1.43; % Sendi S1 default
end
if ~isfield(design,'height_cell_m') || isempty(design.height_cell_m)
    design.height_cell_m = 0.10; % Sendi S1 default
end
if ~isfield(design,'N_cell') || isempty(design.N_cell)
    design.N_cell = 1144; % Sendi S1 default
end

if ~isfield(design,'Ng_ref') || isempty(design.Ng_ref)
    design.Ng_ref = 35; % v10_8 baseline
end
if ~isfield(design,'eps') || isempty(design.eps)
    design.eps = 0.56; % v10_8 baseline
end

if ~isfield(design,'L_flow_fixed_m') || isempty(design.L_flow_fixed_m)
    design.L_flow_fixed_m = 0.0172; % m (v10_8 baseline)
end

% If rb_fixed is still missing, adopt the v10_8 baseline value instead of erroring.
if ~isfield(design,'rb_fixed') || isempty(design.rb_fixed)
    design.rb_fixed = 356.38; % kg/m3 (v10_8 baseline)
end
% -------------------------------------------------------------------

A_cell = design.width_cell_m * design.height_cell_m;
design.A_flow_m2 = max(design.N_cell * A_cell, 1e-12);
design.A_contact_m2 = max(design.A_flow_m2 * design.A_mult, 1e-12);

design.L_flow_m = design.L_flow_fixed_m;
design.Vbed_m3  = design.A_flow_m2 * design.L_flow_m;

design.msorb_kg = design.rb_fixed * design.Vbed_m3;
design.L_contact_m = design.Vbed_m3 / max(design.A_contact_m2, 1e-12);

% Step-2 convention: Vg = Ng_ref * eps*Vbed
design.Vg_m3     = design.Ng_ref * design.eps * design.Vbed_m3;
design.Vg_ref_m3 = design.Vg_m3;
end


%% =========================================================================
% PHYS PARAMS (Sendi SI Note S2)
% =========================================================================
function phys = make_phys_params_SendiNoteS2(design)
phys = struct();

phys.bed = struct();
phys.bed.rb  = design.msorb_kg / max(design.Vbed_m3,1e-12);

phys.eq = struct();
phys.eq.R = 8.314462618;

phys.eq.MW_CO2 = 44e-3;
phys.eq.MW_H2O = 18e-3;
phys.eq.MW_air = 28.97e-3;

phys.eq.Cp_CO2_g = 37.0;
phys.eq.Cp_H2O_g = 34.0;
phys.eq.Cp_air_g = 29.0;
phys.eq.Cp_s     = 1880;

phys.eq.dHads_CO2 = 70e3;
phys.eq.dHads_H2O = 46e3;

% Steam surrogate
phys.eq.r_steam_CO2 = 0.9;
phys.eq.TsatC = 100.0;
phys.eq.Cp_liq = 4180;
phys.eq.hvap  = 2257000;

% COP surrogate
phys.eq.eta2nd  = 0.50;
phys.eq.dTmin_K = 10.0;

% Efficiencies
phys.eq.eta_fan  = 0.60;
phys.eq.eta_comp = 0.60;

% Product CO2 compression
phys.eq.P_CO2_out_bar = 150;
phys.eq.T_comp_K      = 298.15;

% Isotherm params
phys.iso = noteS2_params();

% Optional: CO2 isotherm capacity scale (single scalar, default 1.0)
if isfield(design,'co2_dry_scale') && ~isempty(design.co2_dry_scale)
    phys.iso.mech.dry_scale = design.co2_dry_scale;
end

% Kinetics (Arrhenius + mild pCO2 enhancement)
phys.kin = struct();
phys.kin.kCO2_ref = 3.0e-3;
phys.kin.kH2O_ref = 8.6e-3;
phys.kin.Ea_CO2   = 35e3;
phys.kin.Ea_H2O   = 20e3;
phys.kin.Tref     = 293.15;
phys.kin.p_ref_Pa = 40.0;
phys.kin.beta     = 1.0;
phys.kin.n_exp    = 1.0;

% stage multipliers [vac heat des cool press ads]
phys.kin.mCO2 = [0.2, 0.002, 0.009, 0.001, 0.001, 1.0];
phys.kin.mH2O = [0.20,0.18, 0.10, 0.20, 0.25, 1.00];
end



function iso = noteS2_params()
iso = struct();

% Mechanistic (Sendi SI Note S2)
iso.mech.T0       = 298.15;
iso.mech.qinf0    = 4.86;
iso.mech.chi      = 0.0;
iso.mech.b0       = 2.85e-21;
iso.mech.DH_dry   = 117798;
iso.mech.DH_wet   = 130155;
iso.mech.tau0     = 0.209;
iso.mech.alpha    = 0.523;
iso.mech.phi_max  = 1.000;
iso.mech.fblk_max = 0.433;
iso.mech.k_blk    = 0.795;
iso.mech.phi_dry  = 1.000;
iso.mech.Aexp     = 1.535;
iso.mech.n_blk    = 1.425;

% GAB (Sendi SI Note S2)
iso.gab.qm = 3.63;
iso.gab.C  = 47110;
iso.gab.D  = 0.023744;
iso.gab.F  = 57706;
iso.gab.G  = -47.814;

iso.mech.S_cap = 1.0;
iso.mech.dry_scale = 1.0; % [-] optional CO2 isotherm capacity scale (default OFF)
end


%% =========================================================================
% DRIVER (locked_only | locked_shape_scaled | forced_only | auto)
% =========================================================================
function drv = get_driver_explicit(amb, ctrl, opt, driver_cfg)

% ---- safety: ensure required opt sub-structs exist (wrapper may omit them) ----
opt = ensure_opt_defaults(opt);

mode = lower(strtrim(string(driver_cfg.mode)));

if any(mode == ["locked","lock","locked_only"]), mode = "locked_only"; end
if any(mode == ["scaled","lock_shape","locked_shape","locked_shape_scaled"]), mode = "locked_shape_scaled"; end
if any(mode == ["forced","force","forced_only"]), mode = "forced_only"; end
if any(mode == ["auto"]), mode = "auto"; end

locked_file = "";
if isfield(driver_cfg,'locked_file'), locked_file = string(driver_cfg.locked_file); end

locked_ok = false;
drv_locked = [];
if strlength(locked_file) > 0 && exist(locked_file,'file')==2
    S = load(locked_file);
    if isfield(S,'drv')
        drv_locked = validate_and_finalize_drv(S.drv, opt);
        locked_ok = true;
    end
end

if mode=="locked_only"
    if ~locked_ok
        error('locked_only selected but locked driver not found: %s', locked_file);
    end
    drv = drv_locked;
    drv.meta.source = "locked";
    drv.meta.locked_file = locked_file;
    return;
end

if mode=="locked_shape_scaled"
    if ~locked_ok
        error('locked_shape_scaled selected but locked driver not found: %s', locked_file);
    end
    drv = scale_locked_shape_to_ctrl_amb(drv_locked, amb, ctrl, opt);
    drv = validate_and_finalize_drv(drv, opt);
    drv.meta.source = "locked_shape_scaled";
    drv.meta.locked_file = locked_file;
    return;
end

if mode=="auto" && locked_ok
    try
        drv = validate_and_finalize_drv(drv_locked, opt);
        drv.meta.source = "locked";
        drv.meta.locked_file = locked_file;
        return;
    catch
        drv = scale_locked_shape_to_ctrl_amb(drv_locked, amb, ctrl, opt);
        drv = validate_and_finalize_drv(drv, opt);
        drv.meta.source = "locked_shape_scaled";
        drv.meta.locked_file = locked_file;
        return;
    end
end

% forced_only OR auto fallback
drv = build_forced_profiles(amb, ctrl, opt.dt_s);
drv = validate_and_finalize_drv(drv, opt);
drv.meta.source = "forced";
end


function drv = validate_and_finalize_drv(drv, opt)
req = {'t_s','stage_id','TK_cmd','RH_cmd'};
for i=1:numel(req)
    if ~isfield(drv, req{i}), error('Driver missing field: %s', req{i}); end
end

if isfield(drv,'P_Pa')
    Pvec = drv.P_Pa;
elseif isfield(drv,'Pa')
    Pvec = drv.Pa;
else
    error('Driver missing pressure field: P_Pa (preferred) or Pa');
end

drv.t_s      = drv.t_s(:);
drv.stage_id = uint8(drv.stage_id(:));
drv.TK_cmd   = drv.TK_cmd(:);
drv.RH_cmd   = drv.RH_cmd(:);
drv.P_Pa     = Pvec(:);
drv.Pa       = drv.P_Pa;

N = numel(drv.t_s);
if any([numel(drv.stage_id),numel(drv.TK_cmd),numel(drv.RH_cmd),numel(drv.P_Pa)] ~= N)
    error('Driver vectors must have the same length.');
end

dt = median(diff(drv.t_s));
if abs(dt - opt.dt_s) > 1e-9
    error('Driver dt=%.6g s does not match opt.dt_s=%.6g s.', dt, opt.dt_s);
end

drv.RH_cmd = clamp(drv.RH_cmd, 0, 0.999);
drv.P_Pa   = max(drv.P_Pa, 200.0);
drv.TK_cmd = max(drv.TK_cmd, 200.0);

if ~isfield(drv,'edges_s') || isempty(drv.edges_s)
    drv.edges_s = build_edges_from_stage(drv.t_s, drv.stage_id);
end
if ~isfield(drv,'stage_names') || isempty(drv.stage_names)
    drv.stage_names = {'Vacuum','Heating','Desorption','Cooling','Pressurization','Adsorption'};
end

drv.Pbar  = drv.P_Pa/1e5;
drv.TC_cmd= drv.TK_cmd - 273.15;

if ~isfield(drv,'meta') || isempty(drv.meta), drv.meta = struct(); end
drv.meta.ptraj = opt.ptraj;

% Smooth P edges (optional)
if isfield(opt,'ptraj') && isfield(opt.ptraj,'driver_edge_smooth_enable') && opt.ptraj.driver_edge_smooth_enable
    drv = smooth_driver_pressure_edges(drv, opt.dt_s, opt.ptraj.driver_edge_ramp_s);
end
end

function edges = build_edges_from_stage(t_s, stage_id)
edges = t_s(1);
for k=2:numel(t_s)
    if stage_id(k) ~= stage_id(k-1)
        edges(end+1,1) = t_s(k); %#ok<AGROW>
    end
end
edges(end+1,1) = t_s(end);
end


function drv = build_forced_profiles(amb, ctrl, dt)
[t_s, stage_id, edges_s, stage_names] = build_schedule(ctrl, dt);
Pamb_bar = ambient_pressure_bar(amb.altitude_m);

optp = struct();
optp.T_heat_frac_end = 0.78;
optp.T_heat_shape    = 1.05;
optp.T_des_tau_frac  = 0.55;
optp.T_cool_tau_frac = 0.10;
optp.T_heat_lag_frac = 0.05;

optp.RH_floor_regen  = 0.05;
optp.k_RH_des        = 5.0;

optp.P_blow_tau_frac = 0.20;
optp.P_cool_extra_drop_bar = 0.015;
optp.P_press_tau_frac = 0.35;

[T_C_cmd, RH_cmd, P_bar] = generate_profiles(t_s, stage_id, edges_s, ...
    amb.T_ambC, amb.RH_in, Pamb_bar, ctrl.T_regenC, ctrl.p_vac_bar, optp);

drv = struct();
drv.t_s = t_s(:);
drv.stage_id = stage_id(:);
drv.edges_s = edges_s(:);
drv.stage_names = stage_names;

drv.TC_cmd = T_C_cmd(:);
drv.TK_cmd = T_C_cmd(:) + 273.15;
drv.RH_cmd = clamp(RH_cmd(:), 0, 1);

drv.Pbar = P_bar(:);
drv.P_Pa = P_bar(:) * 1e5;
drv.Pa   = drv.P_Pa;
end

function [t, stage_id, edges, stage_names] = build_schedule(ctrl, dt)
t1 = ctrl.time_vac_s;
t2 = t1 + ctrl.time_heat_s;
t3 = t2 + ctrl.time_des_s;
t4 = t3 + ctrl.time_cool_s;
t5 = t4 + ctrl.time_press_s;
t6 = t5 + ctrl.time_ads_s;

t = (0:dt:t6)';

stage_id = zeros(size(t),'uint8');
for i = 1:numel(t)
    ti = t(i);
    if     ti <  t1, stage_id(i)=1;
    elseif ti <  t2, stage_id(i)=2;
    elseif ti <  t3, stage_id(i)=3;
    elseif ti <  t4, stage_id(i)=4;
    elseif ti <  t5, stage_id(i)=5;
    else,            stage_id(i)=6;
    end
end

edges = [0; t1; t2; t3; t4; t5; t6];
stage_names = {'Vacuum','Heating','Desorption','Cooling','Pressurization','Adsorption'};
end

function [T_C, RH, P_bar] = generate_profiles(t, stage_id, edges, ...
    TambC, RHin, Pamb_bar, TregenC, Pvac_bar, opt)

T_C   = zeros(size(t));
RH    = zeros(size(t));
P_bar = zeros(size(t));

dT = (TregenC - TambC);
T_heat_end = TambC + opt.T_heat_frac_end * dT;

for i = 1:numel(t)
    ti = t(i);
    stg = stage_id(i);

    % Pressure
    if stg == 1
        tau = max(opt.P_blow_tau_frac * (edges(2) - edges(1)), 1e-6);
        P_bar(i) = Pvac_bar + (Pamb_bar - Pvac_bar)*exp(-ti/tau);
    elseif stg == 2 || stg == 3
        P_bar(i) = Pvac_bar;
    elseif stg == 4
        P_bar(i) = max(Pvac_bar - opt.P_cool_extra_drop_bar, 0.001);
    elseif stg == 5
        tau = max(opt.P_press_tau_frac * (edges(6) - edges(5)), 1e-6);
        dtp = max(ti - edges(5), 0);
        Pstart = max(Pvac_bar - opt.P_cool_extra_drop_bar, 0.001);
        P_bar(i) = Pamb_bar - (Pamb_bar - Pstart)*exp(-dtp/tau);
    else
        P_bar(i) = Pamb_bar;
    end

    % Temperature commanded
    if stg == 1
        T_C(i) = TambC;
    elseif stg == 2
        u = clamp01((ti - edges(2)) / max(edges(3) - edges(2), 1e-12));
        u_lag = opt.T_heat_lag_frac;
        if u <= u_lag
            s = smoothstep(u / max(u_lag,1e-12));
            frac = 0.10 * s;
        else
            u2 = (u - u_lag) / max(1 - u_lag, 1e-12);
            frac = 0.10 + 0.90 * (u2.^opt.T_heat_shape);
        end
        T_C(i) = TambC + (T_heat_end - TambC) * frac;
    elseif stg == 3
        u = clamp01((ti - edges(3)) / max(edges(4) - edges(3), 1e-12));
        tau = max(opt.T_des_tau_frac, 1e-3);
        g = (1 - exp(-u/tau)) / (1 - exp(-1/tau));
        T_C(i) = T_heat_end + (TregenC - T_heat_end) * g;
    else
        T3 = TregenC;
        dtc = max(ti - edges(4), 0);
        tau = max(opt.T_cool_tau_frac * (edges(7) - edges(4)), 1e-6);
        T_C(i) = TambC + (T3 - TambC) * exp(-dtc/tau);
    end

    % RH commanded
    if stg == 1
        RH(i) = RHin;
    elseif stg == 2
        u = clamp01((ti - edges(2)) / max(edges(3) - edges(2), 1e-12));
        RH_end_heat = min(0.98, RHin * 1.10);
        RH(i) = (1-u)*RHin + u*RH_end_heat;
    elseif stg == 3
        u = clamp01((ti - edges(3)) / max(edges(4) - edges(3), 1e-12));
        RH0_des  = min(0.98, RHin * 1.10);
        RH_floor = opt.RH_floor_regen;
        RH(i) = RH_floor + (RH0_des - RH_floor) * exp(-opt.k_RH_des * u);
    elseif stg == 4 || stg == 5
        RH(i) = max(opt.RH_floor_regen*0.85, 0.0);
    else
        u = clamp01((ti - edges(6)) / max(edges(7) - edges(6), 1e-12));
        tau = max(0.24, 1e-3);
        p = 1.55;
        rec = 1 - exp(-(u/max(tau,1e-6)).^p);
        RH_start = max(opt.RH_floor_regen*0.85, 0.0);
        RH(i) = RH_start + (RHin - RH_start)*rec;
    end
    RH(i) = clamp(RH(i), 0, 1);
end
end


function drv = scale_locked_shape_to_ctrl_amb(drvL, amb, ctrl, opt)
dt = opt.dt_s;
[tT, sidT, edgesT, stage_names] = build_schedule(ctrl, dt);
N = numel(tT);

TK_cmd_T = zeros(N,1);
RH_cmd_T = zeros(N,1);
P_Pa_T   = zeros(N,1);

PambL_bar = max(drvL.P_Pa(drvL.stage_id==6))/1e5;
PvacL_bar = min(drvL.P_Pa(drvL.stage_id==2 | drvL.stage_id==3))/1e5;
PambL_bar = max(PambL_bar, 0.2);
PvacL_bar = max(min(PvacL_bar, PambL_bar), 0.001);

PambT_bar = ambient_pressure_bar(amb.altitude_m);
PvacT_bar = ctrl.p_vac_bar;

TKL_min = min(drvL.TK_cmd);
TKL_max = max(drvL.TK_cmd);
if (TKL_max - TKL_min) < 1e-6, TKL_max = TKL_min + 1.0; end

TKT_min = amb.T_ambC + 273.15;
TKT_max = ctrl.T_regenC + 273.15;

spanL = TKL_max - TKL_min;
spanT = TKT_max - TKT_min;

idx_adsL = find(drvL.stage_id==6);
if isempty(idx_adsL), RH_adsL = mean(drvL.RH_cmd);
else, RH_adsL = mean(drvL.RH_cmd(idx_adsL));
end
RH_scale = amb.RH_in / max(RH_adsL, 1e-6);

for st = 1:6
    idxL = find(drvL.stage_id==st);
    idxT = find(sidT==st);

    if isempty(idxT), continue; end

    if isempty(idxL)
        TK_cmd_T(idxT) = TKT_min;
        RH_cmd_T(idxT) = amb.RH_in;
        if st==2 || st==3 || st==1 || st==4
            P_Pa_T(idxT) = (PvacT_bar*1e5);
        else
            P_Pa_T(idxT) = (PambT_bar*1e5);
        end
        continue;
    end

    tL = drvL.t_s(idxL);
    uL = (tL - tL(1)) / max(tL(end) - tL(1), 1e-12);

    tS = tT(idxT);
    uT = (tS - tS(1)) / max(tS(end) - tS(1), 1e-12);

    TK_shape = interp1(uL, drvL.TK_cmd(idxL), uT, 'linear', 'extrap');
    RH_shape = interp1(uL, drvL.RH_cmd(idxL), uT, 'linear', 'extrap');
    P_shape  = interp1(uL, drvL.P_Pa(idxL),   uT, 'linear', 'extrap');

    TK_cmd_T(idxT) = TKT_min + spanT * ((TK_shape - TKL_min) / spanL);

    PL_bar = P_shape/1e5;
    denomP = max(PambL_bar - PvacL_bar, 1e-9);
    fracP  = (PL_bar - PvacL_bar) / denomP;
    fracP  = clamp(fracP, -0.25, 1.25);

    PT_bar = PvacT_bar + fracP * (PambT_bar - PvacT_bar);
    P_Pa_T(idxT) = max(PT_bar, 0.001) * 1e5;

    RH_cmd_T(idxT) = clamp(RH_shape * RH_scale, 0, 0.999);
end

drv = struct();
drv.t_s = tT(:);
drv.stage_id = uint8(sidT(:));
drv.edges_s = edgesT(:);
drv.stage_names = stage_names;

drv.TK_cmd = max(TK_cmd_T(:), 200.0);
drv.RH_cmd = clamp(RH_cmd_T(:), 0, 0.999);
drv.P_Pa   = max(P_Pa_T(:), 200.0);
drv.Pa     = drv.P_Pa;
drv.Pbar   = drv.P_Pa/1e5;
drv.TC_cmd = drv.TK_cmd - 273.15;
end


%% =========================================================================
% PERIODIC runner (warm-start override supported)
% =========================================================================
function sim = run_point_periodic(drv, amb, ctrl, design, phys, opt, x_init_override)

% ---- safety: ensure opt.ptraj exists (caller may omit it) ----
opt = ensure_ptraj_defaults(opt);
opt = ensure_aw_defaults(opt);

use_override = (nargin >= 7) && ~isempty(x_init_override);

if use_override
    x0 = x_init_override(:);
    if numel(x0) ~= 6
        error('x_init_override must be 6x1: [qCO2f;qCO2s;qH2Of;qH2Ob;yCO2g;TKbed].');
    end
    x0(1:4) = max(x0(1:4), 0);
    x0(5)   = clamp(x0(5), opt.yco2.y_floor, opt.yco2.y_max);
    x0(6)   = max(x0(6), 200.0);
else
    TKcmd0 = drv.TK_cmd(1);
    P0     = drv.P_Pa(1);
    TKbed0 = TKcmd0;

    [~, Pdry0, aw_gas0] = pH2O_Pdry_aw_chain(RHclamp(drv.RH_cmd(1)), TKcmd0, TKbed0, P0, opt);
    pCO2_in0 = amb.yCO2_air * Pdry0;

    qH2O_eq0 = iso_GAB_H2O_noteS2(TKbed0, aw_gas0, phys.iso, opt.aw.max);

    if opt.h2o_2pool.enable
        fb0 = h2o_fbound_eq(TKbed0, aw_gas0, uint8(drv.stage_id(1)), opt.h2o_2pool.fbound);
        qH2Ob0 = fb0 * qH2O_eq0;
        qH2Of0 = max(qH2O_eq0 - qH2Ob0, 0);
    else
        qH2Of0 = qH2O_eq0;
        qH2Ob0 = 0;
    end
    qH2Otot0 = qH2Of0 + qH2Ob0;

    qCO2tot0 = iso_mech_CO2_noteS2(TKbed0, pCO2_in0, qH2Otot0, phys.iso);
    aw_eff0 = aw_for_co2_effects(uint8(drv.stage_id(1)), TKbed0, aw_gas0, qH2Of0, qH2Ob0, phys.iso, opt);

    if opt.co2_2pool.enable
        fslow0 = co2_fslow_eq(TKbed0, aw_eff0, uint8(drv.stage_id(1)), opt.co2_2pool.fslow);
        qCO2s0 = fslow0 * qCO2tot0;
        qCO2f0 = max(qCO2tot0 - qCO2s0, 0);
    else
        qCO2f0 = qCO2tot0;
        qCO2s0 = 0;
    end

    yCO2g0 = clamp(amb.yCO2_air, opt.yco2.y_floor, opt.yco2.y_max);
    x0 = [qCO2f0; qCO2s0; qH2Of0; qH2Ob0; yCO2g0; TKbed0];
end

relchg = inf; iters = 0;
while (iters < opt.max_iters) && (relchg > opt.tol_rel)
    iters = iters + 1;
    out = run_one_cycle_forward(drv, amb, ctrl, design, phys, opt, x0);

    x1 = out.x_end;
    w = clamp(opt.omega, 0.05, 0.95);
    x1 = w*x1 + (1-w)*x0;

    relchg = norm(x1 - x0,2) / max(1e-12, norm(x0,2));
    x0 = x1;

    fprintf('Periodic iter %02d | relchg = %.3e | Prod=%.2f | Elec=%.3f\n', ...
        iters, relchg, out.ppi.Prod_tCO2_collector_yr, out.ppi.Elec_total_MWh_tCO2);
end

out.periodic = struct('iters',iters,'relchg',relchg);
sim = out;
end


%% =========================================================================
% ONE CYCLE: analytic LDF + yCO2-only gas, ptraj inversion flows
% =========================================================================
function out = run_one_cycle_forward(drv, amb, ctrl, design, phys, opt, x_init)
% Wrapper that optionally performs lightweight early-stop logic by:
%   (1) running one pass with max times to estimate t_ads_eff/t_des_eff
%   (2) rebuilding a forced-only driver with shortened Ads/Des times
%   (3) rerunning once for consistent states/PPIs

opt = ensure_opt_defaults(opt);

if ~isfield(opt,'early_stop') || ~isfield(opt.early_stop,'enable') || ~opt.early_stop.enable
    out = run_one_cycle_forward_base_v11(drv, amb, ctrl, design, phys, opt, x_init);
    out.cyc.t_ads_eff_s = ctrl.time_ads_s;
    out.cyc.t_des_eff_s = ctrl.time_des_s;
    out.cyc.early_stop_changed = false;
    return;
end

out1 = run_one_cycle_forward_base_v11(drv, amb, ctrl, design, phys, opt, x_init);

[ctrl2, did_change, info] = early_stop_adjust_ctrl_from_pass1(out1, ctrl, opt);

if ~did_change
    out = out1;
    out.cyc.t_ads_eff_s = info.t_ads_eff_s;
    out.cyc.t_des_eff_s = info.t_des_eff_s;
    out.cyc.early_stop_changed = false;
    out.cyc.early_stop_info = info;
    return;
end

% Only rebuild drivers when the incoming driver is forced-like.
src = "";
if isfield(out1,'drv') && isfield(out1.drv,'meta') && isfield(out1.drv.meta,'source')
    src = string(out1.drv.meta.source);
end
if ~(contains(lower(src),"forced"))
    out = out1;
    out.cyc.t_ads_eff_s = info.t_ads_eff_s;
    out.cyc.t_des_eff_s = info.t_des_eff_s;
    out.cyc.early_stop_changed = false;
    out.cyc.early_stop_info = info;
    return;
end

drv2 = build_forced_profiles(amb, ctrl2, opt.dt_s);
drv2 = validate_and_finalize_drv(drv2, opt);
drv2.meta.source = "forced_early_stop";

out = run_one_cycle_forward_base_v11(drv2, amb, ctrl2, design, phys, opt, x_init);
out.cyc.t_ads_eff_s = ctrl2.time_ads_s;
out.cyc.t_des_eff_s = ctrl2.time_des_s;
out.cyc.early_stop_changed = true;
out.cyc.early_stop_info = info;
end

function [ctrl2, did_change, info] = early_stop_adjust_ctrl_from_pass1(sim, ctrl_in, opt)
% Determine effective Ads/Des durations from a max-time pass, then shorten
% ctrl.times and rerun once (wrapper handles rerun).
%
% Ads stop (relative): qCO2 >= eta_ads * qCO2_eq
% Des stop (relative + floor): qCO2 <= q_target where
%   q_target = q_floor + (1-eta_des)*(q_des_start - q_floor)
%   q_floor  = max(q_abs_min, f_res*q_ads_end)

ctrl2 = ctrl_in;
did_change = false;
info = struct();

dt = opt.dt_s;

% defaults
eta_ads = 0.95; eta_des = 0.95;
t_ads_min = 0.0; t_des_min = 0.0;
Tgate_frac = 0.95;
Tgate_relC = nan;

if isfield(opt,'early_stop')
    if isfield(opt.early_stop,'eta_ads'), eta_ads = opt.early_stop.eta_ads; end
    if isfield(opt.early_stop,'eta_des'), eta_des = opt.early_stop.eta_des; end
    if isfield(opt.early_stop,'t_ads_min_s'), t_ads_min = opt.early_stop.t_ads_min_s; end
    if isfield(opt.early_stop,'t_des_min_s'), t_des_min = opt.early_stop.t_des_min_s; end
    if isfield(opt.early_stop,'Tgate_frac'), Tgate_frac = opt.early_stop.Tgate_frac; end
    if isfield(opt.early_stop,'Tgate_relC'), Tgate_relC = opt.early_stop.Tgate_relC; end
end

info.eta_ads = eta_ads;
info.eta_des = eta_des;
info.t_ads_min_s = t_ads_min;
info.t_des_min_s = t_des_min;
info.Tgate_frac  = Tgate_frac;
info.Tgate_relC = Tgate_relC;

% ----- Adsorption effective time -----
t_ads_eff = ctrl_in.time_ads_s;
idx_ads = find(sim.drv.stage_id==6);
if ~isempty(idx_ads)
    t0   = sim.drv.t_s(idx_ads(1));
    t_in = sim.drv.t_s(idx_ads) - t0;q = sim.prof.qCO2_total(idx_ads);
% NEW (v14): approach-to-equilibrium criterion (paper/DB-robust)
% approach = (q - q0) / (q_eq - q0)  in [0,1]
if isfield(sim,'prof') && isfield(sim.prof,'qCO2_eq_total')
    qeq = sim.prof.qCO2_eq_total(idx_ads);
else
    qeq = nan(size(q)); % fallback: no early-stop if qeq unavailable
end

q0 = q(1);

% If no meaningful driving force, keep nominal time
den = qeq - q0;
has_drive = any(isfinite(den) & (den > 1e-12));

if isfinite(q0) && has_drive
    mask = (t_in >= t_ads_min) & isfinite(q) & isfinite(qeq) & (den > 1e-12);

    approach = zeros(size(q));
    approach(mask) = (q(mask) - q0) ./ den(mask);
    approach = max(0.0, min(1.0, approach));

    crit = mask & (approach >= eta_ads);

    if any(crit)
        i0 = find(crit,1,'first');
        t_ads_eff = max(t_ads_min, t_in(i0));
        t_ads_eff = min(t_ads_eff, ctrl_in.time_ads_s);
        t_ads_eff = ceil(t_ads_eff/dt)*dt;
        t_ads_eff = min(t_ads_eff, ctrl_in.time_ads_s);
    end
end
end


% ----- Desorption effective time (deterministic crossing) -----
t_des_eff = ctrl_in.time_des_s;

idx_des = find(sim.drv.stage_id==3);
if ~isempty(idx_des)
    t0   = sim.drv.t_s(idx_des(1));
    t_in = sim.drv.t_s(idx_des) - t0;

    q    = sim.prof.qCO2_total(idx_des);
    Tbed = sim.prof.TK_bed(idx_des);

    % Guard: q must have at least one finite sample
    q_start = q(find(isfinite(q), 1, 'first'));
    if isempty(q_start)
        % no usable data -> keep default t_des_eff
    else
        q_ads_end  = sim.cyc.qCO2_ads_end;
        q_abs_min  = opt.co2_2pool.resid.q_abs_min;
        f_res      = opt.co2_2pool.resid.f_res;

        q_floor  = max(q_abs_min, f_res * max(q_ads_end, 1e-12));
        q_target = q_floor + (1-eta_des) * (q_start - q_floor);

                % NEW (v14): temperature gate as a fraction of user-selected T_regenC (Celsius)
        if isfinite(Tgate_relC)
            T_gate_C = Tgate_relC * ctrl_in.T_regenC;
            T_gate   = T_gate_C + 273.15;
        else
            % legacy behavior (kept for backward compatibility)
            T_gate   = (ctrl_in.T_regenC + 273.15) * Tgate_frac;
        end

        % Valid region (apply gates + finite)
        valid = (t_in >= t_des_min) & (Tbed >= T_gate) & isfinite(q);

        % Extract valid sequences
        tv = t_in(valid);
        qv = q(valid);

        if numel(tv) >= 1
            % If already below target at first valid point -> immediate (at first valid time)
            if qv(1) <= q_target
                t_hit = tv(1);

            else
                % Find first crossing from above->below between samples
                above = (qv > q_target);
                below = (qv <= q_target);

                k = find(above(1:end-1) & below(2:end), 1, 'first');

                if ~isempty(k)
                    % Linear interpolation to get deterministic crossing time
                    t1 = tv(k);   t2 = tv(k+1);
                    q1 = qv(k);   q2 = qv(k+1);

                    if q2 ~= q1
                        frac = (q1 - q_target) / (q1 - q2);   % in (0,1]
                        t_hit = t1 + frac * (t2 - t1);
                    else
                        % Flat segment but next is considered below (numerically): take t2
                        t_hit = t2;
                    end

                else
                    % No crossing found in valid window.
                    % Fallback rule (make this explicit):
                    % - if the last valid sample is below target => take last time
                    % - otherwise => treat as "not reached", keep default (full des time)
                    if qv(end) <= q_target
                        t_hit = tv(end);
                    else
                        t_hit = []; % not reached
                    end
                end
            end

            if ~isempty(t_hit)
                % Apply bounds + dt quantization consistently
                t_des_eff = max(t_des_min, t_hit);
                t_des_eff = min(t_des_eff, ctrl_in.time_des_s);

                % Quantize to dt (choose ONE policy; ceil is conservative)
                t_des_eff = ceil(t_des_eff/dt)*dt;
                t_des_eff = min(t_des_eff, ctrl_in.time_des_s);
            end
        end
    end
end

info.t_ads_eff_s = t_ads_eff;
info.t_des_eff_s = t_des_eff;

% --- NaN-safe clamp for effective times (avoid propagation) ---
t_ads_eff = info.t_ads_eff_s;
if ~isfinite(t_ads_eff); t_ads_eff = ctrl_in.time_ads_s; end

t_des_eff = info.t_des_eff_s;
if ~isfinite(t_des_eff); t_des_eff = ctrl_in.time_des_s; end

ctrl2.time_ads_s = max(t_ads_eff, opt.early_stop.t_ads_min_s);
ctrl2.time_des_s = max(t_des_eff, opt.early_stop.t_des_min_s);
assert(isfinite(ctrl2.time_ads_s) && ctrl2.time_ads_s>0, 'early_stop: time_ads_s became NaN/invalid');
assert(isfinite(ctrl2.time_des_s) && ctrl2.time_des_s>0, 'early_stop: time_des_s became NaN/invalid');

end

function out = run_one_cycle_forward_base_v11(drv, amb, ctrl, design, phys, opt, x_init)

dt = opt.dt_s;
t  = drv.t_s(:);
N  = numel(t);


% optional stabilizers (initialized on first use)
pdrv_filt_prev   = NaN;
aw_eff_filt_prev = NaN;

% store CO2 equilibrium (for early-stop logic / diagnostics)
qCO2eqTot_vec = zeros(N,1);
% states
qCO2f = zeros(N,1); qCO2s = zeros(N,1);
qH2Of = zeros(N,1); qH2Ob = zeros(N,1);
yCO2g = zeros(N,1); TKbed = zeros(N,1);

qCO2f(1) = x_init(1); qCO2s(1) = x_init(2);
qH2Of(1) = x_init(3); qH2Ob(1) = x_init(4);
yCO2g(1) = x_init(5); TKbed(1) = x_init(6);

qCO2tot = zeros(N,1); qH2Otot = zeros(N,1);
qCO2tot(1) = qCO2f(1) + qCO2s(1);
qH2Otot(1) = qH2Of(1) + qH2Ob(1);

% ads_end_prev reference for residual floor
q_ads_end_prev_CO2 = qCO2tot(1);

% diags
pCO2_inventory = zeros(N,1);
pCO2_drive     = zeros(N,1);
pH2O_drv = zeros(N,1);
Pdry_v   = zeros(N,1);
aw_gas   = zeros(N,1);
aw_eff   = zeros(N,1);

% kinetics logging for diagnostics
kCO2_fast_vec = nan(N,1);

% peaks for H2O residual floors (keep)
q_peak_H2Ob = qH2Ob(1);
fresid_CO2_active  = false(N,1);
fresid_H2Ob_active = false(N,1);

% adsorption throughput 
% flow_mode = getfield_safe(opt, {'flow_mode'}, 'Aflow');
[Q_ads, ~] = get_Q_ads_and_A(ctrl.u_feed_m_s, design, opt.flow_mode, opt);

% ptraj flow filters
Qout_filt = 0.0; Qin_filt = 0.0;
Qout_vec = zeros(N,1); Qin_vec = zeros(N,1); yin_vec = zeros(N,1);

% dry moles + residual logging
n_dry_vec  = zeros(N,1);
nd_res_all = zeros(N,1);
nd_res_ptr = zeros(N,1);
is_ptraj   = false(N,1);

% adsorption cumulative cap state
ads_cap = struct('active',false,'N_avail_total',0.0,'N_sorb_done',0.0,'N_init_avail',0.0,'N_in_total',0.0, ...
                 'fac_min',1.0,'fac_count',0);

% stage-aware n_dry target LPF state
n_dry_filt_prev = NaN;

% --- stage elapsed time tracker (lightweight) ---
sid_prev = uint8(drv.stage_id(1));
t_stage  = 0.0;
aw_ads0  = NaN;   % adsorption entry aw_eff snapshot


for k=1:N-1
    sid = uint8(drv.stage_id(k));

    % ------------------------------------------------------------------
    % Adsorption-end reference for CO2 residual floor (S2):
    % Update immediately on Ads -> regen transition so the residual floor
    % in the first regen step scales to the current-cycle adsorption end.
    % ------------------------------------------------------------------
    sid_prev0 = uint8(drv.stage_id(max(k-1,1)));
    leaving_ads0 = (sid ~= 6) && (sid_prev0 == 6);
    if leaving_ads0
        q_ads_end_prev_CO2 = qCO2tot(k); % adsorption-end state at stage boundary
    end

    TK_cmd = drv.TK_cmd(k);
    RH_cmd = RHclamp(drv.RH_cmd(k));
    Pk     = drv.P_Pa(k);

    % 1) thermal lag
    tauT = get_tauT_from_id(sid, opt.tauT_s);
    TKbed(k+1) = apply_thermal_lag_step(TKbed(k), TK_cmd, dt, tauT);

    % 2) consistent pH2O/Pdry/aw chain
    [pH2O_k, Pdry_k, aw_gas_k] = pH2O_Pdry_aw_chain(RH_cmd, TK_cmd, TKbed(k), Pk, opt);

    % 3) next-step target dry moles from prescribed P(t), RH(t), Tbed(t)
    Pk1    = drv.P_Pa(k+1);
    RH1    = RHclamp(drv.RH_cmd(k+1));
    TKcmd1 = drv.TK_cmd(k+1);
    [~, Pdry_1, ~] = pH2O_Pdry_aw_chain(RH1, TKcmd1, TKbed(k+1), Pk1, opt);

    n_dry_k = max(Pdry_k * design.Vg_m3 / (phys.eq.R * max(TKbed(k),1.0)), 1e-12);
    n_raw_1 = max(Pdry_1 * design.Vg_m3 / (phys.eq.R * max(TKbed(k+1),1.0)), 1e-12);
    n_dry_vec(k) = n_dry_k;

    % ===== (B) stage-aware n_dry_target filtering =====
    if isnan(n_dry_filt_prev), n_dry_filt_prev = n_dry_k; end
    sid_next = uint8(drv.stage_id(k+1));
    stage_changed = (sid_next ~= sid);
    if stage_changed
    t_stage = 0.0;
    else
    t_stage = t_stage + dt;
    end

% capture adsorption-entry aw_eff "start" value later after aw_eff computed

    if stage_changed && opt.ptraj.enable_smoothing && isfield(opt.ptraj,'reset_flow_filters_on_stage_change') ...
            && opt.ptraj.reset_flow_filters_on_stage_change && (sid ~= 6) && (sid_next ~= 6)
        Qin_filt = 0.0; Qout_filt = 0.0;
    end

    n_dry_raw_next = n_raw_1;

    tau_ndry = opt.ptraj.tau_ndry_s;
    if isstruct(tau_ndry)
        switch double(sid_next)
            case 1, tau_use = tau_ndry.vac;
            case 2, tau_use = tau_ndry.heat;
            case 3, tau_use = tau_ndry.des;
            case 4, tau_use = tau_ndry.cool;
            case 5, tau_use = tau_ndry.press;
            case 6, tau_use = tau_ndry.ads;
            otherwise, tau_use = tau_ndry.des;
        end
    elseif isnumeric(tau_ndry) && numel(tau_ndry)==6
        tau_use = tau_ndry(double(sid_next));
    else
        tau_use = tau_ndry;
    end
    tau_use = max(tau_use, 0);

    if isfield(opt.ptraj,'enable_ndry_filter') && opt.ptraj.enable_ndry_filter
        if isfield(opt.ptraj,'ndry_reset_on_stage_change') && opt.ptraj.ndry_reset_on_stage_change && stage_changed
            n_dry_filt_next = n_dry_raw_next;
            n_dry_filt_prev = n_dry_k;
        else
            if tau_use <= 0
                n_dry_filt_next = n_dry_raw_next;
            else
                n_dry_filt_next = lpf_step(n_dry_filt_prev, n_dry_raw_next, dt, tau_use);
            end
        end
    else
        n_dry_filt_next = n_dry_raw_next;
    end

    ndot_tgt = (n_dry_filt_next - n_dry_filt_prev) / dt;
    n_dry_filt_prev = n_dry_filt_next;
    % ===== END (B) =====

    aw_gas(k) = aw_gas_k;

    % 4) CO2-effects aw
    aw_eff_k = aw_for_co2_effects(sid, TKbed(k), aw_gas_k, qH2Of(k), qH2Ob(k), phys.iso, opt);
    % --- Adsorption-stage aw_eff shaping for CO2 (breakthrough surrogate) ---
if isfield(opt,'ads_shape') && isfield(opt.ads_shape,'enable') && opt.ads_shape.enable && (sid==6)
    % 6=ads
    % choose inlet target: RH_cmd is usually what you want for "forced_only"
    if opt.ads_shape.use_RHcmd
        aw_in = clamp(drv.RH_cmd(k), 0.0, 0.999);
    else
        aw_in = clamp(amb.RH_in, 0.0, 0.999);
    end

    % latch aw_ads0 at adsorption entry
    if stage_changed || isnan(aw_ads0)
        aw_ads0 = aw_eff_k;
    end

    tau = max(opt.ads_shape.tau_aw_ads_s, 1e-6);

    % exponential approach: aw(t)=aw_in - (aw_in-aw0)*exp(-t/tau)
    aw_eff_k = aw_in - (aw_in - aw_ads0)*exp(-t_stage/tau);

    % safety clamp
    aw_eff_k = clamp(aw_eff_k, 0.0, 0.999);
end
    % optional LPF stabilization to suppress stage-boundary spikes
    if isfield(opt,'stab') && isfield(opt.stab,'enable') && opt.stab.enable && ...
            isfield(opt.stab,'aw_lpf_enable') && opt.stab.aw_lpf_enable
        if isnan(aw_eff_filt_prev), aw_eff_filt_prev = aw_eff_k; end
        tau_aw = max(opt.stab.aw_tau_s, 1e-6);
        a_aw = 1 - exp(-dt/tau_aw);
        if stage_changed && isfield(opt.stab,'aw_stagejump_reset') && opt.stab.aw_stagejump_reset
            aw_eff_filt_prev = aw_eff_k;
        else
            aw_eff_filt_prev = aw_eff_filt_prev + a_aw*(aw_eff_k - aw_eff_filt_prev);
        end
        aw_eff_k = clamp(aw_eff_filt_prev, 0.0, 0.999);
    end
    aw_eff(k) = aw_eff_k;
    sid_prev  = sid;

    % 5) pCO2 inventory and drive
    y_k  = clamp(yCO2g(k), opt.yco2.y_floor, opt.yco2.y_max);
    p_inv = y_k * Pdry_k;
    p_inl = amb.yCO2_air * Pdry_k;
    p_drv_raw = select_pCO2_drive(sid, p_inv, p_inl, opt);
    p_drv = p_drv_raw;
    % optional LPF stabilization to suppress stage-boundary spikes
    if isfield(opt,'stab') && isfield(opt.stab,'enable') && opt.stab.enable && ...
            isfield(opt.stab,'pdrv_lpf_enable') && opt.stab.pdrv_lpf_enable
        if isnan(pdrv_filt_prev), pdrv_filt_prev = p_drv; end
        tau_p = max(opt.stab.pdrv_tau_s, 1e-6);
        a_p = 1 - exp(-dt/tau_p);
        if stage_changed && isfield(opt.stab,'pdrv_stagejump_reset') && opt.stab.pdrv_stagejump_reset
            pdrv_filt_prev = p_drv;
        else
            pdrv_filt_prev = pdrv_filt_prev + a_p*(p_drv - pdrv_filt_prev);
        end
        p_drv = max(pdrv_filt_prev, 0.0);
    end
    pCO2_inventory(k) = p_inv;
    pCO2_drive(k)     = p_drv;

    % 6) equilibria
    qH2O_eq_tot = iso_GAB_H2O_noteS2(TKbed(k), aw_gas_k, phys.iso, opt.aw.max);

    % Optional empirical adsorption-only H2O capacity correction (lightweight)
    % - Helps match high-T/high-RH electricity (water inventory -> heat duty & vacuum load)
    % - Only applied in Adsorption stage to avoid "moving the goalposts" during regeneration.
    if sid==6 && isfield(opt,'h2o_T_ads_scale') && isfield(opt.h2o_T_ads_scale,'enable') && opt.h2o_T_ads_scale.enable
        TrefC = 20; if isfield(opt.h2o_T_ads_scale,'TrefC'); TrefC = opt.h2o_T_ads_scale.TrefC; end
        dTspan = 30; if isfield(opt.h2o_T_ads_scale,'dT_spanC'); dTspan = opt.h2o_T_ads_scale.dT_spanC; end
        kgain  = 0.0; if isfield(opt.h2o_T_ads_scale,'k_gain'); kgain = opt.h2o_T_ads_scale.k_gain; end
        Tc = TKbed(k) - 273.15;
        s = (Tc - TrefC)/max(dTspan, 1e-6); s = min(max(s,0),1);
        facH = 1.0 + kgain*s;
        qH2O_eq_tot = qH2O_eq_tot * facH;
    end
    if opt.h2o_2pool.enable
        fb_eq = h2o_fbound_eq(TKbed(k), aw_gas_k, sid, opt.h2o_2pool.fbound);
        qH2OeqB = fb_eq * qH2O_eq_tot;
        qH2OeqF = max(qH2O_eq_tot - qH2OeqB, 0);
    else
        qH2OeqF = qH2O_eq_tot; qH2OeqB = 0;
    end

    qCO2eqTot = iso_mech_CO2_noteS2(TKbed(k), p_drv, qH2Otot(k), phys.iso);
    qCO2eqTot = apply_CO2_dry_scale(qCO2eqTot, aw_eff_k, sid, opt);

    qCO2eqTot_vec(k) = qCO2eqTot;
    if opt.co2_2pool.enable
        fslow_eq = co2_fslow_eq(TKbed(k), aw_eff_k, sid, opt.co2_2pool.fslow);
        qCO2eqS  = fslow_eq * qCO2eqTot;
        qCO2eqF  = max(qCO2eqTot - qCO2eqS, 0);
    else
        qCO2eqF = qCO2eqTot; qCO2eqS = 0;
    end

    % 7) kinetics
    [kCO2_fast, kH2O_fast] = kinetics_LDF(TKbed(k), sid, pCO2_drive(k), aw_eff_k, phys.kin);
    % adsorption-only k shaping (OFF by default): fast early -> relax later
    if isfield(opt,'ads_kinshape') && isfield(opt.ads_kinshape,'enable') && opt.ads_kinshape.enable && (sid==6)
        tauK    = max(opt.ads_kinshape.tau_k_s, 1e-6);
        kboost0 = max(opt.ads_kinshape.kboost0, 0.0);
        fac = 1.0 + kboost0 * exp(-t_stage / tauK);
        % optional clamp to avoid extreme stiffness
        if isfield(opt.ads_kinshape,'fac_max') && ~isempty(opt.ads_kinshape.fac_max)
            fac = min(fac, opt.ads_kinshape.fac_max);
        end
        kCO2_fast = kCO2_fast * fac;
    end
% A_mult -> kinetics scaling (minimal, mild)
kCO2_fast = kCO2_fast * (design.A_mult.^0.5);
kH2O_fast = kH2O_fast * (design.A_mult.^0.5);

    kCO2_fast_vec(k) = kCO2_fast;

    if opt.co2_2pool.enable
        kCO2_slow = clamp(opt.co2_2pool.kslow.ratio * kCO2_fast, opt.co2_2pool.kslow.kmin, opt.co2_2pool.kslow.kmax);
    else
        kCO2_slow = 0;
    end

    if opt.h2o_2pool.enable
        kH2O_bound = clamp(opt.h2o_2pool.kslow.ratio * kH2O_fast, opt.h2o_2pool.kslow.kmin, opt.h2o_2pool.kslow.kmax);
    else
        kH2O_bound = 0;
    end

    % --- Saturation-aware LDF (Adsorption only) ---
alpha_sat = getfield_safe(opt, {'co2_sat','alpha'}, 0.80);  % 0.5~1.0 추천
sat_floor = getfield_safe(opt, {'co2_sat','floor'}, 0.05);  % 과도한 slowdown 방지

kF_eff = kCO2_fast;
kS_eff = kCO2_slow;

if sid == 6
    qeqT = max(qCO2eqTot, 1e-12);           % total eq 기준
    qT   = max(qCO2tot(k), 0.0);            % 현재 total loading
    sat  = 1.0 - min(qT/qeqT, 1.0);         % 1→0
    sat  = max(sat, sat_floor);             % floor
    fac  = sat.^alpha_sat;                  % saturation factor

    kF_eff = kCO2_fast * fac;
    kS_eff = kCO2_slow * fac;               % slow도 같이 곡선화(원치 않으면 이 줄 삭제)
end

aF  = exp(-kF_eff*dt);
aS  = exp(-kS_eff*dt);

qCO2f_next = qCO2eqF + (qCO2f(k) - qCO2eqF)*aF;
qCO2s_next = qCO2eqS + (qCO2s(k) - qCO2eqS)*aS;
    
    aWf = exp(-kH2O_fast*dt);
    aWb = exp(-kH2O_bound*dt);
    
    qH2Of_next = qH2OeqF + (qH2Of(k) - qH2OeqF)*aWf;
    qH2Ob_next = qH2OeqB + (qH2Ob(k) - qH2OeqB)*aWb;

    % Residual floors
    if opt.co2_2pool.enable && opt.co2_2pool.resid.enable
        [qCO2s_next, isActiveC] = apply_CO2_residual_floor_slow( ...
            qCO2s_next, q_ads_end_prev_CO2, TKbed(k), aw_eff_k, sid, ctrl, opt.co2_2pool.resid);
        fresid_CO2_active(k) = isActiveC;
    end
    if opt.h2o_2pool.enable && opt.h2o_2pool.resid.enable
        [qH2Ob_next, isActiveW] = apply_H2O_bound_residual_floor( ...
            qH2Ob_next, q_peak_H2Ob, TKbed(k), sid, ctrl, opt.h2o_2pool.resid);
        fresid_H2Ob_active(k) = isActiveW;
    end

    % clamp
    qCO2f_next = max(qCO2f_next, 0);
    qCO2s_next = max(qCO2s_next, 0);

    qH2Of_next = clamp(qH2Of_next, 0, 80);
    qH2Ob_next = clamp(qH2Ob_next, 0, 80);

    qCO2tot_next = qCO2f_next + qCO2s_next;
    qH2Otot_next = qH2Of_next + qH2Ob_next;

    % enforce qinf cap
    qinf = phys.iso.mech.S_cap * phys.iso.mech.qinf0 * exp(phys.iso.mech.chi*(1 - TKbed(k)/phys.iso.mech.T0));
    if qCO2tot_next > qinf
        excess = qCO2tot_next - qinf;
        qCO2f_next = max(qCO2f_next - excess, 0);
        qCO2tot_next = min(qCO2f_next + qCO2s_next, qinf);
    end

    % Adsorption cumulative supply cap
    entering_ads = (sid==6) && (sid_prev0~=6);
    leaving_ads  = leaving_ads0;

    if isfield(opt,'ads_supply') && isfield(opt.ads_supply,'enable') && opt.ads_supply.enable && entering_ads
        nCO2_g_k = y_k * n_dry_k;
        if isfield(opt,'ads_supply') && isfield(opt.ads_supply,'keep_y_floor_inventory') && opt.ads_supply.keep_y_floor_inventory
            nCO2_min = opt.yco2.y_floor * n_dry_k;
        else
            nCO2_min = 0.0;
        end
        N_init = max(nCO2_g_k - nCO2_min, 0.0);

        ads_cap.active = true;
        ads_cap.N_init_avail  = N_init;
        ads_cap.N_in_total    = 0.0;
        ads_cap.N_avail_total = N_init;
        ads_cap.N_sorb_done   = 0.0;
        ads_cap.fac_min       = 1.0;
        ads_cap.fac_count     = 0;
    end

    if isfield(opt,'ads_supply') && isfield(opt.ads_supply,'enable') && opt.ads_supply.enable && sid==6
        c_in_CO2 = (amb.yCO2_air * Pdry_k) / (phys.eq.R * max(TKbed(k),1.0));
        N_in = Q_ads * c_in_CO2 * dt;
        ads_cap.N_in_total    = ads_cap.N_in_total + N_in;
        ads_cap.N_avail_total = ads_cap.N_avail_total + N_in;

        dn_req = design.msorb_kg * (qCO2tot_next - qCO2tot(k));
        if dn_req > 0
            N_remaining = ads_cap.N_avail_total - ads_cap.N_sorb_done;
            N_remaining = max(N_remaining, 0);

            if dn_req > N_remaining
                fac = N_remaining / max(dn_req, 1e-18);

                df = qCO2f_next - qCO2f(k);
                ds = qCO2s_next - qCO2s(k);

                qCO2f_next = qCO2f(k) + fac*df;
                qCO2s_next = qCO2s(k) + fac*ds;

                qCO2f_next = max(qCO2f_next, 0);
                qCO2s_next = max(qCO2s_next, 0);
                qCO2tot_next = qCO2f_next + qCO2s_next;

                ads_cap.fac_min   = min(ads_cap.fac_min, fac);
                ads_cap.fac_count = ads_cap.fac_count + 1;

                dn_req = fac*dn_req;
            end
            ads_cap.N_sorb_done = ads_cap.N_sorb_done + max(dn_req,0);
        end
    end

    if isfield(opt,'ads_supply') && isfield(opt.ads_supply,'enable') && opt.ads_supply.enable && leaving_ads
        ads_cap.active = false;
    end

    % 8) ptraj inversion flows
    dn_sorb_CO2 = design.msorb_kg * (qCO2tot_next - qCO2tot(k));
    dn_sorb_dot = dn_sorb_CO2 / dt;

    % (C) physical rate limit on required dry-moles rate
    if isfield(opt.ptraj,'enable_ndot_rate_limit') && opt.ptraj.enable_ndot_rate_limit && (sid ~= 6)
        Qlim = opt.ptraj.Qmax_mult_Vg_dt * design.Vg_m3 / dt;

        c_in_dry  = max(Pdry_k / (phys.eq.R * max(TKbed(k),1.0)), 1e-12);
        c_out_dry = max(n_dry_k / max(design.Vg_m3,1e-12),       1e-12);

        if sid == 5
            ndot_req_max = Qlim * c_in_dry;
        else
            ndot_req_max = Qlim * c_out_dry;
        end

        ndot_req = ndot_tgt + dn_sorb_dot;

        lb = -ndot_req_max; ub = +ndot_req_max;
        if sid == 5
            lb = 0;
        else
            ub = 0;
        end
        ndot_req_clip = clamp(ndot_req, lb, ub);
        ndot_tgt = ndot_req_clip - dn_sorb_dot;
    end

    [Qin_raw, Qout_raw, y_in] = flows_from_ptraj( ...
        sid, Q_ads, amb.yCO2_air, n_dry_k, ndot_tgt, dn_sorb_dot, ...
        Pdry_k, TKbed(k), phys.eq.R, design.Vg_m3);

    Qlim = opt.ptraj.Qmax_mult_Vg_dt * design.Vg_m3 / dt;
    Qin_raw  = clamp(Qin_raw,  0, Qlim);
    Qout_raw = clamp(Qout_raw, 0, Qlim);

    if opt.ptraj.enable_smoothing
        a = 1 - exp(-dt/max(opt.ptraj.tau_flow_ramp_s,1e-9));
        Qin_filt  = Qin_filt  + a*(Qin_raw  - Qin_filt);
        Qout_filt = Qout_filt + a*(Qout_raw - Qout_filt);
        Qin  = Qin_filt;
        Qout = Qout_filt;
    else
        Qin  = Qin_raw;
        Qout = Qout_raw;
    end

    Qin_vec(k)  = Qin;
    Qout_vec(k) = Qout;
    yin_vec(k)  = y_in;

    if opt.ptraj.log_ndry_residual
        c_in_dry  = max(Pdry_k / (phys.eq.R * max(TKbed(k),1.0)), 1e-12);
        c_out_dry = max(n_dry_k / max(design.Vg_m3,1e-12), 1e-12);
        bal = Qin*c_in_dry - Qout*c_out_dry - dn_sorb_dot;
        res = ndot_tgt - bal;
        nd_res_all(k) = res;
        if sid ~= 6
            nd_res_ptr(k) = res;
            is_ptraj(k) = true;
        end
    end

    % 9) Update yCO2g
    yCO2g_next = update_yCO2g_from_flows( ...
        yCO2g(k), n_dry_k, n_raw_1, ...
        Qin, Qout, y_in, ...
        Pdry_k, TKbed(k), design.Vg_m3, phys.eq.R, ...
        dn_sorb_CO2, dt, opt.yco2);

    % write states
    qCO2f(k+1)   = qCO2f_next;
    qCO2s(k+1)   = qCO2s_next;
    qCO2tot(k+1) = qCO2tot_next;

    qH2Of(k+1)   = qH2Of_next;
    qH2Ob(k+1)   = qH2Ob_next;
    qH2Otot(k+1) = qH2Otot_next;

    yCO2g(k+1)   = yCO2g_next;

    % peaks
    q_peak_H2Ob = max(q_peak_H2Ob, qH2Ob_next);

    pH2O_drv(k) = pH2O_k;
    Pdry_v(k)   = Pdry_k;
end

% last-point diags
Pk_end = drv.P_Pa(end);
RH_end = RHclamp(drv.RH_cmd(end));
TKcmd_end = drv.TK_cmd(end);
[pH2O_end, Pdry_end, aw_end] = pH2O_Pdry_aw_chain(RH_end, TKcmd_end, TKbed(end), Pk_end, opt);

y_end = clamp(yCO2g(end), opt.yco2.y_floor, opt.yco2.y_max);
pCO2_inventory(end) = y_end * Pdry_end;
pCO2_drive(end)     = select_pCO2_drive(uint8(drv.stage_id(end)), pCO2_inventory(end), amb.yCO2_air*Pdry_end, opt);

pH2O_drv(end) = pH2O_end;
Pdry_v(end)   = Pdry_end;
aw_gas(end)   = aw_end;
aw_eff(end)   = aw_for_co2_effects(uint8(drv.stage_id(end)), TKbed(end), aw_end, qH2Of(end), qH2Ob(end), phys.iso, opt);

n_dry_vec(end)= max(Pdry_end * design.Vg_m3 / (phys.eq.R * max(TKbed(end),1.0)), 1e-12);
kCO2_fast_vec(end) = kCO2_fast_vec(max(end-1,1));

qCO2eqTot_vec(end) = qCO2eqTot_vec(max(end-1,1));
% cycle summary indices
idx_ads_end = find(drv.stage_id==6, 1, 'last'); if isempty(idx_ads_end), idx_ads_end=N; end
idx_des_end = find(drv.stage_id==3, 1, 'last'); if isempty(idx_des_end), idx_des_end=1; end

q_ads_end = qCO2tot(idx_ads_end);
q_des_end = qCO2tot(idx_des_end);
dq_raw = max(q_ads_end - q_des_end, 0);

% adsorption diagnostics
idx_ads = find(drv.stage_id==6);
if ~isempty(idx_ads)
    mean_aw_eff_ads  = mean(aw_eff(idx_ads));
    mean_pCO2drv_ads = mean(pCO2_drive(idx_ads));
else
    mean_aw_eff_ads  = NaN;
    mean_pCO2drv_ads = NaN;
end

% Equilibrium approach at adsorption end (achieved / equilibrium)
TK_end  = TKbed(idx_ads_end);
% For qH2O equilibrium, use gas-phase activity on the bed (aw_gas), not the
% CO2-effective aw (aw_eff) which may be bound-water-masked.
aw_end2 = aw_gas(idx_ads_end);
p_endPa = pCO2_drive(idx_ads_end);

qH2O_eq_end = iso_GAB_H2O_noteS2(TK_end, aw_end2, phys.iso, opt.aw.max);
qCO2_eq_end = iso_mech_CO2_noteS2(TK_end, p_endPa, qH2O_eq_end, phys.iso);
ads_eq_approach = clamp01(q_ads_end / max(qCO2_eq_end, 1e-12));

N_init_avail = 0.0;
N_in_total   = 0.0;

if ~isempty(idx_ads)
    k0 = idx_ads(1);
    P0 = Pdry_v(k0);
    y0 = clamp(yCO2g(k0), opt.yco2.y_floor, opt.yco2.y_max);
    n_dry0 = max(P0 * design.Vg_m3 / (phys.eq.R * max(TKbed(k0),1.0)), 1e-12);
    nCO2_g0 = y0 * n_dry0;

    if isfield(opt,'ads_supply') && isfield(opt.ads_supply,'keep_y_floor_inventory') && opt.ads_supply.keep_y_floor_inventory
        nCO2_min0 = opt.yco2.y_floor * n_dry0;
    else
        nCO2_min0 = 0.0;
    end
    N_init_avail = max(nCO2_g0 - nCO2_min0, 0.0);

    P_ads_vec  = drv.P_Pa(idx_ads);
    RH_ads_vec = RHclamp(drv.RH_cmd(idx_ads));
    TKcmd_ads  = drv.TK_cmd(idx_ads);
    TKbed_ads  = TKbed(idx_ads);

    Pdry_ads = zeros(numel(idx_ads),1);
    for j=1:numel(idx_ads)
        [~, Pdry_ads(j), ~] = pH2O_Pdry_aw_chain(RH_ads_vec(j), TKcmd_ads(j), TKbed_ads(j), P_ads_vec(j), opt);
    end
    pCO2_in_ads = amb.yCO2_air .* Pdry_ads;
    c_in_ads = pCO2_in_ads ./ (phys.eq.R .* max(TKbed_ads,1.0));
    n_in_dot = Q_ads .* c_in_ads;
    N_in_total = trapz(t(idx_ads), n_in_dot);
end

n_work = design.msorb_kg * dq_raw;
n_sup_total = N_init_avail + N_in_total;
n_cap  = min(max(n_work,0), max(n_sup_total,0));
cap_ratio = n_cap / max(n_work, 1e-12);

cyc = struct();
cyc.qCO2_ads_end = q_ads_end;
cyc.qCO2_des_end = q_des_end;
cyc.dqCO2_raw = dq_raw;
cyc.nCO2_working_raw = n_work;
cyc.nCO2_supply_max  = n_sup_total;
cyc.nCO2_captured_eff= n_cap;
cyc.cap_ratio        = cap_ratio;

cyc.mean_aw_eff_ads       = mean_aw_eff_ads;
cyc.mean_pCO2drv_ads_Pa   = mean_pCO2drv_ads;
cyc.ads_eq_approach       = ads_eq_approach;
cyc.qCO2_eq_ads_end_molkg = qCO2_eq_end;

cyc.fresid_CO2_active_pct       = 100*mean(fresid_CO2_active);
cyc.fresid_H2O_bound_resid_pct  = 100*mean(fresid_H2Ob_active);

cyc.ndry_min = min(n_dry_vec);
cyc.ndry_max = max(n_dry_vec);

nd_res_all(end) = nd_res_all(max(end-1,1));
nd_res_ptr(end) = nd_res_ptr(max(end-1,1));

cyc.ndry_res_max_abs_all  = max(abs(nd_res_all));
cyc.ndry_res_mean_abs_all = mean(abs(nd_res_all));
if any(is_ptraj)
    cyc.ndry_res_max_abs_ptraj  = max(abs(nd_res_ptr(is_ptraj)));
    cyc.ndry_res_mean_abs_ptraj = mean(abs(nd_res_ptr(is_ptraj)));
else
    cyc.ndry_res_max_abs_ptraj  = NaN;
    cyc.ndry_res_mean_abs_ptraj = NaN;
end

% PPI
ppi = ppi_calc_sendi_eq23_34_raw(amb, ctrl, phys, design, drv, ...
    qCO2tot, qH2Otot, TKbed, yCO2g, pCO2_inventory, pCO2_drive, ...
    pH2O_drv, Pdry_v, Qin_vec, Qout_vec, yin_vec, cyc, opt, Q_ads);

% pack
out = struct();
out.drv = drv;

out.prof = struct();
out.prof.time_s = t;

out.prof.qCO2_total = qCO2tot;
out.prof.qCO2_eq_total = qCO2eqTot_vec;
out.prof.qCO2_fast  = qCO2f;
out.prof.qCO2_slow  = qCO2s;

out.prof.qH2O_total = qH2Otot;
out.prof.qH2O_fast  = qH2Of;
out.prof.qH2O_bound = qH2Ob;

out.prof.yCO2g = yCO2g;

out.prof.TK_bed = TKbed;
out.prof.TC_bed = TKbed - 273.15;

out.prof.pCO2_inventory_Pa = pCO2_inventory;
out.prof.pCO2_drive_Pa     = pCO2_drive;

out.prof.pH2O_drv_Pa = pH2O_drv;
out.prof.Pdry_Pa     = Pdry_v;

out.prof.aw_gas = aw_gas;
out.prof.aw_eff = aw_eff;
out.prof.RH_cmd = drv.RH_cmd(:);
out.prof.TK_cmd = drv.TK_cmd(:);

out.prof.Qin_m3s  = Qin_vec;
out.prof.Qout_m3s = Qout_vec;
out.prof.yin_CO2  = yin_vec;

out.prof.ndry_res_mol_s_all   = nd_res_all;
out.prof.ndry_res_mol_s_ptraj = nd_res_ptr;
out.prof.is_ptraj_stage       = is_ptraj;

out.prof.kCO2_fast = kCO2_fast_vec;

    % ---- Publication-quality adsorption trace (post-processing only) ----
    % Core state integration is untouched; this only adds "pretty" vectors for figures.
    out.prof.idx_ads = idx_ads;
    out.prof.idx_ads_end = idx_ads_end;

   if isfield(opt,'pub_trace') && isfield(opt.pub_trace,'enable') && opt.pub_trace.enable && ~isempty(idx_ads)
    [qCO2_pub_vec, ads_app_pub_vec] = make_pub_ads_trace( ...
        t, qCO2tot, qCO2eqTot_vec, kCO2_fast_vec, idx_ads, opt.pub_trace);
    out.prof.qCO2_total_pub = qCO2_pub_vec;
    out.prof.ads_eq_approach_pub = ads_app_pub_vec;
    else
    out.prof.qCO2_total_pub = qCO2tot;
    out.prof.ads_eq_approach_pub = nan(size(qCO2tot));
    end


out.cyc = cyc;
out.ppi = ppi;

out.meta = struct();
out.meta.S_cap = phys.iso.mech.S_cap;
out.meta.iso   = phys.iso;

out.x_end = [qCO2f(end); qCO2s(end); qH2Of(end); qH2Ob(end); yCO2g(end); TKbed(end)];
end


%% =========================================================================
% pH2O/Pdry/aw chain helper
% =========================================================================
function [pH2O, Pdry, aw_gas] = pH2O_Pdry_aw_chain(RH_cmd, TK_cmd, TK_bed, P, opt)
RH_cmd = RHclamp(RH_cmd);
P = max(P, 200.0);

mode = lower(strtrim(string(opt.aw.chain_mode)));
switch mode
    case "bed_consistent"
        pH2O = clamp(RH_cmd * psat_H2O_Pa(TK_bed), 0, 0.98*P);
        Pdry = max(P - pH2O, 1.0);
        aw_gas = clamp(pH2O / max(psat_H2O_Pa(TK_bed),1.0), 0, opt.aw.max);

    case "cmd_based"
        pH2O = clamp(RH_cmd * psat_H2O_Pa(TK_cmd), 0, 0.98*P);
        Pdry = max(P - pH2O, 1.0);

        aw_raw = pH2O / max(psat_H2O_Pa(TK_bed), 1.0);

        pol = lower(strtrim(string(opt.aw.cmd_based_cap_policy)));
        switch pol
            case "none"
                aw_gas = aw_raw;
            case "cap_to_rhcmd"
                aw_gas = min(aw_raw, RH_cmd);
            otherwise
                aw_gas = min(aw_raw, 1.0);
        end
        aw_gas = clamp(aw_gas, 0, opt.aw.max);

    otherwise
        error('opt.aw.chain_mode must be bed_consistent or cmd_based.');
end
end

function RH = RHclamp(RH)
RH = clamp(RH, 0, 0.999);
end


%% =========================================================================
% pCO2 drive selector
% =========================================================================
function p = select_pCO2_drive(stage_id, p_inv, p_inl, opt)
% Stage-aware pCO2 drive selector.
% Patch intent:
% - Adsorption (sid==6): inlet-dominant with a hard floor to prevent 0D inventory artifacts
%   from collapsing pCO2 drive in dry cases (Prod under-prediction cluster).
% - Regen: preserve existing behavior (default inventory).

sid = uint8(stage_id);

% ---- Safe defaults if opt.pco2 fields are missing ----
ads_mode = "inlet";          % default adsorption mode (defensive)
reg_mode = "inventory";      % default regen mode
blend_a  = 0.5;
p_floor_frac = 1.0;          % adsorption: enforce p >= p_floor_frac*p_inl (default 1.0 = inlet lock)

if isfield(opt,'pco2')
    if isfield(opt.pco2,'ads_mode')  && strlength(string(opt.pco2.ads_mode))>0
        ads_mode = lower(strtrim(string(opt.pco2.ads_mode)));
    end
    if isfield(opt.pco2,'regen_mode') && strlength(string(opt.pco2.regen_mode))>0
        reg_mode = lower(strtrim(string(opt.pco2.regen_mode)));
    end
    if isfield(opt.pco2,'blend_alpha_inventory')
        blend_a = clamp(opt.pco2.blend_alpha_inventory, 0, 1);
    end
    if isfield(opt.pco2,'ads_p_floor_frac')
        p_floor_frac = max(opt.pco2.ads_p_floor_frac, 0);
    end
end

% ---- Select mode ----
if sid == 6
    mode = ads_mode;
else
    mode = reg_mode;
end

switch mode
    case "inlet"
        p = p_inl;
    case "blend"
        p = blend_a*p_inv + (1-blend_a)*p_inl;
    otherwise % "inventory"
        p = p_inv;
end

% ---- Adsorption-only: enforce inlet dominance floor ----
if sid == 6
    p = max(p, p_floor_frac * p_inl);  % default p_floor_frac=1.0
end

p = max(p, 0);
end

%% =========================================================================
% ptraj inversion: infer Qin/Qout to satisfy n_dry(t)
% =========================================================================
function [Qin, Qout, y_in] = flows_from_ptraj( ...
    sid, Q_ads, yCO2_air, n_dry, ndot_target, dn_sorb_CO2_dot, Pdry, TK, R, Vg)

c_out_dry = max(n_dry / max(Vg,1e-12), 1e-12);
c_in_dry  = max(Pdry / (R * max(TK,1.0)), 1e-12);

Qin = 0.0; Qout = 0.0; y_in = 0.0;

switch double(sid)
    case 6 % Adsorption
        Qin  = Q_ads;
        Qout = Q_ads;
        y_in = yCO2_air;

    case 5 % Pressurization: infer Qin
        Qin  = (ndot_target + dn_sorb_CO2_dot) / c_in_dry;
        Qout = 0.0;
        y_in = yCO2_air;

    otherwise % Regen-side: infer Qout
        Qout = -(ndot_target + dn_sorb_CO2_dot) / c_out_dry;
        Qin  = 0.0;
        y_in = 0.0;
end

Qin  = max(Qin, 0);
Qout = max(Qout, 0);
end


%% =========================================================================
% yCO2 update from flows (semi-implicit outflow)
% =========================================================================
function y_next = update_yCO2g_from_flows( ...
    y_prev, n_dry_k, n_dry_next, Qin, Qout, y_in, Pdry, TK, Vg, R, dn_sorb_CO2, dt, yopt)

y_prev = clamp(y_prev, yopt.y_floor, yopt.y_max);
c_in_CO2 = (y_in * Pdry) / (R * max(TK,1.0));

nCO2_prev = y_prev * max(n_dry_k,1e-12);
rhs = nCO2_prev + (Qin*c_in_CO2)*dt - dn_sorb_CO2;

nCO2_mb = rhs / (1 + (Qout*dt)/max(Vg,1e-12));
nCO2_mb = max(nCO2_mb, yopt.y_floor*max(n_dry_k,1e-12));

nref = max(n_dry_next, 1e-12);
y_mb = clamp(nCO2_mb / nref, yopt.y_floor, yopt.y_max);
y_next = clamp(y_mb, yopt.y_floor, yopt.y_max);
% optional per-step delta-y cap to prevent inventory spikes
if isfield(yopt,'dy_cap') && isfield(yopt.dy_cap,'enable') && yopt.dy_cap.enable
    dymax = yopt.dy_cap.dy_max;
    if ~isempty(dymax) && dymax > 0
        dy = clamp(y_next - y_prev, -dymax, +dymax);
        y_next = clamp(y_prev + dy, yopt.y_floor, yopt.y_max);
    end
end
end


%% =========================================================================
% aw_eff logic for CO2-related effects
% =========================================================================
function aw = aw_for_co2_effects(sid, TK, aw_gas_on_bed, qH2Of, qH2Ob, iso, opt)
% Effective aw used in CO2 equilibrium/kinetics coupling.
% Patch intent:
% - Prevent bound-water (slow pool) from inflating CO2-side aw in dry cases.
% - Minimal physical model: blend gas-phase aw with GAB-inverted aw using free-water fraction.

aw_g = clamp(aw_gas_on_bed, 0, opt.aw.max);

qH2Of = max(qH2Of, 0);
qH2Ob = max(qH2Ob, 0);
qTot  = qH2Of + qH2Ob;

% GAB-inverted aw from total water loading (legacy behavior)
aw_gab = aw_from_GAB_bisect(TK, qTot, iso, opt.aw.max);

% Determine whether to apply bound-water masking (default ON when 2-pool is used)
use_mask = true;
if isfield(opt,'aw') && isfield(opt.aw,'co2_boundmask') && isfield(opt.aw.co2_boundmask,'enable')
    use_mask = logical(opt.aw.co2_boundmask.enable);
end

% Mode handling (keep legacy switch)
mode = "gab";
if isfield(opt,'aw') && isfield(opt.aw,'mode_for_co2_effects') && strlength(string(opt.aw.mode_for_co2_effects))>0
    mode = lower(string(opt.aw.mode_for_co2_effects));
end

switch mode
    case "use_aw_gas"
        aw = aw_g;
    otherwise
        if use_mask && isfield(opt,'h2o_2pool') && isfield(opt.h2o_2pool,'enable') && opt.h2o_2pool.enable
            epsq = 1e-9;
            f_free = qH2Of / max(qTot + epsq, epsq);   % [0..1]
            % Optional sharpening (default 1)
            pwr = 1.0;
            if isfield(opt.aw,'co2_boundmask') && isfield(opt.aw.co2_boundmask,'power')
                pwr = max(opt.aw.co2_boundmask.power, 0.25);
            end
            f_free = f_free^pwr;

            % Blend: aw = aw_gas + f_free*(aw_gab - aw_gas)
            aw = aw_g + f_free*(aw_gab - aw_g);
        else
            % Legacy: use total-water inverted aw
            aw = aw_gab;
        end
end

aw = clamp(aw, 0, opt.aw.max);
aw_g = clamp(aw_gas_on_bed, 0.0, 1.0);   % <-- ADD (gas-phase activity, clamped)
S = stage_ids();

% ======================================================================
% (A) site-level aw floor (independent layer) : ADS only
% ======================================================================
if isfield(opt,'aw_sitefloor') && isfield(opt.aw_sitefloor,'enable') && opt.aw_sitefloor.enable
    aw_g = clamp(aw_gas_on_bed, 0.0, 1.0);
    aw_sf = getfield_safe(opt.aw_sitefloor, {'aw_floor'}, 0.15);
    aw_sf_apply_below = getfield_safe(opt.aw_sitefloor, {'apply_below'}, 0.25);

    if (sid == S.Ads) && (aw_g < aw_sf_apply_below)
        aw = max(aw, aw_sf);
        aw = clamp(aw, 0, opt.aw.max);
    end
end

% ======================================================================
% (B) CO2-side effective activity floor (existing logic preserved) : ADS only
% ======================================================================
if isfield(opt,'co2_awfloor') && isfield(opt.co2_awfloor,'enable') && opt.co2_awfloor.enable
    if sid == S.Ads
        aw_apply_below = 0.25;
        if isfield(opt.co2_awfloor,'apply_below'); aw_apply_below = opt.co2_awfloor.apply_below; end

        if aw_g < aw_apply_below
            if isfield(opt.co2_awfloor,'value')
                aw_floor = opt.co2_awfloor.value;
            elseif isfield(opt.co2_awfloor,'aw_min')
                aw_floor = opt.co2_awfloor.aw_min;
            else
                aw_floor = 0.0;
            end
            aw = max(aw, aw_floor);
            aw = clamp(aw, 0, opt.aw.max);
        end
    end
end

end

function aw = aw_from_GAB_bisect(TK, qH2O_target, iso, aw_max)
qH2O_target = max(qH2O_target, 0);

q_lo = iso_GAB_H2O_noteS2(TK, 0.0, iso, aw_max);
q_hi = iso_GAB_H2O_noteS2(TK, aw_max, iso, aw_max);

if qH2O_target <= q_lo
    aw = 0.0; return;
elseif qH2O_target >= q_hi
    aw = aw_max; return;
end

lo = 0.0; hi = aw_max;
for it=1:18
    mid = 0.5*(lo+hi);
    q_mid = iso_GAB_H2O_noteS2(TK, mid, iso, aw_max);
    if q_mid < qH2O_target
        lo = mid;
    else
        hi = mid;
    end
end
aw = clamp(0.5*(lo+hi), 0, aw_max);
end


%% =========================================================================
% CO2 slow fraction eq (use aw_eff)
% =========================================================================
function f = co2_fslow_eq(~, aw_eff, stage_id, p) %#ok<INUSD>
if stage_id==6
    f0  = p.f0_ads;
    kaw = p.k_aw_ads;
else
    f0  = p.f0_regen;
    kaw = p.k_aw_reg;
end
f = f0 + kaw * aw_eff;
f = clamp(f, p.f_min, p.f_max);
end

%% =========================================================================
% H2O bound fraction eq (use aw_gas_on_bed)
% =========================================================================
function f = h2o_fbound_eq(~, aw_gas_on_bed, stage_id, p) %#ok<INUSD>
if stage_id==6
    f0  = p.f0_ads;
    kaw = p.k_aw_ads;
else
    f0  = p.f0_regen;
    kaw = p.k_aw_reg;
end
f = f0 + kaw * aw_gas_on_bed;
f = clamp(f, p.f_min, p.f_max);
end


%% =========================================================================
% Residual floors
% =========================================================================
function [q_slow_out, isActive] = apply_CO2_residual_floor_slow( ...
    q_slow_in, q_ads_end_prev_total, TK, aw_eff, stage_id, ctrl, p) %#ok<INUSD>

q_slow_out = q_slow_in;
isActive = false;

if ~any(uint8(stage_id) == uint8(p.apply_in_stage_ids(:))), return; end

mode = "ads_scaled";
if isfield(p,'mode') && strlength(string(p.mode))>0
    mode = lower(strtrim(string(p.mode)));
end

if mode == "ads_scaled"
    f0 = p.f_res;

    term_wet = 0.0;
    if isfield(p,'k_wet')
        n_wet = 1.0;
        if isfield(p,'n_wet'), n_wet = p.n_wet; end
        term_wet = p.k_wet * (max(aw_eff,0.0)^n_wet);
    end

    term_T = 0.0;
    if isfield(p,'k_T') && isfield(p,'Tref_K')
        term_T = p.k_T * max((p.Tref_K/max(TK,1.0)) - 1, 0);
    end

    f_resid = f0 + term_wet + term_T;

    fmin = 0.0; fmax = 1.0;
    if isfield(p,'f_min'), fmin = p.f_min; end
    if isfield(p,'f_max'), fmax = p.f_max; end
    f_resid = clamp(f_resid, fmin, fmax);

    q_min = max(p.q_abs_min, f_resid * max(q_ads_end_prev_total, 1e-12));

    if q_slow_out < q_min*(1 - p.tol_frac)
        q_slow_out = q_min;
        isActive = true;
    end
    return;
end

q_min = p.f0 * max(q_ads_end_prev_total,1e-12);
if q_slow_out < q_min*(1 - p.tol_frac)
    q_slow_out = q_min;
    isActive = true;
end
end

function [q_out, isActive] = apply_H2O_bound_residual_floor(q_in, q_peak, TK, sid, ctrl, p)
q_out = q_in;
isActive = false;

if ~any(uint8(sid) == uint8(p.apply_in_stage_ids(:))), return; end

treg = ctrl.time_heat_s + ctrl.time_des_s;

term_T    = p.k_T    * max((p.Tref_K/max(TK,1.0)) - 1, 0);
term_time = p.k_time * max((p.tref_s/max(treg,1.0)) - 1, 0);

f_resid = p.f0 + term_T + term_time;
f_resid = clamp(f_resid, 0, p.f_max);

q_min = max(p.q_abs_min, f_resid * max(q_peak, 1e-12));

if q_out < q_min*(1 - p.tol_frac)
    q_out = q_min;
    isActive = true;
end
end


%% =========================================================================
% ISOTHERMS (GAB + mechanistic)
% =========================================================================
function qH2O = iso_GAB_H2O_noteS2(TK, aw, iso, aw_max)
R = 8.314462618;
aw = clamp(aw, 0, aw_max);

E10p = -44.38*TK + 57220;
E1   = iso.gab.C - exp(iso.gab.D*TK);
E2m9 = iso.gab.F + iso.gab.G*TK;

Cgab = exp((E1   - E10p)/(R*TK));
Kg   = exp((E2m9 - E10p)/(R*TK));

den = (1 - Kg*aw) .* (1 + (Cgab - 1).*Kg.*aw);
den = max(den, 1e-12);

qH2O = iso.gab.qm .* (Kg .* Cgab .* aw) ./ den;
qH2O = max(qH2O, 0);
end

function qCO2 = iso_mech_CO2_noteS2(TK, pCO2_Pa, qH2O, iso)
R = 8.314462618;

qinf = iso.mech.S_cap * iso.mech.qinf0 * exp(iso.mech.chi*(1 - TK/iso.mech.T0));

qH2O_eff = max(qH2O, 1e-9);
w = exp(-iso.mech.Aexp ./ qH2O_eff);

arg_blk = (iso.mech.k_blk * qH2O_eff).^iso.mech.n_blk;
fblk = iso.mech.fblk_max .* (1 - exp(-arg_blk));
fblk = clamp(fblk, 0, iso.mech.fblk_max);

phi_avail = iso.mech.phi_max - fblk;
phi = iso.mech.phi_dry + (phi_avail - iso.mech.phi_dry) .* w;
phi = clamp(phi, 0, iso.mech.phi_max);

DH_ave = (1 - w).*iso.mech.DH_dry + w.*iso.mech.DH_wet;
bT = iso.mech.b0 .* exp(DH_ave./(R*TK));

tauT = iso.mech.tau0 + iso.mech.alpha*(1 - iso.mech.T0/TK);
tauT = max(tauT, 1e-6);

bp = max(bT .* max(pCO2_Pa,0), 0);
den = (1 + bp.^tauT).^(1/tauT);

qCO2 = (phi/iso.mech.phi_dry) .* qinf .* (bp ./ max(den, 1e-12));
qCO2 = clamp(qCO2, 0, qinf);
end


%% =========================================================================
% KINETICS (Arrhenius + mild enhancement vs pCO2)
% =========================================================================
function [kCO2, kH2O] = kinetics_LDF(TK, stage_id, pCO2_Pa, aw_eff, kin)
Rg = 8.314462618; % [J/mol/K] universal gas constant
kCO2 = kin.kCO2_ref .* exp(-kin.Ea_CO2./Rg .* (1./TK - 1./kin.Tref));
kH2O = kin.kH2O_ref .* exp(-kin.Ea_H2O./Rg .* (1./TK - 1./kin.Tref));

% stage multipliers
kCO2 = kCO2 .* kin.mCO2(stage_id);
kH2O = kH2O .* kin.mH2O(stage_id);

% pCO2 enhancement (inlet-based, mild)
kCO2 = kCO2 .* (1 + kin.beta .* (pCO2_Pa/kin.p_ref_Pa).^kin.n_exp);

% extra adsorption-only multiplier (lets us tune early uptake curvature without touching other stages)
% Paper-defense: represents reduced film/micro-pore resistance during adsorption (e.g., higher effective k due to higher driving force / better contacting)
if stage_id==6 && isfield(kin,'kCO2_ads_mult') && isfinite(kin.kCO2_ads_mult) && kin.kCO2_ads_mult>0
    kCO2 = kCO2 .* kin.kCO2_ads_mult;
end

% --- low-aw kCO2 gate (ads-only) ---
prmGate = struct('aw_hi',0.22,'aw_lo',0.14,'fmax',4.0,'power',2);
if isfield(kin,'gate_lowaw'); prmGate = kin.gate_lowaw; end
kCO2 = kCO2_gate_lowaw(kCO2, aw_eff, (stage_id==6), prmGate);
% -----------------------------------

kCO2 = clamp(kCO2, 1e-8, 0.05);
kH2O = clamp(kH2O, 1e-8, 0.10);

end

function kCO2_out = kCO2_gate_lowaw(kCO2_in, aw_eff, isAds, prmGate)
% kCO2_out = kCO2_in * gate(aw_eff)  (adsorption-only)
% - aw_eff >= aw_hi : ~1 (no boost)
% - aw_eff <= aw_lo : fmax (max boost)
% - smooth power-law ramp in between, with clamp

if nargin < 4 || isempty(prmGate)
    prmGate = struct('aw_hi',0.22,'aw_lo',0.14,'fmax',4.0,'power',2);
end

if ~isAds
    kCO2_out = kCO2_in;
    return;
end

aw    = clamp(aw_eff, 0.0, 1.0);
aw_hi = prmGate.aw_hi;
aw_lo = prmGate.aw_lo;
fmax  = prmGate.fmax;
pwr   = prmGate.power;

% guards
if ~(aw_hi > aw_lo) || fmax <= 1.0
    kCO2_out = kCO2_in;
    return;
end

u   = clamp01((aw_hi - aw) / (aw_hi - aw_lo));   % 0@hi, 1@lo
fac = 1.0 + (fmax - 1.0) * (u.^pwr);
fac = clamp(fac, 1.0, fmax);

kCO2_out = kCO2_in .* fac;
end

%% =========================================================================
% THERMAL LAG
% =========================================================================
function tauT = get_tauT_from_id(stage_id, tauT_s)
switch double(stage_id)
    case 1, tauT = tauT_s.vac;
    case 2, tauT = tauT_s.heat;
    case 3, tauT = tauT_s.des;
    case 4, tauT = tauT_s.cool;
    case 5, tauT = tauT_s.press;
    case 6, tauT = tauT_s.ads;
    otherwise, tauT = tauT_s.des;
end
end

function TK_next = apply_thermal_lag_step(TK_prev, TK_cmd, dt, tau)
if tau <= 0
    TK_next = TK_cmd;
else
    a = 1 - exp(-dt/max(tau,1e-9));
    TK_next = TK_prev + a*(TK_cmd - TK_prev);
end
end


%% =========================================================================
% PPI CALC (Sendi Eq23–34 style) + v10_8 corrections/audit
% =========================================================================
function ppi = ppi_calc_sendi_eq23_34_raw(amb, ctrl, phys, design, drv, ...
    qCO2, qH2O, TKbed, yCO2g, pCO2_inv, pCO2_drv, pH2O_drv, Pdry, ...
    Qin_vec, Qout_vec, yin_vec, cyc, opt, Q_ads) %#ok<INUSD>

R   = phys.eq.R;
MWc = phys.eq.MW_CO2;

tcycle = drv.t_s(end);
cycles_per_year = (365*86400) / tcycle;

nCO2 = max(cyc.nCO2_captured_eff, 0);
mCO2_kg = nCO2 * MWc;
tCO2 = max(mCO2_kg/1000, 1e-12);

Prod_tCO2_collector_yr = max((mCO2_kg * cycles_per_year)/1000, 0);
Prod_tCO2_m3bed_yr     = Prod_tCO2_collector_yr / max(design.Vbed_m3, 1e-12);

ppi = struct();
ppi.Prod_tCO2_collector_yr = Prod_tCO2_collector_yr;
ppi.Prod_tCO2_m3bed_yr     = Prod_tCO2_collector_yr / max(design.Vbed_m3, 1e-12);
ppi.breakdown = struct();
ppi.audit = struct();
ppi.COP = NaN;
ppi.Elec_total_MWh_tCO2 = NaN;

if mCO2_kg < 1e-12
    return;
end

is_regen_heatdes = (drv.stage_id==2) | (drv.stage_id==3);
idx = find(is_regen_heatdes);
t  = drv.t_s(:);

if numel(idx) < 3
    return;
end

ti = t(idx);
T  = TKbed(idx);

dTdt    = gradient(T, opt.dt_s);
dqCO2dt = gradient(qCO2(idx), opt.dt_s);

% ---------------- v10_8: H2O 'effective' for ENERGY-ONLY pathway ----------
% Default: no change unless opt.corr.h2o_eff.enable = true
qH2O_energy = qH2O(idx);
f_h2o_eff_vec = ones(size(qH2O_energy));
if isfield(opt,'corr') && isfield(opt.corr,'h2o_eff') && isfield(opt.corr.h2o_eff,'enable') && opt.corr.h2o_eff.enable
    % Gate using ambient RH by default (captures lowT/highRH cases without
    % double-counting regen-side RH trajectory shaping).
    RH_reg  = drv.RH_cmd(idx);
    sid_reg = uint8(drv.stage_id(idx));

    use_amb = true;
    if isfield(opt.corr.h2o_eff,'use_amb_RH_for_gating')
        use_amb = logical(opt.corr.h2o_eff.use_amb_RH_for_gating);
    end
    if use_amb
        RH_reg = amb.RH_in .* ones(size(RH_reg));
    end

    f_h2o_eff_vec = h2o_eff_factor_from_RH(RH_reg, sid_reg, opt.corr.h2o_eff);
    qH2O_energy = qH2O_energy .* f_h2o_eff_vec;
end
dqH2Odt_energy = gradient(qH2O_energy, opt.dt_s);
% -------------------------------------------------------------------------

Ptot = drv.P_Pa(idx);
ctot = Ptot ./ (R .* max(T,1.0));

% gas composition uses inventory-consistent pCO2
yH2O = clamp(pH2O_drv(idx) ./ max(Ptot,1.0), 0, 0.98);
yCO2 = clamp(pCO2_inv(idx) ./ max(Ptot,1.0), 0, 0.95);
yAir = clamp(1 - yH2O - yCO2, 0, 1);

cCO2 = yCO2 .* ctot;
cH2O = yH2O .* ctot;
cAir = yAir .* ctot;

gasCp_vol   = design.eps * (cCO2*phys.eq.Cp_CO2_g + cH2O*phys.eq.Cp_H2O_g + cAir*phys.eq.Cp_air_g);

% --- Dry-sorbent sensible heat correction (regen only) ---
% aw_reg ~ pH2O/psat(T) during regen/heating+des region
aw_reg = clamp(pH2O_drv(idx) ./ max(psat_H2O_Pa(TKbed(idx)), 1.0), 0.0, 0.999);
aw_min_reg = min(aw_reg);

dryCp_gain = getfield_safe(opt, {'energy','dryCp_gain'}, 0.0);  % 0.0~0.25 권장 범위
% dry => (1-aw) 커짐 => solid Cp 약간 증가
solidCp_vol = phys.bed.rb * phys.eq.Cp_s .* (1 + dryCp_gain*(1 - aw_reg));

sens_power_vol = (gasCp_vol + solidCp_vol) .* dTdt;

% v10_8: H2O term uses dqH2Odt_energy (energy-only correction)
rxn_power_vol  = -phys.bed.rb * ( phys.eq.dHads_CO2*dqCO2dt + phys.eq.dHads_H2O*dqH2Odt_energy );

power_vol = sens_power_vol + rxn_power_vol;
Qbed_J = design.Vbed_m3 * trapz(ti, max(power_vol,0));
Qbed_MWh_tCO2 = (Qbed_J/3.6e9) / tCO2;

% Steam surrogate
n_steam = phys.eq.r_steam_CO2 * nCO2;
m_steam = n_steam * phys.eq.MW_H2O;

TsatK = phys.eq.TsatC + 273.15;
TambK = amb.T_ambC + 273.15;

Qsteam_sens_J = m_steam * phys.eq.Cp_liq * max(TsatK - TambK, 0);
Qsteam_gen_J  = m_steam * phys.eq.hvap;

Qsteam_sens_MWh_tCO2 = (Qsteam_sens_J/3.6e9) / tCO2;
Qsteam_gen_MWh_tCO2  = (Qsteam_gen_J /3.6e9) / tCO2;

% COP (robust, weakly sensitive to ambient at extremes)
% - use regen setpoint for TH (reviewer-defensible)
% - clamp TC to avoid unrealistically poor COP at very cold ambient
% - clamp COP bounds for stability

% TH: regen setpoint (if available), else fallback to max(T)
if isfield(phys,'ctrl') && isfield(phys.ctrl,'T_regenC')
    TH = (phys.ctrl.T_regenC + 273.15) + phys.eq.dTmin_K;
else
    TheatingK = max(T);
    TH = TheatingK + phys.eq.dTmin_K;
end

% TC: effective heat sink temperature (clamped)
TC_raw = (amb.T_ambC + 273.15) - phys.eq.dTmin_K;

% clamp: do not allow sink to be colder than ~10C equivalent (tunable, but stable)
TC_min = (10 + 273.15) - phys.eq.dTmin_K;
TC = max(TC_raw, TC_min);

% ensure lift positive
TC = min(TC, TH - 1e-3);

COP_ideal = TH / max(TH - TC, 1e-6);
COP = phys.eq.eta2nd * COP_ideal;

% stability bounds
COP_min = 1.2;
COP_max = 6.0;
COP = min(max(COP, COP_min), COP_max);

E_bed        = Qbed_MWh_tCO2 / COP;
E_steam_sens = Qsteam_sens_MWh_tCO2 / COP;
E_steam_gen  = Qsteam_gen_MWh_tCO2  / COP;

% Vacuum compression work (dry gas)
Patm_tot = ambient_pressure_bar(amb.altitude_m) * 1e5;
pH2O_amb = clamp(amb.RH_in * psat_H2O_Pa(amb.T_ambC + 273.15), 0, 0.98*Patm_tot);
Patm_dry = max(Patm_tot - pH2O_amb, 1.0);

W_vac_J = 0.0;
W_vac_total_J = 0.0; % optional audit
n_dot_dry_sum = 0.0;
n_dot_tot_sum = 0.0;

audit_compare_total = false;
if isfield(opt,'corr') && isfield(opt.corr,'audit') && isfield(opt.corr.audit,'vac_compare_total_enable')
    audit_compare_total = logical(opt.corr.audit.vac_compare_total_enable);
end

for k=1:(numel(t)-1)
    Pdry_k = max(Pdry(k), 1.0);
    if Pdry_k < Patm_dry && Qout_vec(k) > 0
        Tk = max(TKbed(k), 250.0);

        % dry surrogate
        n_dot_dry = (Pdry_k/(R*Tk)) * Qout_vec(k);
        Wmol_dry  = R*Tk*log(Patm_dry/max(Pdry_k,1.0));
        W_vac_J = W_vac_J + (n_dot_dry * Wmol_dry * opt.dt_s);
        n_dot_dry_sum = n_dot_dry_sum + n_dot_dry*opt.dt_s;

        % optional "total-gas" comparison (audit only; OFF by default)
        if audit_compare_total
            Ptot_k = max(drv.P_Pa(k), 1.0);
            n_dot_tot = (Ptot_k/(R*Tk)) * Qout_vec(k);
            Wmol_tot  = R*Tk*log(Patm_tot/max(Ptot_k,1.0));
            W_vac_total_J = W_vac_total_J + (n_dot_tot * Wmol_tot * opt.dt_s);
            n_dot_tot_sum = n_dot_tot_sum + n_dot_tot*opt.dt_s;
        end
    end
end

E_comp_vac = (W_vac_J / max(phys.eq.eta_comp,1e-6)) / 3.6e9 / tCO2;

% --- Vacuum penalty at very dry conditions (targets RH=5~20% under-pred) ---
vac_aw_ref   = getfield_safe(opt, {'energy','vac_aw_ref'}, 0.25);   % 기준 aw
vac_aw_gain  = getfield_safe(opt, {'energy','vac_aw_gain'}, 0.60);  % 0~1 범위 권장
aw_min_reg   = min(aw_reg);                                         % regen 구간 최소 aw
if aw_min_reg < vac_aw_ref
    vac_fac = 1.0 + vac_aw_gain*(vac_aw_ref - aw_min_reg);
    E_comp_vac = E_comp_vac * vac_fac;
end

% Product CO2 compression (captured CO2 only)
Pco2_out = phys.eq.P_CO2_out_bar * 1e5;
TcompK   = phys.eq.T_comp_K;
W_prod_J = nCO2 * R * TcompK * log(Pco2_out / max(Patm_tot,1.0));
E_comp_prod = (W_prod_J / max(phys.eq.eta_comp,1e-6)) / 3.6e9 / tCO2;

E_comp = E_comp_vac + E_comp_prod;

% Fan power (ads stage)
idx_ads0 = find(drv.stage_id==6,1,'first');
if isempty(idx_ads0), idx_ads0=numel(drv.stage_id); end

T_ads = TKbed(idx_ads0);
[mu_air, rho_g] = gas_props_mu_rho(T_ads, max(Pdry(idx_ads0),1.0), phys.eq.MW_air);

u_erg = Q_ads / max(design.A_flow_m2, 1e-12);
Lflow = design.L_flow_m;

dP_cell  = ergun_dp(u_erg, Lflow, rho_g, mu_air, design.eps, design.dp_m);
dP_total = design.N_series_pass * dP_cell + design.dP_misc_Pa;

t_ads_s = stage_duration_s(drv.stage_id, 6, opt.dt_s);
V_inlet = Q_ads * t_ads_s;

E_fan = (V_inlet * dP_total / max(phys.eq.eta_fan,1e-6)) / 3.6e9 / tCO2;

% Totals (pre-correction)
E_heating = (E_bed + E_steam_sens + E_steam_gen);
E_tot = E_heating + E_comp + E_fan;

% v10_7: optional fixed electrical overheads (per collector), [kW]
P_idle_kW_fan = getfield_safe(opt, {'energy','fan_idle_kW'}, 0.0);
P_idle_kW_vac = getfield_safe(opt, {'energy','vac_idle_kW'}, 0.0);
E_idle = ((P_idle_kW_fan + P_idle_kW_vac) .* (tcycle/3600.0) / 1000.0) ./ max(mCO2_kg/1000.0, 1e-12); % [MWh/tCO2]
E_tot = E_tot + E_idle;
% v14: extra per-cycle parasitic electricity (kWh/cycle), helps low-Prod & high T/RH cases
Epar_kWh_cyc = getfield_safe(opt, {'energy','parasitic_kWh_per_cycle'}, 0.0);

% extra cycle-based balance-of-plant (BOP) parasitic load [kW] integrated over full cycle time
% Paper-defense: baseline electrical overhead for controls/valves/heaters/actuation not captured in fan/compressor work
P_bop_kW = getfield_safe(opt, {'energy','parasitic_kW_base'}, 0.0);
if isfinite(P_bop_kW) && P_bop_kW>0
    Epar_kWh_cyc = Epar_kWh_cyc + P_bop_kW * (tcycle/3600.0);
end

% Optional weak dependence on (T,RH): boosts 50C/90% case without exploding others
par_T_slope  = getfield_safe(opt, {'energy','parasitic_T_slope_kWh_per_C'}, 0.0);   % kWh/cycle/°C above Tref
par_RH_slope = getfield_safe(opt, {'energy','parasitic_RH_slope_kWh_per_RH'}, 0.0); % kWh/cycle per RH(0-1) above RHref
TrefC  = getfield_safe(opt, {'energy','parasitic_Tref_C'}, 20.0);
RHref  = getfield_safe(opt, {'energy','parasitic_RHref'}, 0.50);

% Use point ambient conditions (available + cheap + stable)
TmeanC  = getfield_safe(amb, {'T_ambC'}, 20.0);
RHmean  = getfield_safe(amb, {'RH_in'}, 0.50);

par_TRH_slope = getfield_safe(opt, {'energy','parasitic_TRH_slope_kWh_per_C_RH'}, 0.0); % kWh/cycle per (°C * RH) above refs

Tex = max(TmeanC - TrefC, 0.0);
RHex = max(RHmean - RHref, 0.0);

Epar_kWh_cyc = Epar_kWh_cyc ...
    + par_T_slope  * Tex ...
    + par_RH_slope * RHex ...
    + getfield_safe(opt, {'energy','parasitic_TRH_cross_kWh'}, 0.0) * Tex * RHex;

E_par = (Epar_kWh_cyc/1000.0) ./ max(mCO2_kg/1000.0, 1e-12); % [MWh/tCO2]
E_tot = E_tot + E_par;

% ---------------- v10_8: RH-dependent PPI correction layer ----------------
mult_RH = 1.0;
apply_to = 'heating';
if isfield(opt,'corr') && isfield(opt.corr,'ppi_RH') && isfield(opt.corr.ppi_RH,'enable') && opt.corr.ppi_RH.enable
    [mult_RH, apply_to] = ppi_RH_multiplier(amb.RH_in, amb.T_ambC, opt.corr.ppi_RH);

    switch lower(string(apply_to))
        case 'total'
            % scale total AND all components that contribute to total
            E_tot = E_tot * mult_RH;

            E_heating = E_heating * mult_RH;
            E_bed = E_bed * mult_RH;
            E_steam_sens = E_steam_sens * mult_RH;
            E_steam_gen  = E_steam_gen  * mult_RH;

            E_comp = E_comp * mult_RH;
            E_comp_vac  = E_comp_vac  * mult_RH;
            E_comp_prod = E_comp_prod * mult_RH;

            E_fan  = E_fan  * mult_RH;
            E_idle = E_idle * mult_RH;

            % v14: parasitic electricity must be consistent too
            E_par  = E_par  * mult_RH;

        otherwise % 'heating' (default)
            % apply RH correction only to heating-related terms
            E_heating = E_heating * mult_RH;
            E_bed = E_bed * mult_RH;
            E_steam_sens = E_steam_sens * mult_RH;
            E_steam_gen  = E_steam_gen  * mult_RH;

            % recompute totals with unchanged comp/fan/idle/parasitic
            E_tot = E_heating + E_comp + E_fan + E_idle + E_par;
    end
end
% -------------------------------------------------------------------------

ppi.COP = COP;

ppi.breakdown.E_bed        = E_bed;
ppi.breakdown.E_steam_sens = E_steam_sens;
ppi.breakdown.E_steam_gen  = E_steam_gen;
ppi.breakdown.E_heating    = E_heating;
ppi.breakdown.E_comp_vac   = E_comp_vac;
ppi.breakdown.E_comp_prod  = E_comp_prod;
ppi.breakdown.E_comp       = E_comp;
ppi.breakdown.E_fan        = E_fan;
ppi.breakdown.E_idle       = E_idle;

% v14: report parasitic electricity explicitly
ppi.breakdown.E_parasitic  = E_par;

Et = max(E_tot, 1e-12);
ppi.breakdown.frac_heat = ppi.breakdown.E_heating / Et;
ppi.breakdown.frac_comp = ppi.breakdown.E_comp / Et;
ppi.breakdown.frac_fan  = ppi.breakdown.E_fan  / Et;
ppi.breakdown.frac_idle = ppi.breakdown.E_idle / Et;

% v14: add parasitic fraction for clarity (optional but useful)
ppi.breakdown.frac_parasitic = ppi.breakdown.E_parasitic / Et;

ppi.breakdown.corr_ppi_RH_mult  = mult_RH;
ppi.breakdown.corr_ppi_RH_apply = char(apply_to);
ppi.breakdown.corr_h2o_eff_mean = mean(f_h2o_eff_vec);

ppi.Elec_total_MWh_tCO2 = max(E_tot, 0);

% ---------------- v10_8: audit block (consistency checks) -----------------
ppi.audit.Patm_tot_Pa = Patm_tot;
ppi.audit.pH2O_amb_Pa = pH2O_amb;
ppi.audit.Patm_dry_Pa = Patm_dry;

ppi.audit.Pdry_min_Pa = min(Pdry);
ppi.audit.Pdry_max_Pa = max(Pdry);

ppi.audit.mean_yH2O_regen = mean(yH2O);
ppi.audit.mean_RH_regen   = mean(drv.RH_cmd(idx));

ppi.audit.vac_n_dry_mol = n_dot_dry_sum;
ppi.audit.vac_compare_total_enable = audit_compare_total;

if audit_compare_total
    E_comp_vac_total = (W_vac_total_J / max(phys.eq.eta_comp,1e-6)) / 3.6e9 / tCO2;
    ppi.audit.vac_E_comp_total_MWh_tCO2 = E_comp_vac_total;
    ppi.audit.vac_E_comp_dry_MWh_tCO2   = E_comp_vac;
    ppi.audit.vac_E_ratio_total_over_dry = E_comp_vac_total / max(E_comp_vac, 1e-12);
    ppi.audit.vac_n_tot_mol = n_dot_tot_sum;
else
    ppi.audit.vac_E_comp_total_MWh_tCO2 = NaN;
    ppi.audit.vac_E_ratio_total_over_dry = NaN;
    ppi.audit.vac_n_tot_mol = NaN;
end
% -------------------------------------------------------------------------
end


%% =========================================================================
% S_cap baseline calibration (single scalar)
% =========================================================================
function phys = calibrate_Scap_baseline(amb, ctrl, design, phys, opt, drv_base)
bnd = opt.calib.Scap_bounds(:);
Slo = max(bnd(1), 0.05);
Shi = max(bnd(2), Slo + 1e-6);
target = opt.calib.Prod_target_tpy;

x_ws = [];
    function [Prod, sim] = eval_S(S)
        phys2 = phys;
        phys2.iso.mech.S_cap = S;
        sim = run_point_periodic(drv_base, amb, ctrl, design, phys2, opt, x_ws);
        x_ws = sim.x_end;
        Prod = sim.ppi.Prod_tCO2_collector_yr;
    end

if opt.calib.verbose
    fprintf('\n================ S_cap CALIBRATION (baseline only) ================\n');
    fprintf('Target Prod = %.3f tCO2/collector/yr | bounds [%.3f, %.3f]\n', target, Slo, Shi);
end

[Pl, simL] = eval_S(Slo);
[Ph, simH] = eval_S(Shi);

fl = Pl - target;
fh = Ph - target;

bestS = Slo; bestErr = abs(fl); bestSim = simL;
if abs(fh) < bestErr, bestS = Shi; bestErr = abs(fh); bestSim = simH; end

if fl == 0
    bestS = Slo;
elseif fh == 0
    bestS = Shi;
elseif sign(fl) ~= sign(fh)
    for it = 1:opt.calib.max_iter
        Sm = 0.5*(Slo + Shi);
        [Pm, simM] = eval_S(Sm);
        fm = Pm - target;

        if abs(fm) < bestErr
            bestErr = abs(fm);
            bestS = Sm;
            bestSim = simM;
        end

        if opt.calib.verbose
            fprintf('iter %02d | S=%.5f | Prod=%.3f | err=%.3f\n', it, Sm, Pm, fm);
        end

        if abs(fm) <= opt.calib.tol_rel_prod * max(target, 1e-6), break; end
        if sign(fm) == sign(fl), Slo = Sm; fl = fm; else, Shi = Sm; fh = fm; end
    end
else
    if opt.calib.verbose
        fprintf('NOTE: Target not bracketed. Selecting closest bound.\n');
    end
end

phys.iso.mech.S_cap = bestS;

if opt.calib.verbose
    fprintf('CALIB DONE: S_cap=%.5f | |err|=%.3f | Prod=%.3f | Elec=%.3f\n', ...
        bestS, bestErr, bestSim.ppi.Prod_tCO2_collector_yr, bestSim.ppi.Elec_total_MWh_tCO2);
    fprintf('===================================================================\n\n');
end
end


%% =========================================================================
% ERGUN + gas props
% =========================================================================
function [mu, rho] = gas_props_mu_rho(TK, P_Pa, MW)
mu0 = 1.716e-5; T0 = 273.15; S = 111.0;
mu  = mu0 * (TK/T0)^(3/2) * (T0 + S) / (TK + S);

R = 8.314462618;
rho = P_Pa * MW / (R * max(TK,1.0));

mu  = max(mu, 8e-6);
rho = max(rho, 0.2);
end

function dp = ergun_dp(u, L, rho_g, mu_g, eps_bed, d_p)
term1 = 150 * (1 - eps_bed)^2 / eps_bed^3 * mu_g * u / d_p^2;
term2 = 1.75 * (1 - eps_bed)   / eps_bed^3 * rho_g * u^2 / d_p;
dp = (term1 + term2) * L;
end


%% =========================================================================
% psat + ambient pressure
% =========================================================================
function ps = psat_H2O_Pa(TK)
TC = TK - 273.15;
ps_hPa = 6.1121 .* exp((18.678 - (TC./234.5)) .* (TC ./ (257.14 + TC)));
ps = max(ps_hPa .* 100, 1.0);
end

function Pbar = ambient_pressure_bar(alt_m)
P0 = 1.01325; H = 8400;
Pbar = P0 * exp(-alt_m / H);
end


%% =========================================================================
% Reference table (Sendi 12pt)
% =========================================================================
function ref = ref_sendi_12pt()
T = [ 1;  1;  1;  1;  20; 20; 20; 20; 50; 50; 50; 50];
H = [ 5; 20; 50; 90;  5; 20; 50; 90;  5; 20; 50; 90] / 100;

Prod = [61; 64; 62; 57; 54; 57; 56; 52; 22; 32; 36; 33];
Elec = [1.6;1.5;1.7;2.2;1.5;1.4;1.6;2.2;2.4;1.8;1.9;2.85];

t_des = [10000;14000;20000;22000; 10000;14000;20000;22000; 10000;14000;20000;22000];
ref = table(T, H, t_des, Prod, Elec, 'VariableNames', {'T_C','RH','t_des_s','Prod_ref','Elec_ref'});
end


%% =========================================================================
% Utilities
% =========================================================================
function x = clamp(x, xmin, xmax)
x = min(max(x, xmin), xmax);
end

function u = clamp01(u)
u = min(max(u, 0), 1);
end

function s = smoothstep(u)
u = clamp01(u);
s = u.*u.*(3 - 2*u);
end

function dt_stage = stage_duration_s(stage_id, sid, dt)
if numel(stage_id) < 2, dt_stage = 0; return; end
dt_stage = sum(stage_id(1:end-1)==sid) * dt;
end

function [Q_ads, A_used] = get_Q_ads_and_A(u, design, flow_mode, opt)

switch lower(string(flow_mode))
    case 'aflow'
        A_used = design.A_flow_m2;
    case 'acontact'
        A_used = design.A_contact_m2;
    otherwise
        error('flow_mode must be Aflow or Acontact');
end

% legacy.keepQ_with_Acontact:
%  - true  : if flow_mode=Acontact, keep legacy throughput Q=u*A_contact
%  - false : lock throughput to A_flow (paper-defensible default)
keepQ_with_Acontact = getfield_safe(opt, {'legacy','keepQ_with_Acontact'}, false);

if keepQ_with_Acontact && lower(string(flow_mode))=='acontact'
    Q_ads = u * design.A_contact_m2;   % legacy (old behavior) -> Prod ~2.1x back
else
    Q_ads = u * design.A_flow_m2;      % paper-defensible: Q fixed by A_flow
end
end

function x_next = lpf_step(x_prev, x_raw, dt, tau)
if ~(isfinite(tau) && tau > 0)
    x_next = x_raw; return;
end
a = exp(-dt/max(tau,1e-12));
x_next = a*x_prev + (1-a)*x_raw;
end


%% =========================================================================
% Smooth P(t) near stage transitions (finite slope switching)
% =========================================================================
function drv = smooth_driver_pressure_edges(drv, dt, ramp_s)
P = drv.P_Pa(:);
sid = drv.stage_id(:);
N = numel(P);

if ~(isfinite(ramp_s) && ramp_s > 0), return; end

auto_on = false;
rmin = 10; rmax = 90; rgain = 1.8;

if isfield(drv,'meta') && isfield(drv.meta,'ptraj')
    p = drv.meta.ptraj;
    if isfield(p,'driver_edge_ramp_auto_enable'), auto_on = p.driver_edge_ramp_auto_enable; end
    if isfield(p,'driver_edge_ramp_min_s'), rmin = p.driver_edge_ramp_min_s; end
    if isfield(p,'driver_edge_ramp_max_s'), rmax = p.driver_edge_ramp_max_s; end
    if isfield(p,'driver_edge_ramp_gain'),  rgain = p.driver_edge_ramp_gain; end
end

for k = 2:N
    if sid(k) ~= sid(k-1)
        i0 = k-1;

        ramp_s_local = ramp_s;
        if auto_on
            P0 = max(P(i0), 1.0);
            P1 = max(P(k), 1.0);
            jump = abs(P1 - P0) / P0;
            ramp_s_local = ramp_s * (1 + rgain * min(jump, 1.0));
            ramp_s_local = clamp(ramp_s_local, rmin, rmax);
        end

        nR = max(1, round(ramp_s_local / max(dt,1e-12)));
        j2 = min(N, k+nR);

        P0 = P(i0);
        P2 = P(j2);

        m = j2 - i0;
        if m < 2, continue; end

        u = (0:m)'/m;
        s = u.*u.*(3 - 2*u);
        P(i0:j2) = P0 + (P2 - P0)*s;
    end
end

P = max(P, 200.0);
drv.P_Pa = P;
drv.Pa   = P;
drv.Pbar = P/1e5;
end

% ------------------------- small utility -------------------------
function v = getfield_safe(s, field_path, default_v)
v = default_v;
try
    tmp = s;
    for i=1:numel(field_path)
        if ~isstruct(tmp) || ~isfield(tmp, field_path{i})
            return;
        end
        tmp = tmp.(field_path{i});
    end
    v = tmp;
catch
    v = default_v;
end
end

% ------------------------- small utility -------------------------
% Set nested field default if missing OR empty.
% Usage: s = setfield_default(s, {'a','b','c'}, default_value);
function s = setfield_default(s, field_path, default_v)

% Make sure s is a struct
if nargin < 1 || isempty(s) || ~isstruct(s)
    s = struct();
end
if nargin < 2 || isempty(field_path)
    return;
end

f1 = field_path{1};

if numel(field_path) == 1
    if ~isfield(s, f1) || isempty(s.(f1))
        s.(f1) = default_v;
    end
else
    if ~isfield(s, f1) || isempty(s.(f1)) || ~isstruct(s.(f1))
        s.(f1) = struct();
    end
    s.(f1) = setfield_default(s.(f1), field_path(2:end), default_v);
end
end

%% =========================================================================
% Dry-side CO2 equilibrium scaling (single scalar; optional)
% =========================================================================
function q = apply_CO2_dry_scale(q_in, aw_eff, sid, opt)
q = q_in;

if ~isfield(opt,'co2_dry_scale'); return; end
ds = opt.co2_dry_scale;
if ~isfield(ds,'enable') || ~ds.enable; return; end
if ~isfield(ds,'factor'); return; end
fac_hi = ds.factor;
if ~(fac_hi > 1.0); return; end

if sid ~= 6 && sid ~= 5
    return;
end

aw0 = 0.20;
if isfield(ds,'aw_switch'); aw0 = ds.aw_switch; end

% Shape controls (optional):
% - 'p' > 1 concentrates the boost at very low aw, p < 1 spreads it.
p = 1.0;
if isfield(ds,'p'); p = ds.p; end
p = max(p, 0.25);

% Smooth transition width around aw0 (keeps mid/high RH essentially untouched)
w  = 0.03;
if isfield(ds,'aw_smooth'); w  = ds.aw_smooth; end
w = max(w, 1e-6);

% gate g: ~1 when aw_eff<<aw0, ~0 when aw_eff>>aw0
g = 0.5*(1.0 - tanh((aw_eff - aw0)/w));

% power-law intensity toward dry end
x = min(1.0, max(0.0, (aw0 - aw_eff)/max(aw0,1e-9)));
fac_dry = 1.0 + (fac_hi - 1.0) * (x^p);

% blend: only activates when g~1 and fades smoothly around aw0
fac = 1.0 + (fac_dry - 1.0)*g;

q = q_in * fac;
end


%% =========================================================================
% v10_8 CORRECTION HELPERS
% =========================================================================
function [mult, apply_to] = ppi_RH_multiplier(RH_in, T_ambC, p)
% RH-based mild correction (bounded, reviewer-defensible)
% - default: no correction (mult ~ 1)
% - allows 'credit' at low RH (reduced water/latent burden) if k_dry < 0
% - bounded to avoid distorting physics

mult = 1.0;
apply_to = 'heating';

RH = clamp(RH_in, 0, 1);
T  = T_ambC;

RHref  = 0.50;
Tref   = 20;    % [C]
Tscale = 30;    % [C]

% --- knobs (defaults: OFF / very mild) ---
k_dry = 0.0;    % negative => credit at dry (reduce Elec), positive => penalty
k_wet = 0.0;    % optional penalty at wet
k_TRH = 0.0;    % optional cross term (keep 0 for now)

mmin = 0.85;
mmax = 1.20;

% --- read user params (backward compatible) ---
if isfield(p,'apply_to'),    apply_to = string(p.apply_to); end
if isfield(p,'RHref'),       RHref  = p.RHref; end
if isfield(p,'Tref_C'),      Tref   = p.Tref_C; end
if isfield(p,'Tscale_C'),    Tscale = p.Tscale_C; end
if isfield(p,'mult_min'),    mmin   = p.mult_min; end
if isfield(p,'mult_max'),    mmax   = p.mult_max; end

if isfield(p,'k_dry'), k_dry = p.k_dry;
elseif isfield(p,'k'), k_dry = p.k;
end
if isfield(p,'k_wet'), k_wet = p.k_wet; end
if isfield(p,'k_TRH'), k_TRH = p.k_TRH; end

dry = max(RHref - RH, 0.0);
wet = max(RH - RHref, 0.0);
Th  = max((T - Tref)/max(Tscale,1e-9), 0.0);

% mild form
mult = 1 + k_dry*dry + k_wet*wet + k_TRH*Th*wet;
mult = clamp(mult, mmin, mmax);
end

function f = h2o_eff_factor_from_RH(RH_vec, sid_vec, p)
% Returns f<=1 (typ.) applied to qH2O ONLY in energy pathway.
% Default: gate at high RH with smooth transition.
RH = clamp(RH_vec(:), 0, 1);

f_min = 0.75;
f_max = 1.0;
RH_on = 0.70;
RH_width = 0.08;
apply_stage_ids = uint8([2 3]); % default: heat+des only

if isfield(p,'f_min'), f_min = p.f_min; end
if isfield(p,'f_max'), f_max = p.f_max; end
if isfield(p,'RH_on'), RH_on = p.RH_on; end
if isfield(p,'RH_width'), RH_width = p.RH_width; end
if isfield(p,'apply_in_stage_ids'), apply_stage_ids = uint8(p.apply_in_stage_ids(:)); end

RH_width = max(RH_width, 1e-6);

% Smooth gate g ~ 0 below RH_on, g~1 above RH_on
g = 0.5*(1 + tanh((RH - RH_on)/RH_width));
f_raw = f_max - (f_max - f_min)*g; % high RH -> f_min
f_raw = clamp(f_raw, min(f_min,f_max), max(f_min,f_max));

sid_vec = uint8(sid_vec(:));
mask = false(size(RH));
for i=1:numel(apply_stage_ids)
    mask = mask | (sid_vec == apply_stage_ids(i));
end

f = ones(size(RH));
f(mask) = f_raw(mask);
end


%% =========================================================================
% Stage id helper
% =========================================================================
function S = stage_ids()
S = struct('Vac',uint8(1),'Heat',uint8(2),'Des',uint8(3),'Cool',uint8(4),'Press',uint8(5),'Ads',uint8(6));
end


%% =========================================================================
% PTRaj default pack (robustness layer)
% - Several core routines assume opt.ptraj.* exists. Wrappers may omit it.
% - This function defines a conservative default that preserves v10_* behavior:
%   * ptraj inversion is enabled (it is the core closure)
%   * smoothing/filters are OFF by default unless explicitly enabled
% =========================================================================

function opt = ensure_opt_defaults(opt)
% Ensure core-required option fields exist. This allows lightweight wrappers
% to pass partial opt structs without causing runtime errors.
if nargin<1 || isempty(opt), opt = struct(); end

% ---- ptraj (pressure/driver inversion / smoothing) ----
opt = ensure_ptraj_defaults(opt);

% ---- aw / RH handling ----
opt = ensure_aw_defaults(opt);

% ---- flow / geometry coupling ----
if ~isfield(opt,'flow_mode') || isempty(opt.flow_mode)
    % Q = u * A_flow (default). Alternatives: 'Acontact' for u*A_contact
    opt.flow_mode = 'Aflow';
end

% ---- yCO2 guards ----
if ~isfield(opt,'yco2') || isempty(opt.yco2), opt.yco2 = struct(); end
if ~isfield(opt.yco2,'y_floor') || isempty(opt.yco2.y_floor), opt.yco2.y_floor = 1e-12; end
if ~isfield(opt.yco2,'y_max')   || isempty(opt.yco2.y_max),   opt.yco2.y_max   = 0.50;   end

% ---- stage-wise temperature command smoothing time constants ----
if ~isfield(opt,'tauT_s') || isempty(opt.tauT_s)
    opt.tauT_s = struct('vac',120,'heat',180,'des',180,'cool',120,'press',20,'ads',0);
else
    % fill missing fields conservatively
    f = {'vac',120; 'heat',180; 'des',180; 'cool',120; 'press',20; 'ads',0};
    for i=1:size(f,1)
        k=f{i,1};
        if ~isfield(opt.tauT_s,k) || isempty(opt.tauT_s.(k))
            opt.tauT_s.(k)=f{i,2};
        end
    end
end

% ---- CO2 2-pool defaults ----
if ~isfield(opt,'co2_2pool') || isempty(opt.co2_2pool), opt.co2_2pool = struct(); end
if ~isfield(opt.co2_2pool,'enable') || isempty(opt.co2_2pool.enable), opt.co2_2pool.enable = true; end
if ~isfield(opt.co2_2pool,'fslow') || isempty(opt.co2_2pool.fslow)
    opt.co2_2pool.fslow = struct('f0_ads',0.25,'f0_des',0.55,'aw_scale',0.40,'Tmid_K',323.15,'Twidth_K',10,'stage_aware',false);
else
    if ~isfield(opt.co2_2pool.fslow,'f0_ads'), opt.co2_2pool.fslow.f0_ads = 0.25; end
    if ~isfield(opt.co2_2pool.fslow,'f0_des'), opt.co2_2pool.fslow.f0_des = 0.55; end
    if ~isfield(opt.co2_2pool.fslow,'aw_scale'), opt.co2_2pool.fslow.aw_scale = 0.40; end
    if ~isfield(opt.co2_2pool.fslow,'Tmid_K'), opt.co2_2pool.fslow.Tmid_K = 323.15; end
    if ~isfield(opt.co2_2pool.fslow,'Twidth_K'), opt.co2_2pool.fslow.Twidth_K = 10; end
    if ~isfield(opt.co2_2pool.fslow,'stage_aware'), opt.co2_2pool.fslow.stage_aware = false; end
end
if ~isfield(opt.co2_2pool,'kslow') || isempty(opt.co2_2pool.kslow)
    opt.co2_2pool.kslow = struct('ratio',0.10,'aw_pow',0.0,'enable',true,'kmin',2e-6,'kmax',1e-3);
else
    if ~isfield(opt.co2_2pool.kslow,'ratio'), opt.co2_2pool.kslow.ratio = 0.10; end
    if ~isfield(opt.co2_2pool.kslow,'aw_pow'), opt.co2_2pool.kslow.aw_pow = 0.0; end
    if ~isfield(opt.co2_2pool.kslow,'enable'), opt.co2_2pool.kslow.enable = true; end
    if ~isfield(opt.co2_2pool.kslow,'kmin') || isempty(opt.co2_2pool.kslow.kmin), opt.co2_2pool.kslow.kmin = 2e-6; end
    if ~isfield(opt.co2_2pool.kslow,'kmax') || isempty(opt.co2_2pool.kslow.kmax), opt.co2_2pool.kslow.kmax = 1e-3; end
    if opt.co2_2pool.kslow.kmax < opt.co2_2pool.kslow.kmin, tmp=opt.co2_2pool.kslow.kmax; opt.co2_2pool.kslow.kmax=opt.co2_2pool.kslow.kmin; opt.co2_2pool.kslow.kmin=tmp; end
end

% ---- CO2 residual floor (canonical schema) ----
if ~isfield(opt.co2_2pool,'resid') || isempty(opt.co2_2pool.resid)
    opt.co2_2pool.resid = struct();
end

% Alias/mapping from legacy fields -> canonical fields
if isfield(opt.co2_2pool.resid,'stage_ids') && ~isfield(opt.co2_2pool.resid,'apply_in_stage_ids')
    opt.co2_2pool.resid.apply_in_stage_ids = opt.co2_2pool.resid.stage_ids;
end
if isfield(opt.co2_2pool.resid,'frac_floor') && ~isfield(opt.co2_2pool.resid,'f_res')
    opt.co2_2pool.resid.f_res = opt.co2_2pool.resid.frac_floor;
end
if isfield(opt.co2_2pool.resid,'q_min_abs') && ~isfield(opt.co2_2pool.resid,'q_abs_min')
    opt.co2_2pool.resid.q_abs_min = opt.co2_2pool.resid.q_min_abs;
end

% Canonical defaults (match TVSA_user)
if ~isfield(opt.co2_2pool.resid,'enable') || isempty(opt.co2_2pool.resid.enable), opt.co2_2pool.resid.enable = true; end
if ~isfield(opt.co2_2pool.resid,'mode')   || isempty(opt.co2_2pool.resid.mode),   opt.co2_2pool.resid.mode   = 'ads_scaled'; end
if ~isfield(opt.co2_2pool.resid,'apply_in_stage_ids') || isempty(opt.co2_2pool.resid.apply_in_stage_ids), opt.co2_2pool.resid.apply_in_stage_ids = uint8([2 3]); end
if ~isfield(opt.co2_2pool.resid,'f_res')  || isempty(opt.co2_2pool.resid.f_res),  opt.co2_2pool.resid.f_res  = 0.08; end
if ~isfield(opt.co2_2pool.resid,'q_abs_min') || isempty(opt.co2_2pool.resid.q_abs_min), opt.co2_2pool.resid.q_abs_min = 0.05; end
if ~isfield(opt.co2_2pool.resid,'tol_frac')  || isempty(opt.co2_2pool.resid.tol_frac),  opt.co2_2pool.resid.tol_frac  = 1e-4; end

% Recommended optional defaults
if ~isfield(opt.co2_2pool.resid,'k_T')    || isempty(opt.co2_2pool.resid.k_T),    opt.co2_2pool.resid.k_T    = 0.05; end
if ~isfield(opt.co2_2pool.resid,'Tref_K') || isempty(opt.co2_2pool.resid.Tref_K), opt.co2_2pool.resid.Tref_K = 373.15; end
if ~isfield(opt.co2_2pool.resid,'k_wet')  || isempty(opt.co2_2pool.resid.k_wet),  opt.co2_2pool.resid.k_wet  = 0.35; end
if ~isfield(opt.co2_2pool.resid,'n_wet')  || isempty(opt.co2_2pool.resid.n_wet),  opt.co2_2pool.resid.n_wet  = 1.0; end
if ~isfield(opt.co2_2pool.resid,'f_min')  || isempty(opt.co2_2pool.resid.f_min),  opt.co2_2pool.resid.f_min  = 0.02; end
if ~isfield(opt.co2_2pool.resid,'f_max')  || isempty(opt.co2_2pool.resid.f_max),  opt.co2_2pool.resid.f_max  = 0.25; end

% Type normalize
opt.co2_2pool.resid.apply_in_stage_ids = uint8(opt.co2_2pool.resid.apply_in_stage_ids);

% ---- H2O 2-pool defaults ----
if ~isfield(opt,'h2o_2pool') || isempty(opt.h2o_2pool), opt.h2o_2pool = struct(); end
if ~isfield(opt.h2o_2pool,'enable') || isempty(opt.h2o_2pool.enable), opt.h2o_2pool.enable = true; end
if ~isfield(opt.h2o_2pool,'fbound') || isempty(opt.h2o_2pool.fbound)
    opt.h2o_2pool.fbound = struct('f0',0.25,'aw_pow',0.5,'Tmid_K',313.15,'Twidth_K',12,'stage_aware',false,'f0_ads',0.20,'f0_des',0.35);
else
    if ~isfield(opt.h2o_2pool.fbound,'f0'), opt.h2o_2pool.fbound.f0 = 0.25; end
    if ~isfield(opt.h2o_2pool.fbound,'aw_pow'), opt.h2o_2pool.fbound.aw_pow = 0.5; end
    if ~isfield(opt.h2o_2pool.fbound,'Tmid_K'), opt.h2o_2pool.fbound.Tmid_K = 313.15; end
    if ~isfield(opt.h2o_2pool.fbound,'Twidth_K'), opt.h2o_2pool.fbound.Twidth_K = 12; end
    if ~isfield(opt.h2o_2pool.fbound,'stage_aware'), opt.h2o_2pool.fbound.stage_aware = false; end
    if ~isfield(opt.h2o_2pool.fbound,'f0_ads'), opt.h2o_2pool.fbound.f0_ads = 0.20; end
    if ~isfield(opt.h2o_2pool.fbound,'f0_des'), opt.h2o_2pool.fbound.f0_des = 0.35; end
end
if ~isfield(opt.h2o_2pool,'kslow') || isempty(opt.h2o_2pool.kslow)
    opt.h2o_2pool.kslow = struct('ratio',0.15,'aw_pow',0.0,'enable',true,'kmin',5e-7,'kmax',5e-3);
else
    if ~isfield(opt.h2o_2pool.kslow,'ratio'), opt.h2o_2pool.kslow.ratio = 0.15; end
    if ~isfield(opt.h2o_2pool.kslow,'aw_pow'), opt.h2o_2pool.kslow.aw_pow = 0.0; end
    if ~isfield(opt.h2o_2pool.kslow,'enable'), opt.h2o_2pool.kslow.enable = true; end
    if ~isfield(opt.h2o_2pool.kslow,'kmin') || isempty(opt.h2o_2pool.kslow.kmin), opt.h2o_2pool.kslow.kmin = 5e-7; end
    if ~isfield(opt.h2o_2pool.kslow,'kmax') || isempty(opt.h2o_2pool.kslow.kmax), opt.h2o_2pool.kslow.kmax = 5e-3; end
    if opt.h2o_2pool.kslow.kmax < opt.h2o_2pool.kslow.kmin, tmp=opt.h2o_2pool.kslow.kmax; opt.h2o_2pool.kslow.kmax=opt.h2o_2pool.kslow.kmin; opt.h2o_2pool.kslow.kmin=tmp; end
end

% ---- H2O residual floor (canonical schema) ----
if ~isfield(opt.h2o_2pool,'resid') || isempty(opt.h2o_2pool.resid)
    opt.h2o_2pool.resid = struct();
end

% Alias/mapping legacy -> canonical
if isfield(opt.h2o_2pool.resid,'stage_ids') && ~isfield(opt.h2o_2pool.resid,'apply_in_stage_ids')
    opt.h2o_2pool.resid.apply_in_stage_ids = opt.h2o_2pool.resid.stage_ids;
end
if isfield(opt.h2o_2pool.resid,'q_min_abs') && ~isfield(opt.h2o_2pool.resid,'q_abs_min')
    opt.h2o_2pool.resid.q_abs_min = opt.h2o_2pool.resid.q_min_abs;
end

% Canonical defaults (match TVSA_user)
if ~isfield(opt.h2o_2pool.resid,'enable') || isempty(opt.h2o_2pool.resid.enable), opt.h2o_2pool.resid.enable = true; end
if ~isfield(opt.h2o_2pool.resid,'apply_in_stage_ids') || isempty(opt.h2o_2pool.resid.apply_in_stage_ids), opt.h2o_2pool.resid.apply_in_stage_ids = uint8([1 2 3 4 5]); end
if ~isfield(opt.h2o_2pool.resid,'f0')     || isempty(opt.h2o_2pool.resid.f0),     opt.h2o_2pool.resid.f0     = 0.10; end
if ~isfield(opt.h2o_2pool.resid,'f_max')  || isempty(opt.h2o_2pool.resid.f_max),  opt.h2o_2pool.resid.f_max  = 0.55; end
if ~isfield(opt.h2o_2pool.resid,'k_T')    || isempty(opt.h2o_2pool.resid.k_T),    opt.h2o_2pool.resid.k_T    = 0.10; end
if ~isfield(opt.h2o_2pool.resid,'k_time') || isempty(opt.h2o_2pool.resid.k_time), opt.h2o_2pool.resid.k_time = 0.10; end
if ~isfield(opt.h2o_2pool.resid,'Tref_K') || isempty(opt.h2o_2pool.resid.Tref_K), opt.h2o_2pool.resid.Tref_K = 373.15; end
if ~isfield(opt.h2o_2pool.resid,'tref_s') || isempty(opt.h2o_2pool.resid.tref_s), opt.h2o_2pool.resid.tref_s = 20000; end
if ~isfield(opt.h2o_2pool.resid,'q_abs_min') || isempty(opt.h2o_2pool.resid.q_abs_min), opt.h2o_2pool.resid.q_abs_min = 0.05; end
if ~isfield(opt.h2o_2pool.resid,'tol_frac')  || isempty(opt.h2o_2pool.resid.tol_frac),  opt.h2o_2pool.resid.tol_frac  = 1e-4; end

% Type normalize
opt.h2o_2pool.resid.apply_in_stage_ids = uint8(opt.h2o_2pool.resid.apply_in_stage_ids);

% ---- adsorption cumulative supply cap (optional) ----
% Wrapper may not provide this; default OFF.
if ~isfield(opt,'ads_supply') || isempty(opt.ads_supply), opt.ads_supply = struct(); end
if ~isfield(opt.ads_supply,'enable') || isempty(opt.ads_supply.enable), opt.ads_supply.enable = false; end
if ~isfield(opt.ads_supply,'keep_y_floor_inventory') || isempty(opt.ads_supply.keep_y_floor_inventory)
    opt.ads_supply.keep_y_floor_inventory = true;
end
% ---- adsorption shape defaults ----
if ~isfield(opt,'ads_shape') || isempty(opt.ads_shape), opt.ads_shape = struct(); end

% enable: tau를 주면 자동 ON, 아니면 OFF (wrapper/runner 모두 호환)
if ~isfield(opt.ads_shape,'enable') || isempty(opt.ads_shape.enable)
    opt.ads_shape.enable = isfield(opt.ads_shape,'tau_aw_ads_s') && isfinite(opt.ads_shape.tau_aw_ads_s) && (opt.ads_shape.tau_aw_ads_s > 0);
end

if ~isfield(opt.ads_shape,'tau_aw_ads_s') || isempty(opt.ads_shape.tau_aw_ads_s)
    opt.ads_shape.tau_aw_ads_s = 0;   % 0이면 사실상 shaping OFF
end
if ~isfield(opt.ads_shape,'use_RHcmd') || isempty(opt.ads_shape.use_RHcmd)
    opt.ads_shape.use_RHcmd = true;
end

% --- option1: CO2-eq aw reference during adsorption (light + stable) ---
if ~isfield(opt.ads_shape,'aw_CO2_mode') || isempty(opt.ads_shape.aw_CO2_mode)
    % 'ref_ads' = use fixed aw_ref for CO2 equilibrium in adsorption only
    % 'dynamic' = current behavior (aw_eff)
    opt.ads_shape.aw_CO2_mode = 'ref_ads';
end

if ~isfield(opt.ads_shape,'aw_ads_ref') || isempty(opt.ads_shape.aw_ads_ref)
    % if empty, core will compute from adsorption inlet RHcmd (recommended)
    opt.ads_shape.aw_ads_ref = NaN;
end

% ---- adsorption-only k shaping (OFF by default) ----
if ~isfield(opt,'ads_kinshape') || isempty(opt.ads_kinshape), opt.ads_kinshape = struct(); end
if ~isfield(opt.ads_kinshape,'enable') || isempty(opt.ads_kinshape.enable), opt.ads_kinshape.enable = false; end
if ~isfield(opt.ads_kinshape,'tau_k_s') || isempty(opt.ads_kinshape.tau_k_s), opt.ads_kinshape.tau_k_s = 600.0; end
if ~isfield(opt.ads_kinshape,'kboost0') || isempty(opt.ads_kinshape.kboost0), opt.ads_kinshape.kboost0 = 2.0; end


% ---- yCO2 per-step delta cap (inventory spike guard) ----
if ~isfield(opt.yco2,'dy_cap') || isempty(opt.yco2.dy_cap), opt.yco2.dy_cap = struct(); end
if ~isfield(opt.yco2.dy_cap,'enable') || isempty(opt.yco2.dy_cap.enable), opt.yco2.dy_cap.enable = true; end
if ~isfield(opt.yco2.dy_cap,'dy_max') || isempty(opt.yco2.dy_cap.dy_max), opt.yco2.dy_cap.dy_max = 0.02; end

% ---- spike stabilizers (pCO2 drive / aw_eff LPF) ----
if ~isfield(opt,'stab') || isempty(opt.stab), opt.stab = struct(); end
if ~isfield(opt.stab,'enable') || isempty(opt.stab.enable), opt.stab.enable = true; end
if ~isfield(opt.stab,'pdrv_lpf_enable') || isempty(opt.stab.pdrv_lpf_enable), opt.stab.pdrv_lpf_enable = true; end
if ~isfield(opt.stab,'pdrv_tau_s') || isempty(opt.stab.pdrv_tau_s), opt.stab.pdrv_tau_s = 20; end
if ~isfield(opt.stab,'pdrv_stagejump_reset') || isempty(opt.stab.pdrv_stagejump_reset), opt.stab.pdrv_stagejump_reset = true; end
if ~isfield(opt.stab,'aw_lpf_enable') || isempty(opt.stab.aw_lpf_enable), opt.stab.aw_lpf_enable = true; end
if ~isfield(opt.stab,'aw_tau_s') || isempty(opt.stab.aw_tau_s), opt.stab.aw_tau_s = 20; end
if ~isfield(opt.stab,'aw_stagejump_reset') || isempty(opt.stab.aw_stagejump_reset), opt.stab.aw_stagejump_reset = true; end

% ---- early-stop (Ads/Des only) ----
if ~isfield(opt,'early_stop') || isempty(opt.early_stop), opt.early_stop = struct(); end
if ~isfield(opt.early_stop,'enable') || isempty(opt.early_stop.enable), opt.early_stop.enable = false; end
if ~isfield(opt.early_stop,'eta_ads') || isempty(opt.early_stop.eta_ads), opt.early_stop.eta_ads = 0.95; end
if ~isfield(opt.early_stop,'eta_des') || isempty(opt.early_stop.eta_des), opt.early_stop.eta_des = 0.95; end
if ~isfield(opt.early_stop,'t_ads_min_s') || isempty(opt.early_stop.t_ads_min_s), opt.early_stop.t_ads_min_s = 2400; end
if ~isfield(opt.early_stop,'t_des_min_s') || isempty(opt.early_stop.t_des_min_s), opt.early_stop.t_des_min_s = 6000; end
if ~isfield(opt.early_stop,'Tgate_frac') || isempty(opt.early_stop.Tgate_frac), opt.early_stop.Tgate_frac = 0.95; end

% ---- CO2 residual-floor bridging fields (compat with older opt structs) ----
if ~isfield(opt.co2_2pool,'resid') || isempty(opt.co2_2pool.resid), opt.co2_2pool.resid = struct(); end
if ~isfield(opt.co2_2pool.resid,'q_abs_min') && isfield(opt.co2_2pool.resid,'q_min_abs')
    opt.co2_2pool.resid.q_abs_min = opt.co2_2pool.resid.q_min_abs;
end
if ~isfield(opt.co2_2pool.resid,'f_res')
    if isfield(opt.co2_2pool.resid,'frac_floor'), opt.co2_2pool.resid.f_res = opt.co2_2pool.resid.frac_floor; end
end
if ~isfield(opt.co2_2pool.resid,'f_res') || isempty(opt.co2_2pool.resid.f_res), opt.co2_2pool.resid.f_res = 0.35; end
if ~isfield(opt.co2_2pool.resid,'q_abs_min') || isempty(opt.co2_2pool.resid.q_abs_min), opt.co2_2pool.resid.q_abs_min = 0.03; end
% Ensure canonical stage selector + tolerance (for residual-floor functions)
if isfield(opt.co2_2pool.resid,'stage_ids') && ~isfield(opt.co2_2pool.resid,'apply_in_stage_ids')
    opt.co2_2pool.resid.apply_in_stage_ids = opt.co2_2pool.resid.stage_ids;
end
if ~isfield(opt.co2_2pool.resid,'apply_in_stage_ids') || isempty(opt.co2_2pool.resid.apply_in_stage_ids)
    opt.co2_2pool.resid.apply_in_stage_ids = uint8([2 3]);
end
if ~isfield(opt.co2_2pool.resid,'tol_frac') || isempty(opt.co2_2pool.resid.tol_frac)
    opt.co2_2pool.resid.tol_frac = 1e-4;
end
opt.co2_2pool.resid.apply_in_stage_ids = uint8(opt.co2_2pool.resid.apply_in_stage_ids);

% ---- adsorption kinetics-shape clamp ----
if ~isfield(opt,'ads_kinshape') || isempty(opt.ads_kinshape), opt.ads_kinshape = struct(); end
if ~isfield(opt.ads_kinshape,'fac_max') || isempty(opt.ads_kinshape.fac_max), opt.ads_kinshape.fac_max = 3.0; end

% ---- v14 additions: plotting substep + parasitic electricity defaults ----
% Publication-quality adsorption trace (post-processing helper)
if ~isfield(opt,'pub_trace') || isempty(opt.pub_trace), opt.pub_trace = struct(); end
if ~isfield(opt.pub_trace,'dt_max_s') || isempty(opt.pub_trace.dt_max_s), opt.pub_trace.dt_max_s = 10.0; end
if ~isfield(opt.pub_trace,'k_floor')  || isempty(opt.pub_trace.k_floor),  opt.pub_trace.k_floor  = 1e-6; end

% Per-cycle parasitic electricity (kWh/cycle)
if ~isfield(opt,'energy') || isempty(opt.energy), opt.energy = struct(); end
if ~isfield(opt.energy,'parasitic_kWh_per_cycle')     || isempty(opt.energy.parasitic_kWh_per_cycle),     opt.energy.parasitic_kWh_per_cycle = 0.0; end
if ~isfield(opt.energy,'parasitic_T_slope_kWh_per_C') || isempty(opt.energy.parasitic_T_slope_kWh_per_C), opt.energy.parasitic_T_slope_kWh_per_C = 0.0; end
if ~isfield(opt.energy,'parasitic_RH_slope_kWh_per_RH') || isempty(opt.energy.parasitic_RH_slope_kWh_per_RH), opt.energy.parasitic_RH_slope_kWh_per_RH = 0.0; end
if ~isfield(opt.energy,'parasitic_Tref_C')            || isempty(opt.energy.parasitic_Tref_C),            opt.energy.parasitic_Tref_C = 20.0; end
if ~isfield(opt.energy,'parasitic_RHref')             || isempty(opt.energy.parasitic_RHref),             opt.energy.parasitic_RHref  = 0.50; end

end

function opt = ensure_ptraj_defaults(opt)

if ~isfield(opt,'ptraj') || isempty(opt.ptraj)
    opt.ptraj = struct();
end
p = opt.ptraj;

% Core toggles (conservative defaults)
if ~isfield(p,'enable_smoothing'),              p.enable_smoothing = false; end
if ~isfield(p,'tau_flow_ramp_s'),              p.tau_flow_ramp_s = 120.0; end
if ~isfield(p,'reset_flow_filters_on_stage_change'), p.reset_flow_filters_on_stage_change = false; end

% n_dry residual filter knobs
if ~isfield(p,'tau_ndry_s'),                   p.tau_ndry_s = 180.0; end
if ~isfield(p,'enable_ndry_filter'),           p.enable_ndry_filter = true; end
if ~isfield(p,'ndry_reset_on_stage_change'),   p.ndry_reset_on_stage_change = true; end
if ~isfield(p,'log_ndry_residual'),            p.log_ndry_residual = false; end

% Qin/Qout rate-limit (stability guard)
if ~isfield(p,'enable_ndot_rate_limit'),       p.enable_ndot_rate_limit = true; end
if ~isfield(p,'Qmax_mult_Vg_dt'),              p.Qmax_mult_Vg_dt = 3.0; end

% Optional pressure-edge smoothing in driver construction
if ~isfield(p,'driver_edge_smooth_enable'),    p.driver_edge_smooth_enable = false; end
if ~isfield(p,'driver_edge_ramp_s'),           p.driver_edge_ramp_s = 30.0; end

opt.ptraj = p;
end


function opt = ensure_aw_defaults(opt)
% Ensure opt.aw exists with required fields (API-stability for wrappers)
    if ~isfield(opt,'aw') || ~isstruct(opt.aw)
        opt.aw = struct();
    end
    if ~isfield(opt.aw,'max') || isempty(opt.aw.max) || ~isfinite(opt.aw.max)
        opt.aw.max = 0.999; % cap to avoid aw=1 singularities
    end
    if ~isfield(opt.aw,'chain_mode') || strlength(string(opt.aw.chain_mode))==0
        % 'bed_consistent': use bed pH2O/psat for aw_gas
        % 'cmd_based'     : follow commanded RH where applicable (if implemented)
        opt.aw.chain_mode = 'bed_consistent';
    end
    if ~isfield(opt.aw,'cmd_based_cap_policy') || strlength(string(opt.aw.cmd_based_cap_policy))==0
        opt.aw.cmd_based_cap_policy = 'none';
    end
    if ~isfield(opt.aw,'mode_for_co2_effects') || strlength(string(opt.aw.mode_for_co2_effects))==0
        % used by aw_for_co2_effects(): 'gas'|'gab'|'mixed'
        opt.aw.mode_for_co2_effects = 'gas';
    end
    if ~isfield(opt.aw,'co2_boundmask') || ~isstruct(opt.aw.co2_boundmask)
        opt.aw.co2_boundmask = struct();
    end
    if ~isfield(opt.aw.co2_boundmask,'enable') || isempty(opt.aw.co2_boundmask.enable)
        opt.aw.co2_boundmask.enable = false;
    end

    % ---- Publication-quality adsorption trace (post-processing only) ----
    if ~isfield(opt,'pub_trace') || isempty(opt.pub_trace); opt.pub_trace = struct(); end
    if ~isfield(opt.pub_trace,'enable') || isempty(opt.pub_trace.enable)
        opt.pub_trace.enable = 1; % affects ONLY out.prof fields used for plotting/figures
    end
    if ~isfield(opt.pub_trace,'tau_smooth_s') || isempty(opt.pub_trace.tau_smooth_s)
        opt.pub_trace.tau_smooth_s = 180; % small -> light smoothing (seconds)
    end
    if ~isfield(opt.pub_trace,'enforce_monotonic') || isempty(opt.pub_trace.enforce_monotonic)
        opt.pub_trace.enforce_monotonic = 1;
    end

    if ~isfield(opt.aw.co2_boundmask,'power') || isempty(opt.aw.co2_boundmask.power)
        opt.aw.co2_boundmask.power = 1.0;
    end
end


% ========================================================================
% Publication-quality adsorption trace helper (post-processing only)
% ========================================================================
function [q_pub, app_pub] = make_pub_ads_trace(t, q, qeq, kvec, idx_ads, pubopt)
% Publication-quality adsorption trace by exponential reconstruction (LDF).
% - Does NOT change core integration.
% - Only densifies adsorption segment for plotting, using q(t) -> qeq(t) with k(t).

q_pub  = q;
app_pub = nan(size(q));

if isempty(idx_ads), return; end

% Defaults
dt_max = 10.0; % [s] substep max spacing for plotting
if isfield(pubopt,'dt_max_s') && ~isempty(pubopt.dt_max_s)
    dt_max = pubopt.dt_max_s;
end
k_floor = 1e-6;
if isfield(pubopt,'k_floor') && ~isempty(pubopt.k_floor)
    k_floor = pubopt.k_floor;
end

i0 = idx_ads(1);
i1 = idx_ads(end);

% Build densified vectors for adsorption only
t_dense  = t(1:i0);
q_dense  = q(1:i0);
app_dense = app_pub(1:i0);

for k = i0:(i1-1)
    t0 = t(k);  t1 = t(k+1);
    if ~(isfinite(t0) && isfinite(t1)) || t1 <= t0
        continue;
    end

    q0   = q_dense(end);              % start from already reconstructed tail
    qeq1 = qeq(k+1);                  % use next-step equilibrium (better visual)
    kk   = max(kvec(k+1), k_floor);   % use next-step k (more stable)

    nsub = max(1, ceil((t1 - t0)/dt_max));
    dt   = (t1 - t0)/nsub;

    for s = 1:nsub
        ts = t0 + s*dt;
        % LDF analytic step toward qeq1
        q0 = qeq1 + (q0 - qeq1)*exp(-kk*dt);

        t_dense(end+1,1) = ts; %#ok<AGROW>
        q_dense(end+1,1) = q0; %#ok<AGROW>

        denom = max(qeq1 - q(i0), 1e-12);
        app_dense(end+1,1) = (q0 - q(i0)) / denom; %#ok<AGROW>
    end
end

% Append the rest (after adsorption) as-is
t_tail = t(i1+1:end);
q_tail = q(i1+1:end);

t_pub_full = [t_dense; t_tail];
q_pub_full = [q_dense; q_tail];

% Return on original time grid by interpolation (keeps downstream unchanged)
q_pub = interp1(t_pub_full, q_pub_full, t, 'linear', 'extrap');

% Approach on original grid
app_pub = nan(size(t));
denom_all = max(qeq - q(i0), 1e-12);
app_pub(idx_ads) = (q_pub(idx_ads) - q(i0)) ./ denom_all(idx_ads);
end

function [t_hr_dense, q_dense] = build_ads_dense_trace_(t_hr, sid, q, qeq, k, dt_max_hr)
% build_ads_dense_trace_
% - Densify adsorption trace only (sid==6) using substepped integration:
%   dq/dt = k(t) * (qeq(t) - q)
% - k and qeq are linearly interpolated between coarse points.
%
% This is for plotting/diagnostics; it does NOT feed back into the core state.

t_hr = t_hr(:); sid = sid(:); q = q(:); qeq = qeq(:); k = k(:);
N = numel(t_hr);

t_hr_dense = t_hr(1);
q_dense    = q(1);

for i = 1:(N-1)
    dt = t_hr(i+1) - t_hr(i);
    if dt <= 0
        continue;
    end

    isAds = (sid(i) == 6);  % adsorption stage id assumed 6
    if ~isAds
        % keep coarse point only
        t_hr_dense(end+1,1) = t_hr(i+1);
        q_dense(end+1,1)    = q(i+1);
        continue;
    end

    nsub = max(1, ceil(dt / max(1e-9, dt_max_hr)));
    dt_sub = dt / nsub;

    qi   = q_dense(end);
    tcur = t_hr(i);

    for s = 1:nsub
        tnext = tcur + dt_sub;
        % linear interp of k and qeq between endpoints
        a = (tnext - t_hr(i)) / dt;  a = min(max(a,0),1);
        ksub   = (1-a)*k(i)   + a*k(i+1);
        qeqsub = (1-a)*qeq(i) + a*qeq(i+1);

        % one substep (semi-implicit Euler; stable and lightweight)
        qi = (qi + ksub*dt_sub*qeqsub) / (1 + ksub*dt_sub);

        t_hr_dense(end+1,1) = tnext;
        q_dense(end+1,1)    = qi;

        tcur = tnext;
    end
end

% ensure last point matches coarse end (avoid drift)
t_hr_dense(end) = t_hr(end);
q_dense(end)    = q(end);
end

