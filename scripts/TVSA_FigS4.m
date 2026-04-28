function TVSA_make_FigS6_S7_forward_memory_corev2()
% Generate SI Figure S6 and Figure S7 using TVSA_CorePhysics_v2
%
% Figure S6:
%   Effect of adsorption time at 20 C, 50% RH
%   (a) baseline adsorption time: t_ads = 8000 s
%   (b) shortened adsorption time: t_ads = 4000 s
%
% Figure S7:
%   Effect of regeneration severity at 20 C, 50% RH
%   (a) baseline regeneration: T_regen = 100 C, t_des = 20000 s
%   (b) mild regeneration:     T_regen = 90 C,  t_des = 16000 s
%
% No files are exported. Figures are displayed only.

clc; close all;
clear TVSA_CorePhysics_v2
clear functions

core_name = 'TVSA_CorePhysics_v2';
TVSA = @(varargin) feval(core_name, varargin{:});

%% ------------------------------------------------------------------------
% 1. Options
% -------------------------------------------------------------------------
opt = struct();
opt.dt_s      = 5.0;
opt.max_iters = 1;
opt.tol_rel   = 5e-4;
opt.omega     = 1.00;   % irrelevant for max_iters=1; kept explicit

opt.mode = "S2";
opt.driver_baseline = struct('mode','forced_only');

opt.flow_mode = 'Acontact';
opt.legacy.keepQ_with_Acontact = 1;

opt.aw = struct();
opt.aw.max = 0.999;
opt.aw.chain_mode = 'bed_consistent';
opt.aw.cmd_based_cap_policy = 'cap_to_1';
opt.aw.mode_for_co2_effects = 'invert_qH2O';

opt.pco2 = struct();
opt.pco2.ads_mode = 'inlet';
opt.pco2.regen_mode = 'inventory';
opt.pco2.blend_alpha_inventory = 0.20;

opt.yco2 = struct();
opt.yco2.y_floor = 1e-6;
opt.yco2.y_max   = 0.60;
opt.yco2.dy_cap  = struct('enable',1,'dy_max',0.010);

opt.ptraj = struct();
opt.ptraj.enable_smoothing = 1;
opt.ptraj.tau_flow_ramp_s  = 25;
opt.ptraj.Qmax_mult_Vg_dt  = 20;
opt.ptraj.log_ndry_residual = 1;
opt.ptraj.driver_edge_smooth_enable = 1;
opt.ptraj.driver_edge_ramp_s        = 60;   % short physical valve/flow transition smoothing
opt.ptraj.enable_ndry_filter = 1;
opt.ptraj.tau_ndry_s = struct('vac',120,'heat',180,'des',180, ...
                              'cool',120,'press',20,'ads',0);
opt.ptraj.ndry_reset_on_stage_change = 1;
opt.ptraj.reset_flow_filters_on_stage_change = 1;
opt.ptraj.enable_ndot_rate_limit = 1;

opt.stab = struct('enable',1, ...
    'pdrv_lpf_enable',1,'pdrv_tau_s',20,'pdrv_stagejump_reset',1, ...
    'aw_lpf_enable',1,'aw_tau_s',20,'aw_stagejump_reset',1);

opt.tauT_s = struct('vac',0,'heat',900,'des',2000, ...
                    'cool',600,'press',300,'ads',0);

opt.co2_2pool.enable = 1;
opt.co2_2pool.fslow = struct('f0_ads',0.18,'f0_regen',0.18, ...
    'k_aw_ads',0.10,'k_aw_reg',0.00,'f_min',0.05,'f_max',0.45);
opt.co2_2pool.kslow = struct('ratio',0.02,'kmin',2e-6,'kmax',1e-3);
opt.co2_2pool.resid = struct('enable',1,'mode','ads_scaled', ...
    'apply_in_stage_ids',uint8([2 3]), ...
    'f_res',0.08,'q_abs_min',0.05, ...
    'k_dry',0.0,'k_T',0.05,'k_wet',0.35,'n_wet',1.0, ...
    'Tref_K',373.15,'f_min',0.02,'f_max',0.25,'tol_frac',1e-4);

opt.h2o_2pool.enable = 1;
opt.h2o_2pool.fbound = struct('f0_ads',0.15,'f0_regen',0.22, ...
    'k_aw_ads',0.10,'k_aw_reg',0.02,'f_min',0.05,'f_max',0.60);
opt.h2o_2pool.kslow  = struct('ratio',0.03,'kmin',5e-7,'kmax',5e-3);
opt.h2o_2pool.resid  = struct('enable',1, ...
    'apply_in_stage_ids',uint8([1 2 3 4 5]), ...
    'f0',0.10,'f_max',0.55,'k_T',0.10,'k_time',0.10, ...
    'Tref_K',373.15,'tref_s',20000,'q_abs_min',0.05,'tol_frac',1e-4);

% Keep adsorption kinetics finite and stage-independent.
% This avoids artificial entry acceleration after regeneration while preserving
% the same corephysics_v2 equations and cycle-integrated states.
opt.ads_kinshape = struct('enable',0,'kboost0',1.0,'tau_k_s',400,'fac_max',1.0);
opt.co2_awfloor = struct('enable',1,'apply_below',0.10,'value',0.05);
% Mild dry-bed CO2 enhancement for publication-control plots.
% Strong dry-bed enhancement can create an artificial adsorption-entry overshoot
% after low-temperature regeneration, so it is weakened and applied only during
% adsorption, not pressurization.
% Disable dry-bed CO2 enhancement for SI control-variation plots.
% The mechanistic co-adsorption isotherm already accounts for H2O/CO2 coupling;
% keeping an additional dry-bed multiplier can create an artificial adsorption-entry
% overshoot after mild regeneration.
opt.co2_dry_scale = struct('enable',0,'factor',1.00,'aw0',0.15,'p',1.0,'apply_stages',[6]);
opt.h2o_T_ads_scale = struct('enable',1,'TrefC',20,'dT_spanC',30,'k_gain',0.45);

opt.early_stop = struct('enable',0);
opt.pub_trace = struct('enable',0);

%% ------------------------------------------------------------------------
% 2. Baseline control
% -------------------------------------------------------------------------
ctrl = struct();
ctrl.u_feed_m_s   = 0.028;
ctrl.p_vac_bar    = 0.05;
ctrl.T_regenC     = 100;

ctrl.time_vac_s   = 60;
ctrl.time_heat_s  = 2400;
ctrl.time_des_s   = 20000;
ctrl.time_cool_s  = 400;
ctrl.time_press_s = 60;
ctrl.time_ads_s   = 8000;

%% ------------------------------------------------------------------------
% 3. Design
% -------------------------------------------------------------------------
design = struct();
design.width_cell_m  = 1.43;
design.height_cell_m = 0.10;
design.N_cell        = 13*88;

design.L_flow_fixed_m = 0.0172;
design.rb_fixed       = 1000/2.806;
design.eps            = 0.56;
design.dp_m           = 7.50e-4;

design.A_mult         = 2.1;
design.Ng_ref         = 35;
design.dP_misc_Pa     = 120;
design.N_series_pass  = 1;

design = TVSA('finalize_design', design);

%% ------------------------------------------------------------------------
% 4. Physics
% -------------------------------------------------------------------------
phys = TVSA('make_phys', design);
phys.iso.mech.S_cap = 1.1;
phys.kin.kCO2_ref = phys.kin.kCO2_ref * sqrt(design.A_mult);
phys.kin.kCO2_ads_mult = 1.00;

%% ------------------------------------------------------------------------
% 5. Shared ambient condition
% -------------------------------------------------------------------------
amb = struct();
amb.altitude_m = 0;
amb.T_ambC     = 20;
amb.RH_in      = 0.50;
amb.yCO2_air   = 400e-6;

%% ------------------------------------------------------------------------
% 6. Shared y-axis limits
% -------------------------------------------------------------------------
ylims = struct();
ylims.qCO2 = [0 1.6];
ylims.qH2O = [0 14.0];
ylims.T    = [0 105];
ylims.RH   = [0 1.05];
ylims.P    = [0 1.05];

%% ------------------------------------------------------------------------
% 7. Single combined figure: baseline, short adsorption, mild regeneration
% -------------------------------------------------------------------------
ctrl_ads_short = ctrl;
ctrl_ads_short.time_ads_s = 4000;

% Publication-safe mild-regeneration case.
% This represents lower regeneration severity without pushing the 0D model into
% an extreme dry-bed / residual-CO2 mismatch state at adsorption entry.
ctrl_reg_mild = ctrl;
ctrl_reg_mild.T_regenC     = 90;
ctrl_reg_mild.time_des_s   = 16000;
ctrl_reg_mild.time_press_s = 600;

COMBO_cases = struct([]);

COMBO_cases(1).panel = '';
COMBO_cases(1).name  = 'Baseline';
COMBO_cases(1).label = 'T_{regen}=100 ^{\circ}C, t_{ads}=8000 s';
COMBO_cases(1).ctrl  = ctrl;

COMBO_cases(2).panel = '';
COMBO_cases(2).name  = 'Short adsorption';
COMBO_cases(2).label = 't_{ads}=4000 s';
COMBO_cases(2).ctrl  = ctrl_ads_short;

COMBO_cases(3).panel = '';
COMBO_cases(3).name  = 'Mild regeneration';
COMBO_cases(3).label = 'T_{regen}=90 ^{\circ}C, t_{des}=16000 s, t_{press}=600 s';
COMBO_cases(3).ctrl  = ctrl_reg_mild;

COMBO_sims = run_control_cases_forward_memory(TVSA, amb, COMBO_cases, design, phys, opt);

plot_control_comparison_5panel_1x3( ...
    COMBO_sims, COMBO_cases, ylims, ...
    'Single-cycle forward perturbation from baseline end state');

end


%% =========================================================================
function sims = run_control_cases_forward_memory(TVSA, amb, cases, design, phys, opt)
% Run one forward cycle for each control case.
%
% Figure logic:
%   1) No periodic steady-state iteration is imposed (opt.max_iters = 1).
%   2) Baseline is initialized from an adsorption-loaded state and run once.
%   3) Short-adsorption and mild-regeneration cases start from the baseline
%      end-of-cycle state, preserving cycle-to-cycle memory while changing
%      only the control variable of interest.

sims = cell(numel(cases),1);

% --- Baseline: one forward cycle from loaded state ------------------------
ctrl_base = cases(1).ctrl;
opt_base  = opt;
phys_base = phys;

drv_base = TVSA('get_driver', amb, ctrl_base, opt_base, opt_base.driver_baseline);

x0_loaded = make_loaded_initial_state_by_adsorption_preconditioning( ...
    TVSA, drv_base, amb, ctrl_base, design, phys_base, opt_base);

% Cycle starts with vacuum after the previous adsorption step. Keep the solid
% loading memory and reset gas/bed variables to ambient-end conditions.
x0_loaded(5) = amb.yCO2_air;
x0_loaded(6) = amb.T_ambC + 273.15;

sim_base = TVSA('run_point_periodic', ...
    drv_base, amb, ctrl_base, design, phys_base, opt_base, x0_loaded);

sims{1} = sim_base;
print_case_summary(cases(1), sim_base);

% Shared cycle-memory state for control perturbation cases.
x0_memory = sim_base.x_end(:);
x0_memory(1:4) = max(x0_memory(1:4), 0);
x0_memory(5)   = amb.yCO2_air;
x0_memory(6)   = amb.T_ambC + 273.15;

% --- Perturbation cases: one forward cycle from baseline end state --------
for i = 2:numel(cases)

    ctrl_i = cases(i).ctrl;

    opt_i = opt;
    if isfield(cases(i),'opt_override') && ~isempty(cases(i).opt_override)
        opt_i = cases(i).opt_override;
    end
    opt_i.max_iters = 1;

    phys_i = phys;
    if isfield(cases(i),'phys_override') && ~isempty(cases(i).phys_override)
        phys_i = cases(i).phys_override;
    end

    drv_i = TVSA('get_driver', amb, ctrl_i, opt_i, opt_i.driver_baseline);

    sim_i = TVSA('run_point_periodic', ...
        drv_i, amb, ctrl_i, design, phys_i, opt_i, x0_memory);

    sims{i} = sim_i;
    print_case_summary(cases(i), sim_i);
end

end

%% =========================================================================
function print_case_summary(caseinfo, sim)

fprintf('\n%s | %s\n', caseinfo.name, caseinfo.label);
fprintf('  qCO2 start = %.4f | qCO2 end = %.4f mol/kg\n', ...
    sim.prof.qCO2_total(1), sim.prof.qCO2_total(end));
fprintf('  qH2O start = %.4f | qH2O end = %.4f mol/kg\n', ...
    sim.prof.qH2O_total(1), sim.prof.qH2O_total(end));
fprintf('  Tbed start = %.2f C | Tbed end = %.2f C\n', ...
    sim.prof.TC_bed(1), sim.prof.TC_bed(end));

end

%% =========================================================================
function x0_loaded = make_loaded_initial_state_by_adsorption_preconditioning( ...
    TVSA, drv_full, amb, ctrl, design, phys, opt)

drv_ads = extract_single_stage_driver(drv_full, uint8(6), 'Adsorption', opt.dt_s);

x0_clean = [
    0.0;
    0.0;
    0.0;
    0.0;
    amb.yCO2_air;
    amb.T_ambC + 273.15
];

sim_ads = TVSA('run_point_periodic', drv_ads, amb, ctrl, design, phys, opt, x0_clean);

x0_loaded = sim_ads.x_end(:);
x0_loaded(1:4) = max(x0_loaded(1:4), 0);

end

%% =========================================================================
function drv_stage = extract_single_stage_driver(drv, sid_keep, stage_name, dt_s)

idx = find(uint8(drv.stage_id(:)) == sid_keep);

if isempty(idx)
    error('Requested stage id %d not found in driver.', sid_keep);
end

fields_to_copy = fieldnames(drv);
drv_stage = struct();

for i = 1:numel(fields_to_copy)
    fld = fields_to_copy{i};
    val = drv.(fld);

    if isnumeric(val) || islogical(val)
        if isvector(val) && numel(val) == numel(drv.t_s)
            drv_stage.(fld) = val(idx);
        else
            drv_stage.(fld) = val;
        end
    else
        drv_stage.(fld) = val;
    end
end

N = numel(idx);
drv_stage.t_s = (0:N-1)' * dt_s;
drv_stage.stage_id = uint8(6 * ones(N,1));

drv_stage.TK_cmd = drv.TK_cmd(idx);
drv_stage.TC_cmd = drv_stage.TK_cmd - 273.15;
drv_stage.RH_cmd = drv.RH_cmd(idx);
drv_stage.P_Pa   = drv.P_Pa(idx);
drv_stage.Pa     = drv_stage.P_Pa;
drv_stage.Pbar   = drv_stage.P_Pa / 1e5;

drv_stage.edges_s = [0; N*dt_s];
drv_stage.stage_names = {stage_name};

if ~isfield(drv_stage,'meta') || isempty(drv_stage.meta)
    drv_stage.meta = struct();
end
drv_stage.meta.source = "adsorption_preconditioning_driver";

end

%% =========================================================================
%% =========================================================================
function plot_control_comparison_5panel(sims, cases, ylims, fig_title)

fig = figure('Color','w','Units','inches','Position',[0.8 0.8 9.0 7.2]);
tl = tiledlayout(5,2,'TileSpacing','compact','Padding','compact');

for j = 1:2
    sim = sims{j};
    caseinfo = cases(j);
    vars = unpack_sim_vars(sim);

    for r = 1:5
        ax = nexttile((r-1)*2 + j); hold on;

        switch r
            case 1
                plot(vars.t_hr, vars.qCO2, '-', 'Color',[0.00 0.28 0.75], 'LineWidth',1.45);
                ylim(ylims.qCO2);
                if j == 1, ylabel('q_{CO_2} [mol kg^{-1}]'); end
                if j == 2, legend('q_{CO_2}','Location','northeast','Box','off'); end

            case 2
                plot(vars.t_hr, vars.qH2O, '-', 'Color',[0.85 0.28 0.05], 'LineWidth',1.45);
                ylim(ylims.qH2O);
                if j == 1, ylabel('q_{H_2O} [mol kg^{-1}]'); end
                if j == 2, legend('q_{H_2O,total}','Location','northeast','Box','off'); end

            case 3
                plot(vars.t_hr, vars.Tcmd, '--', 'Color',[0.65 0.65 0.65], 'LineWidth',1.1);
                plot(vars.t_hr, vars.Tbed, '-', 'Color',[0.90 0.05 0.02], 'LineWidth',1.45);
                ylim(ylims.T);
                if j == 1, ylabel('T [^{\circ}C]'); end
                if j == 2, legend('T_{cmd}','T_{bed}','Location','southeast','Box','off'); end

            case 4
                plot(vars.t_hr, vars.RHcmd, '-', 'Color',[0.00 0.55 0.15], 'LineWidth',1.35);
                ylim(ylims.RH);
                if j == 1, ylabel('RH [-]'); end
                if j == 3, legend('RH','Location','northeast','Box','off'); end

            case 5
                plot(vars.t_hr, vars.Pbar, '-', 'Color',[0.10 0.10 0.10], 'LineWidth',1.45);
                ylim(ylims.P);
                xlabel('Time [h]');
                if j == 1, ylabel('P [bar]'); end
        end

        xlim([vars.t_hr(1) vars.t_hr(end)]);
        style_axis(ax);
        add_stage_background(ax, vars.edges_hr, vars.stage_names);

        if r == 1
            title(sprintf('%s %s\\newline%s', ...
                caseinfo.panel, caseinfo.name, caseinfo.label), ...
                'FontName','Arial','FontSize',8.5,'FontWeight','normal');
        end
    end
end

sgtitle(fig_title, 'FontName','Arial','FontSize',10,'FontWeight','normal');
set(fig,'Renderer','painters');

end

%% =========================================================================
%% =========================================================================
function plot_control_comparison_5panel_1x3(sims, cases, ylims, fig_title)

ncol = 3;
fig = figure('Color','w','Units','inches','Position',[0.4 0.4 12.0 7.2]);
tiledlayout(5,ncol,'TileSpacing','compact','Padding','compact');

for j = 1:ncol
    sim = sims{j};
    caseinfo = cases(j);
    vars = unpack_sim_vars(sim);

    for r = 1:5
        ax = nexttile((r-1)*ncol + j); hold on;
        switch r
            case 1
                plot(vars.t_hr, vars.qCO2, '-', 'Color',[0.00 0.28 0.75], 'LineWidth',1.45);
                ylim(ylims.qCO2);
                if j == 1, ylabel('q_{CO_2} [mol kg^{-1}]'); end
            case 2
                plot(vars.t_hr, vars.qH2O, '-', 'Color',[0.85 0.28 0.05], 'LineWidth',1.45);
                ylim(ylims.qH2O);
                if j == 1, ylabel('q_{H_2O} [mol kg^{-1}]'); end
            case 3
                plot(vars.t_hr, vars.Tcmd, '--', 'Color',[0.65 0.65 0.65], 'LineWidth',1.1);
                plot(vars.t_hr, vars.Tbed, '-', 'Color',[0.90 0.05 0.02], 'LineWidth',1.45);
                ylim(ylims.T);
                if j == 1, ylabel('T [^{\circ}C]'); end
            case 4
                plot(vars.t_hr, vars.RHcmd, '-', 'Color',[0.00 0.55 0.15], 'LineWidth',1.35);
                ylim(ylims.RH);
                if j == 1, ylabel('RH [-]'); end
            case 5
                plot(vars.t_hr, vars.Pbar, '-', 'Color',[0.10 0.10 0.10], 'LineWidth',1.45);
                ylim(ylims.P);
                xlabel('Time [h]');
                if j == 1, ylabel('P [bar]'); end
        end

        xlim([vars.t_hr(1) vars.t_hr(end)]);
        style_axis(ax);
        add_stage_background(ax, vars.edges_hr, vars.stage_names);

        if r == 1
            title(sprintf('%s\newline%s', caseinfo.name, caseinfo.label), ...
                'FontName','Arial','FontSize',8.5,'FontWeight','normal');
        end
    end
end

sgtitle(fig_title, 'FontName','Arial','FontSize',10,'FontWeight','normal');
set(fig,'Renderer','painters');

end

function vars = unpack_sim_vars(sim)

vars.t_hr = sim.prof.time_s(:) / 3600;

% Use raw qCO2 from corephysics_v2. No plotting-only correction is applied.
vars.qCO2 = sim.prof.qCO2_total(:);

vars.qH2O  = smoothdata(sim.prof.qH2O_total(:), 'movmedian', 9);
vars.Tbed  = sim.prof.TC_bed(:);
vars.Tcmd  = sim.prof.TK_cmd(:) - 273.15;
vars.RHcmd = sim.prof.RH_cmd(:);
vars.Pbar  = sim.drv.P_Pa(:) / 1e5;

vars.edges_hr = sim.drv.edges_s(:) / 3600;
vars.stage_names = sim.drv.stage_names;

end

%% =========================================================================
function add_stage_background(ax, edges_hr, stage_names)

yl = ylim(ax);

colors = [
    1.00 0.92 0.94;  % Vacuum
    1.00 0.96 0.90;  % Heating
    0.98 0.95 0.88;  % Desorption
    0.93 0.97 1.00;  % Cooling
    0.92 0.95 0.98;  % Pressurization
    0.93 0.97 1.00   % Adsorption
];

for s = 1:min(numel(stage_names), numel(edges_hr)-1)

    x1 = edges_hr(s);
    x2 = edges_hr(s+1);
    dur = x2 - x1;

    p = patch(ax, [x1 x2 x2 x1], [yl(1) yl(1) yl(2) yl(2)], ...
        colors(min(s,size(colors,1)),:), ...
        'EdgeColor','none', ...
        'FaceAlpha',0.28, ...
        'HandleVisibility','off');
    uistack(p,'bottom');

    if dur >= 0.20
        text(ax, mean([x1 x2]), yl(2) - 0.045*(yl(2)-yl(1)), stage_names{s}, ...
            'HorizontalAlignment','center', ...
            'VerticalAlignment','top', ...
            'FontName','Arial', ...
            'FontSize',6.8, ...
            'FontWeight','normal', ...
            'Color',[0.05 0.05 0.05], ...
            'HandleVisibility','off');
    end
end

ylim(ax, yl);

end

%% =========================================================================
function style_axis(ax)

set(ax, ...
    'FontName','Arial', ...
    'FontSize',8, ...
    'LineWidth',0.75, ...
    'Box','on', ...
    'TickDir','out', ...
    'XGrid','on', ...
    'YGrid','on', ...
    'GridAlpha',0.12, ...
    'Layer','top');

end
