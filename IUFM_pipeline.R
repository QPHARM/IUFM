# IUFM analysis pipeline, (2026-11-17)
# Run from the folder that contains iufm_input/; results are written to iufm_output/.
# Base R only. Equation, table, figure, and section numbers refer to the article.

IN_DIR <- "iufm_input"; OUT_DIR <- "iufm_output"
dir.create(OUT_DIR, showWarnings = FALSE)
read_in   <- function(f) read.csv(file.path(IN_DIR, f), stringsAsFactors = FALSE)
write_out <- function(x, f) write.csv(x, file.path(OUT_DIR, f), row.names = FALSE)
bind_rows <- function(x) do.call(rbind, x)
last      <- function(x) x[length(x)]

settings <- read_in("in_settings.csv")
cfg <- function(key) { v <- settings$value[settings$key == key]; stopifnot(length(v) == 1); as.numeric(v) }
V0_uL <- cfg("V0_uL"); K_REL <- cfg("k_rel"); STEPS <- cfg("integration_steps")
WINDOW_PCT <- cfg("analysis_window_pct"); CONF <- cfg("conf_level")

dsmap   <- read_in("in_dataset_map.csv")
DSETS   <- dsmap$dataset
LABEL   <- setNames(dsmap$label, DSETS)
MATRIX  <- setNames(dsmap$system, DSETS)
SYSTEMS <- setNames(nm = unique(dsmap$system))
by_ds   <- function(df) split(df, factor(df$dataset, levels = DSETS))

# Weighed masses (mg) to volumes (uL)
density  <- read_in("in_density.csv")
rho      <- function(fr, ds) with(density[density$fraction == fr, ], density_g_per_mL[match(MATRIX[ds], system)])
filtrate <- read_in("in_filtrate.csv")
memdrug  <- read_in("in_membrane_drug.csv")
filtrate$V_filtrate_uL  <- filtrate$V_filtrate_uL  / rho("filtrate",  filtrate$dataset)
filtrate$V_retentate_uL <- filtrate$V_retentate_uL / rho("retentate", filtrate$dataset)
memdrug$pctF            <- memdrug$pctF            / rho("filtrate",  memdrug$dataset)

# ---- Device values (Eqs. 1-3) and condition means (Table S1) ----------------
cond <- lapply(by_ds(filtrate), function(d) {
  m <- function(v, f = mean) as.vector(tapply(v, d$time_min, f))
  data.frame(time_min = sort(unique(d$time_min)), pctF = m(d$V_filtrate_uL) / V0_uL * 100,
             C0 = m(d$C_initial_ng_mL), mean = m(d$C_filtrate_ng_mL), sd = m(d$C_filtrate_ng_mL, sd))
})
C0_ng_mL <- lapply(cond, function(x) mean(x$C0))
pctF     <- lapply(cond, function(x) x$pctF)
HI_PCT   <- sapply(pctF, max)

dev_tab <- bind_rows(lapply(by_ds(filtrate), function(d) {
  d  <- d[order(d$time_min, d$device), ]
  Cf <- d$C_filtrate_ng_mL; Cr <- d$C_retentate_ng_mL; Vf <- d$V_filtrate_uL; Vr <- d$V_retentate_uL
  data.frame(d[c("dataset", "time_min", "device")], FVR_pct = 100 * Vf / V0_uL,
             d[c("C_initial_ng_mL", "C_filtrate_ng_mL", "C_retentate_ng_mL", "V_filtrate_uL", "V_retentate_uL")],
             conventional_fu_pct = 100 * Cf / d$C_initial_ng_mL,
             mass_balance_fu_pct = 100 * Cf * (Vf + Vr) / (Cf * Vf + Cr * Vr),
             recovery_pct        = 100 * (Cf * Vf + Cr * Vr) / (d$C_initial_ng_mL * V0_uL), row.names = NULL)
}))
write_out(dev_tab, "out_device_values.csv")
write_out(bind_rows(lapply(split(dev_tab, list(factor(dev_tab$dataset, DSETS), dev_tab$time_min),
                                 drop = TRUE, lex.order = TRUE), function(g) {
  out <- data.frame(dataset = g$dataset[1], label = LABEL[[g$dataset[1]]], time_min = g$time_min[1], n = nrow(g))
  for (v in names(dev_tab)[-(1:3)]) out[paste0(v, c("_mean", "_sd"))] <- list(mean(g[[v]]), sd(g[[v]]))
  out
})), "out_tableS1.csv")

# ---- Filtrate volume, Eq. 4 (one fit per matrix) ----------------------------
fit_volume <- function(d) {
x <- data.frame(time = d$time_min, Fr = d$V_filtrate_uL / V0_uL)

  as.list(coef(nls(Fr ~ Fmax * (1 - exp(-kf * time)), data = x,
                   start = list(Fmax = max(x$Fr) * 1.1, kf = 0.04))))
}
vol_par <- lapply(SYSTEMS, function(s) fit_volume(filtrate[MATRIX[filtrate$dataset] == s, ]))

# ---- Membrane-associated protein, Eq. 9 -------------------------------------
protein  <- read_in("in_protein.csv")
dep_stat <- lapply(SYSTEMS, function(s) {
  d  <- protein[protein$system == s, ]
  D  <- ave(d$initial_mg, d$compound, d$time_min)
  ok <- !is.na(d$retentate_mg)
  y  <- ((D - d$retentate_mg) / D)[ok]; tm <- d$time_min[ok]
  f  <- function(g) as.vector(tapply(y, tm, g)); n <- f(length)
  data.frame(t = sort(unique(tm)), n = n, cum_mean = f(mean), cum_sd = f(sd),
             cum_sem = sqrt(sum((n - 1) * f(var)) / sum(n - 1)) / sqrt(n))
})

cum_frac_eq9 <- function(t_vec, k_dep, sigma_max, Fmax, kf, n = 2000) sapply(t_vec, function(tt) {
  if (tt <= 0) return(0)
  s  <- seq(0, tt, length.out = n)
  Ff <- Fmax * (1 - exp(-kf * s))
  dF <- Fmax * kf * exp(-kf * s)
  sg <- sigma_max * (1 - exp(-k_dep * s))
  1 - exp(-sum(sg * dF / (1 - Ff)) * (s[2] - s[1]))
})

fit_deposition <- function(obs, Fmax, kf) {
  w    <- function(k, sm) sum(((obs$cum_mean - cum_frac_eq9(obs$t, k, sm, Fmax, kf)) / obs$cum_sem)^2)
  opt  <- function(sm) optimize(function(lk) w(exp(lk), sm), c(log(1e-8), log(50)))
  prof <- function(sm) opt(sm)$objective
  pk   <- function(k) optimize(function(sm) w(k, sm), c(1e-4, 1))$objective
  o    <- optimize(prof, c(0.05, 1), tol = 1e-10)
  sm   <- if (prof(1) <= o$objective) 1 else o$minimum
  k    <- exp(opt(sm)$minimum)
  thr  <- w(k, sm) + qchisq(CONF, 1)
  cut  <- function(f, lo, hi, tol) uniroot(function(x) f(x) - thr, c(lo, hi), tol = tol)$root
  list(sigma_max    = sm,
       sigma_max_lo = if (prof(0.05) < thr) NA else cut(prof, 0.05, sm, 1e-9),
       sigma_max_hi = if (prof(1) < thr) 1 else cut(prof, sm, 1, 1e-9),
       k_dep        = k,
       k_dep_lo     = if (pk(k * 1e-3) < thr) 0 else cut(pk, k * 1e-3, k, 1e-12),
       k_dep_hi     = if (pk(50) < thr) Inf else cut(pk, k, 50, 1e-12),
       pred         = cum_frac_eq9(obs$t, k, sm, Fmax, kf))
}
dep_fit <- lapply(SYSTEMS, function(s) fit_deposition(dep_stat[[s]], vol_par[[s]]$Fmax, vol_par[[s]]$kf))
write_out(bind_rows(lapply(SYSTEMS, function(s)
  data.frame(system = s, dep_stat[[s]][c("t", "n", "cum_mean", "cum_sd")], pred = dep_fit[[s]]$pred))),
  "out_deposition_points.csv")

# ---- IUFM simulation, Eqs. 8-12 ---------------------------------------------
sites <- by_ds(read_in("in_binding_R1.csv"))
free_conc <- function(Ctot, P, KA, n)
  uniroot(function(Cu) Cu + sum(n * P * KA * Cu / (1 + KA * Cu)) - Ctot, c(0, Ctot), tol = 1e-15)$root

simulate_IUFM <- function(phi, k_rel, p, C0, steps = STEPS) {
  dt <- p$T_end / steps
  M_ret <- M_ret0 <- C0 / 1000 * 1e-3 / p$MW * p$V0
  P_tot <- sum(p$P0); frac <- p$P0 / P_tot; Mp0 <- P_tot * p$V0
  M_mem <- M_filt <- V_filt <- F_prev <- 0; V_ret <- p$V0
  out <- matrix(NA_real_, steps, 4, dimnames = list(NULL, c("pctF", "Cfilt", "recovery", "M_mem_ng")))
  for (i in seq_len(steps)) {
    t  <- i * dt
    Ff <- p$Fmax * (1 - exp(-p$kf * t))
    if (Ff * 100 > 60) break
    dV <- p$V0 * (Ff - F_prev); F_prev <- Ff
    Ctot <- M_ret / V_ret
    Cu   <- free_conc(Ctot, Mp0 / V_ret * frac, p$KA, p$n)
    perm <- Cu * dV
    dep  <- p$sigma_max * (1 - exp(-p$k_dep * t)) * (Ctot - Cu) * dV
    to_f <- phi * k_rel * M_mem * dt
    to_r <- (1 - phi) * k_rel * M_mem * dt
    M_ret  <- M_ret - perm - dep + to_r
    M_mem  <- M_mem + dep - to_f - to_r
    M_filt <- M_filt + perm + to_f
    V_ret  <- V_ret - dV; V_filt <- V_filt + dV
    out[i, ] <- c((1 - V_ret / p$V0) * 100, M_filt / V_filt * p$MW * 1e3 * 1000,
                  (M_ret + M_filt) / M_ret0 * 100, M_mem * p$MW * 1e9)
  }
  as.data.frame(out[!is.na(out[, 1]), , drop = FALSE])
}
make_pars <- function(k, sigma_max = dep_fit[[MATRIX[[k]]]]$sigma_max, k_dep = dep_fit[[MATRIX[[k]]]]$k_dep) {
  v <- vol_par[[MATRIX[[k]]]]; b <- sites[[k]]
  list(V0 = V0_uL * 1e-6, MW = b$MW_g_per_mol[1], T_end = 200, Fmax = v$Fmax, kf = v$kf,
       sigma_max = sigma_max, k_dep = k_dep, P0 = b$Cp0_M, KA = b$KA_per_M, n = b$n_sites)
}
simulate <- function(k, phi, k_rel = K_REL, pars = make_pars(k)) simulate_IUFM(phi, k_rel, pars, C0_ng_mL[[k]])
upto     <- function(s, pct) s[s$pctF <= pct, ]

# ---- phi by slope matching and model parameters (Table 1; Section 4.3.5) ----
obs_slope <- function(k) {
  d <- dev_tab[dev_tab$dataset == k, ]
  m <- lm(d$C_filtrate_ng_mL ~ I(ave(d$V_filtrate_uL, d$time_min) / V0_uL * 100))
  c(coef(m)[[2]], confint(m, level = CONF)[2, ])
}
sim_slope <- function(phi, k, k_rel, pars) {
  s <- simulate(k, phi, k_rel, pars); x <- pctF[[k]]
  stopifnot(max(s$pctF) >= max(x))
  coef(lm(approx(s$pctF, s$Cfilt, xout = x)$y ~ x))[[2]]
}
solve_phi <- function(k, target, k_rel = K_REL, pars = make_pars(k)) {
  f <- function(phi) sim_slope(phi, k, k_rel, pars) - target
  if (f(0) >= 0) return(0)
  if (f(1) <= 0) return(1)
  uniroot(f, c(0, 1), tol = 1e-6)$root
}
phi_tab <- bind_rows(lapply(DSETS, function(k) {
  o <- obs_slope(k)
  data.frame(dataset = k, phi = solve_phi(k, o[1]), ci_lo = solve_phi(k, o[2]), ci_hi = solve_phi(k, o[3]),
             slope = o[1])
}))
PHI   <- setNames(phi_tab$phi, DSETS)
SLOPE <- setNames(phi_tab$slope, DSETS)
write_out(bind_rows(lapply(DSETS, function(k) with(phi_tab[phi_tab$dataset == k, ],
  data.frame(dataset = k, label = LABEL[[k]], vol_par[[MATRIX[[k]]]], dep_fit[[MATRIX[[k]]]][1:6],
             k_rel = K_REL, phi, phi_lo = ci_lo, phi_hi = ci_hi)))), "out_parameters.csv")

# ---- Recovery and membrane pool (Section 2.5; Figure 5) ---------------------
S_FIT <- lapply(setNames(nm = DSETS), function(k) simulate(k, PHI[[k]]))
write_out(bind_rows(lapply(DSETS, function(k) {
  s   <- upto(S_FIT[[k]], HI_PCT[[k]])
  rec <- with(dev_tab[dev_tab$dataset == k, ], recovery_pct[time_min == max(time_min)])
  ci  <- t.test(rec, conf.level = CONF)$conf.int
  sim <- last(s$recovery); obs <- mean(rec)
  data.frame(dataset = k, label = LABEL[[k]], pctF = last(s$pctF), recovery_sim = sim, recovery_obs = obs,
             recovery_obs_lo = ci[1], recovery_obs_hi = ci[2], recovery_dev = sim - obs,
             Cfilt_ideal = last(upto(simulate(k, 0, pars = make_pars(k, k_dep = 0)), HI_PCT[[k]])$Cfilt))
})), "out_simulation_summary.csv")

mem_tab <- bind_rows(lapply(DSETS, function(k) {
  d   <- memdrug[memdrug$dataset == k, ]; pts <- sort(unique(d$pctF))
  at  <- function(f, v = d$extracted_ng) sapply(pts, function(x) f(v[d$pctF == x]))
  data.frame(dataset = k, label = LABEL[[k]], time_min = at(min, d$time_min), pctF = pts,
             pool_sim_ng = approx(S_FIT[[k]]$pctF, S_FIT[[k]]$M_mem_ng, xout = pts, rule = 2)$y,
             ext_mean_ng = at(mean), ext_max_ng = at(max))
}))
mem_tab$ratio_sim_ext <- mem_tab$pool_sim_ng / mem_tab$ext_mean_ng
write_out(mem_tab, "out_membrane_check.csv")

# ---- Admissible k_rel range (Section 4.3.5) ---------------------------------
root_or_na <- function(f, lo, hi) if (f(lo) * f(hi) > 0) NA_real_ else uniroot(f, c(lo, hi), tol = 1e-8)$root
krel_tab <- bind_rows(lapply(DSETS, function(k) {
  m   <- mem_tab[mem_tab$dataset == k, ]
  inc <- last(m$ext_max_ng) - m$ext_max_ng[1]
  gap <- function(kr) {
    s    <- simulate(k, solve_phi(k, SLOPE[[k]], kr), kr)
    pool <- approx(s$pctF, s$M_mem_ng, xout = m$pctF, rule = 2)$y
    last(pool) - pool[1] - inc
  }
  data.frame(dataset = k, label = LABEL[[k]],
             k_rel_min = root_or_na(function(kr) sim_slope(1, k, kr, make_pars(k)) - SLOPE[[k]], 1e-4, K_REL),
             k_rel_max = root_or_na(gap, K_REL, 5))
}))
KREL_ADM <- c(max(krel_tab$k_rel_min, na.rm = TRUE), min(krel_tab$k_rel_max, na.rm = TRUE))
write_out(rbind(krel_tab, data.frame(dataset = "all", label = "Admissible range",
                                     k_rel_min = KREL_ADM[1], k_rel_max = KREL_ADM[2])), "out_krel_range.csv")

# ---- Zero-volume extrapolation against FVR^2 (Eq. 14; Section 2.6; Note S2) --
ZV_WINDOWS <- sort(unique(c(25, WINDOW_PCT, 40)))
ZV_KREL    <- sort(unique(c(0.014, K_REL, 0.28)))
stopifnot(min(ZV_KREL) >= KREL_ADM[1], max(ZV_KREL) <= KREL_ADM[2])

write_out(bind_rows(lapply(DSETS[MATRIX[DSETS] == "plasma"], function(k) bind_rows(lapply(ZV_WINDOWS, function(win) {
  g  <- cond[[k]][cond[[k]]$pctF <= win, ]; c0 <- C0_ng_mL[[k]]
  x  <- (g$pctF / 100)^2; m <- lm(g$mean ~ x); ci <- confint(m, level = CONF)[1, ]
  data.frame(dataset = k, label = LABEL[[k]], window_pct = win, n_conditions = nrow(g), C0_ng_mL = c0,
             intercept_ng_mL = coef(m)[[1]], fu_pct = 100 * coef(m)[[1]] / c0,
             fu_lo_pct = 100 * ci[[1]] / c0, fu_hi_pct = 100 * ci[[2]] / c0,
             p_slope = summary(m)$coefficients[2, 4], fu_window_mean_pct = 100 * mean(g$mean) / c0)
})))), "out_zv_extrapolation.csv")

cu0_ng <- function(k) {
  p <- make_pars(k)
  free_conc(C0_ng_mL[[k]] * 1e-6 / p$MW, p$P0, p$KA, p$n) * p$MW * 1e6
}
write_out(bind_rows(lapply(DSETS, function(k) bind_rows(lapply(ZV_KREL, function(kr) {
  ph <- solve_phi(k, SLOPE[[k]], k_rel = kr); s <- simulate(k, ph, k_rel = kr)
  x  <- pctF[[k]][pctF[[k]] <= WINDOW_PCT]; y <- approx(s$pctF, s$Cfilt, xout = x)$y
  b  <- coef(lm(y ~ I((x / 100)^2)))[[1]]; cu <- cu0_ng(k)
  data.frame(dataset = k, label = LABEL[[k]], k_rel = kr, phi = ph, Cu0_ng_mL = cu,
             deviation_pct = (b - cu) / cu * 100)
})))), "out_zv_bias.csv")

# ---- Vertical offset of the simulated profiles in Figure 5A-D (Section 4.3.5)
write_out(bind_rows(lapply(DSETS, function(k) {
  cu <- cu0_ng(k)
  sh <- mean(cond[[k]]$mean) - mean(approx(S_FIT[[k]]$pctF, S_FIT[[k]]$Cfilt, xout = pctF[[k]])$y)
  data.frame(dataset = k, label = LABEL[[k]], Cu0_unshifted_ng_mL = cu, shift_ng_mL = sh, Cu0_shifted_ng_mL = cu + sh)
})), "out_baseline_shift.csv")

writeLines(capture.output(sessionInfo()), file.path(OUT_DIR, "sessionInfo.txt"))
