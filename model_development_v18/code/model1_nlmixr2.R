# Primary joint model re-estimated in nlmixr2 (FOCEI), independently of joint_popk_nlme.py.
#
# The script reads the deposited concentrations itself and writes its own model and likelihood; it shares no
# code with the Python estimator. Each drug follows a 1-compartment model with zero-order infusion at steady
# state (interval 8 h; duration = deposited category, read as hours), written in closed form for estimation.
# The closed form is checked against an ODE solution with steady-state dosing in rxode2 at the individual
# estimates of the primary fit.
#
# nlmixr2 has no syntax for a correlation between the residual errors of 2 endpoints. The model carries one
# (c, between ceftazidime and avibactam in the same sample), so the likelihood is written in its factorized
# form p(caz, avi) = p(caz) p(avi | caz): the avibactam prediction includes b (observed - predicted log
# ceftazidime) and its residual SD is the conditional SD s_avi. This is the same bivariate normal likelihood:
#     sigma_avi = sqrt(s_avi^2 + b^2 sigma_caz^2),   c = b sigma_caz / sigma_avi.
#
# Fits
#   1 primary model with block covariance matrices and with standard-normal deviates (same model, 2
#     parameterizations), each from the same neutral initial values (correlations and b 0.3) and restarted from
#     its own estimates until the OFV changes by less than 0.001
#   2 objective function at the Python estimates (no outer iterations)
#   3 the standardized form with the clearance correlation fixed at 0 and at 0.94, same starts and restarts (the
#     change in objective function against the estimated standardized fit is the likelihood-ratio statistic)
#
# Usage: Rscript model1_nlmixr2.R    writes ../outputs/model1_nlmixr2.csv
#        NLMIXR2_OUTER_OPT=nlminb selects another outer optimizer (default bobyqa, as in nlmixr2est)

suppressPackageStartupMessages(library(nlmixr2))

args <- commandArgs(FALSE)
here <- dirname(normalizePath(sub("^--file=", "", args[grep("^--file=", args)])))
md <- dirname(here)
dd <- file.path(md, "data_external", "dryad_Li2025_CRRT")

# ------------------------------------------------------------------ data ------
caz <- read.csv(file.path(dd, "Ceftazidime_concentration.csv"))
avi <- read.csv(file.path(dd, "Avibactam_concentration.csv"))
med <- read.csv(file.path(dd, "Medication_Information.csv"))
obs <- merge(caz[, c("subjectID", "Time_h", "caz_pre")], avi[, c("subjectID", "Time_h", "avi_pre")],
             by = c("subjectID", "Time_h"))
stopifnot(nrow(obs) == nrow(caz), nrow(obs) == nrow(avi), !anyNA(obs),
          all(med$caz_avi_dose_g == "2+0.5"), all(med$dosage_interval_h == 8),
          all(obs$Time_h >= 0 & obs$Time_h <= 8))
tinf <- setNames(as.numeric(med$infusion_duration_cat), med$subjectID)
stopifnot(all(as.character(obs$subjectID) %in% names(tinf)))
obs$TINF <- tinf[as.character(obs$subjectID)]

mk <- function(cmt, dv) data.frame(ID = obs$subjectID, TIME = obs$Time_h, EVID = 0, CMT = cmt, DV = dv,
                                   TINF = obs$TINF, LCAZ = log(obs$caz_pre))
dat <- rbind(mk("caz", log(obs$caz_pre)), mk("avi", log(obs$avi_pre)))
dat <- dat[order(dat$ID, dat$TIME, dat$CMT == "avi"), ]
rownames(dat) <- NULL
n_obs <- nrow(dat)
cat(sprintf("data: %d patients, %d observations (%d per drug)\n", length(unique(dat$ID)), n_obs, n_obs / 2))

# ----------------------------------------------------------------- models -----
# steady-state concentration within the dosing interval that starts at time 0 (dose 2,000 mg ceftazidime and
# 500 mg avibactam over TINF hours, every 8 hours)
CONC <- "
    k_caz <- cl_caz / v_caz
    k_avi <- cl_avi / v_avi
    e_caz <- exp(-k_caz * TINF)
    e_avi <- exp(-k_avi * TINF)
    c0_caz <- 2000 / TINF / cl_caz * (1 - e_caz) * exp(-k_caz * (8 - TINF)) / (1 - exp(-k_caz * 8))
    c0_avi <- 500 / TINF / cl_avi * (1 - e_avi) * exp(-k_avi * (8 - TINF)) / (1 - exp(-k_avi * 8))
    if (t <= TINF) {
      c_caz <- 2000 / TINF / cl_caz * (1 - exp(-k_caz * t)) + c0_caz * exp(-k_caz * t)
      c_avi <- 500 / TINF / cl_avi * (1 - exp(-k_avi * t)) + c0_avi * exp(-k_avi * t)
    } else {
      c_caz <- (2000 / TINF / cl_caz * (1 - e_caz) + c0_caz * e_caz) * exp(-k_caz * (t - TINF))
      c_avi <- (500 / TINF / cl_avi * (1 - e_avi) + c0_avi * e_avi) * exp(-k_avi * (t - TINF))
    }"
LIK <- "
    lp_caz <- log(c_caz)
    lp_avi <- log(c_avi) + b_res * (LCAZ - lp_caz)
    lp_caz ~ add(s_caz) | caz
    lp_avi ~ add(s_avi) | avi"

primary_model <- function(p) {
  txt <- sprintf("function() {
  ini({
    lcl_caz <- %.10f
    lcl_avi <- %.10f
    lv_caz <- %.10f
    lv_avi <- %.10f
    b_res <- %.10f
    eta_cl_caz + eta_cl_avi ~ c(%.10f, %.10f, %.10f)
    eta_v_caz + eta_v_avi ~ c(%.10f, %.10f, %.10f)
    s_caz <- %.10f
    s_avi <- %.10f
  })
  model({
    cl_caz <- exp(lcl_caz + eta_cl_caz)
    cl_avi <- exp(lcl_avi + eta_cl_avi)
    v_caz <- exp(lv_caz + eta_v_caz)
    v_avi <- exp(lv_avi + eta_v_avi)%s%s
  })
}", log(p$cl[1]), log(p$cl[2]), log(p$v[1]), log(p$v[2]), p$b,
    p$w[1]^2, p$rho * p$w[1] * p$w[2], p$w[2]^2,
    p$w[3]^2, p$rv * p$w[3] * p$w[4], p$w[4]^2, p$s[1], p$s[2], CONC, LIK)
  eval(parse(text = txt))
}

std_model <- function(p, rho_fixed = NA) {
  zr <- if (is.na(rho_fixed)) sprintf("%.10f", atanh(p$rho)) else sprintf("fix(%.10f)", atanh(rho_fixed))
  txt <- sprintf("function() {
  ini({
    lcl_caz <- %.10f
    lcl_avi <- %.10f
    lv_caz <- %.10f
    lv_avi <- %.10f
    lw_cl_caz <- %.10f
    lw_cl_avi <- %.10f
    lw_v_caz <- %.10f
    lw_v_avi <- %.10f
    z_rho <- %s
    z_rv <- %.10f
    b_res <- %.10f
    u_cl_caz ~ fix(1)
    u_cl_avi ~ fix(1)
    u_v_caz ~ fix(1)
    u_v_avi ~ fix(1)
    s_caz <- %.10f
    s_avi <- %.10f
  })
  model({
    r_cl <- tanh(z_rho)
    r_v <- tanh(z_rv)
    cl_caz <- exp(lcl_caz + exp(lw_cl_caz) * u_cl_caz)
    cl_avi <- exp(lcl_avi + exp(lw_cl_avi) * (r_cl * u_cl_caz + sqrt(1 - r_cl^2) * u_cl_avi))
    v_caz <- exp(lv_caz + exp(lw_v_caz) * u_v_caz)
    v_avi <- exp(lv_avi + exp(lw_v_avi) * (r_v * u_v_caz + sqrt(1 - r_v^2) * u_v_avi))%s%s
  })
}", log(p$cl[1]), log(p$cl[2]), log(p$v[1]), log(p$v[2]), log(p$w[1]), log(p$w[2]), log(p$w[3]),
    log(p$w[4]), zr, atanh(p$rv), p$b, p$s[1], p$s[2], CONC, LIK)
  eval(parse(text = txt))
}

# ------------------------------------------------------------ estimation -----
OUTER_OPT <- Sys.getenv("NLMIXR2_OUTER_OPT", "bobyqa")   # the nlmixr2est default, named for the record
run_fit <- function(f, label, outer = TRUE) {
  ctl <- foceiControl(print = 0, covMethod = if (outer) "r,s" else "", outerOpt = OUTER_OPT,
                      maxOuterIterations = if (outer) 5000L else 0L)
  t0 <- Sys.time()
  fit <- suppressWarnings(suppressMessages(nlmixr2(f, dat, est = "focei", control = ctl)))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  th <- fixef(fit)
  if ("lw_cl_caz" %in% names(th)) {
    w <- exp(unname(th[c("lw_cl_caz", "lw_cl_avi", "lw_v_caz", "lw_v_avi")]))
    rho <- tanh(unname(th["z_rho"]))
    rv <- tanh(unname(th["z_rv"]))
  } else {
    om <- fit$omega
    w <- unname(sqrt(diag(om)[c("eta_cl_caz", "eta_cl_avi", "eta_v_caz", "eta_v_avi")]))
    rho <- om["eta_cl_caz", "eta_cl_avi"] / (w[1] * w[2])
    rv <- om["eta_v_caz", "eta_v_avi"] / (w[3] * w[4])
  }
  s_caz <- unname(th["s_caz"])
  b <- unname(th["b_res"])
  sigma_avi <- sqrt(unname(th["s_avi"])^2 + b^2 * s_caz^2)
  cat(sprintf("%-42s OFV %.6f  rho %.6f  (%.0f s)\n", label, fit$objf, rho, secs))
  row <- data.frame(fit = label, ofv = fit$objf, cl_caz = exp(unname(th["lcl_caz"])),
                    cl_avi = exp(unname(th["lcl_avi"])), v_caz = exp(unname(th["lv_caz"])),
                    v_avi = exp(unname(th["lv_avi"])), omega_cl_caz = w[1], omega_cl_avi = w[2],
                    omega_v_caz = w[3], omega_v_avi = w[4], rho = unname(rho), r_v = unname(rv),
                    sigma_caz = s_caz, sigma_avi = sigma_avi, c = b * s_caz / sigma_avi, n_obs = n_obs,
                    seconds = secs, stringsAsFactors = FALSE)
  list(row = row, fit = fit)
}

neutral <- list(cl = c(3, 3), v = c(25, 25), w = c(0.2, 0.2, 0.25, 0.25), rho = 0.3, rv = 0.3, b = 0.3,
                s = c(0.15, 0.15))  # nlmixr2 scales each parameter by its initial value, so none starts at 0
as_start <- function(r) list(cl = c(r$cl_caz, r$cl_avi), v = c(r$v_caz, r$v_avi),
                             w = c(r$omega_cl_caz, r$omega_cl_avi, r$omega_v_caz, r$omega_v_avi),
                             rho = r$rho, rv = r$r_v, b = r$c * r$sigma_avi / r$sigma_caz,
                             s = c(r$sigma_caz, r$sigma_avi * sqrt(1 - r$c^2)))

# Every estimation starts from the neutral values and is restarted from its own estimates until the OFV changes by
# less than 0.001 (at most 4 restarts): the outer optimizer can stop early on the flat surface near the optimum.
fit_converged <- function(make, p0, label) {
  r <- run_fit(make(p0), label)
  first <- r$row$ofv
  n <- 0
  repeat {
    if (n == 4) break
    r2 <- run_fit(make(as_start(r$row)), label)
    n <- n + 1
    done <- abs(r2$row$ofv - r$row$ofv) < 0.001
    if (r2$row$ofv < r$row$ofv) r <- r2
    if (done) break
  }
  r$row$ofv_first_run <- first
  r$row$restarts <- n
  r
}

py <- read.csv(file.path(md, "outputs", "model1_final_parameters.csv"))
P <- setNames(py$estimate, py$parameter)
cr <- P[["corr_residual_caz_avi"]]
python <- list(cl = unname(P[c("CL_caz_L_h", "CL_avi_L_h")]), v = unname(P[c("V_caz_L", "V_avi_L")]),
               w = unname(P[c("omega_CL_caz", "omega_CL_avi", "omega_V_caz", "omega_V_avi")]),
               rho = P[["corr_CL_caz_avi"]], rv = P[["corr_V_caz_avi"]],
               b = cr * P[["sigma_prop_avi"]] / P[["sigma_prop_caz"]],
               s = c(P[["sigma_prop_caz"]], P[["sigma_prop_avi"]] * sqrt(1 - cr^2)))

f1 <- fit_converged(primary_model, neutral, "primary, block covariance")
f3 <- fit_converged(function(p) std_model(p), neutral, "primary, standardized deviates")
f4 <- fit_converged(function(p) std_model(p, 0), neutral, "standardized deviates, rho fixed at 0")
f5 <- fit_converged(function(p) std_model(p, 0.94), neutral, "standardized deviates, rho fixed at 0.94")
f2 <- run_fit(primary_model(python), "primary, at the Python estimates", outer = FALSE)
f2$row$ofv_first_run <- NA
f2$row$restarts <- NA

# ------------------------------------- closed form against an ODE solution -----
eb <- ranef(f1$fit)
th <- fixef(f1$fit)
ind <- data.frame(id = eb$ID, cl_caz = exp(th[["lcl_caz"]] + eb$eta_cl_caz), cl_avi = exp(th[["lcl_avi"]] + eb$eta_cl_avi),
                  v_caz = exp(th[["lv_caz"]] + eb$eta_v_caz), v_avi = exp(th[["lv_avi"]] + eb$eta_v_avi))
ind$TINF <- tinf[as.character(ind$id)]
ev_obs <- merge(data.frame(id = obs$subjectID, time = obs$Time_h, evid = 0, cmt = "caz", amt = 0, rate = 0, ii = 0,
                           ss = 0), ind, by = "id")
# each drug in its own 1-compartment ODE with one steady-state dose (ss = 1) at time 0
ode1 <- rxode2::rxode2("
    d/dt(a) <- -cl / v * a
    oc <- a / v")
ode_conc <- function(cl, v, dose) {
  ev <- do.call(rbind, lapply(seq_len(nrow(ind)), function(i) {
    times <- sort(obs$Time_h[obs$subjectID == ind$id[i]])
    rbind(data.frame(id = ind$id[i], time = 0, evid = 1, amt = dose, rate = dose / ind$TINF[i], ii = 8, ss = 1,
                     cl = cl[i], v = v[i]),
          data.frame(id = ind$id[i], time = times, evid = 0, amt = 0, rate = 0, ii = 0, ss = 0, cl = cl[i], v = v[i]))
  }))
  as.data.frame(rxode2::rxSolve(ode1, events = ev, atol = 1e-12, rtol = 1e-12, ssAtol = 1e-12, ssRtol = 1e-12))$oc
}
oc_caz <- ode_conc(ind$cl_caz, ind$v_caz, 2000)
oc_avi <- ode_conc(ind$cl_avi, ind$v_avi, 500)
cf <- rxode2::rxode2(CONC)
ev_cf <- ev_obs[order(ev_obs$id, ev_obs$time), ]
s_cf <- as.data.frame(rxode2::rxSolve(cf, events = ev_cf))
stopifnot(length(oc_caz) == nrow(s_cf), length(oc_avi) == nrow(s_cf), nrow(s_cf) == nrow(obs))
ode_rel <- max(abs(oc_caz / s_cf$c_caz - 1), abs(oc_avi / s_cf$c_avi - 1))
cat(sprintf("closed form vs ODE at the individual estimates: max relative difference %.3g\n", ode_rel))
stopifnot(ode_rel < 1e-6)

res <- do.call(rbind, list(f1$row, f3$row, f4$row, f5$row, f2$row))
res$dofv_vs_primary <- res$ofv - f3$row$ofv   # the same parameterization as the fixed-rho fits
num <- vapply(res, is.numeric, logical(1))
num[c("n_obs", "restarts")] <- FALSE
res[num] <- lapply(res[num], function(x) ifelse(is.na(x), "", sprintf("%.6f", x)))
res$restarts[is.na(res$restarts)] <- ""
res$ode_check_max_rel_diff <- c(sprintf("%.3g", ode_rel), "", "", "", "")
res$software <- sprintf("R %s; nlmixr2 %s; nlmixr2est %s; rxode2 %s; FOCEI, outer optimizer %s", getRversion(),
                        packageVersion("nlmixr2"), packageVersion("nlmixr2est"), packageVersion("rxode2"), OUTER_OPT)
out <- file.path(md, "outputs", "model1_nlmixr2.csv")
write.csv(res, out, row.names = FALSE, quote = TRUE)
cat("wrote", out, "\n")
