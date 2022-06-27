## ---- include = FALSE---------------------------------------------------------

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  dev = 'png',
  out.width = "90%"
)

## ----setup--------------------------------------------------------------------
library(hydrorecipes)
library(earthtide)
library(ggplot2)
library(splines2)
library(tidyr)

data(transducer, package = "hydrorecipes")
head(transducer)
transducer$datetime_num <- as.numeric(transducer$datetime)

## -----------------------------------------------------------------------------

rec_toll_rasmussen <- recipe(wl~baro+et+datetime_num, transducer) |>
  step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
  step_lead_lag(et, lag = seq(0, 180, 20)) |>
  step_ns(datetime_num, deg_free = 10) |>
  prep()
input_toll_rasmussen <- rec_toll_rasmussen |> bake(new_data = NULL)


## -----------------------------------------------------------------------------
fit_toll_rasmussen <- lm(wl~., input_toll_rasmussen) 

## -----------------------------------------------------------------------------
tidal_freqs <- c(0.89324406, 0.92953571, 0.96644626, 1.00273791, 1.03902956, 
                 1.07594011, 1.86454723, 1.89598197, 1.93227362, 1.96856526, 
                 2.00000000, 2.89841042)
rec_rasmussen_mote <- recipe(wl~baro+datetime_num, transducer) |>
  step_lead_lag(baro, lag = log_lags(100, 86400 * 2 / 120)) |>
  step_harmonic(datetime_num,
                frequency = tidal_freqs,
                cycle_size = 86400,
                keep_original_cols = TRUE) |>
  step_ns(datetime_num, deg_free = 10) |>
  prep()
input_rasmussen_mote <- rec_rasmussen_mote |> bake(new_data = NULL)


## -----------------------------------------------------------------------------
fit_rasmussen_mote <- lm(wl~., input_rasmussen_mote) 

## -----------------------------------------------------------------------------

wave_groups <- earthtide::eterna_wavegroups
wave_groups <- na.omit(wave_groups[wave_groups$time == '1 month', ])
wave_groups <- wave_groups[wave_groups$start > 0.5, ]
latitude  <- 34.0
longitude <- -118.5
rec_dl <- recipe(wl~baro+datetime_num, transducer) |>
  step_distributed_lag(baro, spline_fun = splines2::bSpline,
                       knots = log_lags(15, 86400 * 2 / 120)) |>
  step_earthtide(datetime_num,
                 latitude = latitude,
                 longitude = longitude,
                 astro_update = 1,
                 cutoff = 1e-5,
                 wave_groups = wave_groups) |>
  step_ns(datetime_num, deg_free = 10) |>
  prep()

input_dl <- rec_dl |> bake(new_data = NULL)


## -----------------------------------------------------------------------------
fit_dl <- lm(wl~., input_dl) 

## ----results = 'asis'---------------------------------------------------------

model_results <- rbind(broom::glance(fit_dl),
                       broom::glance(fit_toll_rasmussen),
                       broom::glance(fit_rasmussen_mote)
)
model_results$name <- c('kennel_2020',
                        'toll_rasmussen_2007',
                        'rasmussen_mote_2007'
)
knitr::kable(model_results[,c('name', 'r.squared', 'sigma', 
                              'AIC', 'BIC', 'df.residual', 'nobs')])

## ----fig.height = 11, fig.width = 12, out.width = '90%'-----------------------
pred <- predict_terms(fit = fit_dl, 
                      rec = rec_dl,
                      data = input_dl)
pred <- bind_cols(transducer[, c('datetime', 'wl')], pred)
pred$residuals <- pred$wl - pred$predicted
pred_long <- pivot_longer(pred, cols = !datetime)
levels <-c('intercept', 'ns_datetime_num', 'distributed_lag_baro',
           'earthtide_datetime_num', 'predicted', 'wl', 'residuals')
labels <- c('Intercept', 'Background trend', 'Barometric Component',
            'Earth Tide Component', 'Predicted', 'Water Pressure', 
            'Residuals (obs-mod)')
pred_long$name <- factor(pred_long$name, 
                         levels = levels,
                         labels = labels)
ggplot(pred_long, aes(x = datetime, y = value)) +
  geom_line() + 
  scale_y_continuous(labels = scales::comma) +
  scale_x_datetime(expand = c(0,0)) + 
  ggtitle('Water Level Decomposition Results') + 
  xlab("") + 
  facet_grid(name~., scales = 'free_y') + 
  theme_bw()

## ----fig.height = 5, fig.width = 5, out.width = '47%', fig.show='hold'--------
resp    <- response(fit_dl, rec_dl)
resp_ba <- resp[resp$name == 'cumulative', ]
resp_ba <- resp_ba[resp_ba$term == 'baro', ]

ggplot(resp_ba, aes(x = x * 120 / 3600, y = value)) +
  ggtitle('A: Barometric Loading Response') + 
  xlab('lag (hours)') +
  ylab('Cumulative Response') +
  scale_y_continuous(limits = c(0, 1)) +
  geom_line() + 
  theme_bw()
resp_et <- resp[resp$name %in% c('amplitude', 'phase'), ]
ggplot(resp_et, aes(x = x, xend = x, y = 0, yend = value)) +
  geom_segment() + 
  ggtitle('B: Earthtide Response') +
  xlab('Frequency (cycles per day)') +
  ylab('Phase (radians)   |   Amplitude (dbar)') +
  facet_grid(name~., scales = 'free_y') + 
  theme_bw()

