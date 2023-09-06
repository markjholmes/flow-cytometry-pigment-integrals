# packages
library("tidyverse")
library("fda.usc") # for approximate integration
library("RColorBrewer")
library("viridis")

#% LOAD DATA

# wd of data
wdir <- "data/"

# range of wavelengths to include throughout
lamrange <- 400:900

# load absorption data
abs_dat <- read.csv(paste0(wdir, "pigs-abs.csv")) %>%
    mutate(abs = ifelse(abs < 0, 0, abs)) %>%
    group_by(Pigment) %>%
    nest(dat = lambda:abs) %>%
    rowwise %>%
    mutate(dat = list(as.data.frame(
        approx(dat$lambda, dat$abs, xout = lamrange))
        )) %>%
    unnest(dat) %>%
    rename(Lambda = x, Absorption = y) %>%
    ungroup %>%
    mutate(Absorption = Absorption / max(Absorption, na.rm = TRUE))

# load flourescence data
em_dat <- read.csv(paste0(wdir, "pigs-em.csv")) %>%
    mutate(em = ifelse(em < 0, 0, em)) %>%
    group_by(pigment) %>%
    nest(dat = lambda:em) %>%
    rowwise %>%
    mutate(dat = list(as.data.frame(
        approx(dat$lambda, dat$em, xout = lamrange))
        )) %>%
    unnest(dat) %>%
    rename(Lambda = x, Flourescence = y, Pigment = pigment) %>%
    ungroup %>%
    mutate(Flourescence = Flourescence / max(Flourescence, na.rm = TRUE))

# merge
absem_dat <- full_join(abs_dat, em_dat) %>%
    pivot_longer(Absorption:Flourescence,
        names_to = "Type", values_to = "Value") %>%
    na.omit

# lasers
lasers_dat <- read.csv(paste0(wdir, "lasers.csv"))

# cytometer channels
channels_dat <- read.csv(paste0(wdir, "channels.csv")) %>%
    rowwise %>%
    mutate(ID = paste(Receptor, substring(Laser, 1, 1),
        collapse = "", sep = "-"))

# define colour palettes

# pigments
pig_cols <- setNames(brewer.pal(3, "Dark2"), unique(abs_dat$Pigment))

# channels
channels_cols <- setNames(viridis(nrow(channels_dat)), sort(channels_dat$ID))

#% PLOT PIGMENTS, LASERS, AND CHANNELS

ggplot() +
    # rectangle for channel bandwidths
    geom_rect(data = channels_dat,
        aes(xmin = LambdaCentre - LambdaPM,
            xmax = LambdaCentre + LambdaPM,
            fill = ID, ymin = 0, ymax = 2),
        alpha = .2, color = "grey50", lty = 2) +
    # vertical lines for channel lasers
    geom_vline(data = lasers_dat,
        aes(xintercept = Lambda),
        col = c("blue", "red"), lwd = 1) +
    # lines for different pigments abs and em
    geom_line(data = absem_dat,
        aes(x = Lambda, y = Value, lty = Type, col = Pigment),
        lwd = 1) +
    geom_label(data = lasers_dat, aes(x = Lambda, y = 1.1, label = Laser)) +
    coord_cartesian(ylim = c(0, 1.2), expand = 0) +
    scale_fill_manual(values = channels_cols) +
    scale_color_manual(values = pig_cols) +
    theme_bw() +
    labs(y = "Absorption/flourescence",
        fill = "Cytometer channel",
        lty = "Type") +
    scale_y_continuous(breaks = seq(0, 1, len = 5))

ggsave("outputs/lasers_channels_pigments.pdf",
    width = 8, height = 5, units = "in")

#% FOR EACH CYTO CHANNEL, COMPUTE THE EXCITATION OF EACH PIGMENT

channels_lasers <- left_join(channels_dat, lasers_dat) %>%
    left_join(abs_dat, relationship = "many-to-many") %>%
    dplyr::select(-LambdaCentre, -LambdaPM) %>%
    # bc this is a relative scale, it's useful now to scale by the maximum value
    ungroup %>%
    mutate(Absorption = Absorption / max(Absorption))

#% FOR EACH CHANNEL, COMPUTE THE ABSORPTION FOR EACH PIGMENT

# get the bandwidth as a logical vector
channels_bw <- channels_dat %>%
    mutate(
        lambda_min = LambdaCentre - LambdaPM,
        lambda_max = LambdaCentre + LambdaPM,
        lambda_range = list(lamrange),
        in_bandwidth = list(lambda_min:lambda_max),
        logical_lambda = list(lambda_range %in% in_bandwidth),
        Lambda = list(data.frame(
            Lambda = lambda_range, BW = logical_lambda))
        ) %>%
    dplyr::select(ID, Laser, Receptor, Lambda) %>%
    unnest(Lambda)

# combine that with the emission spectra
channels_abs <- full_join(channels_bw, em_dat,
        relationship = "many-to-many") %>%
    arrange(ID, Pigment, Lambda) %>%
    dplyr::filter(BW) %>%
    na.omit %>%
    # we don't care about the actual maximum flourescence, only the flourescence by channel
    group_by(ID) #%>%
    # mutate(Flourescence = Flourescence / max(Flourescence))

#% PLOTTING 

ggplot(channels_abs) +
    aes(x = Lambda, y = Flourescence, col = ID, lty = Pigment) +
    geom_line(linewidth = 1) +
    facet_wrap(.~Laser, nrow = 2) +
    geom_point(data = channels_abs %>%
        group_by(Pigment, ID) %>%
        dplyr::filter(Lambda == min(Lambda) | Lambda == max(Lambda)),
        pch = 1, size = 5) +
    theme_bw() +
    labs(title = "Proportion of maximum flourescence captured by the channel bands",
        caption = "Circles indicate the ends of the channel bands")

# approx integral of each channel
channel_pig_abs <- channels_abs %>%
    group_by(ID, Laser, Receptor, BW, Pigment) %>%
    summarise(int_ind = int.simpson2(Lambda, Flourescence)) %>%
    group_by(ID) %>%
    mutate(int_ind = int_ind / max(int_ind))

ggplot(channel_pig_abs) +
    aes(x = ID, y = int_ind, fill = Pigment) +
    geom_col(width = 0.25,
        position = position_dodge(width = 0.5, preserve = "single")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    labs(title = "Proportion of integral that is due to each pigment",
        caption = "This does not take differences in excitation into account")

#% SUMMARISE, FOR EACH CHANNEL, THE TOTAL ABSORPTION FOR EACH PIGMENT

channel_pig_totabs <- channel_pig_abs %>%
    group_by(ID, Laser, Receptor, BW) %>%
    summarise(int_tot = sum(int_ind)) %>%
    ungroup %>%
    mutate(int_tot = int_tot / max(int_tot))

#% FOR EACH CHANNEL, GET THE RELATIVE EXCITATION * ABSORPTION FOR EACH PIGMENT

# rescale
channel_pig_relabs <- left_join(channel_pig_abs, channel_pig_totabs) %>%
    rowwise %>%
    mutate(int_rel = int_ind * int_tot)

ggplot(channel_pig_relabs) +
    aes(x = ID, y = int_rel, fill = Pigment) +
    geom_col(width = 0.25, 
        position = position_dodge(width = 0.5, preserve = "single")) +
    scale_y_sqrt(expand = expansion(mult = c(0, 0.05))) +
    theme_bw()

# % CONSTRUCT FUNCTION

pigmentation <- function(channels...) {

}

# % TEST DATA TO SEE 