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

channels_lasers <- left_join(lasers_dat, abs_dat)

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
    labs(title = "Proportion of max. flourescence captured by the channels",
        caption = "Circles indicate the ends of the channel bands")

# approx integral of each channel
channel_pig_abs <- channels_abs %>%
    group_by(ID, Laser, Receptor, BW, Pigment) %>%
    summarise(int_ind = int.simpson2(Lambda, Flourescence))# %>%
    # group_by(ID) %>%
    # mutate(int_ind = int_ind / sum(int_ind))

ggplot(channel_pig_abs) +
    aes(x = ID, y = int_ind, fill = Pigment) +
    geom_col(width = 0.25,
        position = position_dodge(width = 0.5, preserve = "single")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    labs(title = "Proportion of integral that is due to each pigment",
        caption = "This does not take differences in excitation into account")

#% FOR EACH CHANNEL, GET THE RELATIVE EXCITATION * ABSORPTION FOR EACH PIGMENT

# rescale
channel_pig_relabs <- left_join(channel_pig_abs, channels_lasers) %>%
    rename(Excitation = Absorption, Flourescence = int_ind) %>%
    rowwise %>%
    mutate(Relative_Flourescence = Flourescence * Excitation)

# this is disgusting
RA_max <- max(channel_pig_relabs$Relative_Flourescence)

channel_pig_relabs$Relative_Flourescence <- 
    channel_pig_relabs$Relative_Flourescence / RA_max

ggplot(channel_pig_relabs) +
    aes(x = ID, y = Relative_Flourescence, fill = Pigment) +
    geom_col(width = 0.5,
        position = position_dodge(width = 0.5, preserve = "single")) +
    theme_bw() +
    scale_y_sqrt(expand = expansion(mult = c(0, 0.05)))

ggsave("outputs/flourescence_scaled_by_excitation.pdf",
    width = 8, height = 5, units = "in")

# basically we're going to use this to compute weighted means
    # e.g., chlorophyll is the mean chlorophyll across all channels weighted by
    # how much it flouresces in each channel

# final data for use in model
weights <- channel_pig_relabs %>%
    rename(channel = ID, pigment = Pigment, weight = Relative_Flourescence) %>%
    ungroup %>%
    dplyr::select(channel, pigment, weight) %>%
    tidyr::complete(channel, pigment, fill = list(weight = 0)) %>%
    arrange(channel, pigment) %>%
    pivot_wider(names_from = pigment, values_from = weight)

write.csv(weights, "data/weights.csv", row.names = FALSE)

pdf("outputs/heatmap.pdf", width = 9, height = 7)

column_to_rownames(.data = weights, "channel") %>%
    as.matrix %>% 
    t %>%
    heatmap(margins = c(10, 10))

dev.off()

# % DEMONSTRATION

# data to test on
dat <- readRDS("data/final-pops-mon.RData") %>%
    dplyr::filter(date < 11) %>%
    dplyr::select(date.time, date, strain, species,
            treat, repl, pop, contains("av.")) %>%
        rename_with(~ gsub("av.", "", .x, fixed = TRUE),
            starts_with("av")) %>%
    dplyr::select(date.time, date, strain, species,
        treat, repl, pop, FSC, SSC, sort(colnames(.))) %>%
    # do i need to standardise across the channels?
    rowwise %>%
    mutate(Chlorophyll = weighted.mean(
            c(GRN.B, NIR.B, NIR.R, RED.B, RED.R, YEL.B),
            weights$Chlorophyll),
        Phycocyanin = weighted.mean(
            c(GRN.B, NIR.B, NIR.R, RED.B, RED.R, YEL.B),
            weights$Phycocyanin),
        Phycoerythrin = weighted.mean(
            c(GRN.B, NIR.B, NIR.R, RED.B, RED.R, YEL.B),
            weights$Phycoerythrin)
        )

ggplot(dat %>%
    dplyr::select(date.time:pop, Chlorophyll:Phycoerythrin) %>%
    mutate(pop = log10(pop)) %>%
    pivot_longer(pop:Phycoerythrin) %>%
    mutate(name = factor(name,
        levels = c("pop", "Chlorophyll", "Phycocyanin", "Phycoerythrin"),
        labels = c("Population", "Chlorophyll", "Phycocyanin", "Phycoerythrin"))
        )) +
    scale_x_continuous() +
    aes(x = date, y = value, col = strain) +
    facet_grid(name ~ treat, scales = "free") +
    stat_summary(geom = "pointrange", fun.data = mean_se,
        position = position_dodge(width = 0.33)) +
    theme_bw()
