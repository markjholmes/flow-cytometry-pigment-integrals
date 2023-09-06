library("tidyverse")
library("fda.usc")
library("RColorBrewer")
library("viridis")

# wd of data ====
wdir <- "data/"

# load absorption data ====
abs.dat <- read.csv(paste0(wdir, "pigs-abs.csv")) %>%
  mutate(abs = ifelse(abs < 0, 0, abs)) %>%
  group_by(Pigment) %>%
  summarise(interp = as.data.frame(approx(lambda, abs, xout = 400:700))) %>%
  unnest(interp) %>%
  rename(Lambda = x, Absorption = y) %>%
  ungroup %>%
  mutate(Absorption = Absorption/max(Absorption, na.rm = TRUE)) %>%
  as.data.frame

# define a colour palette for the pigments for later
pig.cols <- setNames(brewer.pal(3, "Dark2"), unique(abs.dat$Pigment))

# load flourescence data ====
flo.dat <- read.csv(paste0(wdir, "pigs-emm.csv")) %>%
  mutate(emm = ifelse(emm < 0, 0, emm)) %>%
  group_by(pigment) %>%
  summarise(interp = as.data.frame(approx(lambda, emm, xout = 400:900))) %>%
  unnest(interp) %>%
  rename(Lambda = x, Flourescence = y, Pigment = pigment) %>%
  ungroup %>%
  mutate(Flourescence = Flourescence/max(Flourescence, na.rm = TRUE)) %>%
  as.data.frame

# merge ====
absflo.dat <- full_join(abs.dat, flo.dat) %>%
  pivot_longer(Absorption:Flourescence, names_to = "Type", values_to = "Value") %>% 
  na.omit

absflo.annot <- absflo.dat %>%
  dplyr::select(Type) %>% 
  distinct %>% 
  mutate(Label = LETTERS[1:2])

## lasers ====
lasers.dat <- read.csv(paste0(wdir, "lasers.csv"))

## cytometer channels ====
channels.dat <- read.csv(paste0(wdir, "channels.csv")) %>%
  rowwise %>% 
  mutate(ID = paste(Receptor, substring(Laser, 1, 1), collapse = "", sep = "-")) 

### channel plot colours =====
channels.cols <- setNames(viridis(nrow(channels.dat)), sort(channels.dat$ID))

#### plot ====
ggplot() +
  geom_rect(data = channels.dat,
            aes(xmin = LambdaCentre - LambdaPM,
              xmax = LambdaCentre + LambdaPM,
                fill = ID, ymin = 0, ymax = 2),
                alpha = .2, color = "grey50", lty = 2) +
  geom_vline(data = lasers.dat,
    aes(xintercept = Lambda),
    col = c("blue", "red"), lwd = 1) +
  geom_line(data = absflo.dat,
            aes(x = Lambda, y = Value, lty = Type, col = Pigment), lwd = 1) +
  geom_label(data = lasers.dat, aes(x = Lambda, y = 1.1, label = Laser)) +
  coord_cartesian(ylim = c(0, 1.2), expand = 0) +
  scale_fill_manual(values = channels.cols) +
  scale_color_manual(values = pig.cols) +
  theme_bw() +
  labs(y = "Absorption/flourescence",
       fill = "Cytometer channel",
       lty = "Type") +
  scale_y_continuous(breaks = seq(0, 1, len = 5))

ggsave("outputs/lasers_channels_pigments.pdf", width = 8, height = 5, units = "in")

# APPROXIMATE INTEGRATION USING SIMPSONS METHOD



pigments <- function() {

}

# load data ====
abs.dat <- read.csv(paste0(wdir, "pigs-abs.csv"))

## define a colour palette for the pigments for later
pig.cols <- setNames(brewer.pal(3, "Dark2"), unique(abs.dat$Pigment))

## fix data ====
abs.out <- abs.dat %>%
  mutate(abs = ifelse(abs < 0, 0, abs)) %>%
  group_by(Pigment) %>%
  summarise(interp = as.data.frame(approx(lambda, abs, xout = 400:700))) %>%
  unnest(interp) %>%
  rename(Lambda = x, Absorption = y) %>%
  as.data.frame

# load data ====
flo.dat <- read.csv(paste0(wdir, "pigs-emm.csv"))

## fix data ====
flo.out <- flo.dat %>%
  mutate(emm = ifelse(emm < 0, 0, emm)) %>%
  group_by(pigment) %>%
  summarise(interp = as.data.frame(approx(lambda, emm, xout = 400:900))) %>%
  unnest(interp) %>%
  rename(Lambda = x, Flourescence = y, Pigment = pigment) %>%
  as.data.frame

lasers.dat <- read.csv(paste0(wdir, "lasers.csv"))

## cytometer channels ====
channels.dat <- read.csv(paste0(wdir, "channels.csv")) %>% 
  rowwise %>% 
  mutate(ID = paste(Receptor, substring(Laser, 1, 1), collapse = "", sep = "-"))

### finter abs spec by laser wavelength ====
lasers.dat.full <- abs.out %>%
  right_join(lasers.dat, by = "Lambda") %>% 
  as_tibble

#### plot ====
ggplot(lasers.dat.full, aes(x = Pigment, y = Absorption, fill = Laser)) + 
  geom_col(position = position_dodge(), col = "black") +
  scale_fill_manual(values = c("Red" = "Red", "Blue" = "Navy")) + 
  theme_bw() +
  coord_cartesian(ylim = c(0,1), xlim = c(.5,3.5), expand = 0) + 
  labs(y = "Proportion of maximum excitation",
       fill = "Cytometer\nlaser") +
  theme(legend.position = c(0.87, 0.73), # bottom right
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = "white", color = "black"))

# create range of detection by the cytometer ====
channels.dat.full <- channels.dat %>% 
  group_by(Laser, Receptor, ID) %>%
  mutate(LMin = LambdaCentre - LambdaPM,
         Lmax = LambdaCentre + LambdaPM) %>%
  summarise(Lambda = list(LMin:Lmax)) %>%
  unnest(Lambda) %>%
  group_by(ID) %>%
  group_split()

# integrate total detection ====
integrated.totals <- lapply(channels.dat.full, function(i) {
  lamb <- i$Lambda
  laser <- unique(i$Laser)
  receptor <- unique(i$Receptor)
  id <- unique(i$ID)
  
  out <- flo.out %>%
    group_by(Pigment) %>%
    mutate(Flourescence = Flourescence / sum(Flourescence, na.rm = T)) %>%
    dplyr::filter(Lambda %in% lamb) %>%
    na.omit %>%
    summarise(int.tot = int.simpson2(Lambda, Flourescence)) %>%
    mutate(Laser = laser, Receptor = receptor, ID = id)
}) %>%
  bind_rows()

# plot ====
ggplot(data = integrated.totals %>% 
         dplyr::filter(Pigment %in% names(pig.cols),
                       ID %in% names(channels.cols))) +
  aes(x = Pigment, y = int.tot, fill = ID) +
  geom_col(position = position_dodge(), col = "black") +
  scale_fill_manual(values = channels.cols) +
  theme_bw() +
  coord_cartesian(ylim = c(0,1.5), xlim = c(.5,3.5), expand = 0) + 
  labs(y = "Proportion of detected flourescence",
       fill = "Cytometer\nchannel") +
  theme(legend.position = c(0.5, 0.8),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = "white", color = "black"))

## combined ====
combined.totals <- full_join(lasers.dat.full, integrated.totals) %>%
  mutate(Detection = Absorption * int.tot) 

### plot ====
ggplot(combined.totals %>% 
         dplyr::filter(Pigment %in% names(pig.cols),
                       ID %in% names(channels.cols))) +
  aes(x = Pigment, y = Detection, fill = ID) +
  geom_col(position = position_dodge(), col = "black") +
  theme_bw() +
  scale_fill_manual(values = channels.cols) +
  coord_cartesian(ylim = c(0,1), xlim = c(.5,3.5), expand = 0) + 
  labs(y = "Detection x excitation",
       fill = "Cytometer\nchannel") +
  theme(legend.position = c(0.5, 0.8),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = "white", color = "black"))

# load data ====
spec.dat <- read.csv(paste0(wdir, "species-pigments.csv")) %>%
  pivot_longer(cols = Chlorophyll:Phycoerythrin, names_to = "Pigment", values_to = "Conc") %>%
  mutate(strain = as.factor(ifelse(PigmentType == "VIII", "2434", "2524")),
         strain2 = as.factor(ifelse(PigmentType == "VIII", "2383", "2375"))) 

# plot amount of each pigment by species ====
ggplot(spec.dat) +
  aes(x = Pigment, y = Conc, fill = strain, pattern_fill = strain2) + 
  ggpattern::geom_col_pattern(
    position = position_dodge(), col = 1, show.legend = FALSE,
    pattern_size = 0, pattern_density = 0.5) + 
  theme_bw() + 
  scale_fill_manual(values = s_pal) +
  ggpattern::scale_pattern_fill_manual(values = s_pal) +
  coord_cartesian(xlim = c(.5, 3.5), ylim = c(0,10), expand = 0) +
  labs(fill = NULL, pattern_fill = NULL,
       y = expression("Approximate pigment concentration (pg"%.%"cell"^-1*")")) +
  theme(legend.position = c(0.5, 0.9),
        legend.direction = "horizontal",
        legend.background = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = "white", color = "black"))
