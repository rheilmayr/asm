# Import packages
library(bife)
library(tidyverse)

# Load primary dataset
df <- read.csv("data\\long.csv") # note: working directory needs to be root directory of project
df <- df %>%
  filter((dist_amb>-300) & (dist_aml<300) & (year>2001) & (legal_amazon==1)) %>%
  mutate(soy_suit = (suit>0 & GAEZsuit>40),
         post_2005 = year>2005,
         state_year = paste0(state, "_", year))

reg_df <- df %>%
  select(mb2_vdefor, soy_suit, post_2005, biome, state_year, state, year, temp, trmm, roaddist, urbandist, pa, set) %>%
  drop_na()

## Main regression results - building up to full model (T3)
bc_APE <- function(bife_model){
  apes <- get_APEs(bias_corr(bife_model))
  return(apes)
}
ddd_reg <- bife(mb2_vdefor ~ (soy_suit : biome : post_2005) +
                  (soy_suit : biome) + (biome : post_2005) + (soy_suit : post_2005) + 
                  soy_suit + biome + temp + trmm + roaddist + urbandist + pa + set | state_year, 
                data = reg_df)
ddd <- ddd_reg %>% bc_APE()


## Modify stata tex file with results
coef <- format(summary(ddd)[12,1], digits = 3, nsmall = 3)
if (summary(ddd)[12,4]<0.01) {
  stars = "***"
} else if (summary(ddd)[12,4]<0.05) {
  stars = "**"
} else if (summary(ddd)[12,4]<0.1) {
  stars = "*"
}
coef <- paste0(coef, stars)
se <- format(summary(ddd)[12,2], digits = 3, nsmall = 3)

tex <- "data\\tables\\ts4_robustness_wcox.tex"
new_tex <- "data\\tables\\ts4_robustness_amend.tex"
table <- readChar(tex, file.info(tex)$size)
table <- table %>% str_replace("0.583", coef)
table <- table %>% str_replace("0.390", se)
table <- table %>% str_replace("467", "279360")
table <- table %>% str_replace("27", "563")
table %>% write(file = new_tex)
