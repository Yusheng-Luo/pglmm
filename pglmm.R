library(phyr)
library(openxlsx)
library(picante)
library(tidyr)
library(dplyr)
library(MuMIn)
library(lme4)
library(glmm.hp)
library(rdacca.hp)
library(hier.part)
library(lmerTest)
library(sjstats)
library(learnasreml)
comm <- read.csv("community.csv",2)
row.names(comm) <- comm$site
envi <- read.xlsx("envi.xlsx", 3)
dat <- tidyr::gather(comm, key = "sp", value = "freq", -site) %>% left_join(envi, by = "site")
dat$pa = as.numeric(dat$freq > 0)

effects <- ranef(test2)
print(effects)



phy <- read.tree("tree1.newick")
LM_1 <- lmer(freq ~ N + P + SOM + K + PH + (1|zone), data = dat)
LM_2 <- lmer(freq ~ N + (N|zone) + (1|zone), data = dat)
summary(LM_2)


summary(LM_1)
glmm.hp(LM_1)
performance::r2(LM_1)
anova(LM_1, type = "I")
  test1 = phyr::pglmm(freq ~ 1 + N + P + SOM + K + (1|sp__) + (1|site) + (1|sp__@site),
                    data = dat, family = "gaussian", REML = FALSE,
                    cov_ranef = list(sp = phy))

test2 = phyr::pglmm(pa ~ 1 + Ele + (1|sp__) + (1|site) + (1|sp__@site), 
                    data = dat, family = "binomial", REML = FALSE,
                    cov_ranef = list(sp = phy))
summary(test2)

summary(test1)
anova(test1, test2)

ggplot(data = dat, aes(Ele, freq)) + geom_point() + geom_smooth(method = "lm")
# 
# test3 = phyr::pglmm(freq ~ 1 + SOM + (1|sp__) + (1|site) + (1|sp__@site), 
#                     data = dat, family = "gaussian", REML = FALSE,
#                     cov_ranef = list(sp = phy))

library(rr2)

install.packages("rr2")
rr2::R2_pred(test1)
p_value(test1)
r2(test1)

r.squaredLR(test1)


library(performance)
r2(test1)







