
library(spaa)
library(reshape2)
library(openxlsx)
library(phyr)
library(V.PhyloMaker2)
library(V.PhyloMaker)
library(plantlist)
library(LPSC)
library(patchwork)
library(picante)
library(phytools)

dat <- read.xlsx("胶州.xlsx", 5)
dat_t <- acast(dat, formula = site ~ species,
               value.var = "abun",
               fill = 0)
dat_t <- as.data.frame(dat_t)

tree <- phylo.maker(dat, tree=GBOTB.extended.TPL, nodes=nodes.info.1.TPL, scenarios="S3")
write.tree(tree$scenario.3, "GBOTB.extended.LCVP.newick")

tree_TPL <- phylo.maker(dat_tree, tree=GBOTB.extended.TPL, nodes=nodes.info.1.TPL, scenarios="S3")
write.tree(tree_TPL$scenario.3, "GBOTB.extended.TPL.tre")


write.tree(test$scenario.3, "test.newick")


test_tree <- read.tree("test.newick")

plot(tree_1)


tree_1 <- read.tree("GBOTB.extended.LCVP.newick")
tree_2 <- read.tree("GBOTB.extended.TPL.tre")

comm.pd <- pd(dat_1, tree_2)

phyr::psd(dat_t, tree_1)
phyr::ps
dat_t <- as.data.frame(dat_t)
comm <- t(dat_t)
comm.pd <- pd(comm, tree_2)
head(comm.pd)

comm.bc.dist <- vegdist(dat_t, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
# plot cluster diagram
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")


comm.bc.mds <- metaMDS(dat_t, dist = "bray")
stressplot(comm.bc.mds)
ordiplot(comm.bc.mds, display = "sites", type = "text")

phy.dist <- cophenetic(tree_1)
# calculate ses.mpd
comm.sesmpd <- ses.mpd(dat_t, phy.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
head(comm.sesmpd)

phyr::psv(dat_t, tree_1)


picante::pd(dat_t, tree_1)
ses.pd(dat_t, tree_1)
ses.pd(samp, tree, null.model = c("taxa.labels", "richness", "frequency",
                                  "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
       runs = 999, iterations = 1000, include.root=TRUE)
install.packages("codyn")
library(codyn)
data(package="codyn")
data("knz_001d")
data("collins08")
data("pplots")

library(qgraph)
data(package = "qgraph")


# 李颖数据探索 ------------------------------------------------------------------
library(openxlsx)
library(qgraph)
dat_hz <- read.xlsx("Data-Leaf trait network.xlsx", 1)
big5Graph <- qgraph(cor(dat),alpha = 0.01,borders=FALSE)

library("psych")


qgraph(big5Graph,layout="spring")
install.packages("igraph")
library(igraph)
igraph::

  
library(tidyverse)
library(igraph)
library(psych)



### 2.1 计算相关性系数，建议使用spearman系数。
cor <- corr.test(dat_hz, use = "pairwise",method="pearson",adjust="holm", alpha=.05,ci=FALSE) # ci=FALSE,不进行置信区间计算，数据量较大时，可以加快计算速度。
cor.r <- data.frame(cor$r) # 提取R值
cor.p <- data.frame(cor$p) # 提取p值
colnames(cor.r) = rownames(cor.r)
colnames(cor.p) = rownames(cor.p) # 变量名称中存在特殊字符，为了防止矩阵行名与列名不一致，必须运行此代码。
write.csv(cor.r,"cor.r.csv",quote = FALSE,col.names = NA,row.names = TRUE) # 保存结果到本地
write.csv(cor.p,"cor.p.csv",quote = FALSE,col.names = NA,row.names = TRUE)

# knitr::kable(
#   head(cor.r),
#   caption = "cor.r"
# )
# head(cor.p)

cor.r[abs(cor.r) < 0.2 | cor.p > 0.05] = 0
cor.r = as.matrix(cor.r)
g = graph_from_adjacency_matrix(cor.r,mode = "undirected",weighted = TRUE,diag = FALSE)
cor.r$node1 = rownames(cor.r) 
cor.p$node1 = rownames(cor.p)

r = cor.r %>% 
  gather(key = "node2", value = "r", -node1) %>%
  data.frame()

p = cor.p %>% 
  gather(key = "node2", value = "p", -node1) %>%
  data.frame()
head(r)
head(p)

#### 将r和p值合并为一个数据表
cor.data <- merge(r,p,by=c("node1","node2"))
cor.data

or.data <- cor.data %>%
  filter(abs(r) >= 0.2, p <= 0.05, node1 != node2) %>%
  mutate(
    linetype = ifelse(r > 0,"positive","negative"), # 设置链接线属性，可用于设置线型和颜色。
    linesize = abs(r) # 设置链接线宽度。
  ) # 此输出仍有重复链接，后面需进一步去除。
head(cor.data)

c(as.character(cor.data$node1),as.character(cor.data$node2)) %>%
  as_tibble() %>%
  group_by(value) %>%
  summarize(n=n()) -> vertices
colnames(vertices) <- c("node", "n")
head(vertices)
g <- graph_from_data_frame(cor.data, vertices = vertices, directed = FALSE )
g
is.simple(g)
E(g)$weight <- 1
g <- igraph::simplify(g,
                      remove.multiple = TRUE,
                      remove.loops = TRUE,
                      edge.attr.comb = "first")

#g <- delete.vertices(g,which(degree(g) == 0)) # 删除孤立点
is.simple(g)
E(g)$weight <- 1
is.weighted(g)
vcount(g) # 节点数目：35
ecount(g) # 链接数:585


### 3.4 计算节点链接数
V(g)$degree <- degree(g)
#vertex.attributes(g)
#edge.attributes(g) 

### 4.1 准备网络图布局数据
#?layout_in_circle # 帮助信息中，有其它布局函数。
layout1 <- layout_in_circle(g) # 径向布局适合节点较少的数据。
layout2 <- layout_with_fr(g) # fr布局。
layout3 <- layout_on_grid(g) # grid布局。
layout4 <- layout_with_mds(g)
layout5 <- layout_with_fr(g)
layout6 <- layout_with_gem(g)
head(layout1)


color <- c(rgb(65,179,194,maxColorValue = 255),
           rgb(255,255,0,maxColorValue = 255),
           rgb(201,216,197,maxColorValue = 255))

# names(color) <- unique(V(g)$type) # 将颜色以节点分类属性命名
# V(g)$point.col <- color[match(V(g)$type,names(color))] # 设置节点颜色。
# #names(color2) <- unique(V(g)$type) # 如果想要节点颜色与背景颜色不一致，则可以为节点单独设置一个颜色集。
# #V(g)$point.col <- color2[match(V(g)$type,names(color2))]
# 
# #### 边颜色按照相关性正负设置
# #E(g)$color <- ifelse(E(g)$linetype == "positive",rgb(255,215,0,maxColorValue = 255),"gray50")
# E(g)$color <- ifelse(E(g)$linetype == "positive","red",rgb(0,147,0,maxColorValue = 255))


plot.igraph(g, layout=layout6)



install.packages("herblabel", repos="http://R-Forge.R-project.org")
install.packages("Rcpp")
library(herblabel)

install.packages("forestmangr")
library(forestmangr)
data("exfm16")
head(exfm16)
remotes::install_github(repo = "ForestManagementGeodatabase/fmgr")


install.packages("MtreeRing")
library(MtreeRing)
library(forestSAS)
library(spatstat)


library(ape) # for phylogeny
library(plyr)
library(dplyr, quietly = TRUE)
library(pez) # for commuityPGLMM function
library(parallel) #
install.packages("pez")

load("d_li_data.RData")
install.packages("remotes")
library(remotes)
remotes::install_github("helixcn/fgeo.habitat")
library(fgeo)
??fgeo
remotes::install_github("helixcn/fgeo.habitat")
remotes::install_github("forestgeo")
fgeo.analyze::fgeo_habitat()

library(plantlist)
CTPL("淡竹")

# 崂山样地 --------------------------------------------------------------------
library(spatstat)
library(openxlsx)
library(ggplot2)
library(dplyr)
library(LPSC)
library(plantlist)
library(reshape2)
library(spaa)
library(vegan)
library(rdacca.hp)
library(forestSAS)
library(ggpubr)
library (MASS)
dat_ls <- read.xlsx("崂山木本.xlsx",5)
dat_spe <- read.xlsx("species.xlsx", 2)

setwd("F:/博士研究计划2024/崂山英文想法/冗余分析/")
dat_plot <- read.xlsx("spe.xlsx", 3)
raw_abun <- acast(dat_plot, formula = subplot~spe,
                  value.var = "abun",
                  fill = 0)
# 冗余分析 --------------------------------------------------------------------

windowsFonts(A=windowsFont("Times New Roman"), B=windowsFont("Arial"))
dat_plot <- read.xlsx("崂山木本.xlsx", 15)
raw_abun <- acast(dat_plot, formula = subplot~spe,
                  value.var = "abun",
                  fill = 0)

raw_abun <- as.data.frame(raw_abun)
write.csv(raw_abun, "raw_abun.csv")
raw_soil <- read.xlsx("崂山木本.xlsx", 7)


row.names(raw_soil) <- raw_soil[,1]

raw_soil <-  raw_soil[,-1]

dca <- decorana(veg = t(raw_abun))
dca
raw_abun <- decostand(raw_abun, method = "hellinger")
res <- rda(raw_abun, raw_soil, scale = FALSE)
plot(res)

res_1 <- rda(raw_abun, raw_soil, scale = TRUE)
plot(res_1)
res_1 <- rda(raw_abun~., raw_soil, scale = TRUE)
anova(res_1, permutations = how(nperm = 999))
summary(res_1)
RsquareAdj(res_1)


anova(res, permutations = how(nperm = 999))
summary(res)
res <- rda(raw_abun~., raw_soil, scale = FALSE)
anova(res, permutations = how(nperm = 999))
summary(res)
raw_abun <- as.matrix(raw_abun)
raw_soil <- as.matrix(raw_soil)

RDA.Perm = permutest(res, permutations = 999)
RDA.Perm
envfit <- envfit(res, raw_soil, permutations  = 999)
envfit
RsquareAdj(res)


zh <- summary(res_1)$concont$importance
zh <- round(zh, 4)

pdat <- res_1$CCA
samples <-data.frame(sample = row.names(pdat$u),RDA1 = pdat$u[,1],RDA2 = pdat$u[,2])
species<-data.frame(spece = row.names(pdat$v),RDA1 = pdat$v[,1],RDA2 = pdat$v[,2])
envi<-data.frame(en = row.names(pdat$biplot),RDA1 = pdat$biplot[,1],RDA2 = pdat$biplot[,2])
species <- read.xlsx("RDA_sp.xlsx")
p <- ggplot() +
  geom_point(data = species, aes(x = RDA1, y = RDA2),size = 2,  shape = 16) +
  # stat_ellipse(data = species,aes(x = RDA1, y = RDA2), level = 0.95,linetype = "dashed", show.legend = F, size = 1.5) +
  # annotate('text', label = 'Evergreen', x = 0.25, y = 0, size = 5, colour = 'black') +
  # annotate('text', label = 'Deciduous', x = -0.22, y = -0.25, size = 3, colour = 'black') +
  geom_hline(aes(yintercept = 0), colour="gray88", linetype="dashed") +
  geom_vline(aes(xintercept = 0), colour="gray88", linetype="dashed")  +
  #stat_ellipse(level = 0.95, show.legend = F, size = 1.5, color = species$group.ad) +
  #geom_point(data = species, aes(x=RDA1, y=RDA2, size = 5,color = group$ad)) +
  #scale_color_gradient(low = "#f47720",high = "#88c4e8")+
  #geom_point(aes(clour =factor(1)), size = 5)  +
  #geom_point(color = "green", size = 1.5)
  #scale_colour_manual(values =c("green","yellow"), ) +
  #geom_segment(data = species, aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), color = 'blue') +
  geom_text(data = species,aes(x = RDA1*1.1, y = RDA2*1.1, label = spece), size = 5, color = '#f47720') +
  geom_segment(data = envi,aes(x=0, xend= RDA1, y=0, yend= RDA2), size = 1,arrow = arrow(length = unit(0.3, "cm")), colour = '#88c4e8') +
  geom_text(data = envi,aes(x = RDA1*1.1, y = RDA2*1.1, label = en), size = 3, colour = 'black', check_overlap = FALSE) +
  #geom_point(data = samples, aes(x=RDA1, y=RDA2),size = 3, shape = 25, color = "black", fill = "black")  +
  theme_bw() + theme(axis.ticks.length = unit(-5, "pt"), panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"), panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12,face = 'bold'), axis.title.x=element_text(vjust=1, size=13), axis.title.y=element_text(size=13)) +
  theme(legend.position = 'none') +
  xlab(paste('RDA1 (', zh[2,1]*100, '%)', sep ='')) + ylab(paste('RDA2 (', zh[2,2]*100, '%)', sep =''))
p

write.xlsx(species, "RDA_sp.xlsx")


dist_matrix <- vegdist(raw_abun, method = "bray")
rda_hp <- rdacca.hp(dist_matrix, raw_soil, method = "dbRDA", type = "R2")
p_2 <- plot(rda_hp)
p_2
p_6 <- ggarrange(p, p_2, 
                 labels = c("A", "B"),
                 ncol = 2, nrow = 1)


p_6

ggsave("rdahp_LS_1.png", width = 14, height = 6.5, dpi = 500)

save(dat_ls, dat_plot, dca, envfit, envi,pdat, res, samples, raw_abun, raw_soil,species,zh, file = "ls_RDA.RData")


ggsave("rdahp_LS_1.svg", width = 14, height = 6.5, dpi = 300)


# 菌根分组



p <- ggplot() +
  geom_point(data = species, aes(x = RDA1, y = RDA2, color = type), size = 5) +
  # stat_ellipse(data = species,aes(x = RDA1, y = RDA2, color = type), level = 0.95,linetype = "dashed", show.legend = F, size = 1.5) +
  # annotate('text', label = 'AM', x = 0.25, y = 0, size = 5, colour = 'black') +
  # annotate('text', label = 'EM', x = -0.22, y = -0.25, size = 3, colour = 'black') +
  # annotate('text', label = 'AEM', x = -0.22, y = -0.25, size = 3, colour = 'black') +
  geom_hline(aes(yintercept = 0), colour="gray88", linetype="dashed") +
  geom_vline(aes(xintercept = 0), colour="gray88", linetype="dashed")  +
  #stat_ellipse(level = 0.95, show.legend = F, size = 1.5, color = species$group.ad) +
  #geom_point(data = species, aes(x=RDA1, y=RDA2, size = 5,color = group$ad)) +
  # scale_color_gradient(low = "#f47720",high = "#88c4e8")+
  #geom_point(aes(clour =factor(1)), size = 5)  +
  #geom_point(color = "green", size = 1.5)
  #scale_colour_manual(values =c("green","yellow"), ) +
  #geom_segment(data = species, aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), color = 'blue') +
  geom_text(data = species,aes(x = RDA1*1.1, y = RDA2*1.1, label = spece, color = type), size = 3) +
  geom_segment(data = envi,aes(x=0, xend= RDA1, y=0, yend= RDA2), size = 1,arrow = arrow(length = unit(0.3, "cm")), colour = 'black') +
  geom_text(data = envi,aes(x = RDA1*1.1, y = RDA2*1.1, label = en), size = 3, colour = 'black', check_overlap = FALSE) +
  #geom_point(data = samples, aes(x=RDA1, y=RDA2),size = 3, shape = 25, color = "black", fill = "black")  +
  theme_bw() + theme(axis.ticks.length = unit(-5, "pt"), panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"), panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12,face = 'bold'), axis.title.x=element_text(vjust=1, size=13), axis.title.y=element_text(size=13)) +
  theme(legend.position = 'bottom') +
  xlab(paste('RDA1 (', zh[2,1]*100, '%)', sep ='')) + ylab(paste('RDA2 (', zh[2,2]*100, '%)', sep =''))
p


ggsave("ls_jg.png", width = 7, height = 6.5, dpi = 300)









# 种间联结计算 ------------------------------------------------------------------
setwd("F:/博士研究计划2024/大样地数据/种间联结计算/")
windowsFonts(A=windowsFont("Times New Roman"), B=windowsFont("Arial"))
library(spaa)
spaa::sp.assoc(dat_t)
# sub.sp.matrix(raw_abun, freq = 0.5, common = NULL)
dat_abunance <- read.xlsx("dat_i.xlsx", 3, rowNames = TRUE)
niche_width <- spaa::niche.width(dat_abunance, method = "shannon")
write.xlsx(niche_width, "shannon.xlsx")
niche_width_1 <- spaa::niche.width(dat_abunance, method = "levins")
write.xlsx(niche_width_1, "levins.xlsx")
spaa::niche.width()
raw_abun <- as.matrix(raw_abun)

pianka <- spaa::niche.overlap(raw_abun, method = "pianka")
pianka <- as.data.frame(pianka)
spaa::sp.assoc(raw_abun)
pianka <- spaa::niche.overlap.pair(raw_abun, method = "pianka")

pianka <-niche.overlap(raw_abun, method = c("pianka"))
pianka <- as.matrix(pianka)
pianka[upper.tri(pianka)] = NA
pianka <- reshape2::melt(pianka, na.rm = TRUE)

cor <- cor(raw_abun, method = "pearson")
cor <- as.matrix(cor)
cor[upper.tri(cor)] = NA
cor <- reshape2::melt(cor, na.rm = TRUE)

write.xlsx(pianka, "F:/博士研究计划2024/崂山英文想法/冗余分析/pianka.xlsx")
write.xlsx(CC, "F:/博士研究计划2024/崂山英文想法/冗余分析/CC.xlsx")
write.xlsx(cor, "F:/博士研究计划2024/崂山英文想法/冗余分析/cor.xlsx")
# 重要值与生态位
library(ggplot2)
library(ggsci) 
library(ggExtra) 
library(ape)
library(caper)
library(phangorn)
library(phylolm)
library(picante)
library(ggpmisc)        
library(palmerpenguins)
library(openxlsx)
dat_niche <- read.xlsx("F:/博士研究计划2024/大样地数据/系统发育指标/levins.xlsx", 2)
dat_niche_1 <- read.xlsx("F:/博士研究计划2024/大样地数据/系统发育指标/levins.xlsx", 3,rowNames = TRUE)
windowsFonts(A=windowsFont("Times New Roman"), B=windowsFont("Arial"))
library(nlme)
library(geiger)
library(cowplot)


apply(dat_niche_1, 2, Kcalc, phy)

multiPhylosignal(dat_niche_1, multi2di(phy))

combined <- match.phylo.data(phy, dat_niche_1)


root.gls <- gls(levins~iv, data = dat_niche_1)
anova(root.gls)
summary(root.gls)
root.pgls <- gls(levins~iv,
                 correlation = corBrownian(value = 1, phy),
                 data = dat_niche_1)
anova(root.pgls)
summary(root.pgls)



library(caper)




# 绘图
plot(levins~iv, data = dat_niche_1, 
     xlab = "SRL (specific root length)",
     ylab = "Root tissue density")
abline(coef(root.gls), lwd = 2, col = "black")
abline(coef(root.pgls), lwd = 2, col = "red")
legend("bottomleft", 
       legend = c("GLS fit", "Phylogenetic GLS fit"), 
       lwd = 2, col = c("black", "red"))





p1 = ggplot(dat_niche_1,aes(x=iv,y=levins))+
  geom_point(size=2.6,shape = 21,color="gray2",alpha=1)+
  # geom_smooth(method = 'lm', se = T, level=0.95,size=1.0) +
  scale_fill_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_color_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_x_continuous(limit = c(0,0.16))+
  scale_y_continuous(limit = c(0,15))+
  # abline(coef(root.gls), lwd = 2, col = "black") +
  # abline(coef(root.pgls), lwd = 2, col = "red") +
  theme_bw()+labs(x="Important value",y="Niche width(Levins)")+
  # annotate('text', label = 'R2 = 0.39, p = 0.001', x =-0.20, y = 2, size =3.9,color="#DD5F60")+
  # annotate('text', label = 'R2 = 0.48, p = 0.001', x =-0.20, y = 1.9, size =3.9,color="#9BCD9B")+
  theme(axis.text=element_text(colour='black',size=9))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +theme(legend.position=c(0.8, 0.9)) +
  geom_abline(slope= 116.12789 , intercept= 1.57514 , lwd = 1, col = "#9BCD9B") +
  geom_abline(slope= 117.33058 , intercept= 1.46854, lwd = 1, col = "#DD5F60") 

p1 




root.gls <- gls(shannon~iv, data = dat_niche_1)
anova(root.gls)
summary(root.gls)
root.pgls <- gls(shannon~iv,
                 correlation = corBrownian(value = 1, phy),
                 data = dat_niche_1)
anova(root.pgls)
summary(root.pgls)


p2 = ggplot(dat_niche_1,aes(x=iv,y=shannon))+
  geom_point(size=2.6,shape = 21,color="gray2",alpha=1)+
  # geom_smooth(method = 'lm', se = T, level=0.95,size=1.0) +
  # scale_fill_manual(values=c("#DD5F60","#9BCD9B"))+
  # scale_color_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_x_continuous(limit = c(0,0.16))+
  scale_y_continuous(limit = c(0,3))+
  # abline(coef(root.gls), lwd = 2, col = "black") +
  # abline(coef(root.pgls), lwd = 2, col = "red") +
  theme_bw()+labs(x="Important value",y="Niche width(Shannon)")+
  # annotate('text', label = 'R2 = 0.39, p = 0.001', x =-0.20, y = 2, size =3.9,color="#DD5F60")+
  # annotate('text', label = 'R2 = 0.48, p = 0.001', x =-0.20, y = 1.9, size =3.9,color="#9BCD9B")+
  # theme(axis.text=element_text(colour='black',size=9))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())  +
  geom_abline(slope= 27.924189 , intercept= 0.539952 , lwd = 1, col = "#9BCD9B") +
  geom_abline(slope= 31.84385 , intercept= 0.32351, lwd = 1, col = "#DD5F60") 

p2 


cowplot::plot_grid(p1,p2,nrow=1, rel_widths = c(1, 1), labels=LETTERS[1:2])

ggsave("hgpgls_1.png", width = 12, height = 7,dpi = 500)
ggsave("F:/博士研究计划2024/崂山英文想法/图_1/hgpgls_1.pdf", width = 12, height = 7,dpi = 500)

ggsave("F:/博士研究计划2024/崂山英文想法/图_1/hgpgls_1.png", width = 12, height = 7,dpi = 500)
getwd()



















# 生态位/系统发育信号检验
dat_n <- read.xlsx("F:/博士研究计划2024/大样地数据/phyr/dat_phyr.xlsx", 5)
# phy_tree研究系统发育树
phylod.niche <- phylo.d(data=dat_n, phy=phy_tree, names = tip)
phylod.bird
plot(phylod.bird)

# 物种间相关性计算
setwd("F:/博士研究计划2024/大样地数据/种间联结计算/")
library(spaa)
library(openxlsx)
library(corrplot)
library(ggcorrplot)
library(Cairo)
dat_cor <- read.xlsx("dat_i.xlsx",5)
dat_corr <- read.xlsx("dat_i.xlsx",5,rowNames = TRUE)
spaa::sp.assoc(dat_cor)
spaa::plotnetwork(dat_corr)
dat_ <- spaa::sp.pair(dat_corr)
dat_$Pearson

dat_pic <- cor(dat_corr, method = "pearson") 
p_1 <- corrplot(dat_pic, method = "circle")
p_2 <- corrplot(dat_pic, method = "square")

# pearson相关性分析
cor_Matrix <- round(cor(dat_corr, method = "pearson"),2)
# 查看相关矩阵输出结果
cor_Matrix
# 提取p值矩阵
p.mat_Matrix <-round(cor_pmat(dat_corr,method = "pearson"),2)
col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"),alpha = TRUE)
# 先绘制相关系数数字样式的热图


Cairo::CairoPNG( 
  filename = "name.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


corrplot(cor_Matrix, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set(100))
# 绘制下三角方块样式的热图并使用add参数把下三角热图添加到系数数字热图中(组合)
corrplot(cor_Matrix, method = "square", type="lower", order = "hclust", add=TRUE, diag=FALSE, 
         col = col_set(100), tl.cex = 0.7, tl.col = "black", tl.pos = "n", cl.pos="n", 
         cl.length=5, addgrid.col = "grey70", outline = "grey60", p.mat = p.mat_Matrix, 
         insig = "label_sig", sig.level = c(.01, .05), pch.cex = 1, pch.col = "black")

dev.off()
dat_pearson <- cor(dat_corr, method = "pearson")

dat_pearson[upper.tri(dat_pearson)]= NA
dat_pearson_matrix <- reshape2::melt(dat_pearson, na.rm = T)
write.csv(dat_pearson_matrix, "dat_pearson_matrix.csv")


# spearman相关性分析
cor_Matrix <- round(cor(dat_corr, method = "spearman"),2)
# 查看相关矩阵输出结果
cor_Matrix
# 提取p值矩阵
p.mat_Matrix <-round(cor_pmat(dat_corr,method = "spearman"),2)
col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"),alpha = TRUE)
# 先绘制相关系数数字样式的热图


Cairo::CairoPNG( 
  filename = "spearman.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


corrplot(cor_Matrix, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set(100))
# 绘制下三角方块样式的热图并使用add参数把下三角热图添加到系数数字热图中(组合)
corrplot(cor_Matrix, method = "square", type="lower", order = "hclust", add=TRUE, diag=FALSE, 
         col = col_set(100), tl.cex = 0.7, tl.col = "black", tl.pos = "n", cl.pos="n", 
         cl.length=5, addgrid.col = "grey70", outline = "grey60", p.mat = p.mat_Matrix, 
         insig = "label_sig", sig.level = c(.01, .05), pch.cex = 1, pch.col = "black")

dev.off()


dat_spearman <- cor(dat_corr, method = "spearman")

dat_spearman[upper.tri(dat_spearman)]= NA
dat_spearman_matrix <- reshape2::melt(dat_spearman, na.rm = T)
write.csv(dat_spearman_matrix, "dat_spearman_matrix.csv")
# AC联结系数
Cairo::CairoPNG( 
  filename = "AC.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率
# ac_toal <- spaa::sp.pair(dat_corr)
# AC <- ac_toal$AC
# AC <- as.matrix(AC)

corrplot(AC, method = "square",type = "lower")
# p.mat_Matrix <-round(AC,2)

p.mat_Matrix <- read.csv("AC.csv", row.names = 1)
p.mat_Matrix <- as.matrix(p.mat_Matrix)
col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"),alpha = TRUE)
corrplot(p.mat_Matrix, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set(100))
dev.off()


# PC联结系数
ac_toal$Jaccard
PC <- ac_toal$Jaccard
PC <- as.matrix(PC)
PC <- round(PC,2)
write.csv(PC, "PC.csv")
p.mat_Matrix <- read.csv("PC.csv", row.names = 1)
p.mat_Matrix <- as.matrix(p.mat_Matrix)

Cairo::CairoPNG( 
  filename = "PC.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"),alpha = TRUE)
corrplot(p.mat_Matrix, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set(100))
dev.off()

# 卡方联结系数
ac_toal$chisq
KF <- ac_toal$chisq
KF <- as.matrix(KF)
# KF <- round(KF,2)
write.csv(KF, "KF.csv")
p.mat_Matrix <- read.csv("KF.csv", row.names = 1)
p.mat_Matrix <- as.matrix(p.mat_Matrix)

Cairo::CairoPNG( 
  filename = "KF.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


col_set <- colorRampPalette(c("#77C034","white" ,"lightskyblue"),alpha = TRUE)
corrplot(p.mat_Matrix, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set(100))
dev.off()

# 生态位重叠系数
# install.packages("devtools")
devtools::install_github("Hy4m/linkET", force = TRUE)
packageVersion("linkET")
library(linkET)
niche.overlap()
niche_levins <- spaa::niche.overlap(dat_corr, method = "levins")
View(niche_levins)
niche_levins <- as.matrix(niche_levins)

niche_levins <- round(niche_levins,2)




Cairo::CairoPNG( 
  filename = "niche_levins.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


col_set <- colorRampPalette(c("#77C034" ,"lightskyblue"),alpha = TRUE)
corrplot(niche_levins, method = "square", type="lower",order = "hclust", tl.cex = 0.7,
         tl.col = "black",col = col_set())
dev.off()
corrplot(niche_levins, method = "square")
??qcorrplot
qcorrplot(niche_levins)


Cairo::CairoPNG( 
  filename = "niche_levins.png", # 文件名称
  width = 7,           # 宽
  height = 7,          # 高
  units = "in",        # 单位
  dpi = 300)           # 分辨率


set_corrplot_style()
qcorrplot(niche_levins, type = "lower", diag = FALSE, size = 0.1) + geom_square() + scale_fill_gradientn(colours = c("#77C034" ,"lightskyblue")) +
  geom_couple(aes(niche_levins),label.colour = "black",
                                        label.fontface=2,
              label.size =4,
              drop = T)
dev.off()







qcorrplot(niche_levins, type="lower", drop = TRUE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#77C034" ,"white","lightskyblue"),
                       name = NULL) +
  theme(legend.key.height = unit(0.5, 'cm'),legend.key.width = unit(3.1, 'cm'),legend.position = "bottom")

ggsave("niche_levins.png", width = 7, height = 7, dpi = 300)
niche_levins[upper.tri(niche_levins)]= NA
niche_levins_matrix <- reshape2::melt(niche_levins, na.rm = T)
write.csv(niche_levins_matrix, "niche_levins_matrix.csv")

# niche_schoener

niche_schoener <- spaa::niche.overlap(dat_corr, method = "schoener")

niche_schoener <- as.matrix(niche_schoener)
niche_schoener[upper.tri(niche_schoener)]= NA
niche_schoener_matrix <- reshape2::melt(niche_schoener, na.rm = T)
write.csv(niche_schoener_matrix, "niche_schoener_matrix.csv")

# niche_pianka

niche_pianka <- spaa::niche.overlap(dat_corr, method = "pianka")
niche_pianka <- as.matrix(niche_pianka)
niche_pianka[upper.tri(niche_pianka)]= NA
niche_pianka_matrix <- reshape2::melt(niche_pianka, na.rm = T)
write.csv(niche_pianka_matrix, "niche_pianka_matrix.csv")

qcorrplot(niche_pianka, type="lower", drop = TRUE) +
  geom_square() +
  scale_fill_gradientn(colours = c("#77C034" ,"white","lightskyblue"),
                       name = NULL) +
  theme(legend.key.height = unit(0.5, 'cm'),legend.key.width = unit(3.1, 'cm'),legend.position = "bottom")

ggsave("niche_pianka.png", width = 7, height = 7, dpi = 300)







# 谱系距离计算
library(picante)
lsyd_tree <- read.tree("tree.newick")


phylodist <- as.matrix(as.dist(cophenetic.phylo(lsyd_tree)))
phylodist <- scale(phylodist)
phylodist[upper.tri(phylodist)]= NA
phylodist_matrix <- reshape2::melt(phylodist, na.rm = T)
write.xlsx(phylodist_matrix, "phylodist_matrixx.xlsx")

# 回归
library(ggplot2)
library(ggpubr)
library(car)
dat_niche_hg <- read.xlsx("phylodist_matrix.xlsx", 3)
data_n_hg_p <- read.xlsx("phylodist_matrix.xlsx", 5)
data_n_hg_s <- read.xlsx("phylodist_matrix.xlsx", 6)
fit_lm<-lm(pearson~niche_levins, data = dat_niche_hg)
summary(fit_lm)
fit_lm<-lm(pearson~niche_pianka, data = dat_niche_hg)
summary(fit_lm)
p1 = ggplot(data_n_hg_p,aes(x=pearson,y=niche))+
  geom_point(size=2.6,aes(fill=group),shape = 21,color="gray2",alpha=1)+
  geom_smooth(aes(color=group), method = 'lm', se = T, level=0.95,size=1.0) +
  scale_fill_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_color_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_x_continuous(limit = c(-0.5,1))+
  scale_y_continuous(limit = c(0,2))+
  theme_bw()+labs(x="Pearson correlation coefficient",y="Niche overlap")+
  annotate('text', label = 'R2 = 0.38, p < 0.001', x =-0.20, y = 2, size =3.9,color="#DD5F60")+
  annotate('text', label = 'R2 = 0.48, p < 0.001', x =-0.20, y = 1.9, size =3.9,color="#9BCD9B")+
  theme(axis.text=element_text(colour='black',size=9))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +theme(legend.position=c(0.8, 0.9))
p1 

p2 = ggplot(data_n_hg_s,aes(x=spearman,y=niche))+
  geom_point(size=2.6,aes(fill=group),shape = 21,color="gray2",alpha=1)+
  geom_smooth(aes(color=group), method = 'lm', se = T, level=0.95,size=1.0) +
  scale_fill_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_color_manual(values=c("#DD5F60","#9BCD9B"))+
  scale_x_continuous(limit = c(-0.5,1))+
  scale_y_continuous(limit = c(0,2))+
  theme_bw()+labs(x="Spearman correlation coefficient",y="Niche overlap")+
  annotate('text', label = 'R2 = 0.28, p = 0.001', x =-0.20, y = 2, size =3.9,color="#DD5F60")+
  annotate('text', label = 'R2 = 0.37, p = 0.001', x =-0.20, y = 1.9, size =3.9,color="#9BCD9B")+
  theme(axis.text=element_text(colour='black',size=9))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank()) +theme(legend.position=c(0.8, 0.9))
p2 

library(cowplot)
cowplot::plot_grid(p1,p2,nrow=1, rel_widths = c(1, 1), labels=LETTERS[1:2])

ggsave("hg_1.png", width = 12, height = 7,dpi = 500)
ggsave("hg_1.pdf", width = 12, height = 7,dpi = 500)

fit_lm<-lm(spearman~niche_levins, data = dat_niche_hg)
summary(fit_lm)
fit_lm<-lm(spearman~niche_pianka, data = dat_niche_hg)
summary(fit_lm)



AC <- read.csv("F:/博士研究计划2024/大样地数据/种间联结计算/AC.csv")
colnames(AC) <- AC[,1]









# 构建系统发育树 -----------------------------------------------------------------
library(spaa)
library(reshape2)
library(openxlsx)
library(phyr)
library(V.PhyloMaker2)
library(V.PhyloMaker)
library(plantlist)
library(LPSC)
library(patchwork)
library(picante)
library(phytools)
dat_tree <- read.xlsx("RDA_sp.xlsx")
dat_sp_ls <- CTPL(dat_tree$spece) 
dat <- dat_sp_ls[, c(1,3,5,8)]
write.xlsx(dat, "dat_sp_latin.xlsx")

dat_c_tree <- read.xlsx("dat_sp_latin.xlsx", 4)

tree_TPL <- phylo.maker(dat_c_tree, tree=GBOTB.extended.TPL, nodes=nodes.info.1.TPL, scenarios="S3")
write.tree(tree_TPL$scenario.3, "tree.newick")
plot(tree_TPL$scenario.3)

dat_abun <- read.xlsx("崂山木本.xlsx", 13)
dat_t <- acast(dat_abun, formula = site ~ spe,
               value.var = "abun",
               fill = 0)
dat_t <- read.xlsx("dat_i.xlsx")
write.xlsx(dat_t, "dat_i.xlsx")
tree <- read.tree("tree.newick")
phyr::psd(dat_t, tree)


comm.pd <- pd(dat_t, tree)
head(comm.pd)

comm.bc.dist <- vegdist(dat_t, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
# plot cluster diagram
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")


comm.bc.mds <- metaMDS(dat_t, dist = "bray")
stressplot(comm.bc.mds)
ordiplot(comm.bc.mds, display = "sites", type = "text")

phy.dist <- cophenetic(tree)
# calculate ses.mpd
comm.sesmpd <- ses.mpd(dat_t, phy.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
head(comm.sesmpd)
write.xlsx(comm.sesmpd, "comm.sesmpd.xlsx")
write.xlsx(comm.pd, "comm.pd.xlsx")
comm.ses.mntd <- ses.mntd(dat_t, phy.dist, null.model = "richness", abundance.weighted = FALSE, runs = 999)
write.xlsx(comm.ses.mntd, "comm.ses.mntd.xlsx")
phyr::psv(dat_t, tree)
write.xlsx(phyr::psv(dat_t, tree), "psv.xlsx")

comm.ses.pd <- ses.pd(dat_t, tree, runs = 99)
write.xlsx(comm.ses.pd, "ses.pd.xlsx")

picante::pd(dat_t, tree)
ses.pd(dat_t, tree)
ses.pd(samp, tree, null.model = c("taxa.labels", "richness", "frequency",
                                  "sample.pool", "phylogeny.pool", "independentswap", "trialswap"),
       runs = 999, iterations = 1000, include.root=TRUE)



psr <- phyr::psr(dat_t, tree)
write.xlsx(psr,"psr.xlsx")


dat <- read.xlsx("崂山木本.xlsx", 13)
dat_t <- acast(dat, formula = site ~ spe,
               value.var = "abun",
               fill = 0)


# 谱系结构指数箱线图
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(ggbreak)


library(magrittr)
library(ggsci)
library(ggpmisc)
dat_nri <- read.xlsx("comm.sesmpd.xlsx", 3)

p_box <- ggplot(dat_nri, aes(x = type, y = value))+
  ##'@绘制柱状图
  geom_bar(aes(fill = type), stat = "summary", position = position_dodge(1),
           color = "black",
           fun = mean, size = 0.5)+
  geom_jitter(aes(fill = type), szie = 2.5, alpha = 0.8, shape = 21) +
  scale_fill_manual(values = c("#386cb0","#1b9e77")) + 
  geom_hline(yintercept=1.96, linetype='dotted', col = '#7D4444')+
  annotate("text", x = 0.7, y = 2.0, label = "1.96", vjust = -0.5) + 
  geom_hline(yintercept=-1.96, linetype='dotted', col = '#7D4444')+
  annotate("text", x = 0.7, y = -2.0, label = "-1.96", vjust = -0.5) + 
  # annotate("text", x = 1.5, y = 3, label = "*") + 

  
  ##'@添加误差线
  stat_summary(fun.data = "mean_sd", geom = "errorbar",
               width = 0.2, size = 1)+
  ##'@添加显著性
  geom_signif(comparisons = list(c("NRI","NTI")),
              map_signif_level= F,  ##'@T:显示*号，F显示数字
              tip_length=0, 
              size=1, 
              test = "t.test")+  ##'@t.test, wilcox.test 
  ##'@X轴和Y轴坐标
  labs(x = "", y = "Degree",title = NULL)+
  ##'@设置颜色
  scale_fill_manual(values = c("#386cb0","#1b9e77"))+
  theme_classic()+
  theme(axis.line = element_line(size = 1),  ## 粗细
        text=element_text(family = "sans",colour ="black",size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.ticks = element_line(size = 1,colour = "black"),
        strip.text = element_text(color = "black",size = 16),
        axis.title = element_text(color = "black",size = 18),
        legend.position = "none",
        strip.background = element_blank()
  )

p_box

dat_hg <- read.xlsx("comm.ses.mntd.xlsx", 3)
my.formula <- NRI ~ NTI
p_line <- ggplot(dat_hg, aes(NRI,NTI)) +
  geom_smooth(method = "lm") + 
  geom_point() +
  theme_classic()+
  stat_poly_eq(formula = my.formula, 
               aes(label = paste(..eq.label.., 
                                 ..rr.label.., 
                                 sep = "~~~")), 
               parse = TRUE) +         
  theme(axis.line = element_line(size = 1),  ## 粗细
        text=element_text(family = "sans",colour ="black",size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.ticks = element_line(size = 1,colour = "black"),
        strip.text = element_text(color = "black",size = 16),
        axis.title = element_text(color = "black",size = 18),
        legend.position = "none",
        strip.background = element_blank()
  )
p_line



t_nk <- ggscatter(dat_hg, x = "NRI", y = "NTI",
                              size = 1.5,
                              add = "reg.line",  # 添加回归线
                              add.params = list(color = "#386cb0", fill = "#1b9e77", size = 1),  # 自定义回归线的颜色
                              conf.int = TRUE  # 添加置信区间
                    ) +
    stat_cor(method = "pearson", label.x = -1.2, label.y = 1.5, label.sep = "\n") +
  labs(x = "NRI", y = "NTI",title = NULL)+
  theme(axis.text.y=element_text(size=4),  
  axis.text.x=element_text(size=4)) +
  theme(axis.line = element_line(size = 1),  ## 粗细
        text=element_text(family = "sans",colour ="black",size = 12),
        axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black",size = 12),
        axis.ticks = element_line(size = 1,colour = "black"),
        strip.text = element_text(color = "black",size = 16),
        axis.title = element_text(color = "black",size = 18),
        legend.position = "none",
        strip.background = element_blank()
  )

t_nk 
p_1 <- ggarrange(p_box, t_nk, 
                 labels = c("a", "b"),
                 ncol = 2, nrow = 1)

p_1 
ggsave("p_box.png", width = 7, height = 6.5, dpi = 300)

# NRI与NTI作图
dat_nti_nri <- read.xlsx("F:/博士研究计划2024/大样地数据/系统发育指标/comm.sesmpd.xlsx", 4)
LM_results<- lm(NRI~NTI,data =dat_nti_nri)
summary(LM_results)


R2_expression <- expression(paste(" ", R^2 , "= ", .29))

#add text to plot
text(x =-0.5, y =1.5, label = R2_expression)


p <- ggplot(dat_nti_nri, aes(x=NTI, y=NRI, colour = 21))+
  geom_jitter( position=position_jitter(0.17), size=2)+
  geom_smooth(formula =" y ~x",method = 'lm',level=0.95, se = T, color = "#1b9e77",size=1.4) +
  theme_bw()+theme(axis.text=element_text(colour='black',size=9))+
  theme(axis.text.x=element_text(vjust=1,size=20,face = "bold")) +
  theme(axis.text.y=element_text(vjust=1,size=20,face = "bold")) +

  annotate('text', label = R2_expression, x =-0.9, y =1.3, size =8) +
  annotate('text', label = "p = 0.004*", x =-0.3, y =1.3, size =8) +
  annotate('text', label = 'y = 0.60x - 0.39', x =-0.5, y =1.6, size =8) +
  theme(panel.grid=element_blank()) +
  theme(legend.position = 'none') +
  theme(axis.title = element_text(color = "black",size = 18),
        legend.position = "none",
        strip.background = element_blank()
  ) +
  theme(panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid")) +

theme_classic()+
  theme(axis.line = element_line(size = 1),  ## 粗细
        text=element_text(family = "sans",colour ="black",size = 15),
        axis.text.x = element_text(color = "black", size = 14),
        axis.text.y = element_text(color = "black",size = 14),
        axis.ticks = element_line(size = 1,colour = "black"),
        strip.text = element_text(color = "black",size = 16),
        axis.title = element_text(color = "black",size = 20),
        legend.position = "none",
        strip.background = element_blank()
  )
p
ggsave("nri_hg.png",height = 10, width = 11, dpi = 500)


getwd()


# 发育树与物种多度做图

lsyd_tree <- read.tree("tree.newick")
lsyd_comm_1 <- read.xlsx("raw_abun.xlsx")
lsyd_comm_1 <- as.matrix(dat_t)
write.csv(lsyd_comm_1, "dat_t.csv")
lsyd_comm_1 <- read.xlsx("F:/博士研究计划2024/大样地数据/系统发育指标/dat_i.xlsx", 3)

dat_1 <- tidyr::gather(lsyd_comm_1, key = "sp", value = "freq", -site) 



dat.veg2.wide_1 = dcast(dat_1, site~sp, value.var = "freq")
row.names(dat.veg2.wide_1) = dat.veg2.wide_1$site; dat.veg2.wide_1$site = NULL

pdf(file = "yd_abu_0607.pdf", width = 15, height = 12)
par(mar = c(0.1,0.1,1.1,0.0))
plot(lsyd_tree,
     show.tip.label=T, cex=1, x.lim=1220,
     label.offset=4, edge.width=1, direction = "rightwards")
for(i in 1:25){
  tiplabels(tip = which(lsyd_tree$tip.label %in%
                          names(dat.veg2.wide_1)[
                            dat.veg2.wide_1[i, ] > 0]), pch = 20,
            cex=dat.veg2.wide_1[i, ][dat.veg2.wide_1[i, ]>0]/15,
            adj = c(970-i*30, 0.5),
            col = "blue")
}
# mtext(text = "1958", side = 2, line = 1, at = 1000, col = "blue")
mtext(text = "Quadrats (20 × 20m)", side = 3, line = 0, at = 800,col = "black")
dev.off()


H# 01数据

lsyd_tree <- read.tree("tree.newick")
lsyd_comm <- read.xlsx("raw_abun.xlsx", 4)

dat <- tidyr::gather(lsyd_comm, key = "sp", value = "freq", -site) 



dat.veg2.wide = dcast(dat, site~sp, value.var = "freq")
row.names(dat.veg2.wide) = dat.veg2.wide$site; dat.veg2.wide$site = NULL

pdf(file = "dune_phylo_2.pdf", width = 10, height = 12)
par(mar = c(0.1,0.1,1.1,0.0))
plot(lsyd_tree,
     show.tip.label=T, cex=1, x.lim=1200,
     label.offset=4, edge.width=1, direction = "rightwards")
for(i in 1:25){
  tiplabels(tip = which(lsyd_tree$tip.label %in%
                          names(dat.veg2.wide)[
                            dat.veg2.wide[i, ] > 0]), pch = 20,
            cex=dat.veg2.wide[i, ][dat.veg2.wide[i, ]>0]/3,
            adj = c(950-i*30, 0.5),
            col = "blue")
}
# mtext(text = "1958", side = 2, line = 1, at = 1000, col = "blue")
mtext(text = "Quadrats (20 × 20m)", side = 3, line = 0, at = 600, col = "black")
dev.off()

# 谱系主坐标分析 -----------------------------------------------------------------
setwd("F:/博士研究计划2024/大样地数据/谱系主坐标分析")
library(PCPS)
library(ggplot2)
library(openxlsx)
library(magrittr)
library(ggnewscale)
library(scales)
library(viridis) 
library(ggpubr)
library(vegan)
library(picante)
library(phytools)
lsyd_tree <- read.tree("tree.newick")
phylodist <- as.matrix(as.dist(cophenetic.phylo(lsyd_tree)))

# comm <- read.csv("W.csv", row.names = 1)
# p.dist <- read.csv("phylodist.csv", row.names = 1)
# comm <- read.xlsx("pcps.xlsx", rowNames = TRUE)
# res <- pcps(comm,p.dist)
phylodist <- as.data.frame(phylodist)
comm <- read.xlsx("raw_abun.xlsx", 5,rowNames = TRUE)
envi <- read.xlsx("sumpcps.xlsx", 3, rowNames = TRUE)
comm <- t(comm)
res_phytool <- phyl.pca(lsyd_tree, comm) 
plot(res_phytool)
biplot(res_phytool)
plot(res_phytool$V)

summary(res_phytool)

res <- pcps(comm,phylodist)
tbpca_ef <- envfit(dat_pcps[,2:3]~., data = envi, perm = 999, choices = c(1,2), display = 'sites')
tbpca_ef
tbpca_ef_adj <- tbpca_ef
tbpca_ef_adj$vectors$pvals <- p.adjust(tbpca_ef_adj$vectors$pvals, method = 'bonferroni')
tbpca_ef_adj
# 环境箭头
tbpca_env <- data.frame(cbind(tbpca_ef_adj$vectors$arrows, tbpca_ef_adj$vectors$r, tbpca_ef_adj$vectors$pvals))
names(tbpca_env) <- c('PCPS1', 'PCPS2', 'r2', 'p.adj')
write.csv(tbpca_env, 'pca_env.csv', quote = FALSE)

tbpca_env_cor <- cor(envi[ ,c('TN', 'TP')], dat_pcps[,2:3], method = 'pearson')

#或者由排序坐标和 r2 反向推算出环境变量与样方得分的相关性
arrow_heads <- tbpca_ef_adj$vectors$arrows       #提取 env 向量（坐标）矩阵
r2 <- tbpca_ef_adj$vectors$r #提取 env r2
arrow_heads * sqrt(r2)   #即可得相关性，结果和上述“tbpca_env”是一致的




res_p <- res$P
res$P
windowsFonts(A=windowsFont("Times New Roman"), B=windowsFont("Arial"))
sum_pcps <- summary(res, choices = c(1, 2))$scores$scores.sites
sum_pcps <- as.data.frame(sum_pcps)
write.xlsx(sum_pcps, "sumpcps.xlsx")


plot(res)
plot(res, display = "text", groups = c(rep("Clade-A", 1), rep("Clade-B", 1),
                                       rep("Clade-C", 4),rep("Clade-D", 24),rep("Clade-E", 26)),cex=0.6)

class(res)

data_1 <- list(summary(res)$scores)

data_site <- data.frame("site" = data_1[[1]]$scores.sites)
data_species <- as.data.frame(data_1[[1]]$scores.species)
colnames(data_site) <- c("pcps1", "pcps2")
colnames(data_species) <- c("PCPS1", "PCPS2")
clade <- read.xlsx("clade.xlsx",2)
nri <- read.xlsx("NRI.xlsx", 2)
nri <- round(nri, 2)
data_species <- cbind(data_species, clade)
data_site <- cbind(data_site, nri)
centroid <- aggregate(data_species[,1:2], by = list(Clade = data_species$Clade),
                      data = data_species,
                      FUN = mean)

pca.results1 <- dplyr::left_join(data_species,centroid, by = "Clade")

library(RColorBrewer)
library(scales)
show_col(brewer.pal(5, "Set1"))
cols = brewer.pal(5, "Set1")
cols
fix(data_site)

round(data_site[,3],2)
write.xlsx(data_site, "site.xlsx")
SITE <- read.xlsx("site.xlsx", 2)
fix(SITE)
p <- ggplot(pca.results1,aes(x=PCPS1.x,y=PCPS2.x))+
  geom_point() +
  geom_vline(xintercept = 0,lty="dashed")+
  geom_hline(yintercept = 0,lty="dashed")+
  geom_point(aes(color=Clade))+
  geom_segment(aes(xend=PCPS1.y,yend=PCPS2.y,color=Clade),
               show.legend = F) +
  geom_label(data = centroid, 
             aes(x = PCPS1, y = PCPS2, label = Clade, fill = Clade), size = 5, 
             show.legend = FALSE,
             color="white") +
  scale_fill_manual(values = cols)+
  theme(legend.position = "top") + 
  geom_point(data = SITE, aes(x=pcps1, y=pcps2)) 

p


nri <- read.xlsx("NRI.xlsx")
shapiro.test(nri$NRI)
t.test(nri$NRI, mu = 0)

# 

dca <- decorana(veg = t(res_p))
dca
res_p <- decostand(res_p, method = "hellinger")
res_rda <- rda(res_p, raw_soil, scale = FALSE)
plot(res_rda)

res_1 <- rda(res_p, raw_soil, scale = TRUE)
plot(res_1)
res_1 <- rda(res_p~., raw_soil, scale = TRUE)
anova(res_1, permutations = how(nperm = 999))
summary(res_1)
RsquareAdj(res_1)


anova(res, permutations = how(nperm = 999))
summary(res)
res <- rda(raw_abun~., raw_soil, scale = FALSE)
anova(res, permutations = how(nperm = 999))
summary(res)
raw_abun <- as.matrix(raw_abun)
raw_soil <- as.matrix(raw_soil)

RDA.Perm = permutest(res, permutations = 999)
RDA.Perm
envfit <- envfit(res, raw_soil, permutations  = 999)
envfit
RsquareAdj(res)


zh <- summary(res_1)$concont$importance
zh <- round(zh, 4)

pdat <- res_1$CCA
samples <-data.frame(sample = row.names(pdat$u),RDA1 = pdat$u[,1],RDA2 = pdat$u[,2])
species<-data.frame(spece = row.names(pdat$v),RDA1 = pdat$v[,1],RDA2 = pdat$v[,2])
envi<-data.frame(en = row.names(pdat$biplot),RDA1 = pdat$biplot[,1],RDA2 = pdat$biplot[,2])
species <- read.xlsx("RDA_sp.xlsx")

# 崂山论文PCPS绘图
dat_pcps <- read.xlsx("sumpcps.xlsx", 2)

p <- ggplot() +
  geom_point(data = dat_pcps, aes(x = PCPS1, y = PCPS2),size = 2,  shape = 16) +
  # stat_ellipse(data = species,aes(x = RDA1, y = RDA2), level = 0.95,linetype = "dashed", show.legend = F, size = 1.5) +
  # annotate('text', label = 'TN*', x = -0.5, y = -0.39, size = 5, colour = 'black') +
  # annotate('text', label = 'TP*', x = -0.5, y = -0.40, size = 3, colour = 'black') +
  geom_hline(aes(yintercept = 0), colour="gray88", linetype="dashed") +
  geom_vline(aes(xintercept = 0), colour="gray88", linetype="dashed")  +
  #stat_ellipse(level = 0.95, show.legend = F, size = 1.5, color = species$group.ad) +
  #geom_point(data = species, aes(x=RDA1, y=RDA2, size = 5,color = group$ad)) +
  #scale_color_gradient(low = "#f47720",high = "#88c4e8")+
  #geom_point(aes(clour =factor(1)), size = 5)  +
  #geom_point(color = "green", size = 1.5)
  #scale_colour_manual(values =c("green","yellow"), ) +
  #geom_segment(data = species, aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), color = 'blue') +
  geom_text(data = dat_pcps,aes(x = PCPS1*1.1, y = PCPS2*1.1, label = site), size = 5, color = '#f47720') +
  geom_segment(data = tbpca_env_cor,aes(x=0, xend= PCPS1, y=0, yend= PCPS2), size = 1,arrow = arrow(length = unit(0.3, "cm")), colour = '#88c4e8') +
  annotate('text', label = 'TN*', x = -0.5, y = -0.37, size = 5, colour = 'black') +
  annotate('text', label = 'TP*', x = -0.5, y = -0.43, size = 5, colour = 'black') +
  # geom_text(data = tbpca_env_cor,aes(x = PCPS1*1.1, y = PCPS2*1.1, label = row.names(tbpca_env_cor)), size = 3, colour = 'black', check_overlap = FALSE) +
  #geom_point(data = samples, aes(x=RDA1, y=RDA2),size = 3, shape = 25, color = "black", fill = "black")  +
  theme_bw() + theme(axis.ticks.length = unit(-5, "pt"), panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"), panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12,face = 'bold'), axis.title.x=element_text(vjust=1, size=13), axis.title.y=element_text(size=13)) +
  theme(legend.position = 'none') +
  xlab(paste('PCPS1 (', 62.59, '%)', sep ='')) + ylab(paste('PCPS2 (', 15.12, '%)', sep =''))
p
ggsave("PCPS_2.png", height = 8, width = 9,dpi = 500)



??aov

# 回归
library(MASS)
library(glmm.hp)
sum_pcps <- read.xlsx("sumpcps.xlsx")
lm_1 <- lm(p1~pH + EC + TN + TP + SOC + M + W+ U, data = sum_pcps)
summary(lm_1)
step.model <- stepAIC(lm_1, direction = "both", trace = FALSE)
summary(step.model)

glmm.hp(step.model, type = "R2")



lm_2 <- lm(p2~pH + EC + TN + TP + SOC + M + W+ U, data = sum_pcps)
summary(lm_2)
step.model_1 <- stepAIC(lm_2, direction = "both", trace = FALSE)
summary(step.model_1)












# 主坐标分析
setwd("F:/博士研究计划2024/大样地数据/主坐标分析")
library(openxlsx)
library(vegan)
library(ggplot2)
dat_comm <- read.xlsx("raw_abun.xlsx", rowNames = TRUE) 
group_comm <- read.xlsx("group.xlsx", rowNames = TRUE)


library(ggplot2)
library(ade4)   # 用于计算PcoA
library(vegan)  # 用于计算距离




# PCoA计算
df.dist = vegdist(dat_comm,method='bray')    #基于euclidean距离
pcoa =  dudi.pco(df.dist,
                 scannf = F,   # 一种逻辑值，指示是否应该显示特征值条形图
                 nf=2)         # 保留几个维度的坐标信息



# 置换多元方差分析
dim <- vegdist(dat_comm, method = "bray", diag = T, upper = T)
dim
dim <- as.matrix(dim)



Adonis <- adonis2(dim ~ sites$group,
                  distance = "bray",
                  permutations = 999)
Adonis





adonis <- paste0("adonis R2: ",round(Adonis$R2,2), 
                 "; P-value: ", Adonis$`Pr(>F)`)

# 整理绘图所需的数据
data = pcoa$li
data$name = rownames(data)
data$group = group_comm$group

# 绘图
ggplot(data,aes(x = A1,
                y = A2,
                color = group,
                group = group,
                fill = group
))+
  geom_point()+
  theme_classic()+
  geom_vline(xintercept = 0, color = 'gray', size = 0.4) +   # 在0处添加垂直线条
  geom_hline(yintercept = 0, color = 'gray', size = 0.4) +
  stat_ellipse(aes(x=A1,    # 添加置信区间圈
                   y=A2,
  ),
  geom = "polygon",
  level = 0.95,
  alpha=0.4)+
  geom_text(                # 添加文本标签
    aes(label=name),   
    vjust=1.5,            
    size=2,
    color = "black"
  )+
  labs(  # 更改x与y轴坐标为pcoa$eig/sum(pcoa$eig)
    x = paste0("PCoA1 (",as.character(round(pcoa$eig[1] / sum(pcoa$eig) * 100,2)),"%)"),
    y = paste0("PCoA2 (",as.character(round(pcoa$eig[2] / sum(pcoa$eig) * 100,2)),"%)")
  )


raw_env <- read.xlsx("envi.xlsx", rowNames = TRUE)

envfit <- envfit(pcoa$li, raw_env, permutations  = 999)
env <- envfit$vectors$arrows
env <- as.data.frame(env)

p.adjust(envfit$vectors$pvals, method = 'bonferroni')

sites <- data.frame(row.names(group_comm),PCoA1 = pcoa$li[,1],PCoA2 = pcoa$li[,2], group_comm$group)
envi <- data.frame(row.names(env), PCoA1 = env$A1,PCoA2 = env$A2)

colnames(sites)[1:4] <- c("site","PCoA1","PCoA2", "group")
colnames(envi)[1] <- "site"
windowsFonts(A=windowsFont("Times New Roman"), B=windowsFont("Arial"))

p <- ggplot() +
  geom_point(data = sites, aes(x = PCoA1, y = PCoA2, size = 5, color = group, shape = group),show.legend = FALSE) +
  stat_ellipse(data = sites,aes(x = PCoA1, y = PCoA2, group = group, fill = group),level = 0.95, geom = "polygon", alpha=0.4)+
  #stat_ellipse(data = species,aes(x = RDA1, y = RDA2, group = group.ad), level = 0.95,linetype = "dashed", show.legend = F, size = 1.5, color = species$group.ad) +
  # annotate('text', label = 'Evergreen', x = 0.25, y = 0, size = 5, colour = 'black') +
  # annotate('text', label = 'Deciduous', x = -0.22, y = -0.25, size = 3, colour = 'black') +
  geom_hline(aes(yintercept = 0), colour="gray88", linetype="dashed") +
  geom_vline(aes(xintercept = 0), colour="gray88", linetype="dashed")  +
  #stat_ellipse(level = 0.95, show.legend = F, size = 1.5, color = species$group.ad) +
  #geom_point(data = species, aes(x=RDA1, y=RDA2, size = 5,color = group$ad)) +
  # scale_color_gradient(low = "#f47720",high = "#88c4e8")+
  #geom_point(aes(clour =factor(1)), size = 5)  +
  #geom_point(color = "green", size = 1.5)
  #scale_colour_manual(values =c("green","yellow"), ) +
  #geom_segment(data = species, aes(x=0, xend= RDA1, y=0, yend= RDA2 ), arrow = arrow(length = unit(0.3, "cm")), color = 'blue') +
  geom_text(data = sites,aes(x = PCoA1*1.1, y = PCoA2*1.1, label = site), size = 5, color = 'black', check_overlap = FALSE,show.legend = FALSE) +
  geom_segment(data = envi,aes(x=0, xend= PCoA1, y=0, yend= PCoA2), size = 1,arrow = arrow(length = unit(0.3, "cm")), colour = 'black') +
  geom_text(data = envi,aes(x = PCoA1*1.1, y = PCoA2*1.1, label = site), size = 5, colour = 'black', check_overlap = TRUE) +
  #geom_point(data = samples, aes(x=RDA1, y=RDA2),size = 3, shape = 25, color = "black", fill = "black")  +
  theme_bw() + theme(axis.ticks.length = unit(-5, "pt"), panel.border = element_rect(fill=NA,color="black", size=1.5, linetype="solid"), panel.grid.major=element_line(colour=NA), panel.grid.minor = element_blank()) +
  theme(axis.text = element_text(size=12,face = 'bold'), axis.title.x=element_text(vjust=1, size=13), axis.title.y=element_text(size=13)) +
  xlab(paste('PCoA1 (', 20.83, '%)', sep ='')) + ylab(paste('PCoA2 (', 17.91, '%)', sep =''))+
  theme(axis.title.x = element_text(size = 18)) +
  theme(axis.title.y = element_text(size = 18)) +
  theme(legend.position = c(0.9, 0.8)) +
  ggtitle(label = adonis) +
  theme(plot.title = element_text(size = 10,hjust = 0.9 , vjust = - 10 )) 
  

p


# theme(plot.subtitle=element_text(face="bold.italic", #字体
#                                  color="black", #颜色
#                                  size=9,  #大小
#                                  hjust=0.5, #位置
#                                  vjust=0,
#                                  angle=360), # 角度
#       plot.caption=element_text(face="bold",color="black",size=10))


ggsave("pcoa.svg", height = 10, width = 10,dpi = 300)

ggsave("pcoa_1.png", height = 8, width = 9,dpi = 500)

# 指示物种分析
install.packages("indicspecies")
library(indicspecies)

otu_transposed <- dat_comm
group_labels <- group_comm$group

# 计算指示物种
set.seed(123) # 为了结果可重现
indicator_results <- multipatt(otu_transposed, # OTU表
                               group_labels, # 分组标签
                               func = "r.g", # 选择指示物种方法
                               control = how(nperm = 999)) # 根据需要调整nperm

# 汇总指示物种结果
summary(indicator_results, # 指示物种结果
        alpha = 1, # 显著性水平
        indvalcBBp = TRUE) # 是否计算指示物种的p值

# 选择并提取显著指示物种
significant_indicators <- indicator_results$sign
significant_matrix <- as.matrix(significant_indicators[which(significant_indicators$p.value < 0.05), ])





network_data <- list()

# 定义处理列表
# 需要根据你自己的处理名称进行调整，但是前面需要加上"s."前缀
treatments <- c("s.A", "s.B", "s.C")

# 对每个处理进行遍历
for(treatment in treatments) {
  # 选择当前处理下指示值为1的OTU
  species <- rownames(significant_indicators)[significant_indicators[[treatment]] == 1]
  
  # 检查是否有指示物种，如果有，则添加到网络数据中
  if(length(species) > 0) {
    for(specie in species) {
      network_data[[length(network_data) + 1]] <- data.frame(
        Source = treatment,
        Target = specie,
        Weight = significant_indicators[specie, "stat"]
      )
    }
  }
}

# 合并所有边到一个数据框中
network_df <- do.call(rbind, network_data)

# 如果network_data为空，这将避免错误
if(length(network_data) == 0) {
  network_df <- data.frame(Source = character(), Target = character(), Weight = numeric())
}








library(igraph) # 加载igraph包

g <- graph_from_data_frame(d = network_df, directed = FALSE) # 从数据框创建网络

# 可视化网络
p_indicator <- plot(g, vertex.size=10, # 设置节点大小
     vertex.label.cex=0.8, # 设置节点标签大小
     edge.arrow.size=.5,  # 设置边箭头大小
     layout=layout_with_fr(g)) # 使用Fruchterman-Reingold布局
p_indicator

ggsave("indicator.png", height = 10, width = 9, dpi = 500)





rm(list = ls())
# phyr --------------------------------------------------------------------
setwd("F:/博士研究计划2024/大样地数据/phyr")
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
library(devtools)
library(learnasreml)
library(lattice)
library(glmmTMB)
library(ggplot2)


dat_trait <- read.xlsx("F:/博士研究计划2024/大样地数据/种间联结计算/每木调查_崂山.xlsx", 5)
phy_tree_ <- read.tree("tree.newick")
plot(phy_tree)
comm_phyr <- read.xlsx("raw_abun.xlsx", 5)
row.names(comm_phyr) <- comm_phyr$site
envi_phyr <- read.xlsx("envi-2.xlsx")
dat_phyr <- tidyr::gather(comm_phyr, key = "sp", value = "freq", -site) %>% left_join(envi_phyr, by = "site")


dat_phyr$pa = as.numeric(dat_phyr$freq > 0)

# 亲缘关系相近的物种是否会出现在同一个样方？
mod_1 <- communityPGLMM(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
                        data = dat_phyr, cov_ranef = list(sp = phy_tree))

mod_1_r <- communityPGLMM(freq ~ 1 + (1|sp__) + (1|site) , 
                        data = dat_phyr, cov_ranef = list(sp = phy_tree))
       
# P-value of the underdispersion can be calculated with likelihood ratio test
pchisq(2*(mod_1$logLik - mod_1_r$logLik), df = 1, lower.tail = F)/2

# 0.4903404

# pvalue <- pchisq(2*(mod_1$logLik - mod_1_r$logLik), df = 1, lower.tail = FALSE)
# pvalue


summary(mod_1)
summary(mod_1_r)
# 
# 
# Test phylogenetic overdispersion
mod_2 <- communityPGLMM(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), 
                        data = dat_phyr, cov_ranef = list(sp = phy_tree), repulsion = TRUE)

mod_2_r <- communityPGLMM(freq ~ 1 + (1|sp__) + (1|site) , 
                          data = dat_phyr, cov_ranef = list(sp = phy_tree))

pchisq(2*(mod_2$logLik - mod_2_r$logLik),
       df = 1, lower.tail = F)/2

# 5.882197e-08


summary(mod_2)
summary(mod_2_r)

rr2::R2_pred(mod_2)
# 0.2701217
rr2::R2_pred(mod_2_r)
# 0.2917726

#LOG
dat_phyr$Y <- log(dat_phyr$freq + 1)

mod_1 <- communityPGLMM(Y ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
                        data = dat_phyr, cov_ranef = list(sp = phy_tree))

mod_1_r <- communityPGLMM(Y ~ 1 + (1|sp__) + (1|site) , 
                          data = dat_phyr, cov_ranef = list(sp = phy_tree))




# test_1 = phyr::pglmm(freq ~ 1 + N + P + SOM + K + (1|sp__) + (1|site) + (1|sp__@site),
#                     data = dat, family = "gaussian", REML = FALSE,
#                     cov_ranef = list(sp = phy))
# 
# 
# 
# 
# 
# summary(LM_1)
# glmm.hp(LM_1)
# performance::r2(LM_1)
# anova(LM_1, type = "I")
# test1 = phyr::pglmm(freq ~ 1 + N + P + SOM + K + (1|sp__) + (1|site) + (1|sp__@site),
#                     data = dat, family = "gaussian", REML = FALSE,
#                     cov_ranef = list(sp = phy))
# 
# test2 = phyr::pglmm(pa ~ 1 + Ele + (1|sp__) + (1|site) + (1|sp__@site), 
#                     data = dat, family = "binomial", REML = FALSE,
#                     cov_ranef = list(sp = phy))
# summary(test2)



# possion distribution

mod_1 <- pglmm(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
                        data = dat_phyr,family = "poisson", cov_ranef = list(sp = phy_tree_),add.obs.re= TRUE)
summary(mod_1)

mod_2 <- pglmm(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), 
                        data = dat_phyr, family = "poisson",cov_ranef = list(sp = phy_tree_), repulsion = TRUE,add.obs.re= TRUE)
summary(mod_2)

pchisq(2*(mod_1$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2

mod_2_r <- pglmm(freq ~ 1 + (1|sp__) + (1|site) , 
                          data = dat_phyr, family = "poisson", cov_ranef = list(sp = phy_tree))

summary(mod_2_r)

pchisq(2*(mod_2$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2

# binomial 检验群落系统发育结构1104

mod_1 <- pglmm(pa ~ 1 + (1|site) + (1|sp__@site),
               data = dat_phyr,family = "binomial", cov_ranef = list(sp = phy_tree_),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_1)

mod_2 <- pglmm(pa ~ 1  + (1|site) + (1|sp__@site), 
               data = dat_phyr, family = "binomial",cov_ranef = list(sp = phy_tree_), repulsion = TRUE,s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_2)

pchisq(2*(mod_1$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2 #0.5

mod_2_r <- pglmm(pa ~ 1 + (1|sp__) + (1|site) , 
                 data = dat_phyr, family = "binomial", cov_ranef = list(sp = phy_tree_),s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_2_r)

pchisq(2*(mod_2$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2#0.5


# possion 检验群落系统发育结构1104

mod_1 <- pglmm(freq ~ 1 + (1|site) + (1|sp__) + (1|sp__@site),
               data = dat_phyr,family = "poisson", cov_ranef = list(sp = phy_tree_),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_1)

mod_2 <- pglmm(pa ~ 1  + (1|site) + (1|sp__@site), 
               data = dat_phyr, family = "poisson",cov_ranef = list(sp = phy_tree_), repulsion = TRUE,s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_2)

pchisq(2*(mod_1$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2 #0.5

mod_2_r <- pglmm(pa ~ 1 + (1|sp__) + (1|site) , 
                 data = dat_phyr, family = "poisson", cov_ranef = list(sp = phy_tree_),s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_2_r)

pchisq(2*(mod_2$logLik - mod_2_r$logLik), df = 1, lower.tail = F)/2#0.5






# 
test_1 = pglmm(freq ~ U + (1|sp__) + (1|site) + (1|sp__@site),
                     data = dat_phyr, family = "gaussian", REML = FALSE,
                     cov_ranef = list(sp = phy_tree), repulsion = TRUE)


summary(test_1)

test_2 <- pglmm(freq ~ U + M + W + pH+ SOC + TN + TP + EC + (1|sp__) + (1|site) + (1|sp__@site),
                     data = dat_phyr, family = "gaussian", REML = FALSE,
                     cov_ranef = list(sp = phy_tree), repulsion = TRUE)

summary(test_2)

test_3 = phyr::pglmm(freq ~ 1 + U + M + W + TP+ SOC + TN + (1|sp__) + (1|site) + (1|sp__@site),
                     data = dat_phyr, family = "gaussian", REML = FALSE,
                     cov_ranef = list(sp = phy_tree), repulsion = TRUE)

summary(test_3)


test_4 = phyr::pglmm(freq ~ 1 + U + M + W + TP+ SOC + TN + (1|sp__) + (1|site),
                     data = dat_phyr, family = "gaussian", REML = FALSE,
                     cov_ranef = list(sp = phy_tree))

summary(test_4)

# start.model = lmer(pa ~ 1 + (1|sp) + (1|site), data = dat_phyr, REML = FALSE)
# SigAIC(start.model)
# AIC(start.model)


data("comm_a")
??comm_a

save(mod_1, mod_1_r, mod_2, dat_phyr, phy_tree,file = "lsyd.RData")

load("lsyd.RData")


# phylo_pattern
library(vegan)
library(stringr) # work with strings
library(dplyr) # data manipulation
library(brranching) # for pyhlomatic, to get the phylogeny
library(tidyr) # data manipulation
library(reshape2) # long to wide table
library(ape) # plot phylogeny
library(pez) # for communityPGLMM
library(lme4)
library(phylolm) # to test phylo signal of traits
library(picante)
library(phylocomr)
##################### analyses ---------------------
#### Patterns (Table 1) ----


dat_phyr$freq <- log(dat_phyr$freq + 1)

# #########亲缘关系相近的物种是否会出现在同一个样方？
mod_1 <- pglmm(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site),
                        data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))

mod_1_r <- pglmm(freq ~ 1 + (1|sp__) + (1|site) , 
                          data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))

# P-value of the underdispersion can be calculated with likelihood ratio test
pchisq(2*(mod_1$logLik - mod_1_r$logLik), df = 1, lower.tail = F)/2

# 0.1762258
summary(mod_1)
summary(mod_1_r)
# 
# 
# Test phylogenetic overdispersion
mod_2 <- pglmm(freq ~ 1 + (1|sp__) + (1|site) + (1|sp__@site), 
                        data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

mod_2_r <- pglmm(freq ~ 1 + (1|sp__) + (1|site) , 
                          data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))

pchisq(2*(mod_2$logLik - mod_2_r$logLik),
       df = 1, lower.tail = F)/2

# 1.375661e-14
summary(mod_2)
summary(mod_2_r)


mod_2__ <- pglmm(freq ~ 1 + pH + EC + TN + TP + SOC + W + M + U + (1|sp__) + (1|site) + (1|sp__@site), 
               data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))


mod_2__ <- pglmm(freq ~ 1 + TN + TP + SOC  + (1|sp__) + (1|site) + (1|sp__@site), 
                 data = dat_phyr, family = "gaussian", cov_ranef = list(sp = phy_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_2__)


dat_d <- read.xlsx("F:/博士研究计划2024/崂山英文想法/功能性状/dat_phyr.xlsx",3)

mod_2___ <- pglmm(Y ~ 1 + iv + niches + nichel + (1|sp__) + (1|site) + (1|sp__@site), 
                 data = dat_d, family = "gaussian", cov_ranef = list(sp = phy_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_2___)

mod_3___ <- pglmm(Y ~ 1 + iv + niches + nichel+ (1|iv)+ (1|niches)+ (1|nichel) + (1|sp__) + (1|site) + (1|sp__@site), 
                  data = dat_d, family = "gaussian", cov_ranef = list(sp = phy_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_3___)



# 物种分布对环境的响应 -------------------------------------------------------------------
library(openxlsx)
phy_tree_ls <- read.tree("tree.newick")
comm_phyr_ls <- read.xlsx("raw_abun.xlsx", 5)
row.names(comm_phyr_ls) <- comm_phyr_ls$site
envi_phyr_ls <- read.xlsx("envi-2.xlsx")
dat_phyr_ls <- tidyr::gather(comm_phyr_ls, key = "sp", value = "freq", -site) %>% left_join(envi_phyr_ls, by = "site")
dat_phyr_ls$freq <- log(dat_phyr_ls$freq + 1)


phylo_signal_slopes_envi(dat = dat_phyr_ls, phylo = phy_tree_ls)

# this functional can be used to search environmental variables that closely related
# species have similar responses to. This may provide insight about future traits
# to measure.
# dat: long format for both veg and env data
# veg.long: used only if dat is NULL
# envi: environmental variables as columns, sites as rows

phylo_signal_slopes_envi = function(dat = NULL, veg.long = NULL, phylo, envi = NULL,  
                                    binary = FALSE){
  # remove singleton sp and include envi variables and scale them.
  if(is.null(dat)){
    veg.long = filter(veg.long, sp %in% phylo$tip.label)
    dat = left_join(veg.long, envi, by = "site")
  } else{
    dat = filter(dat, sp %in% phylo$tip.label)
  }
  
  dat$sp <- as.factor(dat$sp)
  dat$site <- as.factor(dat$site)
  
  nspp <- nlevels(dat$sp)
  nsite <- nlevels(dat$site)
  
  # transform frequency data
  dat$presence <- as.numeric(dat$freq > 0)
  dat$Y <- log(dat$freq + 1)
  
  # the covar matrix to be standardized to have determinant  of 1
  phy <- drop.tip(phylo, tip=phylo$tip.label[
    !phylo$tip.label %in% unique(dat$sp)])
  Vphy <- vcv(phylo)
  Vphy <- Vphy[order(phy$tip.label), order(phy$tip.label)]
  Vphy <- Vphy/max(Vphy)
  Vphy <- Vphy/det(Vphy)^(1/nspp)
  
  show(c(nspp, Ntip(phy)))
  if(nspp != Ntip(phy)){
    stop("The vegetation data and the phylogeny have different number of species")
  }
  
  # output
  envir <- data.frame(var=names(dat)[4:(dim(dat)[2] - 2)], 
                      lmer.logLik=NA, lmer.logLik0=NA, lmer.X.Pr=NA, 
                      PGLMM.intercept.star.s2=NA, PGLMM.intercept.phy.s2=NA, 
                      PGLMM.slope.star.s2=NA, PGLMM.slope.phy.s2=NA, 
                      PGLMM.resid.s2=NA, PGLMM.slopephy.logLik=NA,
                      PGLMM.slopephy.logLik0=NA, PGLMM.slope.phy.Pr=NA, 
                      PGLMM.slope.star.logLik=NA, PGLMM.slope.star.logLik0=NA, 
                      PGLMM.slope.star.Pr=NA, PGLMM.intercept.phy.Pr = NA,
                      PGLMM.B = NA, PGLMM.B.Pr = NA, PGLMM.B.se = NA)
  levels(envir$var) <- names(dat)[4:(dim(dat)[2] - 2)] # presence, Y at the end
  
  # start with the fourth column
  for (i in 4:(dim(dat)[2] - 2)) {
    print(i)
    if(binary == FALSE){
      # lmer
      z.1 <- lmer(Y ~ 1 + (1 | sp) + as.matrix(dat[, i]) + (0 + as.matrix(dat[, i]) | sp), 
                  data = dat, REML=F)
      z.0 <- lmer(Y ~ 1 + (1 | sp) + as.matrix(dat[, i]), data = dat, REML=F)
      envir[i-3,1] <- names(dat)[i]
      envir[i-3,2:4] <- c(logLik(z.1), logLik(z.0), 
                          pchisq(2 * (logLik(z.1) - logLik(z.0)), df=1, lower.tail=F)/2)
    }
    
    if(binary == TRUE){
      # lmer
      z.1 <- glmer(presence ~ 1 + (1 | sp) + as.matrix(dat[, i]) + (0 + as.matrix(dat[, i])| sp), 
                   data = dat, family = binomial)
      z.0 <- glmer(presence ~ 1 + (1 | sp) + as.matrix(dat[, i]), 
                   data = dat, family = binomial)
      envir[i-3,1] <- names(dat)[i]
      envir[i-3,2:4] <- c(logLik(z.1), logLik(z.0), 
                          pchisq(2 * (logLik(z.1) - logLik(z.0)), df=1, lower.tail=F)/2)
    }
    # re.1.site <- list(1, site = dat$site, covar = diag(nsite))
    re.1.star <- list(1, sp = dat$sp, covar = diag(nspp))
    re.1.phy <- list(1, sp = dat$sp, covar = Vphy)
    
    re.X.star <- list(as.matrix(dat[, i]), sp = dat$sp, covar = diag(nspp))
    re.X.phy <- list(as.matrix(dat[, i]), sp = dat$sp, covar = Vphy)
    
    # the goal is to test the significance of re.X.phy (phylo signal of slopes in abund~envi)
    if(binary == FALSE){
      z <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian",
                          sp = dat$sp, site = dat$site, 
                          random.effects = list(re.1.star, re.1.phy, re.X.star, re.X.phy), 
                          REML = F, verbose = F, s2.init=.1)
      print(z$convcode)
      
      z0 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = list(re.1.star, re.1.phy, re.X.star), 
                           REML = F, verbose = F, s2.init=z$ss[c(1,2,3,5)]^2)
      print(z0$convcode)
      
      z00 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                            sp = dat$sp, site = dat$site, 
                            random.effects = list(re.1.star, re.1.phy), 
                            REML = F, verbose = F, s2.init=z$ss[c(1,2,3)]^2)
      print(z00$convcode)
      
      z000 <- communityPGLMM(Y ~ 1 + as.matrix(dat[, i]), data = dat, family = "gaussian", 
                             sp = dat$sp, site = dat$site, 
                             random.effects = list(re.1.star), 
                             REML = F, verbose = F, s2.init=z$ss[c(1,3)]^2)
      print(z000$convcode)
      
      envir[i-3,5:19] <- c(z$ss^2, 
                           z$logLik, 
                           z0$logLik, 
                           pchisq(2 * (z$logLik - z0$logLik), df=1, lower.tail=F)/2, 
                           z0$logLik, 
                           z00$logLik, 
                           pchisq(2 * (z0$logLik - z00$logLik), df=1, lower.tail=F)/2,
                           pchisq(2 * (z00$logLik - z000$logLik), df=1, lower.tail=F)/2,
                           z$B[2], z$B.pvalue[2], z$B.se[2])  
    }
    
    if(binary == TRUE){
      z <- communityPGLMM(presence ~ 1 + as.matrix(dat[, i]), data = dat, family = "binomial", 
                          sp = dat$sp, site = dat$site, 
                          random.effects = list(re.1.star, re.1.phy, re.X.star, re.X.phy), 
                          REML = F, verbose = F, s2.init = .1)
      print(z$convcode)
      w0 = communityPGLMM.binary.LRT(z, re.number = 4)
      z0 <- communityPGLMM(presence ~ 1 + as.matrix(dat[, i]), data = dat, family = "binomial", 
                           sp = dat$sp, site = dat$site, 
                           random.effects = list(re.1.star, re.1.phy, re.X.star), 
                           REML = F, verbose = F, s2.init = .1)
      print(z0$convcode)
      w00 = communityPGLMM.binary.LRT(z0, re.number = 3)
      #         z00 <- communityPGLMM(presence ~ 1 + dat[, i], data = dat, family = "binomial", 
      #                               sp = dat$sp, site = dat$site, 
      #                               random.effects = list(re.1.star, re.X.star), 
      #                               REML = F, verbose = F, s2.init=z$ss[c(1,3,5)]^2)
      envir[i-3,5:19] <- c(z$ss^2, NA, # no s2.resid
                           NA, # no logLik output
                           w0$LR, 
                           w0$Pr, 
                           NA, # no logLik output
                           w00$LR, 
                           w00$Pr,
                           NA, # intercept.phy not really interesting
                           z$B[2], z$B.pvalue[2], z$B.se[2])
      names(envir)[names(envir) == "PGLMM.slopephy.logLik0"] = "slopephy.logLik.minus.loglik0"
      names(envir)[names(envir) == "PGLMM.s2phy.logLik0"] = "s2phy.logLik.minus.loglik0"
    }
    write.xlsx(envir, file = "envir_ls.xlsx")
  }
  
  envir$lmer.X.Pr = round(envir$lmer.X.Pr, 5)
  envir$PGLMM.slope.star.Pr = round(envir$PGLMM.slope.star.Pr, 5)
  envir$PGLMM.slope.phy.Pr = round(envir$PGLMM.slope.phy.Pr, 5)
  envir$PGLMM.intercept.phy.Pr = round(envir$PGLMM.intercept.phy.Pr, 5)
  envir$PGLMM.B.Pr = round(envir$PGLMM.B.Pr, 5)
  select(envir, var, lmer.X.Pr, PGLMM.slope.star.Pr, PGLMM.slope.phy.Pr, PGLMM.B, PGLMM.B.Pr) %>% 
    rename(slope.var.lmer.Pr = lmer.X.Pr,
           slope.var.indep.pglmm.Pr = PGLMM.slope.star.Pr, 
           slope.var.phylo.pglmm.Pr = PGLMM.slope.phy.Pr, 
           slope.estimated.pglmm = PGLMM.B, 
           slope.estimated.pglmm.Pr = PGLMM.B.Pr)
}



### 本代码有GEB原文附件提供，由二傻统计略作注
### 该案例不需额外加载数据，安装必要的包后，直接复制粘贴即可运行模型和出图
set.seed(22)
library(tidyverse); library(ape); library(brms)
#Calulcate the correlation between plant height and range area
df_ <- data.frame(
  SPECIES = c('Abies alba', 'Abies borisii-regis', 'Abies cephalonica', 'Abies cilicica', 'Abies nordmanniana', 'Abies pinsapo', 'Acer campestre', 'Acer heldreichii', 'Acer monspessulanum', 'Acer platanoides', 'Acer pseudoplatanus', 'Aesculus hippocastanum', 'Alnus cordata', 'Alnus glutinosa', 'Alnus incana', 'Alnus alnobetula', 'Arbutus unedo', 'Betula pendula', 'Betula pubescens', 'Buxus balearica', 'Buxus sempervirens', 'Carpinus betulus', 'Carpinus orientalis', 'Castanea sativa', 'Cedrus libani', 'Celtis australis', 'Cornus mas', 'Cornus sanguinea', 'Corylus avellana', 'Euonymus europaeus', 'Fagus sylvatica', 'Frangula alnus', 'Fraxinus angustifolia', 'Fraxinus excelsior', 'Fraxinus ornus', 'Ilex aquifolium', 'Juglans regia', 'Juniperus communis', 'Juniperus drupacea', 'Juniperus excelsa', 'Juniperus foetidissima', 'Juniperus oxycedrus', 'Juniperus phoenicea', 'Juniperus thurifera', 'Larix decidua', 'Liquidambar orientalis', 'Olea europaea', 'Ostrya carpinifolia', 'Picea abies', 'Picea omorika', 'Picea orientalis', 'Pinus brutia', 'Pinus cembra', 'Pinus halepensis', 'Pinus heldreichii', 'Pinus mugo', 'Pinus nigra', 'Pinus peuce', 'Pinus pinaster', 'Pinus pinea', 'Pinus sylvestris', 'Pistacia atlantica', 'Pistacia lentiscus', 'Pistacia terebinthus', 'Platanus orientalis', 'Populus alba', 'Populus nigra', 'Populus tremula', 'Prunus avium', 'Prunus padus', 'Prunus spinosa', 'Quercus cerris', 'Quercus coccifera', 'Quercus faginea', 'Quercus frainetto', 'Quercus ilex', 'Quercus petraea', 'Quercus pubescens', 'Quercus pyrenaica', 'Quercus robur', 'Quercus suber', 'Quercus trojana', 'Salix alba', 'Salix caprea', 'Salix eleagnos', 'Sambucus nigra', 'Sorbus aria', 'Sorbus aucuparia', 'Sorbus domestica', 'Sorbus torminalis', 'Taxus baccata', 'Tetraclinis articulata', 'Tilia cordata', 'Tilia platyphyllos', 'Tilia tomentosa', 'Ulmus glabra', 'Ulmus laevis', 'Ulmus minor'),
  PLANT.HEIGHT = c(6.465,5.477,6.325,5.477,6.259,5.477,3.731,3.873,3.254,4.782,5.266,4.732,4.032,4.381,3.716,1.688,2.944,4.64,4.381,1.257,1.849,4.201,2.573,5.231,6.481,3.6,2.476,1.929,2.492,1.995,5.789,1.905,4.776,5.126,2.823,3.48,4.412,2.458,3.464,4.472,3.873,2.921,2.452,3.622,5.98,5.477,2.728,3.192,5.997,5.362,6.708,4.189,4.385,3.277,3.121,1.766,5.623,5,5.02,4.179,5.306,1.3,1.936,2,4.852,4.61,5.112,4.528,3.719,3.072,1.536,5.313,1.565,4.472,5.014,4.458,5.403,3.756,4.472,5.257,3.972,2.926,4.587,3.038,2.717,2.17,3.536,3.045,3.767,4.095,3.979,2.646,4.751,5.465,4.763,5.174,4.816,4.905),
  RANGE.AREA = c(58.609,19.206,19.396,27.555,38.932,14.251,82.144,43.094,68.2,86.193,71.636,20.345,26.835,91.61,100.216,101.997,64.611,112.959,104.989,33.151,49.425,77.608,54.926,48.84,33.512,64.68,71.724,84.369,88.933,82.225,75.277,97.49,81.475,87.753,61.964,69.003,84.791,111.544,29.588,66.075,43.14,72.279,62.965,41.492,41.693,20.847,69.018,57.519,80.593,11.533,33.434,43.89,38.024,50.975,20.381,46.85,52.579,15.987,48.655,40.825,102.269,78.063,62.349,69.644,54.872,100.962,91.811,114.642,78.551,108.543,89.709,63.645,61.917,45.999,54.088,62.754,77.281,70.533,48.224,87.163,49.997,41.121,97.16,100.37,60.889,85.877,85.195,110.686,66.53,72.476,66.336,45.731,90.209,72.794,52.633,88.311,84.59,85.689)
)
head(df_)
#N.B.: plant height (= unit in meters) is square-root transformed to improve normality; range area was Box-Cox transformed (see SI online for the complete dataset).
#Plant height was retrieved from FloraVeg.EU (https://floraveg.eu/taxon/); 
#Check the manuscript for additional details on how trait data and range attributes were compiled

#Phylogenetic tree:
#N.B.: the phylogenetic tree was obtained from the 鈥榩hylo.maker()鈥? function from the V.PhyloMaker R package (Jin & Qian, 2019; https://doi.org/10.1111/ecog.04434)
treetxt <- '((((((((Sambucus_nigra:102.692841,Ilex_aquifolium:102.692841)campanulids:4.04857,(Olea_europaea:17.739301,((Fraxinus_angustifolia:5.330423,Fraxinus_excelsior:5.330423):1.081153,Fraxinus_ornus:6.411576):11.327725):89.00211)mrcaott248ott320:5.599318,Arbutus_unedo:112.340729)mrcaott248ott650:2.225571,(Cornus_sanguinea:47.710427,Cornus_mas:47.710427):66.855873)mrcaott248ott27233:9.167937,((((((((Sorbus_domestica:2.706285,Sorbus_aucuparia:2.706285):5.893129,(Sorbus_torminalis:3.837911,Sorbus_aria:3.837911):4.761503):42.133854,((Prunus_spinosa:18.078495,Prunus_avium:18.078495):8.851155,Prunus_padus:26.92965):23.803618):48.228528,((Celtis_australis:79.170641,((Ulmus_minor:0.536514,Ulmus_glabra:0.536514):5.899536,Ulmus_laevis:6.43605):72.734592):6.322774,Frangula_alnus:85.493416):13.46838)Rosales.rn.d8s.tre:12.186212,(((((((Carpinus_betulus:8.226376,Carpinus_orientalis:8.226376):4.624557,Ostrya_carpinifolia:12.850933):19.198365,Corylus_avellana:32.049298):32.8096,(Betula_pendula:29.352098,Betula_pubescens:29.352098):35.5068):3.416551,(((Alnus_glutinosa:15.974943,Alnus_incana:15.974943):4.606542,Alnus_cordata:20.581485):18.979252,Alnus_alnobetula:39.560737):28.714712):22.373581,Juglans_regia:90.64903):7.230533,(((Quercus_robur:11.776698,Quercus_pubescens:11.776698,Quercus_ilex:11.776698,(((((Quercus_frainetto:3.966982,Quercus_faginea:3.966982):0.778652,Quercus_petraea:4.745634):0.424375,Quercus_pyrenaica:5.170009):3.507011,((Quercus_cerris:1.05467,Quercus_trojana:1.05467):0.52343,Quercus_suber:1.5781):7.09892):2.417941,Quercus_coccifera:11.094961):0.681738):0.532457,Castanea_sativa:12.309156):27.922252,Fagus_sylvatica:40.231407):57.648155):13.268445)mrcaott371ott2511:4.637557,((((Salix_caprea:0.308197,Salix_eleagnos:0.308197):7.034235,Salix_alba:7.342432):30.69623,((Populus_alba:1.574183,Populus_tremula:1.574183):1.537515,Populus_nigra:3.111698):34.926964):73.981924,Euonymus_europaeus:112.020586)mrcaott2ott1479:3.764979)mrcaott2ott371:2.793039,((Tilia_cordata:2.967151,(Tilia_platyphyllos:0.801068,Tilia_tomentosa:0.801068):2.166083):101.290908,(((((Acer_campestre:3.268696,Acer_platanoides:3.268696):9.197756,(Acer_heldreichii:6.709758,Acer_monspessulanum:6.709758):5.756694):1.629005,Acer_pseudoplatanus:14.095457):31.041674,Aesculus_hippocastanum:45.137131):34.78121,((Pistacia_atlantica:2.730903,Pistacia_terebinthus:2.730903):9.157292,Pistacia_lentiscus:11.888195):68.030146):24.339718)mrcaott96ott378:14.320545)mrcaott2ott96:3.826488,Liquidambar_orientalis:122.405092)mrcaott2ott2464:1.329145)Pentapetalae:4.81969,(Buxus_sempervirens:9.812122,Buxus_balearica:9.812122):118.741805)mrcaott2ott8379:1.770604,Platanus_orientalis:130.324531)mrcaott2ott969:194.725497,((((Juniperus_foetidissima:5.036479,((Juniperus_thurifera:1.964459,Juniperus_excelsa:1.964459):2.205496,Juniperus_phoenicea:4.169955):0.866525,((Juniperus_oxycedrus:3.777148,Juniperus_communis:3.777148):1.0901,Juniperus_drupacea:4.867248):0.169232):15.340157,Tetraclinis_articulata:20.376637):47.444059,Taxus_baccata:67.820696):59.002062,((((((((Pinus_mugo:4.876341,Pinus_sylvestris:4.876341):0.701058,Pinus_nigra:5.577399):3.316714,Pinus_heldreichii:8.894113):0.291606,((Pinus_pinea:4.318439,Pinus_pinaster:4.318439):2.865478,(Pinus_halepensis:2.227706,Pinus_brutia:2.227706):4.956211):2.001802):11.673147,(Pinus_cembra:5.599878,Pinus_peuce:5.599878):15.258988):14.157556,((Picea_abies:1.57471,Picea_omorika:1.57471):0.729596,Picea_orientalis:2.304306):32.712116):3.781157,Larix_decidua:38.797579):2.897999,((Abies_borisii-regis:4.548573,((((Abies_cilicica:0.692976,Abies_nordmanniana:0.692976):0.094111,Abies_cephalonica:0.787087):0.172634,Abies_pinsapo:0.959721):0.240809,Abies_alba:1.20053):3.348044):28.575012,Cedrus_libani:33.123586):8.571992):85.12718):198.22727)Spermatophyta;'
phy.tree <- read.tree(text = treetxt)

#Set names of phylogenetic tree and arrange the dataset accordingly
rownames(df) <- str_replace(df$SPECIES, ' ', '_') 


df <- read.xlsx("F:/博士研究计划2024/大样地数据/种间联结计算/每木调查_崂山.xlsx", 5)

df <- df[match(phy_tree$tip.label, df[,1]),]


#2. Run MR-PMMs ####

#2.1 Set the VCV matrix
C <- vcv.phylo(phy_tree, corr = T) 

#2.2 Run phylogenetic models
# rename traits:(dbh/h)
names(df)[names(df) %in% 'dbh'] <- 'dbh'
names(df)[names(df) %in% 'h'] <- 'h'
yi1 <- df[['iv']]
yi2 <- df[['levins']]

# fit MR-PMM:
#!!! N.B.: It may takes up to ~7 minutes to run !!!#
st <- Sys.time()

df$sp <- df[,1]

model.formula <- bf(mvbind(levins, iv) ~ (1|a|gr(sp, cov = C))) + set_rescor(TRUE)

MRPMM <- brm(model.formula,
             family = gaussian(),
             data   = df, 
             data2  = list(C = C),
             control = list(adapt_delta = 0.99),
             cores  = 4)
et <- Sys.time()-st
message(paste0('DONE. Total elapsed time: ',paste(round(et,2)),' ',attr(et,'units')))


#3. Get model posterior stats####
### 分别提取两因变量的截距和残差的相关性
#3.1 get 95% Credible intervals
median.ci95 <- as.data.frame(posterior_summary(MRPMM, probs = c(0.025, 0.975)))[c(5,8),]
median.ci95$dim <- ifelse(str_detect(rownames(median.ci95),'rescor'), 'ind', 'phy')

#3.2 get 50% Credible intervals
median.ci50 <- as.data.frame(posterior_summary(MRPMM, probs = c(0.25, 0.70)))[c(5,8),] %>% dplyr::select(-Estimate, -Est.Error)
median.ci50$dim <- ifelse(str_detect(rownames(median.ci50),'rescor'), 'ind', 'phy')

#3.3 inspect the results
res <- median.ci95 %>% left_join(median.ci50,'dim') %>%
  dplyr::rename(median=Estimate) %>%
  dplyr::select(median, dim, Q2.5, Q97.5, Q25, Q70)
res

### 对两种相关性（系统发育造成的截距相关性，独立于系统发育的残差相关性）及其置信区间作图
#3.4 visualize the results
ggplot() +
  geom_vline(xintercept = 0, col = 'black', lty = 'longdash') +
  geom_linerange(data = res, aes(y=dim, xmin=Q2.5, xmax=Q97.5, col=dim),linewidth =.7) +
  geom_pointrange(data = res, aes(y= dim,  x=median, xmin=Q25, xmax=Q70, col=dim), 
                  size=.8, shape=21, linewidth =2, fill='black') +
  theme_bw() + 
  labs(y='', col='Corr.', x='Estimated correlation') + 
  lims(x=c(-1,1))


### 检验数据与模型的匹配性
#4. Posterior Predictive Checks####
pp_check(MRPMM, resp= 'yi1', ndraws = 100) + labs(x='PLANT HEIGHT')
pp_check(MRPMM, resp= 'yi2', ndraws = 100) + labs(x='RANGE AREA')

##########示例数据####################
########################################
data <- data.frame(
  Trait = factor(rep(c("Plant lifespan", "Plant height", "Seed mass", "Stem specific density", 
                       "Xylem conduct diameter", "Xylem vulnerability"), each = 2), 
                 levels = c("Plant lifespan", "Plant height", "Seed mass", "Stem specific density", 
                            "Xylem conduct diameter", "Xylem vulnerability")),
  Type = rep(c("CTC", "IND"), times = 6),
  Estimate = c(-0.2, -0.1, -0.3, -0.15, 0.1, 0.05, -0.25, -0.2, 0.2, 0.15, -0.3, -0.25),  # Example data
  CI_low = c(-0.3, -0.2, -0.4, -0.25, 0, -0.05, -0.35, -0.3, 0.1, 0.05, -0.4, -0.35),    # Example data
  CI_high = c(-0.1, 0, -0.2, -0.05, 0.2, 0.15, -0.15, -0.1, 0.3, 0.25, -0.2, -0.15)     # Example data
)

# 准备阴影数据
shading_data <- data.frame(
  ymin = seq(0.5, 5.5, by = 2),
  ymax = seq(1.5, 6.5, by = 2),
  xmin = -1.1,
  xmax = 1.1
)

ggplot(data, aes(x = Estimate, y = Trait, color = Type)) +
  # 添加阴影
  geom_rect(data = shading_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE, fill = "gray90") +
  geom_point(position = position_dodge(width = 0.5), size = 3) +
  geom_errorbar(aes(xmin = CI_low, xmax = CI_high), width = 0.2, 
                position = position_dodge(width = 0.5)) +
  theme_minimal() +
  scale_color_manual(values = c("CTC" = "blue", "IND" = "red")) +
  labs(title = "Whole-plant / stem traits",
       x = "Correlation estimate",
       y = "",
       color = "Type:") +
  theme(axis.text.y = element_text(size = 10),
        strip.text = element_text(size = 12),
        legend.position = "bottom") +
  annotate("text", x = -1.1, y = 1:6, label = "⊖", size = 4) +
  annotate("text", x = -1.1, y = c(5, 6), label = "⊕", size = 4) +
  geom_vline(xintercept = 0, linetype = "dashed")





# A non-significant institution effect
fit1 <- coxph(Surv(time, status) ~ ph.ecog + age, data=lung,
              subset=(!is.na(inst)))
fit2 <- coxme(Surv(time, status) ~ ph.ecog + age + (1|inst), lung)
anova(fit1, fit2)
lung
# Shrinkage effects (equivalent to ridge regression)
temp <- with(lung, scale(cbind(age, wt.loss, meal.cal)))
rfit <- coxme(Surv(time, status) ~ ph.ecog + (temp | 1), data=lung)




# 竞争指数计算 ------------------------------------------------------------------
setwd("F:/博士研究计划2024/大样地数据/竞争指数计算")
library(ForestStatTool)
library(openxlsx)
library(parallel)

dat_compete <- read.xlsx("竞争计算.xlsx",2)

dat_compete <- Coord_Move(Data = dat_compete, Plot = dat_compete$Plot, X=dat_compete$X, Y=dat_compete$Y, Range_xy = c(18,18))


rawdata0 <- subset(dat_compete, !is.na(dat_compete$D))
newdata2 <- Coord_Remove(Data = rawdata0, Plot = rawdata0$Plot, X = rawdata0$X, Y = rawdata0$Y, D = rawdata0$D)
view(newdata2)
write.xlsx(newdata2, "newdata2.xlsx")
Coord_Remove(Data = newdata2, Plot = newdata2$Plot, X = newdata2$X, Y = newdata2$Y, D = newdata2$D)





b_1 <- CI.dist.base(Data = newdata2, Shape = "rectangle", Neighbor = "number",k = 4, Correct = "translation", Origin = c(0,0), Range_xy = c(100,100))
b_1
write.xlsx(b_1, "CI.dist.xlsx")
# 与距离无关
a_1 <- CI.nondist(Plot = dat_compete$Plot, Tag =dat_compete$Tag, D = dat_compete$D, SP = dat_compete$SP, S=1)

write.xlsx(a_1, "CI.nodist.xlsx")


# 箱线图NRI
library(ggplot2)
dat_ <- read.xlsx("CI.dist.xlsx", 11)
ggplot(dat_, aes(x = Group, y = Value)) +
  geom_jitter(aes(color = Group),position = position_jitter(0.4), alpha = 0.5, size = 2.5)+
  geom_boxplot(fill = "#ffe788") +
  theme_classic() 
# geom_hline(yintercept = 1.96, linetype = "dotted", col = "black")+
# geom_hline(yintercept = 0, linetype = "dotted", col = "black")+
# geom_hline(yintercept = -1.96, linetype = "dotted", col = "black")
# 箱线图NRIorNTI
wilcox.test(Value~Group, data = dat_)
ggsave("box_2.pdf", width = 5,height = 4,dpi = 300)

dat_ <- read.xlsx("CI.dist.xlsx", 11)
ggplot(dat_, aes(x = Group, y = Value)) +
  geom_jitter(aes(color = Group),position = position_jitter(0.4), alpha = 0.5, size = 2.5)+
  geom_boxplot(fill = "#ffe788") +
  theme_classic() 
  # geom_hline(yintercept = 1.96, linetype = "dotted", col = "black")+
  # geom_hline(yintercept = 0, linetype = "dotted", col = "black")+
  # geom_hline(yintercept = -1.96, linetype = "dotted", col = "black")


ggsave("box_1.pdf", width = 5,height = 4,dpi = 300)


dat <- read.xlsx("CI.dist.xlsx", 10)
shapiro.test(dat$CI15_intra_sq)
shapiro.test(dat$CI15_inter_sq)
wilcox.test(Value~Group, data = dat_)
wilcox.test(Value~Group, data = dat_)
# 空间点格局分析 -----------------------------------------------------------------
library(spatstat)
library(openxlsx)
library(ggplot2)
library(ForestStatTool)
library(devtools)
library(onpoint)
library(tidyverse)
library(dplyr)
library(idar)

library(spatgraphs)
library(spatialsegregation)

dat_point_ls <- read.xlsx("F:/博士研究计划2024/大样地数据/重要值计算/point.xlsx",6)

dat_point_ls$SP<-as.factor(dat_point_ls$SP)

dat_point_ls_ls <- dplyr::filter(dat_point_ls, spe%in% c("枹栎", "栓皮栎"))
dat_bl_hs <- dplyr::filter(dat_point_ls, spe%in% c("枹栎", "黑松"))
sp_ <- ppp(dat_point_ls$X,dat_point_ls$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_point_ls$SP)
sp_ls<-ppp(dat_point_ls_ls$X,dat_point_ls_ls$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_point_ls_ls$spe)

sp_bl <- dplyr::filter(dat_point_ls, spe%in% "枹栎")
sp_bl$spe<-as.factor(sp_bl$spe)
sp_bl <- ppp(sp_bl$X,sp_bl$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = sp_bl$spe)

# the following are equivalent



# # calculate L-function
# l_function <- Lest(sp_, correction = "Ripley")
# plot(l_function)
# 

# # center L-function to zero
# # center_l_function <- center_l_function(l_function)
# l_function_centered <- center_l_function(sp_, correction = "Ripley")
# 
# o_ring <- estimate_o_ring(sp_)
# 
# oring_envelope <- envelope(sp_, fun = estimate_o_ring, nsim = 199, verbose = FALSE)

# plot_quantums(oring_envelope, xlab = "r(m)",ylab = "g(r)")
# ggsave("oring.png", width = 14, height = 6.5, dpi = 500)

# heterogenous_pattern
null_model_hetero <- simulate_heterogenous_pattern(sp_, nsim = 199)

hetero <- envelope(sp_, fun = pcf, 
                   funargs = list(correction = "Ripley", divisor = "d"),
                   simulate = null_model_hetero, nsim = 199, 
                   verbose = FALSE)

plot_quantums(hetero, xlab = "r(m)",ylab = "g(r)" ,title = "Woody plant species")
ggsave("hetero0923.png", width = 14, height = 6.5, dpi = 500)
??envelope

oring_envelope <- envelope(sp_, fun = estimate_o_ring, nsim = 199, verbose = FALSE)

# marks(spab) <- ifelse(marks(spab) > 5, yes = "adult", no = "seedling")

null_model_antecedent <- simulate_antecedent_conditions(sp_ls, i = "枹栎", j = "栓皮栎",nsim = 199)

antecedent <- envelope(sp_ls, fun = pcf, 
                       funargs = list(correction = "Ripley", divisor = "d"),
                       simulate = null_model_antecedent, nsim = 199, 
                       verbose = FALSE)

plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus serrata VS Quercus variabilis")))
ggsave("antecedent.png", width = 14, height = 6.5, dpi = 500)




# 枹栎VS黑松
dat_bl_hs
sp_bl_hs <- ppp(dat_bl_hs$X,dat_bl_hs$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_bl_hs$spe)
null_model_antecedent <- simulate_antecedent_conditions(sp_bl_hs, i = "枹栎", j = "黑松",nsim = 199)

antecedent <- envelope(sp_ls, fun = pcf, 
                       funargs = list(correction = "Ripley", divisor = "d"),
                       simulate = null_model_antecedent, nsim = 199, 
                       verbose = FALSE)

plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus serrata VS Pinus thunbergii")))
ggsave("antecedent_hs_bl.png", width = 14, height = 6.5, dpi = 500)


# 栓皮栎VS黑松

dat_spl_hs <- dplyr::filter(dat_point_ls, spe%in% c("栓皮栎", "黑松"))
sp_spl_hs <- ppp(dat_spl_hs$X,dat_spl_hs$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_spl_hs$spe)
null_model_antecedent <- simulate_antecedent_conditions(sp_spl_hs, i = "栓皮栎", j = "黑松",nsim = 199)

antecedent <- envelope(sp_spl_hs, fun = pcf, 
                       funargs = list(correction = "Ripley", divisor = "d"),
                       simulate = null_model_antecedent, nsim = 199, 
                       verbose = FALSE)

plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus variabilis VS Pinus thunbergii")))
ggsave("antecedent_hs_spl.png", width = 14, height = 6.5, dpi = 500)

# 栓皮栎VS朴树
dat_spl_ps <- dplyr::filter(dat_point_ls, spe%in% c("栓皮栎", "朴树"))
sp_spl_ps <- ppp(dat_spl_ps$X,dat_spl_ps$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_spl_ps$spe)
null_model_antecedent <- simulate_antecedent_conditions(sp_spl_hs, i = "栓皮栎", j = "朴树",nsim = 199)

antecedent <- envelope(sp_spl_hs, fun = pcf, 
                       funargs = list(correction = "Ripley", divisor = "d"),
                       simulate = null_model_antecedent, nsim = 199, 
                       verbose = FALSE)

plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus variabilis VS Celtis sinensis")))
ggsave("antecedent_hs_spl.png", width = 14, height = 6.5, dpi = 500)


# 枹栎

dat_bl <- read.xlsx("F:/博士研究计划2024/大样地数据/重要值计算/point.xlsx",5)
dat_bl_adult <- dplyr::filter(dat_bl, dat_bl$dbh > 10)
sp_bl_adult <- ppp(dat_bl_adult$X,dat_bl_adult$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_bl_adult$spe)

# null_model_antecedent <- simulate_antecedent_conditions(sp_bl, i = "栓皮栎", j = "朴树",nsim = 199)
# 
# antecedent <- envelope(sp_spl_hs, fun = pcf, 
#                        funargs = list(correction = "Ripley", divisor = "d"),
#                        simulate = null_model_antecedent, nsim = 199, 
# #                        verbose = FALSE)
# 
# plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus variabilis VS Celtis sinensis")))
# ggsave("antecedent_hs_spl.png", width = 14, height = 6.5, dpi = 500)

null_model_hetero <- simulate_heterogenous_pattern(sp_bl_adult, nsim = 199)

hetero <- envelope(sp_bl_adult, fun = pcf, 
                   funargs = list(correction = "Ripley", divisor = "d"),
                   simulate = null_model_hetero, nsim = 199, 
                   verbose = FALSE)

plot_quantums(hetero, xlab = "r(m)",ylab = "g(r)" ,title = "All woody plants")
ggsave("hetero.png", width = 14, height = 6.5, dpi = 500)

# 栓皮栎
dat_spl <- dplyr::filter(dat_point_ls, spe%in% c("栓皮栎"))
sp_spl <- ppp(dat_spl$X,dat_spl$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_spl$spe)
sp_bl$marks
# null_model_antecedent <- simulate_antecedent_conditions(sp_bl, i = "栓皮栎", j = "朴树",nsim = 199)
# 
# antecedent <- envelope(sp_spl_hs, fun = pcf, 
#                        funargs = list(correction = "Ripley", divisor = "d"),
#                        simulate = null_model_antecedent, nsim = 199, 
# #                        verbose = FALSE)
# 
# plot_quantums(antecedent,xlab = "r(m)", ylab = "g(r)", title = expression(italic("Quercus variabilis VS Celtis sinensis")))
# ggsave("antecedent_hs_spl.png", width = 14, height = 6.5, dpi = 500)

null_model_hetero <- simulate_heterogenous_pattern(sp_spl, nsim = 199)

hetero <- envelope(sp_spl, fun = pcf, 
                   funargs = list(correction = "Ripley", divisor = "d"),
                   simulate = null_model_hetero, nsim = 199, 
                   verbose = FALSE)

plot_quantums(hetero, xlab = "r(m)",ylab = "g(r)" ,title = "All woody plants")
ggsave("hetero.png", width = 14, height = 6.5, dpi = 500)

plot(sp_spl$x, sp_spl$y)
plot(sp_bl$x, sp_bl$y)

fit0<-ppm(sp_spl ~ 1) # fit the stationary Poisson process
dclf.test(fit0, fun = Lest)# CSR model, 99 simulation

g11.enve.CSR=function(p){enve=envelope(p,pcf,r=seq(0,100,1))enve} 
sp.res1<-g11.enve.CSR(sp)
plot(sp.res1)

plot(simulate(fit0))


library(ForestStatTool)
ForestStatTool::Coord_Remove()
install.packages("forestcompetition")
dat <- read.xlsx("F:/博士研究计划2024/大样地数据/重要值计算/point.xlsx",3)
dat_compete <- filter(dat, dat$D > 5)
dat_compete_1 <- Coord_Move(Data = dat_compete, Plot = dat_compete$Plot, X=dat_compete$X, Y=dat_compete$Y, Range_xy = c(2,2))
write.xlsx(dat_compete_1, "dat_compete_1.xlsx")
dat_compete_2 <- Coord_Remove(Data = dat_compete, Plot = dat_compete$Plot, X=dat_compete$X, Y=dat_compete$Y, D = dat_compete$D)
write.xlsx(dat_compete_2, "dat_compete_2.xlsx")
b1 <- CI.dist.base(dat_compete_2, Shape = "rectangle", Neighbor = "number", k = 4, Correct = "translation", Origin = c(0,0), Range_xy = c(100,100))
write.xlsx(b1, "b1.xlsx")


# isarF(sp_, 20, target = "枹栎")
# mciF(sp_, 20, target = "枹栎")
# 
# bl <- segregationFun(sp_, fun = "isar", 20)
# 
# 
# 
# 
# 
# 
# C <- vcv.phylo(lsyd_tree, corr = T)
# C[upper.tri(C)] = NA
# C
# CC <- reshape2::melt(C, na.rm = TRUE)



# onpoint -----------------------------------------------------------------
library(onpoint)
library(spatstat)
data(spruces)
dat_point_ls <- read.xlsx("F:/博士研究计划2024/大样地数据/重要值计算/point.xlsx",5)

dat_point_ls$spe<-as.factor(dat_point_ls$spe)

# dat_point_ls_ls <- dplyr::filter(dat_point_ls, spe%in% c("枹栎", "栓皮栎"))
# dat_bl_hs <- dplyr::filter(dat_point_ls, spe%in% c("枹栎", "黑松"))
sp_ <- ppp(dat_point_ls$X,dat_point_ls$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_point_ls$spe)
# sp_ls<-ppp(dat_point_ls_ls$X,dat_point_ls_ls$Y,window=owin(xrange=c(0,100),yrange=c(0,100)), marks = dat_point_ls_ls$spe)

plot_quantums(antecedent, ylab = "g(r)")



null_model_antecedent <- simulate_antecedent_conditions(sp_, 
                                                        i = "seedling", j = "adult", nsim = 199)
antecedent <- envelope(sp_, fun = pcf, 
                       funargs = list(correction = "Ripley", divisor = "d"),
                       simulate = null_model_antecedent, nsim = 199, 
                       verbose = FALSE)




#

# 贝叶斯PGLMM ----------------------------------------------------------------
########################################
########################################
# Supporting Information:
# Appendix S1 to: Phylogenetic generalized linear mixed modeling presents 
# novel opportunities for eco-evolutionary synthesis
# Phylogenetic generalized linear mixed modeling presents novel opportunities for eco-evolutionary synthesis
########################################
########################################

# Some general notes:

# * Ignore 'divergent transitions' and 'rejecting initial value'
#   warnings for these examples. With empirical dataset, as with all
#   Bayesian models, you should perform model diagnostic checks (e.g.,
#   posterior predictive checks) before using the output of a
#   model. See Gelman et al. (2013; Bayesian Data Analysis) for more
#   details; if you are concerned in real datasets, consider
#   increasing the 'adapt_delta' argument, reducing the 'stepsize',
#   or increasing length/number of chains.
# * We use 'set.seed' to make the simulations below completely
#   reproducible. Since this makes all random numbers in your R
#   session the same, you may wish to restart R if doing further
#   simulation/analysis work.
# * To make the Bayesian posteriors reproducible, we also specify a
#   'seed' in each stan call. This is not something that would require
#   you to re-start your computer inbetween runs, but is probably
#   something that you wouldn't want to use in your own analyses.
# * You may still get slightly different (but qualitatively the same)
#   results if using an old verison of R (<3.5; released in July
#   2019). Later versions of R randomly permute data slightly
#   differently to address a known bug ("significant user-visible
#   changes in 3.6.0";
#   https://cran.r-project.org/doc/manuals/r-devel/NEWS.html).
# * To make these models fit in a reasonable timeframe on a laptop, we
#   work with fairly simply datasets (in particular for Question
#   2.2). By all means use this code as a starting point for more
#   complex simulations and analyses. Recently, we have been exploring
#   using variation Bayes methods (rstan::vb) for bigger problems,
#   and while this is (at the time of writing) an experimental
#   approach in rstan, we have found it to work well.
# * Running "options(mc.cores=XXX)", where XXX is the number of
#   threads your computer can run (2 or 4 on most laptops) will make
#   the model fitting *much* faster
# * The code below will install packages if they are missing from your
#   computer, which means it will take longer to run the first time.

########################################
########################################
# R packages and functions #############
########################################
########################################
# Load / install packages
if(!require("geiger")){
  install.packages("geiger", dependencies=TRUE); library(geiger)}
if(!require("mvtnorm")){
  install.packages("mvtnorm", dependencies=TRUE); library(mvtnorm)}
if(!require("pez")){
  install.packages("pez", dependencies=TRUE); library(pez)}
if(!require("lme4")){
  install.packages("lme4", dependencies=TRUE); library(lme4)}
if(!require("rstan")){
  install.packages("rstan", dependencies=TRUE); library(rstan)}
library(geiger)
library(mvtnorm)
library(pez)
library(lme4)
library(rstan)
########################################
# Community simulation functions #######
########################################
# Inverse logit function
inv.logit <- function(x)
  exp(x) / (exp(x)+1)
# Environmental filtering
make.env.comm <- function(trait, env, occupancy, tree, rescale=TRUE){
  data <- expand.grid(species=names(trait), site=names(env), stringsAsFactors=FALSE)
  data$trait <- trait[data$species]
  data$env <- env[data$site]
  data$occupancy <- occupancy[data$species]
  data$num_species <- as.numeric(gsub("s","",data$species))
  data$num_sites <- as.numeric(factor(data$site))
  data$prob <- data$trait * data$env + data$occupancy
  if(rescale)
    data$prob <- inv.logit(data$prob)
  data$presence <- rbinom(nrow(data), 1, data$prob)
  return(data)
}
# Scramble competition
make.comp.comm <- function(comp.trait, n.spp, env, tree, scale=2){
  # Make overdispersed VCV from trait
  vcv <- as.matrix(dist(comp.trait))
  vcv <- max(vcv) - vcv
  vcv <- diag(length(env)) %x% solve(vcv)
  
  # Make overdispersed VCV from phylogeny
  # vcv <- solve(vcv(tree))/scale
  # vcv <- diag(length(env)) %x% vcv
  
  # Draw presences/absences
  pres.abs <- matrix(inv.logit(rmvnorm(1, sigma=vcv)), nrow=length(env), byrow=TRUE)
  for(i in seq_len(nrow(pres.abs)))
    pres.abs[i,] <- ifelse(pres.abs[i,]>=sort(pres.abs[i,], decreasing=TRUE)[sample(n.spp,1)], 1, 0)
  comm <- matrix(as.numeric(t(pres.abs)), nrow=length(env), ncol=length(comp.trait), dimnames=list(names(env), tree$tip.label), byrow=TRUE)
  
  # Format data and return
  data <- comparative.comm(tree, comm, traits=data.frame(trait=comp.trait), env=data.frame(env=env))
  data <- as.data.frame(data)
  data$site <- factor(data$site)
  data$species <- factor(data$species)
  data$num_species <- as.numeric(gsub("s","",data$species))
  data$num_sites <- as.numeric(factor(data$site))
  return(data)
}
# Association simulation
make.assoc <- function(base.data, base.tree, base.trait, partner.trait){
  # Make prob. of association matrix
  assoc <- 1 - plnorm(abs(outer(base.trait, partner.trait, FUN=`-`)))
  
  # Add in partner interactions on basis of probabilities
  partner.comm <- matrix(0, nrow=nrow(base.data), ncol=length(partner.trait))
  for(i in seq_len(nrow(base.data))){
    if(base.data$presence[i]==1)
      partner.comm[i,] <- as.numeric(runif(length(partner.trait)) <= assoc[base.data$num_species[i],])
  }
  
  # Reformat data and return
  data <- base.data[rep(seq_len(nrow(base.data)), length(partner.trait)),]
  names(data)[c(6,8,9)] <- c("num_base_species","base_prob","base_presence")
  data$num_partner_species <- rep(seq_along(partner.trait), each=nrow(partner.comm))
  data$partner_presence <- as.numeric(partner.comm)
  return(data)
}

########################################
########################################
# Simulate data used throughout ########
########################################
########################################
# Basics
set.seed(1234)
n.spp <- 50
tree <- sim.bdtree(n=n.spp, seed=1234)
n.sites <- 50
site.names <- apply(as.matrix(expand.grid(letters,letters,stringsAsFactors=FALSE)), 1, paste, collapse="")[seq_len(n.sites)]
env.sites <- setNames(rnorm(n.sites), site.names)

# Overall OCCupancies (intercepts) and traits (slopes) under three
#   different models of evolution: CONStrained (phy.signal), LABile
#   (phylogenetically random), and COMpetition (repulsion)
cons.occ <- sim.char(tree, .05, root=0)[,,1]
lab.occ <- setNames(sample(cons.occ), names(cons.occ))
comp.occ <- setNames(rmvnorm(1, rep(0,n.spp), solve(vcv(tree))/10)[1,], tree$tip.label)
cons.trait <- sim.char(tree, .05, root=0)[,,1]
lab.trait <- setNames(sample(cons.trait), names(cons.trait))
comp.trait <- setNames(rmvnorm(1, rep(0,n.spp), solve(vcv(tree))/10)[1,], tree$tip.label)

# Make datasets based around traits to exemplify biologically plausible scenarios
cons.env <- make.env.comm(cons.trait, env.sites, lab.occ, tree)
# - Constrained evolution of species' environmental responses,
#   environmental assembly on the basis of those traits, and
#   phylogenetically neutral overall occupancies
ghost.comp <- make.env.comm(comp.trait, env.sites, lab.occ, tree)
# - Evolution to minimise competition (the 'ghost of competition
#   past') with environmental assembly on the basis of those traits,
#   and phylogenetically neutral overall occupancies
s.tree <- drop.tip(tree, 16:50)
cons.comp <- make.comp.comm(cons.occ[1:15], c(2,5), env.sites[1:15], s.tree)
# - Constrained evolution of species' environmental responses, and
#   competitive (scramble) assembly on the basis of those same traits
# NOTE: we're now working with fewer sites, because these competition
#   models are quite computationally intensive and we would like these
#   examples to run in a reasonable timeframe


# Make association dataset
base.n.spp <- 10; partner.n.spp <- 10
base.tree <- sim.bdtree(n=base.n.spp, seed=2345); partner.tree <- sim.bdtree(n=partner.n.spp, seed=3456)
assoc.n.sites <- 10; assoc.env.sites <- setNames(rnorm(assoc.n.sites), letters[seq_len(assoc.n.sites)])
base.env.trait <- sim.char(base.tree, .05, root=0)[,,1]; base.occ <- setNames(sample(sim.char(base.tree, .05, root=0)[,,1]), names(base.env.trait))
base.assoc.trait <- sim.char(base.tree, .05, root=0)[,,1]; partner.assoc.trait <- sim.char(partner.tree, .05, root=0)[,,1]
# - First, make a "dummy" dataset with a 'base' species and a
#   'partner' associated species, but have them be independent. This
#   is to provide a test of what happens when two species are independent
base.env <- make.env.comm(base.env.trait, assoc.env.sites, base.occ, base.tree)
# - Now simulate a "base" species and a "partner" species, where the
#   partner is dependent on the base species (i.e., is found in
#   association with that base species)
assoc.data <- make.assoc(base.env, base.tree, base.assoc.trait, partner.assoc.trait)


########################################
########################################
# Q1: Environmental responses  #########
########################################
########################################

########################################
# Model 1.1: presence ~ trait * env ####
one.one <- glmer(presence ~ env * trait + (1|species), family=binomial, data=cons.env)
summary(one.one)
# - Significant interaction between environment and trait. Note the random
#   effect accounts for species-level pseudoreplication; if desired
#   you can add in site-level terms too but for this toy example we do
#   not think they are necessary

########################################
# Model 1.2: presence ~ trait * env ####
#                 WITHOUT phylogeny ####
# - First, let's fit it using frequentist methods
one.two <- glmer(presence ~ (env|species), family=binomial, data=cons.env)
#NOTE: Sometimes, you may get a boundary fit warning/error. If so, try
#      the fit again (perhaps with a different random number seed or
#      optimiser). This is a difficult model to numerically fit.
summary(one.two)

# Extract the latent trait (the coefficients of the interaction term),
# then re-order to match the simulated trait and plot
latent.traits <- setNames(ranef(one.two)$species$env, rownames(ranef(one.two)$species))
latent.traits <- latent.traits[names(cons.trait)]
plot(latent.traits ~ cons.trait)
abline(0,1, col="red", lwd=2)
cor.test(latent.traits, cons.trait)
# ... Not a perfect reconstruction (remember, there is error added in
#     the simulation) but not bad either

# - Second, let's do this using Bayesian methods. Go through the
#   comments below *slowly* because everything else builds on this

one.two.stan.code <- "
// Tell rstan what our data are
data{
int Ntotal;              // Number of rows in our dataset
int Nspp;                // Number of species
int presence[Ntotal];    // The response variable (presence/absence)
real env[Ntotal];        // An explanatory variable (environment)
int spp[Ntotal];         // An explanatory variable (species)
}
// The parameters/coefficients we want to estimate and report back
parameters{
vector[Nspp] spp_int;  // Each species' overall occupancy
vector[Nspp] spp_slp;      // Each species' environmental response
}
// Some calculations rstan is going to perform internally to help model-fitting
transformed parameters{
vector[Ntotal] predictions;   // Make a variable to hold our model predictions
// Loop over all our input data and specify our model, which is:
// a species' intercept (overall occupancy) +  env response x the environment
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i];
}
// Fit our model to our predictions
model{
// Species' occupancies and responses are drawn from uninformative priors
spp_int ~ normal(0, 0.5); 
spp_slp ~ normal(0, 0.5);
// The model itself: our presences are drawn from our predictions
presence ~ bernoulli_logit(predictions);
}
"
one.two.stan.model <- stan(model_code=one.two.stan.code,
                           data=list(Ntotal=nrow(cons.env), Nspp=n.spp,
                                     presence=cons.env$presence,
                                     env=cons.env$env, spp=cons.env$num_species),
                           iter=1000, chains=4, seed=123456)

# Now let's match the output up (in much the same way as before, only
#   now using 'extract' to get the information from the posterior
#   distribution of our rstan model fit)
latent.traits <- apply(extract(one.two.stan.model)$spp_slp, 2, median)
names(latent.traits) <- paste0("s", seq_along(latent.traits))
latent.traits <- latent.traits[names(cons.trait)]
plot(latent.traits, cons.trait)
cor.test(latent.traits, cons.trait)
# ... and, again, we've reconstructed the trait

########################################
# Model 1.3: presence ~ trait * env ####
#                 WITH phylogeny    ####
# Below is a modified Bayesian version of the model, with the changes
#   from model 1.3 highlighted
# NOTE: We supply a function to apply a lambda transformation
#       called 'lambda_vcv' inside the stan code; feel free to
#       focus only on the comments or skip over that function definition
#       on your first read-through
one.three.stan.code <- "
// Function to transform a phylogenetic VCV according to Pagel's Lambda
// - and multiply through by overall sigma
functions {
matrix lambda_vcv(matrix vcv, real lambda, real sigma){
matrix[rows(vcv),cols(vcv)] local_vcv; // Make a local copy of the VCV to modify
local_vcv = vcv * lambda;              // Lambda transforms are just a multiplier
for(i in 1:rows(local_vcv))            // ... but we do have to loop over the matrix
local_vcv[i,i] = vcv[i,i];           // ... and make the diagonal the same as it was before
return(local_vcv * sigma);             // Return the transformed matrix x sigma (overall variation)
}
}
data {
int Ntotal;
int Nspp;
int presence[Ntotal];
real env[Ntotal];
int spp[Ntotal];
//    
matrix[Nspp,Nspp]Vphy;     // Give the phylogeny as data
}
parameters{
vector[Nspp] spp_int;
vector[Nspp] spp_slp;
// Lambda transforms for the intercepts and slopes
real<lower=0> lam_int;     // (priors --> cannot be negative) 
real<lower=0> lam_slp;
// Coefficients for the NON-phylogenetically-derived variance in model terms
real<lower=0> null_int;    // (priors --> cannot be negative)
real<lower=0> null_slp;
// Coefficients for the mean intercepts/slopes
real mean_int;
real mean_slp;
}
transformed parameters{
vector[Ntotal] predictions;
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i];
}
model{
// Specify priors
mean_int ~ normal(0,0.5);
mean_slp ~ normal(0,0.5);
lam_int ~ normal(0,0.5);
lam_slp ~ normal(0,0.5);
null_int ~ normal(0,0.5);
null_slp ~ normal(0,0.5);

// Now we draw our species coefficients, incorporating the lambda-transformed phylogeny
spp_int ~ multi_normal(rep_vector(mean_int,Nspp), lambda_vcv(Vphy,lam_int,null_int));
spp_slp ~ multi_normal(rep_vector(mean_slp,Nspp), lambda_vcv(Vphy,lam_slp,null_slp)); 
//
presence ~ bernoulli_logit(predictions);
}
"
one.three.stan.model <- stan(model_code=one.three.stan.code,
                             data=list(Ntotal=nrow(cons.env), Nspp=n.spp,
                                       presence=cons.env$presence,
                                       env=cons.env$env, spp=cons.env$num_species, Vphy=vcv(tree)),
                             iter=1000, chains=2, seed=123456)

# Summarise the model fit, focusing on the parameters we care about
summary(one.three.stan.model, pars=c("lam_int","lam_slp"))$summary
# ... A lambda value of 0 suggests a lack of evolutionary constraint on present-day
#     ecology, while a value of 1 suggests strong constraint. Here we have strong 
#     phylogenetic signal on the phylogenetically conserved environmental 
#     responses (lam_slp)

# We can even get histograms of the posterior distributions of Lambda
# if we like
hist(unlist(extract(one.three.stan.model, "lam_slp")))

# (An aside: the latent trait we estimated in Q1.2 does have
# phylogenetic signal too, albeit when estimated using a different
# method)
fitContinuous(tree, latent.traits, model="lambda")
# (Look at the 'lambda' parameter in the output above)

# We can also get estimates for the latent traits and occupancies (we
# could have done this before, but there was no need). We're going to
# use a wrapper function for this now, to make this code easier to
# read and generalise.
.grab.latent.trait <- function(model, name, known.trait){
  latent <- apply(extract(model)[[name]], 2, mean)
  names(latent) <- paste0("s", seq_along(latent))
  return(latent[names(known.trait)])
}
cor.test(.grab.latent.trait(one.three.stan.model, "spp_slp", cons.trait), cons.trait)
cor.test(.grab.latent.trait(one.three.stan.model, "spp_int", lab.occ), lab.occ)
# ... great, nice positive correlations (i.e., we accurately estimated the traits)

########################################
########################################
# Q2: Competition  #####################
########################################
########################################

########################################
# Model 2.1: the ghost of competition past
two.one.stan.code <- "
data{
int Ntotal;
int Nspp;
int presence[Ntotal];
real env[Ntotal];
int spp[Ntotal];
matrix[Nspp,Nspp] Vphy;
matrix[Nspp,Nspp] inv_Vphy; // The inverse of the phylogeny
//                             (used to measure repulsive evolution)
}
parameters{
vector[Nspp] spp_int;
vector[Nspp] spp_slp;
real<lower=0> phy_int;
real<lower=0> phy_slp;
real<lower=0> null_int;
real<lower=0> null_slp;
// Coefficients for the inverse-phylogenetically-derived variance in model terms
real<lower=0> invphy_int;       // (with priors specified too) 
real<lower=0> invphy_slp;       //
real mean_int;
real mean_slp;
}
transformed parameters{
vector[Ntotal] predictions;
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i];
}
model{
// Specify priors
mean_int ~ normal(0,0.5);
mean_slp ~ normal(0,0.5);
phy_int ~ normal(0,0.5);
phy_slp ~ normal(0,0.5);
invphy_int ~ normal(0,0.5);
invphy_slp ~ normal(0,0.5);
null_int ~ normal(0,0.5);
null_slp ~ normal(0,0.5);

// Now we draw our species coefficients to measure the importance of phylogeny and its inverse
// - notice how we're drawing from a covariance matrix with null (non-phylogenetic), phylogenetic, and the inverse of the phylogenetic, components
spp_int ~ multi_normal(rep_vector(mean_int,Nspp), diag_matrix(rep_vector(null_int,Nspp)) + phy_int*Vphy + invphy_int*inv_Vphy); 
spp_slp ~ multi_normal(rep_vector(mean_slp,Nspp), diag_matrix(rep_vector(null_slp,Nspp)) + phy_slp*Vphy + invphy_slp*inv_Vphy);
presence ~ bernoulli_logit(predictions);
}
"
two.one.stan.model <- stan(model_code=two.one.stan.code,
                           data=list(Ntotal=nrow(ghost.comp),
                                     Nspp=max(ghost.comp$num_species),
                                     presence=ghost.comp$presence,
                                     env=ghost.comp$env, spp=ghost.comp$num_species,
                                     Vphy=vcv(tree), inv_Vphy=solve(vcv(tree))),
                           iter=1000, chains=2, seed=123456)

# Let's compare the relative support for species' traits evolving
# under Brownian motion (phy_slp), under repulsion (invphy_slp -
# the inverse of the phylogenetic covariance matrix), or null (no
# phylogenetic signal at all) using a wrapper function
.rel.support <- function(model, ...){
  terms <- list(...)
  vars <- setNames(numeric(length(terms)), names(terms))
  for(i in seq_along(terms))
    vars[i] <- median(extract(model)[[terms[[i]]]])
  return(setNames(vars / sum(vars), unlist(terms)))
}
.rel.support(two.one.stan.model, "phy_slp", "invphy_slp", "null_slp")
# ... repulsion is most likely (these numbers sum to one: ~65% support
#     for repulsion, ~30% for null, but definitely not Lambda (~2%))

.rel.support(two.one.stan.model, "phy_int", "invphy_int", "null_int")
# ... contrasting with the intercepts, which were labile (~75%
#     support), and show no evidence of repulsion.

# NOTE: These models are limited by the simulated data; 15 species in an assemblage 
#       is not a lot to fit with a phylogeny, and results in an unusual distribution 
#       of entries for the inverse of the phylogenetic covariance matrix (and
#       subsequent model-fitting and simulation issues). We've simulated only
#       15 species here because these models are computationally intensive, 
#       and we want to ensure these examples can be run on a laptop. Playing around 
#       with larger datasets will result in better performance, and also 
#       (if you like) you can change to simulating using phylogeny (and not a 
#       two-step process of phylogeny --> trait --> assembly) to increase 
#       statistical power.


# Model 2.2: competition in the present
two.two.stan.code <- "
data{
int Ntotal;
int presence[Ntotal];
matrix[Ntotal,Ntotal] InvPhySite;   // Matrix of phylogeny x sites to detect competition
}
parameters{
real<lower=0.0001, upper=10> null_sites;       // Null variance of sites
real<lower=0.0001, upper=10> inv_sites;        // Phylogenetic repulsion of sites
vector[Ntotal] predictions;
}
model{
predictions ~ multi_normal(rep_vector(0,Ntotal), diag_matrix(rep_vector(fabs(null_sites),Ntotal)) + fabs(inv_sites)*InvPhySite);
presence ~ bernoulli_logit(predictions);
}
"

two.two.stan.model <- stan(model_code=two.two.stan.code,
                           data=list(Ntotal=nrow(cons.comp),
                                     Nspp=max(cons.comp$num_species),
                                     Nsite=max(cons.comp$num_site),
                                     presence=sample(cons.comp$presence),
                                     spp=cons.comp$num_species,
                                     InvPhySite=((solve(vcv(s.tree))/10) %x% diag(max(cons.comp$num_site)))),
                           iter=2000, chains=4, seed=123456)
# NOTE: We are working with a matrix with as many rows and columns as
#       we have data entries - this is a computationally difficult
#       model to fit (be patient and careful with its output)
# Let's see the support for phylogenetic repulsion in these data
inv.sites <- extract(two.two.stan.model)[["inv_sites"]]
sum(inv.sites > .2) / length(inv.sites)
#...If this value is greater than 0, we have support for overdispersion.
#   In this case, ~95% of the posterior density is above 0.2, which we would 
#   argue is reasonably strong support for phylogenetic overdispersion. Note
#   that there are many alternative ways of testing for this (this
#   approach is similar to that in Ives and Helmus 2011), and we have
#   used an unhelpful (uniform) prior here.

# NOTE: Again, we are using very little data here (only 15 species) because
#       these models are computationally intensive. Playing around with larger 
#       datasets will result in better performance, and simulating using phylogeny 
#       (and not a two-step process of phylogeny --> trait --> assembly) to increase 
#       statistical power.


########################################
########################################
# Q3: Species' associations  ###########
########################################
########################################

########################################
# Model 3.1: presence ~ env * phylogeny + other.species
# (Contrast this model 1.3)
three.one.stan.code <- "
functions {
matrix lambda_vcv(matrix vcv, real lambda, real sigma){
matrix[rows(vcv),cols(vcv)] local_vcv; // Make a local copy of the VCV to modify
local_vcv = vcv * lambda;              // Lambda transforms are just a multiplier
for(i in 1:rows(local_vcv))            // ... but we do have to loop over the matrix
local_vcv[i,i] = vcv[i,i];           // ... and make the diagonal the same as it was before
return(local_vcv * sigma);             // Return the transformed matrix x sigma (overall variation)
}
}
data{
int Ntotal;
int Nspp;
int presence[Ntotal];
real env[Ntotal];
int spp[Ntotal];
matrix[Nspp,Nspp]Vphy;
//
int Npartner;                              // # of partner species (c.f., Nspp)
int partner[Ntotal];                       // Give partner presence as data
matrix[Nspp,Nspp]VpartnerPhy;              // Partner phylogeny
}
parameters{
vector[Nspp] spp_int;
vector[Nspp] spp_slp;
vector[Nspp] partner_int;                    // Allow for partner presence/absence to predict base species' presence/absence
real<lower=0> lam_int;
real<lower=0> lam_slp;
real<lower=0> null_int;
real<lower=0> null_slp;
real mean_int;
real mean_slp;
// Coefficients for the partner-derived variance in model terms
real<lower=0> lam_partner_int;     // (phylogenetic component)
real<lower=0> null_partner_int;    // (non-phylogenetic)
real mean_partner_int;             // (and the overall mean here)
}
transformed parameters{
vector[Ntotal] predictions;
for (i in 1:Ntotal)
predictions[i] = spp_int[spp[i]] + spp_slp[spp[i]]*env[i] + partner_int[partner[i]]*partner[i];
// Notice extra term added at the end for the partner (analogous to spp_intercepts)
}
model{
// Priors
lam_int ~ normal(0, 0.5);
lam_slp ~ normal(0, 0.5);
null_int ~ normal(0, 0.5);
null_slp ~ normal(0, 0.5);
lam_partner_int ~ normal(0, 0.5);
null_partner_int ~ normal(0, 0.5);
mean_int ~ normal(0, 0.5);
mean_slp ~ normal(0, 0.5);
mean_partner_int ~ normal(0, 0.5);
spp_int ~ multi_normal(rep_vector(mean_int,Nspp), lambda_vcv(Vphy,lam_int,null_int));
spp_slp ~ multi_normal(rep_vector(mean_slp,Nspp), lambda_vcv(Vphy,lam_slp,null_slp));
partner_int ~ multi_normal(rep_vector(mean_partner_int,Npartner), lambda_vcv(VpartnerPhy,lam_partner_int,null_partner_int));
// Notice the intercepts term is the same as the above, but now we address the partner species
presence ~ bernoulli_logit(predictions);
}
"

# First fit: as a null, let's model the 'base' species that doesn't
# depend on the partner species, to see what finding no association
# looks like
three.one.stan.model.base <- stan(model_code=three.one.stan.code,
                                  data=list(Ntotal=nrow(assoc.data), Nspp=base.n.spp,
                                            presence=assoc.data$base_presence,
                                            env=assoc.data$env, spp=assoc.data$num_base_species, Vphy=vcv(base.tree),
                                            Npartner=partner.n.spp, partner=assoc.data$num_partner_species,
                                            VpartnerPhy=vcv(partner.tree)),
                                  iter=1000, chains=2, seed=123456)

# Test the relative importance of the partner on the base species' presences
.rel.support(three.one.stan.model.base, "null_int", "lam_int", "null_partner_int", "lam_partner_int")
# ... very low values for the partner terms - this species' presences
#     isn't determined by the other species

########################################
# Now let's repeat the model fitting, but for the *other* species
#   ('partner') that is driven by the presence of the 'base' species
three.one.stan.model.partner <- stan(model_code=three.one.stan.code,
                                     data=list(Ntotal=nrow(assoc.data), Nspp=partner.n.spp,
                                               presence=assoc.data$partner_presence,
                                               env=assoc.data$env, spp=assoc.data$num_partner_species, Vphy=vcv(partner.tree),
                                               Npartner=base.n.spp, partner=assoc.data$num_base_species,
                                               VpartnerPhy=vcv(base.tree)),                
                                     iter=1000, chains=2, seed=123456)
# -  notice that, above, we've swapped "partner" and "base" around,
#    otherwise the data input is the same

.rel.support(three.one.stan.model.partner, "null_int", "lam_int", "null_partner_int", "lam_partner_int")
# ... now we can see the impact of partner on species' distributions,
#     explaining >50% of the variance (70% + 25%)

summary(three.one.stan.model.partner, par="lam_partner_int")$summary
# ... and we can see that it has some phylogenetic signal (as it
#     should, because we simulated it to) - it is remarkable we can
#     detect this given we have only 10 species (in each clade) in
#     this dataset!




# 草本与木本联系 -----------------------------------------------------------------
# #########草本与木本植物的联系
setwd("F:/博士研究计划2024/崂山英文想法/林下植被")
library(devtools)
library(plantlist)
library(openxlsx)
library(V.PhyloMaker2)
library(V.PhyloMaker)
library(ape)
library("devtools")
library(lme4)
dat_ls_herb <- read.xlsx("herb.xlsx",9)
herb_latin__ <- CTPL(dat_ls_herb$spe)
write.xlsx(herb_latin__, "herb_latin__.xlsx")


# 草本植物系统发育树构建
library(spaa)
library(reshape2)
library(openxlsx)
library(phyr)
library(V.PhyloMaker2)
library(V.PhyloMaker)
library(plantlist)
library(patchwork)
library(picante)
library(phytools)
library(rr2)
library(phylolm)
# 构建系统发育树
dat_tree <- read.xlsx("herb_latin__.xlsx", 2)
tree_TPL <- phylo.maker(dat_tree, tree=GBOTB.extended.TPL, nodes=nodes.info.1.TPL, scenarios="S3")
write.tree(tree_TPL$scenario.3, "herb.newick")
herb_tree <- read.tree("herb.newick")
woody_tree <- read.tree("tree.newick")
plot(herb_tree)
# 计算系统发育指数

herb_abun <- read.xlsx("herb.xlsx",9) 
herb_abun[is.na(herb_abun)] <- 0
dat_phyr <- tidyr::gather(herb_abun, key = "sp", value = "freq", -site)
write.xlsx(dat_phyr, "herb_dat_phyr.xlsx")
dat_herb_phyr <- read.xlsx("herb_dat_phyr.xlsx",2)

# binomial
# ############系统发育相近的物种是否会出现在相同的样地中--------
rr2::R2_pred(mod_2)
herb_mod_1 <- pglmm(pa_herb ~ 1 + (1|herb__) + (1|site) + (1|herb__@site),
             data = dat_herb_phyr,family = "binomial", cov_ranef = list(herb = herb_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(herb_mod_1)
# rr2::R2_pred(mod_2)
pchisq(2*(herb_mod_1$logLik - herb_mod_2_r$logLik), df = 1, lower.tail = F)/2


herb_mod_2 <- pglmm(pa_herb ~ 1 + (1|herb__) + (1|site) + (1|herb__@site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(herb = herb_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065),repulsion = TRUE)

summary(herb_mod_2)

pchisq(2*(herb_mod_2$logLik - herb_mod_2_r$logLik), df = 1, lower.tail = F)/2

herb_mod_2_r <- pglmm(pa_herb ~ 1 + (1|herb__) + (1|site) , 
                 data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree))
summary(herb_mod_2_r)

rr2::R2_pred(herb_mod_1, herb_mod_2_r)






Model <- pglmm(pa ~ 1 + (1|sp__) + (1|sp__@herb) , 
                      data = dat_herb_phyr, family = "binomial", cov_ranef = list(sp = woody_tree, herb = herb_tree))
summary(Model)
# 分解R2
rr2::R2_pred(Model)


#############草本的分布是否收到了木本植物的影响------------------------
Model_ <- pglmm(pa_herb ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|herb__) + (1|herb__@sp__) + (1|site) , 
               data = dat_herb_phyr, family = "binomial", cov_ranef = list(sp = woody_tree, herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model_)
rr2::R2_pred(Model_)
rr2::R2_lik(Model_)
pchisq(2*(Model_ $logLik -Model__$logLik), df = 1, lower.tail = F)/2 # 0.2974474

Model_1 <- pglmm(pa_herb ~ pH + EC + TN + TP + SOC + U + M + W + (1|herb__) + (1|herb__@sp) + (1|site), 
                data = dat_herb_phyr, family = "binomial", cov_ranef = list(sp = woody_tree, herb = herb_tree),repulsion = TRUE,s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model_1)
rr2::R2_pred(Model_1)
rr2::R2_lik(Model_1)
rr2::R2_lik(Model_, Model_1)
pchisq(2*(Model_$logLik -Model_1$logLik), df = 1, lower.tail = F)/2 # 0.3652864

Model__ <- pglmm(pa_herb ~ pH + EC + TN + TP + SOC+ U + M + W+ (1|herb__)+ (1|site), 
                 data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model__)
rr2::R2_pred(Model__)
rr2::R2_lik(Model__)
rr2::R2_lik(Model_, Model__)


Model__1<- pglmm(pa_herb ~ pH + EC + TN + TP + SOC+ U + M + W+ (1|herb__) + (1|sp)+ (1|site), 
                data = dat_herb_phyr, family = "binomial", cov_ranef = list(sp = woody_tree, herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model__1)
rr2::R2_pred(Model__1)
rr2::R2_lik(Model__1)
rr2::R2_lik(Model__1, Model__)


#############相近的草本是否分布相同的样地------------------------######“没有显著的系统发育结构”
Model_4 <- pglmm(pa_herb ~ 1 +(1|herb__) + (1|site) + (1|herb__@site) , 
                data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model_4)
rr2::R2_pred(Model_4)
rr2::R2_lik(Model_4)
pchisq(2*(Model_4$logLik -Model___$logLik), df = 1, lower.tail = F)/2 # 0.5

Model_44 <- pglmm(pa_herb ~ 1 + (1|herb__) + (1|site) + (1|herb__@sp) , 
                 data = dat_herb_phyr, family = "binomial", cov_ranef = list(sp = woody_tree, herb = herb_tree),repulsion = TRUE,s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model_44)
rr2::R2_pred(Model_44)
rr2::R2_lik(Model_44)
pchisq(2*(Model_44$logLik -Model___$logLik), df = 1, lower.tail = F)/2 # 0.2641454

Model___ <- pglmm(pa_herb ~ 1 + (1|herb__)+ (1|site), 
                 data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Model__)
rr2::R2_pred(Model__)
rr2::R2_lik(Model__)

#######################################探讨系统发育树、林分空间结构、土壤理化性质对草本植物分布的影响---#################
# 分解R2 Rlik2
###############################系统发育树的影响######################
Mod_ <- pglmm(pa_herb ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|herb__)+ (1|site), 
                data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))#全模型
summary(Mod_)
rr2::R2_pred(Mod_)
rr2::R2_lik(Mod_)  # 0.1161513

Mod_1 <- pglmm(pa_herb ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|herb)+ (1|site), 
               data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
rr2::R2_lik(Mod_1) # 0.1140276
rr2::R2_lik(Mod_, Mod_1) # 0.002397038
pchisq(2*(Mod_$logLik - Mod_1$logLik), df = 1, lower.tail = F)/2#显著性检验# 0.01818315

###############################林分空间结构的影响######################
Mod_11 <- pglmm(pa_herb ~ 1 + pH + EC + TN + TP + SOC  + (1|herb__)+ (1|site), 
               data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Mod_11)
rr2::R2_lik(Mod_11) # 0.1141113
rr2::R2_lik(Mod_, Mod_11) # 0.002302777
pchisq(2*(Mod_$logLik - Mod_11$logLik), df = 1, lower.tail = F)/2 ##0.02012383
###############################土壤理化性质的影响######################
Mod_111 <- pglmm(pa_herb ~ 1 +  U + M + W + (1|herb__)+ (1|site), 
                data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))

rr2::R2_lik(Mod_111) #0.1131988
rr2::R2_lik(Mod_, Mod_111) #0.003329367
pchisq(2*(Mod_$logLik - Mod_111$logLik), df = 1, lower.tail = F)/2 # 0.00681197

#保存数据
save(Mod_, Mod_1, Mod_11, Mod_111, dat_herb_phyr, woody_tree, herb_tree, file = "herb.RData")
load("herb.RData")


######################################################################################################################################################

# 系统发育树R2
rr2::R2_lik(Mod_, Mod_1) # 0.00249593
pchisq(2*(Mod_$logLik - Mod_1$logLik), df = 1, lower.tail = F)/2 # 0.01635621
R2(Mod_, Mod_1)


Mod__1 <- pglmm(pa_herb ~ 1 + pH + EC + TN + TP + SOC + (1|herb__)+ (1|site), 
                 data = dat_herb_phyr, family = "binomial",cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Mod__1)

#林分空间结构R2 
rr2::R2_lik(Mod_, Mod__1) # 0.002762319
pchisq(2*(Mod_$logLik - Mod_1$logLik), df = 1, lower.tail = F)/2 # 0.01635621

Mod___1 <- pglmm(pa_herb ~  U + M + W+ (1|herb__)+ (1|site), 
                 data = dat_herb_phyr, family = "binomial", cov_ranef = list(herb = herb_tree),s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Mod___1)
# 土壤理化性质R2
rr2::R2_lik(Mod_, Mod___1) #0.003363021
pchisq(2*(Mod_$logLik - Mod___1$logLik), df = 1, lower.tail = F)/2 #0.006578574


Mod____1 <- pglmm(pa_herb ~ pH + EC + U + M + W + (1|herb)+ (1|site), 
               data = dat_herb_phyr, family = "binomial",s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(Mod___)

rr2::R2_lik(Mod_, Mod___)
pchisq(2*(Mod__$logLik - Mod___$logLik), df = 1, lower.tail = F)/2
rr2::R2_pred(Mod___)

M <- glm(pa_herb ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|site), 
              data = dat_herb_phyr, family = "binomial")
R2_lik(Mod_, M)
summary(M)

##########################探讨系统发育树、林分空间结构、土壤理化性质对木本植物分布的影响---#################
#木本植物
#系发育树对木本植物分布的R2###
mod_2 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp__) + (1|site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_2)
R2_lik(mod_2) #0.157834


mod_2_ <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp) + (1|site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_2_)
R2_lik(mod_2_) #0.1606897
R2_lik(mod_2, mod_2_) #-0.003402465

pchisq(2*(mod_2$logLik - mod_2_$logLik), df = 1, lower.tail = F)/2



mod__2 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp) + (1|site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree))

summary(mod__2)
R2_lik(mod__2) #0.1605479
pchisq(2*(mod_2$logLik - mod__2$logLik), df = 1, lower.tail = F)/2 
R2_lik(mod_2,mod__2)

#############################

#检验系统竞争排除的作用
mod_3 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|sp__@site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_3)
R2_lik(mod_3) #  0.03647256

mod__3 <- glm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W, 
                data = dat_herb_phyr, family = binomial)

mod__3 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W, 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod__3)
R2_lik(mod__3) # 0.02069404
pchisq(2*(mod_3$logLik - mod__3$logLik), df = 1, lower.tail = F)/2 ##0.5 
pchisq(2*(-915.535 - -930.3568), df = 1, lower.tail = F)/2 #2.596189e-08 

R2_lik(mod_3,mod__3) # 0.01611194
mod_3$logLik
logLik(mod__3)

#检验土壤的作用
mod__4 <- glm(pa ~ 1 + U + M + W , 
                data = dat_herb_phyr, family = binomial)
R2_lik(mod__4)# 0.01415083
logLik(mod__4)
pchisq(2*(-915.535 - -936.4334), df = 1, lower.tail = F)/2 #5.063451e-11
R2_lik(mod_3,mod__4) #0.02264213


mod__4 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W + (1|site) + (1|sp__@site), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod__4) #0.04198417

# 空间结构
mod__5 <- glm(pa ~ 1 + pH + EC + TN + TP + SOC , 
                data = dat_herb_phyr, family = binomial)
summary(mod__5)
R2_lik(mod_3, mod__5) #0.03359654
logLik(mod__5)
pchisq(2*(-915.535 - -946.7186), df = 1, lower.tail = F)/2 #1.425156e-15
R2_lik(mod__5) #0.002976003


??glm
# GLM
library(glmm.hp)
library(lme4)
mod5 <- glm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site) , 
               data = dat_herb_phyr, family = binomial)
AIC(mod5)
summary(mod5)
glmm.hp(mod5)


write.xlsx(dat_phyr, "dat_phyr.xlsx")

dat <- data.frame(dat_phyr[,c(1:3,12:13)], scale(dat_phyr[,c(4:11)]))



mod5_gs <- glm(freq ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site) , 
            data = dat_herb_phyr, family = gaussian)
AIC(mod5, mod5_gs)

mod_5_gs_ <- pglmm(Y ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site)  + (1|sp__@site), 
                  data = dat, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_)#0.02703203

summary(mod_5_gs_)

mod_5_gs_limit_ <- pglmm(Y ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site), 
                        data = dat, family = "gaussian",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_limit) # 0.02950642










# 高斯分布------
mod_5_gs <- pglmm(freq ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site)  + (1|sp__@site), 
               data = dat_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs)#0.06275911
#####################limit_gaussian
mod_5_gs_limit <- pglmm(freq ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site), 
                  data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_limit) # 0.02950642
R2_lik(mod_5_gs,mod_5_gs_limit) # 0.03426369
pchisq(2*(mod_5_gs$logLik -mod_5_gs_limit$logLik), df = 1, lower.tail = F)/2 #7.515291e-16
#####################pH_gaussian
mod_5_gs_ph <- pglmm(freq ~ 1 + EC + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                        data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_ph) # 0.05973478
R2_lik(mod_5_gs,mod_5_gs_ph) # 0.003216466
pchisq(2*(mod_5_gs$logLik -mod_5_gs_ph$logLik), df = 1, lower.tail = F)/2 #  0.007659037

#####################EC_gaussian
mod_5_gs_EC <- pglmm(freq ~ 1 + pH + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                  data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_EC)# 0.05874028
R2_lik(mod_5_gs,mod_5_gs_EC) #0.00426963
pchisq(2*(mod_5_gs$logLik -mod_5_gs_EC$logLik), df = 1, lower.tail = F)/2 # 0.002599682
#####################Tn_gaussian
mod_5_gs_TN <- pglmm(freq ~ 1 + pH + EC + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                     data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_TN)#0.003229133
R2_lik(mod_5_gs,mod_5_gs_TN) #0.007558831
pchisq(2*(mod_5_gs$logLik -mod_5_gs_TN$logLik), df = 1, lower.tail = F)/2 # 
#####################SOC_gaussian
mod_5_gs_SOC <- pglmm(freq ~ 1 + pH + EC + TP + TN + U + M+ W +(1|site) + (1|sp__@site), 
                     data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_SOC)# 0.05943984
R2_lik(mod_5_gs,mod_5_gs_SOC) #0.003529039
pchisq(2*(mod_5_gs$logLik -mod_5_gs_SOC$logLik), df = 1, lower.tail = F)/2 # 0.005541726
#####################TP_gaussian
mod_5_gs_TP <- pglmm(freq ~ 1 + pH + EC + SOC + TN + U + M+ W +(1|site) + (1|sp__@site), 
                      data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_TP)#0.05832123
R2_lik(mod_5_gs,mod_5_gs_TP) #0.004712735
pchisq(2*(mod_5_gs$logLik -mod_5_gs_TP$logLik), df = 1, lower.tail = F)/2 # 0.001661481
#####################W_gaussian
mod_5_gs_W <- pglmm(freq ~ 1 + pH + EC + SOC + TN + U + M+ TP +(1|site) + (1|sp__@site), 
                     data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_W)#
R2_lik(mod_5_gs,mod_5_gs_W) #
pchisq(2*(mod_5_gs$logLik -mod_5_gs_TP$logLik), df = 1, lower.tail = F)/2 # 
#####################m_gaussian
mod_5_gs_M <- pglmm(freq ~ 1 + pH + EC + SOC + TN + U + W+ TP +(1|site) + (1|sp__@site), 
                    data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_M)#0.06007386
R2_lik(mod_5_gs,mod_5_gs_M) #0.002856881
pchisq(2*(mod_5_gs$logLik -mod_5_gs_M$logLik), df = 1, lower.tail = F)/2 # 0.01115613
#####################U_gaussian
mod_5_gs_U <- pglmm(freq ~ 1 + pH + EC + SOC + TN + M + W+ TP +(1|site) + (1|sp__@site), 
                    data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_5_gs_U)#0.05956943
R2_lik(mod_5_gs,mod_5_gs_U) #0.003391724
pchisq(2*(mod_5_gs$logLik -mod_5_gs_U$logLik), df = 1, lower.tail = F)/2 # 0.006386007

# 绘图表达-------------------------（高斯分布）
df_gs <- data.frame(matrix(ncol = 3, nrow = 9))
name_col <- c("Limit", "pH", "EC", "TN", "TP", "SOC", "U", "M", "W")
name_row <- c("variable", "value", "p-value")
rownames(df_gs) <- name_col
colnames(df_gs) <- name_row
df_gs[,1] <- name_col
# save(mod__3, mod__4, mod_5, mod_5_,woody_tree, herb_tree, file = "woody.RData")
partialR2 <- c(R2_lik(mod_5_gs,mod_5_gs_limit), R2_lik(mod_5_gs,mod_5_gs_ph),R2_lik(mod_5_gs,mod_5_gs_EC),R2_lik(mod_5_gs,mod_5_gs_TN),R2_lik(mod_5_gs,mod_5_gs_TP),R2_lik(mod_5_gs,mod_5_gs_SOC),R2_lik(mod_5_gs,mod_5_gs_U),R2_lik(mod_5_gs,mod_5_gs_M),R2_lik(mod_5_gs,mod_5_gs_W))
df_gs[,2] <- partialR2 
pvalue <- c(pchisq(2*(mod_5_gs$logLik -mod_5_gs_limit$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_ph$logLik), df = 1, lower.tail = F)/2 ,pchisq(2*(mod_5_gs$logLik -mod_5_gs_EC$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_TN$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_TP$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_SOC$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_U$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_M$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5_gs$logLik -mod_5_gs_TP$logLik), df = 1, lower.tail = F)/2)
df_gs[,3] <- pvalue

label <- c("***", "***", "***", "***", "***", "***", "***", "***", "***")
library(ggplot2)
ggplot(data = df_gs) + 
  geom_histogram(aes(x = reorder(variable,-value), y = value, fill = variable), stat = "identity") +
  # geom_text(aes(x = reorder(variable,-value), y = value,label = label)) +
  ylab("Partial R^2lik") + xlab("Variables")+scale_y_continuous(expand = c(0,0),limits = c(0,0.04)) +
  theme(#legend.position = "none",
    axis.text.x=element_text(angle=0,hjust=0.5),
    text = element_text(size=10,face="bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=10,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))


ggsave("limit_soil_stand_gs.pdf", width= 7 , height= 4)
ggsave("limit_soil_stand_gs.png", width= 7 , height= 4, dpi = 500)


model_logit <- glm(pa ~ 1 + pH + TN + TP + SOC + U + M+ W, data = dat_herb_phyr,family = binomial(link = "logit"))

summary(model_logit)
vif(model_logit)


vif(glm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M+ W, data = dat_herb_phyr))
####################林分空间结构对木本植物分布的影响#####

mod_55 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_55)

# (1|sp__@site)
mod_555 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
R2_lik(mod_555)
R2_lik(mod_55, mod_555)
pchisq(2*(mod_55$logLik - mod_555$logLik), df = 1, lower.tail = F)/2 # 1.361654e-08


mod_5 <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_5)
R2_lik(mod_5) # 0.04198417
mod_5$AIC
R2_lik(mod__3,mod_5) #0.01678644
pchisq(2*(mod__3$logLik - mod_5$logLik), df = 1, lower.tail = F)/2 # 1.361654e-08

#####################limit
mod_5_ <- pglmm(pa ~ 1 +pH + EC + TN + TP + SOC + U + M + W + (1|site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_5_)
R2_lik(mod_5_) #0.02525773
R2_lik(mod_5,mod_5_) #0.01715986
pchisq(2*(mod_5$logLik - mod_5_$logLik), df = 1, lower.tail = F)/2 # 9.527503e-09

#####################空间结构
mod_5__ <- pglmm(pa ~ 1 +pH + EC + TN + TP + SOC + (1|site), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_5__)
R2_lik(mod_5__) #0.01868753
R2_lik(mod_5,mod_5__) #0.02374028
pchisq(2*(mod_5$logLik - mod_5__$logLik), df = 1, lower.tail = F)/2 # 1.773933e-11

#####################土壤
mod_5___ <- pglmm(pa ~ 1 + U + M + W + (1|site), 
                 data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_5___)
R2_lik(mod_5___) #0.02040705
R2_lik(mod_5,mod_5___) #0.02202661
pchisq(2*(mod_5$logLik - mod_5___$logLik), df = 1, lower.tail = F)/2 # 9.114158e-11

# 逐个因子
# pH
mod_5_ph <- pglmm(pa ~ 1 +  EC + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
               data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_5_ph)
R2_lik(mod_5_ph) # 0.04129525

R2_lik(mod_5,mod_5_ph) #0.00071859
pchisq(2*(mod_5$logLik - mod_5_ph$logLik), df = 1, lower.tail = F)/2 # 0.1260256

# EC
mod_5_EC <- pglmm(pa ~ 1 +  pH + TN + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                  data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_EC) # 0.04112237

R2_lik(mod_5,mod_5_EC) #0.0008987506
pchisq(2*(mod_5$logLik - mod_5_EC$logLik), df = 1, lower.tail = F)/2 # 0.1000971

# TN
mod_5_TN <- pglmm(pa ~ 1 +  pH + EC + TP + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                  data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_TN) # 0.04106696

R2_lik(mod_5,mod_5_TN) #0.0009564905
pchisq(2*(mod_5$logLik - mod_5_TN$logLik), df = 1, lower.tail = F)/2 # 0.09316305
# TP
mod_5_TP <- pglmm(pa ~ 1 +  pH + TN + EC + SOC + U + M+ W +(1|site) + (1|sp__@site), 
                  data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_TP) # 0.04128123

R2_lik(mod_5,mod_5_TP) #0.0007332042
pchisq(2*(mod_5$logLik - mod_5_TP$logLik), df = 1, lower.tail = F)/2 # 0.1236415
# SOC
mod_5_SOC <- pglmm(pa ~ 1 +  pH + TN + EC + TP + U + M+ W +(1|site) + (1|sp__@site), 
                  data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_SOC) # 0.03986809

R2_lik(mod_5,mod_5_SOC) #0.002203939
pchisq(2*(mod_5$logLik - mod_5_SOC$logLik), df = 1, lower.tail = F)/2 # 0.0223937

# U
mod_5_U <- pglmm(pa ~ 1 +  pH + TN + EC + TP + SOC + M+ W +(1|site) + (1|sp__@site), 
                   data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_U) # 0.04105642

R2_lik(mod_5,mod_5_U) #0.0009674662
pchisq(2*(mod_5$logLik - mod_5_U$logLik), df = 1, lower.tail = F)/2 # 0.09190912

# M
mod_5_M <- pglmm(pa ~ 1 +  pH + TN + EC + TP + SOC + U+ W +(1|site) + (1|sp__@site), 
                 data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_M) # 0.04009618

R2_lik(mod_5,mod_5_M) #0.001966843
pchisq(2*(mod_5$logLik - mod_5_M$logLik), df = 1, lower.tail = F)/2 # 0.02901129

# W
mod_5_W <- pglmm(pa ~ 1 +  pH + TN + EC + TP + SOC + U+ M +(1|site) + (1|sp__@site), 
                 data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

R2_lik(mod_5_W) #  0.03785441

R2_lik(mod_5,mod_5_W) #0.004292236
pchisq(2*(mod_5$logLik - mod_5_W$logLik), df = 1, lower.tail = F)/2 # 0.002540768

save(mod_5, mod_5_ph, mod_5_EC, mod_5_TN, mod_5_TP, mod_5_SOC, mod_5_U, mod_5_M,mod_5_W,df,dat_herb_phyr,phy_tree, file = "stand_soil_ls_woody.RData")
load("stand_soil_ls_woody.RData")
# 绘图表达-------------------------（二项分布）
df <- data.frame(matrix(ncol = 3, nrow = 9))
name_col <- c("Limit", "pH", "EC", "TN", "TP", "SOC", "U", "M", "W")
name_row <- c("variable", "value", "p-value")
rownames(df) <- name_col
colnames(df) <- name_row
df[,1] <- name_col
save(mod__3, mod__4, mod_5, mod_5_,woody_tree, herb_tree, file = "woody.RData")
partialR2 <- c(R2_lik(mod_5,mod_5_), R2_lik(mod_5,mod_5_ph),R2_lik(mod_5,mod_5_EC),R2_lik(mod_5,mod_5_TN),R2_lik(mod_5,mod_5_TP),R2_lik(mod_5,mod_5_SOC),R2_lik(mod_5,mod_5_U),R2_lik(mod_5,mod_5_M),R2_lik(mod_5,mod_5_W))
df[,2] <- partialR2 
pvalue <- c(pchisq(2*(mod_5$logLik - mod_5_$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_ph$logLik), df = 1, lower.tail = F)/2 ,pchisq(2*(mod_5$logLik - mod_5_EC$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_TN$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_TP$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_SOC$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_U$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_M$logLik), df = 1, lower.tail = F)/2,pchisq(2*(mod_5$logLik - mod_5_W$logLik), df = 1, lower.tail = F)/2)
df[,3] <- pvalue

label <- c("***", "", "", "", "", "*", "", "*", "**")

ggplot(data = df) + 
  geom_histogram(aes(x = reorder(variable,-value), y = value, fill = variable), stat = "identity") +
  geom_text(aes(x = reorder(variable,-value), y = value,label = label)) +
  ylab("Partial R^2") + xlab("Variables")+scale_y_continuous(expand = c(0,0),limits = c(0,0.018)) +
  theme(#legend.position = "none",
    axis.text.x=element_text(angle=0,hjust=0.5),
    text = element_text(size=10,face="bold"),
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    strip.text.y = element_blank(), 
    strip.text = element_text(size=10,face="bold"),
    strip.background = element_blank(),
    axis.line = element_line(color = 'black'))
  
ggsave("limit_soil_stand.pdf", width= 7 , height= 4)
ggsave("limit_soil_stand.png", width= 7 , height= 4, dpi = 500)
getwd()
write.xlsx(df, "limit_soil_stand.xlsx")
# library(dplyr)
# library(ForestStatTool)
# dat_newdata2 <- read.xlsx("newdata2.xlsx",5)
# fix(dat_newdata2)
# dat_newdata2_stru <- Unit_Struc(Data = dat_newdata2, k = 4,
#            Shape = "rectangle", Correct = "single", Origin = c(0,0), Range_xy =
#              c(20,20),  Buf_dist = 2, Bind = TRUE)
# 
# dat_newdata1 <- Coord_Move(Data = dat_newdata2, Plot = dat_newdata2$Plot, X=dat_newdata2$X,
#                        Y=dat_newdata2$Y, Range_xy = c(2,2))
# 
# newdata2 <- Coord_Remove(Data = dat_newdata2, Plot = dat_newdata2$Plot, X =
#                            dat_newdata2$X, Y = dat_newdata2$Y, D = dat_newdata2$D)
# write.xlsx(newdata2, "newdata22.xlsx")
# dat_newdata2_stru <- Unit_Struc(Data = dat_newdata2, k = 4,
#                                 Shape = "rectangle", Correct = "single", Origin = c(0,0), Range_xy =
#                                   c(20,20),  Buf_dist = 2, Bind = TRUE)
# 
# write.xlsx(dat_newdata2_stru$tree_value, "tree_value.xlsx")
# write.xlsx(dat_newdata2_stru$stand_value, "stand_value.xlsx")





#############################
summary(mod__2)
R2_lik(mod__2)
pchisq(2*(mod_2$logLik - mod__2$logLik), df = 1, lower.tail = F)/2
R2_lik(mod_2,mod__2)


mod_2$AIC
mod__2$AIC

VarCorr(mod_2)


mod_2222 <- glmer(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp), 
               data = dat_herb_phyr, family = "binomial")
summary(mod_2222)
partialR2(mod_2222)

##############################
mod_22 <- pglmm(freq ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp__) + (1|site), 
               data = dat_herb_phyr, family = "poisson",cov_ranef = list(sp = woody_tree))
summary(mod_22)
R2_lik(mod_22)

mod__22 <- pglmm(freq ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp), 
                data = dat_herb_phyr, family = "gaussian",cov_ranef = list(sp = woody_tree))

summary(mod__22)
R2_lik(mod__22)
pchisq(2*(mod_22$logLik - mod__22$logLik), df = 1, lower.tail = F)/2
R2_lik(mod_22,mod__22)




hist(dat_herb_phyr$U)
hist(dat_herb_phyr$W)
hist(dat_herb_phyr$M)




trait_try <- read.table("F:/博士研究计划2024/崂山英文想法/trait/36676.txt", header = FALSE, sep = "\t")
install.packages("data.table")
library(data.table)
dt <- fread("F:/博士研究计划2024/崂山英文想法/trait/36676.txt", nrows = 10000, header = TRUE)







mod_2_ <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp) + (1|sp__@site), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_2_)
R2_lik(mod_2_)


mod_22_ <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp), 
                data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree),  s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_22_)
R2_lik(mod_22_)
R2_lik(mod_2, mod_2_)

pchisq(2*(mod_2_$logLik - mod_22_$logLik), df = 1, lower.tail = F)/2


mod_222_ <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp), 
                 data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree), repulsion = TRUE, s2.init = c(0.98, 0.001, 0.001, 0.0065))
summary(mod_222_)
R2_lik(mod_222_)


mod_22222_ <- pglmm(pa ~ 1 + pH + EC + TN + TP + SOC + U + M + W +(1|sp__), 
                   data = dat_herb_phyr, family = "binomial",cov_ranef = list(sp = woody_tree),  s2.init = c(0.98, 0.001, 0.001, 0.0065))

summary(mod_22222_)
R2_lik(mod_22222_)

















