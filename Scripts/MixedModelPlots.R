# libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(afex)
library(emmeans)

source("/Volumes/Data/zoocon/Rev_Hippocampgoal/Scripts/LMM_plot.R")

data = read.csv("/Volumes/Data/zoocon/Rev_Hippocampgoal/Data/long_cue_period_RSA_results.csv")

## Same Sequence Same Context 
## Figure 2B
target = c("sseq scon no sub", "dseq scon no sub", "sseq dcon no sub", "dseq dcon no sub")

# subset data
long_data_AOV= data %>% 
  filter(contrast %in% target) %>%
  mutate(same_context = (contrast %in% c("sseq scon no sub", "dseq scon no sub")),
         same_sequence = (contrast %in% c("sseq scon no sub","sseq dcon no sub")))

# change to factors 
long_data_AOV$same_context = as.factor(long_data_AOV$same_context)
long_data_AOV$same_sequence = as.factor(long_data_AOV$same_sequence)

## HPC
HIPP_MERGE_BL.intercept = mixed(PS ~ same_sequence*same_context + (1|subject), data = filter(long_data_AOV, ROI == "HIPP_MERGE_BL"), method = "LRT")
HIPP_MERGE_BL.intercept.F = mixed(PS ~ same_sequence*same_context + (1|subject), data = filter(long_data_AOV, ROI == "HIPP_MERGE_BL"))
HIPP_MERGE_BL.inercept.lmer = lmer(PS ~ same_sequence*same_context + (1|subject), data = filter(long_data_AOV, ROI == "HIPP_MERGE_BL"))
print(HIPP_MERGE_BL.intercept)
p=mixed_model_plot(HIPP_MERGE_BL.inercept.lmer,'seq*con')
ggsave('/Volumes/Data/zoocon/Rev_Hippocampgoal/Figures/Fig2B_HPC_BL_EMM_seq_con_bars.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

## For reviewers
tmp_emm = emmeans(HIPP_MERGE_BL.inercept.lmer, ~same_sequence*same_context)
tmp_df = as.data.frame(tmp_emm)
tmp_indiv_df = HIPP_MERGE_BL.inercept.lmer@frame
write.csv(tmp_df, "~/Downloads/seq_con_means.csv")
write.csv(tmp_indiv_df, "~/Downloads/seq_con_indivs.csv")

## Converging vs. Diverging 
target = c("dstart send diff context no sub", "sstart dend diff context no sub", 
           "dstart dend diff context no sub","dstart send same context no sub", 
           "sstart dend same context no sub", "dstart dend same context no sub",
          "sseq scon no sub", "sseq dcon no sub")
long_data_con_div_AOV = data %>% 
  filter(contrast %in% target)
  
## add overlap and context
long_data_con_div_AOV = long_data_con_div_AOV %>% 
  mutate(same_context = !(contrast %in% c("dstart send diff context no sub", "sstart dend diff context no sub", "dstart dend diff context no sub", 
                                          "sseq dcon no sub")),
         converge = (contrast %in% c("dstart send same context no sub",  "dstart send diff context no sub")),
         diverge = (contrast %in% c("sstart dend same context no sub",  "sstart dend diff context no sub")),
         no_overlap = (contrast %in% c("dstart dend same context no sub",  "dstart dend diff context no sub")),
         full_overlap = (contrast %in% c("sseq scon no sub", "sseq dcon no sub"))) 

# make a new row and combine levels
long_data_con_div_AOV$overlap = NaN
long_data_con_div_AOV[long_data_con_div_AOV$converge,"overlap"] = "converge"
long_data_con_div_AOV[long_data_con_div_AOV$diverge,"overlap"] = "diverge"
long_data_con_div_AOV[long_data_con_div_AOV$no_overlap,"overlap"] = "no_overlap" 
long_data_con_div_AOV[long_data_con_div_AOV$full_overlap,"overlap"] = "full_overlap" 

long_data_con_div_AOV$overlap <- as.factor(long_data_con_div_AOV$overlap)
long_data_con_div_AOV$same_context <- as.factor(long_data_con_div_AOV$same_context)
long_data_con_div_AOV = droplevels(long_data_con_div_AOV) # drop un used levels

## HPC 
con_div_intercept.HIPP_MERGE_BL = mixed(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "HIPP_MERGE_BL"), method = "LRT")
con_div_intercept.HIPP_MERGE_BL.F = mixed(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "HIPP_MERGE_BL"))
con_div_intercept.HIPP_MERGE_BL.lmer = lmer(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "HIPP_MERGE_BL"))
print(con_div_intercept.HIPP_MERGE_BL)

## For reviewers
tmp_emm = emmeans(con_div_intercept.HIPP_MERGE_BL.lmer, ~overlap*same_context)
tmp_df = as.data.frame(tmp_emm)
tmp_indiv_df = con_div_intercept.HIPP_MERGE_BL.lmer@frame
write.csv(tmp_df, "~/Downloads/win_con_means.csv")
write.csv(tmp_indiv_df, "~/Downloads/win_con_indivs.csv")

context_effect = emmeans(con_div_intercept.HIPP_MERGE_BL.lmer, pairwise ~ same_context | overlap, lmer.df = "asymptotic", adjust = "none" )
tmp_df = as.data.frame(context_effect$emmeans)


# within context 
p = main_effect_plot(con_div_intercept.HIPP_MERGE_BL.lmer, "overlap") 
print(p)
ggsave('/Volumes/Data/zoocon/Rev_Hippocampgoal/Figures/Fig2C_HPC_BL_EMM_win_context.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

# between context
p = main_effect_plot(con_div_intercept.HIPP_MERGE_BL.lmer, "overlap*context" )
print(p)
ggsave('/Volumes/Data/zoocon/Rev_Hippocampgoal/Figures/Fig2D_HPC_BL_EMM_btwn_context.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

## pairwise comparisons
temp=emmeans(con_div_intercept.HIPP_MERGE_BL.lmer,pairwise  ~ same_context*overlap,lmer.df = "asymptotic", adjust = "none") # stats 
# contrast(temp, "consec", simple = "each", combine = T, adjust = "none", interaction = T)
a = contrast(temp,simple = "each",interaction = "pairwise",adjust = "none") # 13 is Converge > Diverge, Converge > Diverge  22 is converge >full overlap 24,  282 
a$contrasts$`simple contrasts for contrast`[13]


# hard coded 
# overlap - estimate = 0.00906, se = 0.00362, z = 2.5, p = 0.012
# diverge - estimate = 0.00679, se = 0.00362, z = 1.878, p = 0.0604
# full overlap - estimate = 0.00939, se = 0.00362, z = 2.597, p = 0.0094
# no overlap - estimate = 0.00155, se = 0.00362, z = 0.430, p = 0.6675

# same_context = TRUE:
# overlap_pairwise          estimate      SE  df z.ratio p.value
# converge - diverge         0.00791 0.00362 Inf  2.188  0.0286 
# converge - full_overlap   -0.00470 0.00362 Inf -1.299  0.1941 
# converge - no_overlap      0.00532 0.00362 Inf  1.471  0.1414 
# diverge - full_overlap    -0.01261 0.00362 Inf -3.487  0.0005 
# diverge - no_overlap      -0.00260 0.00362 Inf -0.718  0.4730 
# full_overlap - no_overlap  0.01002 0.00362 Inf  2.769  0.0056 

# for the between context interaction plot
# [287] "(FALSE diverge - TRUE diverge) - (FALSE no_overlap - TRUE no_overlap)" / Diverging btwn contexts vs. no overlap btwn contexts
# contrast_pairwise                                                     estimate      SE  df z.ratio p.value
# (FALSE diverge - TRUE diverge) - (FALSE no_overlap - TRUE no_overlap)  0.00524 0.00512 Inf 1.024   0.3058 

# [27] "(FALSE converge - TRUE converge) - (FALSE no_overlap - TRUE no_overlap)" / Converging btwn contexts vs. no overlap btwn contexts       
# contrast_pairwise                                                       estimate      SE  df z.ratio p.value
# (FALSE converge - TRUE converge) - (FALSE no_overlap - TRUE no_overlap)  -0.0106 0.00512 Inf -2.075  0.0380 

# contrast                                estimate      SE  df z.ratio p.value
# FALSE converge - TRUE converge         -0.009060 0.00362 Inf -2.505  0.0122 
# FALSE diverge - TRUE diverge            0.006792 0.00362 Inf  1.878  0.0604 
# FALSE full_overlap - TRUE full_overlap -0.009393 0.00362 Inf -2.597  0.0094
# FALSE no_overlap - TRUE no_overlap      0.001553 0.00362 Inf  0.430  0.6675 


## VISUAL ROIs
con_div_intercept.V12_BL = mixed(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "V12_BL"), method = "LRT")
con_div_intercept.V12_BL.lmer = lmer(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "V12_BL"))
print(con_div_intercept.V12_BL)
p=mixed_model_plot(con_div_intercept.V12_BL.lmer,'overlap*con')
ggsave('/Volumes/Data/zoocon/Rev_Hippocampgoal/Figures/FigS2F_V12_BL_EMM_con_div.eps', device = "eps", width = 10, height = 8, units = 'in', dpi = 600)

## For reviewers
tmp_emm = emmeans(con_div_intercept.V12_BL.lmer, ~overlap*same_context)
tmp_df = as.data.frame(tmp_emm)
tmp_indiv_df = con_div_intercept.V12_BL.lmer@frame
write.csv(tmp_df, "~/Downloads/win_con_means_V12.csv")
write.csv(tmp_indiv_df, "~/Downloads/win_con_indivs_V12.csv")

## Move analysis
# grab only relevant contrasts
target = c("same mov dcon no sub", "share mov dcon no sub", "no share mov dcon no sub",
           "same mov scon no sub", "share mov scon no sub", "no share mov scon no sub")
long_data_moves_AOV = data %>% 
  filter(contrast %in% target)

## add overlap and context
long_data_moves_AOV = long_data_moves_AOV %>% 
  mutate(same_context = !(contrast %in% c("same mov dcon no sub", "share mov dcon no sub", "no share mov dcon no sub")),
         same_move = (contrast %in% c("same mov dcon no sub", "same mov scon no sub")) ,
         share_move = (contrast %in% c("share mov dcon no sub",  "share mov scon no sub")) ,
         no_move = (contrast %in% c("no share mov dcon no sub",  "no share mov scon no sub"))) 

# make a new row and combine levels
long_data_moves_AOV$move = NaN
long_data_moves_AOV[long_data_moves_AOV$same_move,"move"] = "same_move"
long_data_moves_AOV[long_data_moves_AOV$share_move,"move"] = "share_move"
long_data_moves_AOV[long_data_moves_AOV$no_move,"move"] = "no_move" 

long_data_moves_AOV$move <- as.factor(long_data_moves_AOV$move)
long_data_moves_AOV$same_context <- as.factor(long_data_moves_AOV$same_context)
long_data_moves_AOV = droplevels(long_data_moves_AOV) # drop un used levels

## HPC 
move_intercept.HIPP_MERGE_BL = mixed(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "HIPP_MERGE_BL"), method = "LRT")
move_intercept.HIPP_MERGE_BL.lmer = lmer(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "HIPP_MERGE_BL"))
print(move_intercept.HIPP_MERGE_BL)
p=mixed_model_plot(move_intercept.HIPP_MERGE_BL.lmer,'motor')
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/FigS2B_HPC_BL_EMM_moves.eps', device = "eps", width = 9, height = 7, units = 'in', dpi = 600)

# pairwise comparisons for reviewer
temp=emmeans(move_intercept.HIPP_MERGE_BL.lmer,pairwise  ~ same_context*move,lmer.df = "asymptotic", adjust = "none") # stats 

## For reviewers
tmp_emm = emmeans(move_intercept.HIPP_MERGE_BL.lmer, ~same_context*move)
tmp_df = as.data.frame(tmp_emm)
tmp_indiv_df = move_intercept.HIPP_MERGE_BL.lmer@frame
write.csv(tmp_df, "~/Downloads/moves_means_hipp.csv")
write.csv(tmp_indiv_df, "~/Downloads/moves_indivs_hipp.csv")

## BA4ap
move_intercept.BA4ap = mixed(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "BA4ap"), method = "LRT")
move_intercept.BA4ap.lmer = lmer(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "BA4ap"))
print(move_intercept.BA4ap)
p=mixed_model_plot(move_intercept.BA4ap.lmer,'motor')
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/FigS2D_BA4ap_BL_EMM_moves.eps', device = "eps", width = 9, height = 7, units = 'in', dpi = 600)

# pairwise comparisons for reviewer
temp=emmeans(move_intercept.BA4ap.lmer,pairwise  ~ same_context*move,lmer.df = "asymptotic", adjust = "none") # stats 

## For reviewers
tmp_emm = emmeans(move_intercept.BA4ap.lmer, ~same_context*move)
tmp_df = as.data.frame(tmp_emm)
tmp_indiv_df = move_intercept.BA4ap.lmer@frame
write.csv(tmp_df, "~/Downloads/moves_means_ba4ap.csv")
write.csv(tmp_indiv_df, "~/Downloads/moves_indivs_ba4ap.csv")
