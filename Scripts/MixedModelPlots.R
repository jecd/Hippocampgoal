# libraries
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(afex)
library(emmeans)

source("/Volumes/Data/zoocon/Hippocampgoal/Scripts/LMM_plot.R")

data = read.csv("/Volumes/Data/zoocon/Hippocampgoal/Data/long_cue_period_RSA_results.csv")

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
HIPP_MERGE_BL.inercept.lmer = lmer(PS ~ same_sequence*same_context + (1|subject), data = filter(long_data_AOV, ROI == "HIPP_MERGE_BL"))
print(HIPP_MERGE_BL.intercept)
p=mixed_model_plot(HIPP_MERGE_BL.inercept.lmer,'seq*con')
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/Fig2B_HPC_BL_EMM_seq_con_bars.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

## Converging vs. Diverging 
target = c("dstart send diff context no sub", "sstart dend diff context no sub", "dstart dend diff context no sub",
          "dstart send same context no sub", "sstart dend same context no sub", "dstart dend same context no sub",
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
con_div_intercept.HIPP_MERGE_BL = mixed(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "HIPP_MERGE_BL"))
con_div_intercept.HIPP_MERGE_BL.lmer = lmer(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "HIPP_MERGE_BL"))
print(con_div_intercept.HIPP_MERGE_BL)

# within context 
p = main_effect_plot(con_div_intercept.HIPP_MERGE_BL.lmer, "overlap") 
print(p)
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/Fig2C_HPC_BL_EMM_win_context.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

# between context
p = main_effect_plot(con_div_intercept.HIPP_MERGE_BL.lmer, "overlap*context" )
print(p)
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/Fig2D_HPC_BL_EMM_btwn_context.eps',device = "eps", width = 8, height = 6, units = 'in', dpi = 600)

## VISUAL ROIs
con_div_intercept.V12_BL = mixed(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "V12_BL"))
con_div_intercept.V12_BL.lmer = lmer(PS ~ same_context*overlap + (1|subject), data =  filter(long_data_con_div_AOV, ROI == "V12_BL"))
print(con_div_intercept.V12_BL)
p=mixed_model_plot(con_div_intercept.V12_BL.lmer,'overlap*con')
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/FigS2F_V12_BL_EMM_con_div.eps', device = "eps", width = 10, height = 8, units = 'in', dpi = 600)


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

## BA4ap
move_intercept.BA4ap = mixed(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "BA4ap"), method = "LRT")
move_intercept.BA4ap.lmer = lmer(PS ~ same_context*move + (1|subject), data =  filter(long_data_moves_AOV, ROI == "BA4ap"))
print(move_intercept.BA4ap)
p=mixed_model_plot(move_intercept.BA4ap.lmer,'motor')
ggsave('/Volumes/Data/zoocon/Hippocampgoal/Figures/FigS2D_BA4ap_BL_EMM_moves.eps', device = "eps", width = 9, height = 7, units = 'in', dpi = 600)
