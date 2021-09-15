# libraries
require(ggplot2)
require(plyr)
require(dplyr)
require(tidyr)
require(afex)
require(emmeans)

# plot EMMs
mixed_model_plot = function(model_2_plot, analysis_type){
  
  tmp_plot_name = as.character(model_2_plot@call)[3] # name of data 
  
  if (analysis_type == "seq*con"){
    # save model 
    tmp_emm = emmeans(model_2_plot, ~same_sequence*same_context)
    tmp_df = as.data.frame(tmp_emm)
    
    # organizing factors
    tmp_df$same_sequence=mapvalues(tmp_df$same_sequence, from = c(TRUE,FALSE), to = c("Same Sequence", "Diff Sequence"))
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_sequence = mapvalues(indiv_df$same_sequence, from = c(TRUE,FALSE), to = c("Same Sequence", "Diff Sequence"))
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    
    # MAKE PLOT
    p=tmp_df %>% 
      ggplot(aes(x = same_sequence, y = emmean, fill = same_context)) + 
      geom_bar(stat = "identity", aes(color = same_context, fill = same_context), width = 1, position = position_dodge2(width =  1)) +
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2, position = position_dodge(width = 1)) +
      scale_color_manual(values = c("Same Context" = "dodgerblue4", "Diff Context" = "red3")) + 
      scale_fill_manual(values = c("Same Context" = "dodgerblue4", "Diff Context" = "red3")) + 
      geom_jitter(aes(y = PS, color = NULL), position = position_jitterdodge(dodge.width = 1), color = "black", cex = 2 , data = indiv_df) +
      scale_x_discrete(limits = c("Same Sequence", "Diff Sequence")) +
      labs(title = paste(tmp_plot_name), 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Sequence Pair") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank()) + 
      theme(panel.grid.minor = element_blank()) + 
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold"))
    
  } else if (analysis_type == "overlap*con") {
    
    # get emms
    tmp_emm = emmeans(model_2_plot, ~overlap*same_context)
    tmp_df = as.data.frame(tmp_emm)

    # organizing factors
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    tmp_df$overlap = mapvalues(tmp_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    tmp_df$overlap = factor(tmp_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    indiv_df$overlap = mapvalues(indiv_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    indiv_df$overlap = factor(indiv_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    # plot
    p = tmp_df %>% 
      ggplot(aes(x = same_context, y = emmean, fill = overlap)) + 
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2) +
      geom_point(aes(color = overlap), size = 5) +
      scale_x_discrete(labels = c("Same Context" = "Same \n Context", "Diff Context" = "Diff \n Context")) +
      facet_grid(~overlap) + 
      geom_jitter(aes(y = PS, color = NULL), color = "black", width = 0.1, cex = 2 , data = indiv_df) + 
      labs(title = tmp_plot_name, 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Context Pair") + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank()) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=18, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold")) +
      scale_fill_manual("legend", values = c("Converging" = "#56B4E9", "Same Sequence" = "#0072B2", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      scale_color_manual("legend", values = c("Converging" = "#56B4E9", "Same Sequence" = "#0072B2", "Diverging" = "#009E73", "No Overlap" = "#D55E00"))
    
  } else if (analysis_type == "motor") {
    
    # get emms
    tmp_emm = emmeans(model_2_plot, ~move*same_context)
    tmp_df = as.data.frame(tmp_emm)
    
    # organizing factors
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    tmp_df$move = mapvalues(tmp_df$move, from = c("same_move", "share_move", "no_move"), to = c("Same Moves", "Shared Moves", "No Moves"))
    tmp_df$move = factor(tmp_df$move, levels = c("Same Moves", "Shared Moves", "No Moves"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    indiv_df$move = mapvalues(indiv_df$move, from = c("same_move", "share_move", "no_move"), to = c("Same Moves", "Shared Moves", "No Moves"))
    indiv_df$move = factor(indiv_df$move, levels = c("Same Moves", "Shared Moves", "No Moves"))
    
    # plot
    p = tmp_df %>% 
      ggplot(aes(x = same_context, y = emmean, fill = move)) + 
      scale_color_manual(values = c("Same Moves" = "#0072B2", "Shared Moves" = "#009E73", "No Moves" = "#D55E00")) + 
      scale_fill_manual(values = c("Same Moves" = "#0072B2", "Shared Moves" = "#009E73", "No Moves" = "#D55E00")) + 
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2) +
      geom_point(aes(color = move), size = 5) +
      scale_x_discrete(labels = c("Same Context" = "Same \n Context", "Diff Context" = "Diff \n Context")) +
      facet_grid(~move) + 
      geom_jitter(aes(y = PS, color = NULL), color = "black", width = 0.1, cex = 2 , data = indiv_df) + 
      labs(title = tmp_plot_name, 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Context Pair") + 
      theme_bw() + 
      theme(panel.grid.minor = element_blank()) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold")) 
    
  } else {
    error('incorrect analysis type supplied')
  }
  
  # display and output plot
  print(p)
  return(p)
  
}

main_effect_plot <- function(model_2_plot, analysis_type) {
  tmp_plot_name = as.character(model_2_plot@call)[3] # name of data 
  
  if (analysis_type == "overlap") {
    # save model
    context_effect = emmeans(con_div_intercept.HIPP_MERGE_BL.lmer, pairwise ~ same_context | overlap, lmer.df = "asymptotic", adjust = "none" )
    tmp_df = as.data.frame(context_effect$emmeans)
    
    # organizing factors
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    tmp_df$overlap = mapvalues(tmp_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    tmp_df$overlap = factor(tmp_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    indiv_df$overlap = mapvalues(indiv_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    indiv_df$overlap = factor(indiv_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    overlap_effect = tmp_df[tmp_df$same_context == "Same Context", ]
    overlap_effect_indiv = indiv_df[indiv_df$same_context == "Same Context", ]
    # plot
    p = overlap_effect %>% 
      ggplot(aes(x = overlap, y = emmean, fill = overlap)) + 
      geom_bar(stat = "identity", aes(color = overlap, fill = overlap), width = 0.9) +
      scale_color_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      scale_fill_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2, position = position_dodge(width = 1)) +
      geom_jitter(aes(y = PS, color = NULL), color = "black", width = 0.1, cex = 2 , data = overlap_effect_indiv) +
      scale_x_discrete(labels = c("Same Sequence" = "Same \n Sequence")) +
      labs(title = tmp_plot_name, 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Overlap") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank()) + 
      theme(panel.grid.minor = element_blank()) + 
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +       
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold"))
    
  } else if(analysis_type == "overlap*context") {
    
    # save model
    context_effect = emmeans(con_div_intercept.HIPP_MERGE_BL.lmer, pairwise ~ same_context | overlap, lmer.df = "asymptotic", adjust = "none" )
    tmp_df = as.data.frame(context_effect$emmeans)
    
    # organizing factors
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    tmp_df$overlap = mapvalues(tmp_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    tmp_df$overlap = factor(tmp_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    indiv_df$overlap = mapvalues(indiv_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    indiv_df$overlap = factor(indiv_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    
    # subset same and different contexts. subset to original DF into same context. Take out the subset data and rename
    tmp_df[tmp_df$same_context == "Same Context", "emmean" ] = tmp_df[tmp_df$same_context == "Same Context", "emmean"] - 
      tmp_df[tmp_df$same_context == "Diff Context","emmean"]
    cx_effect = tmp_df[tmp_df$same_context == "Same Context", ] 
    
    diffs = indiv_df[indiv_df$same_context == "Same Context", 'PS'] - indiv_df[indiv_df$same_context == "Diff Context",'PS'] 
    indiv_df[indiv_df$same_context == "Same Context","PS"] = diffs
    indiv_df = filter(indiv_df, same_context == "Same Context")
    
    # plot
    p = cx_effect %>% 
      ggplot(aes(x = overlap, y = emmean, fill = overlap)) + 
      geom_bar(stat = "identity") +
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2) +
      geom_jitter(aes(y = PS, color = NULL), color = "black", width = 0.1, cex = 2 , data = indiv_df) +
      theme_bw() + 
      theme(panel.grid.minor = element_blank()) + 
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold")) +
      scale_color_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      scale_fill_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      scale_x_discrete(labels = c("Same Sequence" = "Same \n Sequence")) +
      labs(title = tmp_plot_name, 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Overlap") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank()) + 
      theme(panel.grid.minor = element_blank()) + 
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +       
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold")) 
    
  } else if(analysis_type == "overlap_diff_context") { 
    # save model
    context_effect = emmeans(con_div_intercept.HIPP_MERGE_BL.lmer, pairwise ~ same_context | overlap, lmer.df = "asymptotic", adjust = "none" )
    tmp_df = as.data.frame(context_effect$emmeans)
    
    # organizing factors
    tmp_df$same_context=mapvalues(tmp_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    tmp_df$same_context=factor(tmp_df$same_context, levels = c("Same Context", "Diff Context"))
    tmp_df$overlap = mapvalues(tmp_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    tmp_df$overlap = factor(tmp_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    # for the indiv people
    indiv_df=model_2_plot@frame
    indiv_df$same_context = mapvalues(indiv_df$same_context, from = c(TRUE,FALSE), to = c("Same Context", "Diff Context"))
    indiv_df$same_context=factor(indiv_df$same_context, levels = c("Same Context", "Diff Context"))
    indiv_df$overlap = mapvalues(indiv_df$overlap, from = c("converge", "full_overlap", "diverge", "no_overlap"), to = c("Converging", "Same Sequence", "Diverging", "No Overlap"))
    indiv_df$overlap = factor(indiv_df$overlap, levels = c("Same Sequence", "Converging", "Diverging", "No Overlap"))
    
    overlap_effect = tmp_df[tmp_df$same_context == "Diff Context", ]
    overlap_effect_indiv = indiv_df[indiv_df$same_context == "Diff Context", ]
    
    # plot
    p = overlap_effect %>% 
      ggplot(aes(x = overlap, y = emmean, fill = overlap)) + 
      geom_bar(stat = "identity", aes(color = overlap, fill = overlap), width = 0.9) +
      scale_color_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      scale_fill_manual(values = c("Converging" = "#56B4E9", "Same Sequence" = "dodgerblue4", "Diverging" = "#009E73", "No Overlap" = "#D55E00")) + 
      geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0, size = 2, position = position_dodge(width = 1)) +
      geom_jitter(aes(y = PS, color = NULL), color = "black", width = 0.1, cex = 2 , data = overlap_effect_indiv) +
      scale_x_discrete(labels = c("Same Sequence" = "Same \n Sequence")) +
      labs(title = tmp_plot_name, 
           y = "Pattern Similarity \n (Estimated Marginal Mean)", 
           x = "Overlap") + 
      theme_bw() + 
      theme(panel.grid.major = element_blank()) + 
      theme(panel.grid.minor = element_blank()) + 
      theme(panel.border = element_blank()) +
      theme(axis.line = element_line(color = 'black')) +       
      theme(plot.title = element_text(hjust = 0.5)) + 
      theme(axis.text=element_text(size=20, face="bold"),
            axis.title=element_text(size=25,face="bold"),
            title=element_text(size=17.5,face="bold")) + 
      theme(legend.position = "none")  + 
      theme(strip.text.x = element_text(size = 17.5, face = "bold")) 
  }
}
