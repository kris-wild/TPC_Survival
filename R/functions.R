########
# functions used for figures: 

library(lme4)
library(lmerTest)

# Define the function
extract_anova_dataframe <- function(model) {
  # Perform ANOVA on the model
  anova_summary <- anova(model)
  
  # Initialize a dataframe to store the results
  anova_df <- data.frame(
    Effect = rownames(anova_summary),
    df1 = anova_summary$NumDF,
    df2 = anova_summary$DenDF,
    F_value = anova_summary$`F value`,
    p_value = ifelse(anova_summary$`Pr(>F)` < 0.01, "<0.01", sprintf("%.3f", anova_summary$`Pr(>F)`)),
    stringsAsFactors = FALSE
  )
  
  return(anova_df)
}
