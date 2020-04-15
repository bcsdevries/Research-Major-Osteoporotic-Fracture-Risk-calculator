# T-test and corresponding p-values for ANOVA plot
t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}
Cox_RSF_MICE <- t.test2(0.6251563, 0.5936189, sqrt(0.001053383), sqrt(0.0008693254), 10, 10)
Cox_RSF_regular <- t.test2(0.6251563, 0.5929356, sqrt(0.001053383), 0.0262929, 10, 10)
Cox_ANN <- t.test2(0.6251563, 0.5880895,sqrt(0.001053383),sqrt(0.00177484),10,10)
RSF_MICE_regular <- t.test2(0.5929356, 0.5936189, 0.0262929, sqrt(0.0008693254), 10, 10)
RSF_MICE_ANN <- t.test2(0.5880895, 0.5936189, sqrt(0.00177484), sqrt(0.0008693254), 10, 10)
RSF_regular_ANN <- t.test2(0.5929356, 0.5880895, 0.0262929, sqrt(0.00177484), 10, 10)

ANOVA_outcome <- cbind(Cox_RSF_MICE, Cox_RSF_regular, Cox_ANN,
                       RSF_MICE_regular, RSF_MICE_ANN, RSF_regular_ANN)
