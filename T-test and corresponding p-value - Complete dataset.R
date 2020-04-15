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


Cox_RSF_MICE <- t.test2(0.6966725, 0.6878885, sqrt(0.0002840391), sqrt(0.0003295021), 10, 10)
Cox_RSF_regular <- t.test2(0.6966725, 0.6874454, sqrt(0.0002840391), 0.01352235, 10, 10)
Cox_ANN <- t.test2(0.6966725, 0.6696976,sqrt(0.0002840391),sqrt(0.001562243),10,10)
RSF_MICE_regular <- t.test2(0.6874454, 0.6878885, 0.01352235, sqrt(0.0003295021), 10, 10)
RSF_MICE_ANN <- t.test2(0.6696976, 0.6878885, sqrt(0.001562243), sqrt(0.0003295021), 10, 10)
RSF_regular_ANN <- t.test2(0.6874454, 0.6696976, 0.01352235, sqrt(0.001562243), 10, 10)

ANOVA_outcome <- cbind(Cox_RSF_MICE, Cox_RSF_regular, Cox_ANN,
                       RSF_MICE_regular, RSF_MICE_ANN, RSF_regular_ANN)
