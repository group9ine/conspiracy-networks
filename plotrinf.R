rinf <- function(x, a, b, p){
  return(a/b*(1-x) - a/b*log(1-x) + p - a/b - p*(1-x) - x)
}
x <- seq(0,1,0.01)
library(ggplot2)
ggplot() + geom_hline(aes(x=x, yintercept=0), linewidth=sz_small, linetype = 'dashed') +
  geom_line(aes(x=x, y=rinf(x, 0.15, 0.3, 0.1)), color=sir_pal[3], linewidth=sz_small) +
  theme_sir()
