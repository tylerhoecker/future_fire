rescale01 <- function(x){(x-min(x))/(max(x)-min(x))}

reference = rnorm(1000, .1, 0.1)
small = rnorm(1000, 1, 0.1)
large = rnorm(1000, 5, 0.1)

df <- data.frame('values' = c(reference, small, large),
                 'period' = c(rep('reference',1000),rep('small',1000),rep('large',1000))) %>% 
  mutate(across(values, ~ rescale01(.x)*100, .names = 'z_{.col}'))

ggplot(df) +
  geom_histogram(aes(x = z_values, fill = period), color = 'black', position = 'identity') +
  theme_bw(base_size = 16) +
  labs(x = 'Fire severity', y = 'Count') 


ref_dens <- density(filter(df, period == 'reference')$values, bw = 0.1) 
small_dens <- density(filter(df, period == 'small')$values, bw = 0.1) 
large_dens <- density(filter(df, period == 'large')$values, bw = 0.1) 


ref_mat <- cbind(ref_dens$y, ref_dens$x)
small_mat <- cbind(small_dens$y, small_dens$x)
large_mat <- cbind(large_dens$y, large_dens$x)

emdist::emdr(ref_mat, ref_mat, max.iter = 10000)
emdist::emdr(ref_mat, small_mat, max.iter = 10000)
emdist::emdr(ref_mat, large_mat, max.iter = 10000)

mkd1 <- Compositional::mkde(as.matrix(iris[1:25, 1:4]), thumb = "scott" )
df1 <- as.matrix(cbind('dens' = mkd1, iris[1:25,1:4]))

mkd2 <- Compositional::mkde(as.matrix(iris[100:150, 1:4]), thumb = "scott" )
df2 <- as.matrix(cbind('dens' = mkd2, iris[100:150,1:4]))
emdist::emd(df1, df2, max.iter = 10000)


ggplot(df1) +
  geom_tile(aes(x = Sepal.Length, y = Petal.Length, fill = dens))
