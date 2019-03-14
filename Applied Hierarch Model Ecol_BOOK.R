############ Applied hierarchical modeling ecology ###############
#################### Chapter 8  ########################

#### Book section 8.2.5 ####
# Define function to compute cell probs for binned distance sampling
cp.ri <- function(radius1, radius2, sigma){
  Pi <- 3.141593
  a <- Pi*radius2^2 - Pi*radius1^2
  integrate(function(x, s = sigma) exp(-x^2/(2+s^2))*x,
            radius1,
            radius2)$value*(2*Pi/a)
}

# Define distance intervals and compute multinomial probabilities
delta <- 0.5 # Width of distance bins
B <- 3 # Max count distance
dist.breaks <- seq(0, B, delta) # Make the interval cut points
nD <- length(dist.breaks) - 1
sigma <- 1
p.x <- rep(NA, nD) # Conditional detection probabilities
for(i in 1:nD){
  p.x[i] <- cp.ri(dist.breaks[i], dist.breaks[i + 1], sigma = 1)
}

area <- 3.141593 * dist.breaks[-1]^2
ring.area <- diff(c(0, area))

# Pr(detection| in ring)*Pr(in ring)
cp <- p.x*ring.area/sum(ring.area)

#### Book section 8.2.5.1 ####
sim.pdata <- function(N = 1000,
                      sigma = 1,
                      B = 3,
                      keep.all = FALSE){
# Function simulate coordinates of individuals on a square
# Square is [0,2*B] x [0,2*B], with a count location on the center
# point(B, B)
# Function arguments:
#      N: total population size in the square
#      Sigma: scale of half-normal detection function
#      B: circle radius
#      keep.all: return the data for y = 0 individuals or not
  
  
# Plot the detection function
par(mfrow = c(1, 2))
curve(exp(-x^2/(2*sigma^2)),
      0,
      B,
      xlab = "Distance (x)",
      ylab = "Detection prob.",
      lwd = 2,
      main = "Detection function",
      ylim = c(0, 1))

# Simulated and plot simulated data
library(plotrix)
u1 <- runif(N, 0, 2*B) # (u1, u2) coordinates of N individuals
u2 <- runif(N, 0, 2*B)
d <- sqrt((u1 - B)^2 + (u2 - B)^2) # Distance to center point of square
plot(u1,
     u2,
     asp = 1,
     pch = 1,
     main = "Point transect")
N.real <- sum(d <= B) # Population size inside of count circle

# Can only count individuals in the circle, so set to zero detection probability of individuals in the corners (thereby truncating them)
p <- ifelse(d < B, 1, 0) * exp(-d*d/(2*(sigma^2)))

# Now we decide whether each individual is detected or not
y <- rbinom(N, 1, p)
points(u1[d <= B], u2[d <= B], pch = 16, col = "black")
points(u1[y == 1], u2[y == 1], pch = 16, col = "blue")
points(B, B, pch = "+", cex = 3, col = "red")
draw.circle(B, B, B)

# Put all the data in the matrix 
#     Note: we don't care about y, u, or v naomally

if(!keep.all){
  u1 <- u1[y == 1]
  u2 <- u2[y == 1]
  d <- d[y == 1]
}
return(list(N = N, sigma = sigma, B = B, u1 = u1, u2 = u2, d = d, y = y, N.real = N.real))
}

# Obtain a data set by distance sampling population of N = 1000 out to a distance of B = 3
set.seed(1234)
tmp <- sim.pdata(N = 1000,
                 sigma = 1,
                 keep.all = FALSE,
                 B = 3)


#### Package AHMBook ####

install.packages("AHMbook")
library(AHMbook)
ls("package:AHMbook")
?sim.pdata

set.seed(12345)
simu.data <- sim.pdata(N = 1000,
                       sigma = 1,
                       B = 3,
                       keep.all = TRUE,
                       show.plot = TRUE)

tmp <- sim.pdata(N = 1000,
                 sigma = 1,
                 keep.all = FALSE,
                 show.plot = FALSE)
summary(tmp)
str(tmp)

# Bin the data and tabulate the bin frequencies. Be sure to pad the 0s !
delta <- 0.5 # Width of distance bins
B <- 3 # MAx count distance
dist.breaks <- seq(0, B, delta) # Make the interval cut points
dclass <- tmp$d %/% delta + 1 # %/% division entiere
nD <- length(dist.breaks) - 1 # How many intervals do you have
y.obs <- table(dclass)
y.padded <- rep(0, nD)
names(y.padded) <- 1:nD
y.padded[names(y.obs)] <- y.obs
y.obs <- y.padded
