Oblig 3
================
Isak Edvardsen
3 april 2020

``` r
# Simulation from Gamma(alpha,1) distribution using accept-reject 
# and a Gamma(a,b) proposal distribution.

gammasim <- function(n=10000, alfa=3.7) 
{
  a=3
  b=a/alfa
  e = NULL
  x = NULL
  
  for (i in 1:n) {
    
    y = rgamma(1,a,b)
    u = runif(1)
    
    k = b^(-a)*y^(alfa-a)*exp(-(1-b)*y)
      
    if (u < k) x = c(x,y)
  }
  return(x)
}
p = gammasim()
hist(p, freq = FALSE)
curve(dgamma(x,3.7,1), add=TRUE,col='red')
```

![](Oblig3_files/figure-gfm/Gamma%20simulation,%20might%20be%20useful-1.png)<!-- -->

``` r
# Monte carlo integration

mci_int <- function(n) 
{
  n <- 10000 
  u <- rgamma(n, 3.7, 1)
  seq <- 1:(n)
  
  # The function to be integrated:
  
  mci.ex <- function(x) x*log(x)
  
  par(mfrow=c(1,3))
  
  curve(mci.ex(x),xlab="x",ylab="w(x)", from = 0, to = 2)
  
  # The monte carlo sum:
  
  w_x <- mci.ex(u)
  
  hplot <- cumsum(w_x)/seq(1:n)
  
  # stdh <- sqrt( cumsum(w_x^2)/seq - (w_x)^2) Nan produced???
  
  par(new=F)
  
  hist(w_x, xlab="Generated Values of w(X) where X~Gamma(3.7, 1)",
       freq=F,col="green",breaks=30,ylab="",main="")
  
  par(new=F)
  
  plot(hplot,type="l",col="blue",xlab="Iteration",
       ylab="Estimate", xlim=c(1,n), ylim=c(1,7))
  
  # par(new=T)
  
  # plot(hplot+stdh/sqrt(seq),type="l",col="green",xlab="",
  #      ylab="",lty=2)
  
  # par(new=T)
  
  # plot(hplot-stdh/sqrt(seq),type="l",col="red",
  #     xlab="",ylab="",lty=2)
  
  return(1/n*sum(w_x))  
}

mci_int(10000)
```

![](Oblig3_files/figure-gfm/Monte%20carlo%20integration-1.png)<!-- -->

    ## [1] 5.344904

``` r
# Computes the MCI approximation to E[Xlog(X)] for X~Gamma(3.7, 1).

mci.gamma = function(n=10000)
{   
    x = rgamma(n,3.7, 1)
    return(mean(x*log(x)))
}

mci.gamma()
```

    ## [1] 5.255851

``` r
riemann.gamma = function(n=6000, alfa = 3.7, lambda =1)
{   
    x = rgamma(n, alfa, lambda)
    x = sort(x)
    
    f <- function(x, alfa = 3.7, lambda =1) lambda^alfa * 1/gamma(alfa)*x^(alfa-1)*exp(-lambda*x)
      
    g <- function(x) x*log(x)
    

    k = NULL
    for (i in 1:n-1) k[i] <- ( g(x[i]) * f(x[i]) * (x[i+1]-x[i]) )
    
    plot(cumsum(k), type = 'l', lwd=2)
    
    # georg: k <- ( g(x[-n]) * f(x[-n]) * (x[-1]-x[-n]) )

    #Z = x
    #p <- function(x) cumsum( g(Z[x]) * f(Z[x]) * (Z[x+1]-Z[x])/x )
    #curve((p(x)), from = 1, to=n-1 )
    
    
    
    return(sum(k))
}

k = riemann.gamma(n=5000)
```

![](Oblig3_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
k
```

    ## [1] 5.319834

``` r
h <- function(x, alfa=3.7) alfa*log(x) + log(log(x)) - x 
  
h_d <- function(x, alfa=3.7) alfa/x + 1/(x*log(x)) - 1

h_dd <-function(x, alfa=3.7) return(-(alfa*(log(x))^2 + log(x)+1)/(x^2*(log(x))^2))

sd <-function(x, n) sqrt(-1/(n*h_dd(x0)))

x0 = 4.3773076886

sigma = sd(x0, 10000)


I = exp(h(x0)) * sqrt(-2*pi/h_dd(x0))*(pnorm(sqrt(-h_dd(x0))*(Inf-x0)) - pnorm(sqrt(-h_dd(x0))*(0-x0)))

I/gamma(3.7)
```

    ## [1] 5.156945
