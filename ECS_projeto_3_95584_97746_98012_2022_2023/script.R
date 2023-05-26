##################################################
##################projeto 3
##################################################
##################################################
#Goncalo Freitas - 98012
##############################
#Tiago Alvim - 95584
##############################
#Vasco Costa - 97746
##############################


###########################################################################
#exercicio 1
###########################################################################
set.seed(21)
install.packages("coda")
library(coda)

x<-c(4.0,4.1,3.9,4.4,3.2,4.0,3.7,4.2,4.5,4.3,3.6,1.9,3.3,1.9,2.9,2.7,2.4,2.9,
3.8,3.5,2.7,3.9,2.8,3.3,2.9,3.8,4.4,5.1,5.2,7.2,6.2,4.8,4.0,2.7,4.4,3.4,4.2,
4.8,5.3,4.5,4.1,4.0,2.9,0.8,5.2,7.3,5.1,5.3,7.1,8.1,7.8,6.9,7.5,6.0,5.0,5.3,
4.8,4.3,5.8,4.6,4.5,4.1,4.6,6.4,6.3,6.2,6.2,6.8,7.5,7.4,7.0,6.7,7.5,6.1,5.7,
5.4,5.3,4.0,3.7,2.5,0.8,1.3,3.3,4.1,5.7,4.3,3.5,3.8,2.0,3.8,4.1,1.8,3.0,4.7,
6.2,6.0,5.3,4.4,3.4,4.7,4.5,3.7,4.3,1.6,2.9,3.6,3.7,3.9,4.6,5.0,5.3,4.7,6.5, 
5.7,5.8,8.0,7.4,6.1,7.6)

# histograma dos dados
hist(x)

#valores iniciais dos parâmetros
n<-length(x)
N<-11000   #11000    7000       
M<-matrix(0,N,3)   

M[1,1]<- alfa0 <- 0.25
M[1,2]<- gama0 <- 0.05
M[1,3]<- rho0 <- 0.25 

#S0
x0<-x[2:n]
S0<-sum(x0)
#S1
x1<-x[1:(n-1)]
S1<-sum(x1)

f <- function(u,b,d){
    	((log(  1+(exp(b*d) - 1 ) * u )) / b)*(((log(  1+(exp(b*d) - 1 ) * u )) / b)>0)*
(((log(  1+(exp(b*d) - 1 ) * u )) / b)<d)  
}

for(i in 2:N){ #numero de simulacoes
    alfa = rgamma(1,n,S0-M[i-1,3]-(n-1)*M[i-1,2])
    M[i,1] = alfa
    
    gama_star = min(x[2:n]-M[i-1,3]*x[1:n-1]) #e o que esta no enunciado
    gama=NaN
    while(is.nan(gama)){ # pode dar nan ao usar a exponencial nas restricoes por isso um ciclo ate se obter um valor valido
        u = runif(1,0,1)
        gama = f(u,(n-1)*alfa,gama_star) #usa se o alfa ja deste ciclo
    }
    M[i,2] = gama

    rho_star = min(1,(x[2:n]-gama)/(x[1:n-1])) #usa-se o gama ja deste ciclo tambem
    rho = NaN
    while(is.nan(rho)){ # pode dar nan ao usar a exponencial nas restricoes 
        rho = f(u,alfa*S1,rho_star)
    }
    M[i,3] = rho
}

M

colMeans(M, dims = 1)
colMeans(M[3000:N,], dims = 1)
summary(M[3000:N,])
sd(M[3000:N,1])
sd(M[3000:N,2])
sd(M[3000:N,3])
var(M[3000:N,])
#mcmc(M)
plot(mcmc(M))
#---------- plot dos percentis (evolution)
cumuplot(mcmc(M))
geweke.diag(mcmc(M))#valor do teste (compare means)
geweke.plot(mcmc(M),frac1=0.1,frac2=0.5,auto.layout = TRUE)

###########################################################################
#exercicio 2
###########################################################################
#c
########################

# Nota: como os valores de x são fixos não faz sentido repetir o processo algumas vezes
	
# Instalação do pacote stats4, necessário
install.packages("stats4")
library(stats4)

# Ficheiro 'caudais.txt'
x = c(23.8,8.5,44.8,51.0,4.8,19.1,47.8,25.4,14.3,66.2,37.0,48.6,28.3,21.0,72.5,28.3,7.9,47.1,64.4,10.7,19.6,19.0,16.9,22.7,
            65.3,33.3,31.0,14.2,7.3,73.4,44.8,50.2,40.4,57.6,32.6,24.0,31.0,73.4,33.9,84.0,17.0,23.7,74.7,31.2,33.1)

# histograma dos dados
hist(x)


# Estimativa inicial de sigma e mu
sigma_i = sqrt(6)/pi * sd(x)                    #valor inicial de sigma para estimativa
sigma_i
mu_i = mean(x) - 0.5772*sigma_i                 #valor inicial de mu para estimativa
mu_i

#################### Uniroot ####################

obj = function(sigma,n,x){
mu = sigma*log(n/sum(exp(-x/sigma)))
return (-sum( -mu*exp((mu-x)/sigma)/sigma^2+x*exp((mu-x)/sigma)/sigma^2-(x-mu)/sigma^2  )-n/sigma )
}

n=length(x)
u = uniroot(obj,lower=5,upper=100,n=n,x=x)
iter = u$iter
# Cálculo das estimativas
sigma.hat = u$root
sigma.hat

mu.hat = sigma.hat*log(n/sum(exp(-x/sigma.hat)))
mu.hat

# Cálculo dos Erros relativos
err.sigma = abs(sigma_i - sigma.hat)/sigma_i * 100
err.sigma

err.mu = abs(mu_i - mu.hat)/mu_i * 100
err.mu



#################### Optimize ####################

# como é 2 dimensões não se usa a função optimize, como nas aulas, mas sim a optim - Ref : https://www.graduatetutor.com/uncategorized/optimization-using-r-programming/
# optim command minimizes a given function, so we need to consider -L* (slide 13)

logL = function(par){# - log-verosimilhança , ie, -L*
	sigma = par[1]
	mu = par[2]
    var = (x-mu)/sigma
	Loglik = ( -length(x)*log(sigma) - sum( var + exp(- var ) ) )
	return( - Loglik)
}

optim(c(sigma_i,mu_i), logL)

# Cálculo das estimativas
sigma.hat = optim(c(sigma_i,mu_i), logL)$par[1]
sigma.hat

mu.hat = optim(c(sigma_i,mu_i), logL)$par[2]
mu.hat

# Cálculo dos Erros relativos
err.sigma = abs(sigma_i - sigma.hat)/sigma_i * 100
err.sigma

err.mu = abs(mu_i - mu.hat)/mu_i * 100
err.mu

# Cálculo dos intervalos de confiança (95%) - Ref : https://rstudio-pubs-static.s3.amazonaws.com/107801_2785d7d7a49744539ef21eaaebe4fe5a.html
parameter.fits = optim(c(sigma_i,mu_i), logL, hessian = T)
hessian <- parameter.fits$hessian
hessian.inv <- solve(hessian)
parameter.se <- sqrt(diag(hessian.inv))
parameter.se

# intervalo - sigma
sigma.low = sigma.hat - 1.96 * parameter.se[1]
sigma.upp = sigma.hat + 1.96 * parameter.se[1]

sigma.low
sigma.upp

# intervalo - mu
mu.low = mu.hat - 1.96 * parameter.se[2]
mu.upp = mu.hat + 1.96 * parameter.se[2]

mu.low
mu.upp

########################
#d
########################

caudais = c(23.8,8.5,44.8,51.0,4.8,19.1,47.8,25.4,14.3,66.2,37.0,48.6,28.3,21.0,72.5,28.3,7.9,47.1,64.4,10.7,19.6,19.0,16.9,22.7,
            65.3,33.3,31.0,14.2,7.3,73.4,44.8,50.2,40.4,57.6,32.6,24.0,31.0,73.4,33.9,84.0,17.0,23.7,74.7,31.2,33.1)

sigma_i = sqrt(6)/pi * sd(caudais)

m<-20000   #chain length
aquecimento<-m/4 #burn-in
set.seed(21)            # seed para obter sempre os  mesmos valores

verosimilhanca = function(caudais,sigma,mu){
	return ((1/sigma)^length(caudais)*prod(exp(-(caudais-mu)/sigma-exp(-(caudais-mu)/sigma))))
}


x = numeric(m)
x[1]<-runif(1,0,100)
x[1]
u<-runif(m,0,1)
k=0
for (i in 2:m){
	xt<-x[i-1]
    	y<-rchisq(1,df=xt) # proponent distribution 
	
 	num<-verosimilhanca(caudais,sigma_i,y)*dchisq(xt,df=y)

	den<-verosimilhanca(caudais,sigma_i,xt)*dchisq(y,df=xt)

    if (u[i]<=min(num/den,1)){
   		x[i]<-y
    }
    else{
        x[i]<-xt  #y os rejected
	  k<-k+1
	}
}
print( (k/m)*100 )	# % de rejeições
mean(x[aquecimento:m]) #estimativa


# Cálculo do Erro relativo
mu_i = mean(caudais) - 0.5772*sigma_i  # valor que se pretende estimar
mu_i

mu.hat = mean(x[aquecimento:m])    #estimativa

err = abs(mu_i - mu.hat)/mu_i * 100
err