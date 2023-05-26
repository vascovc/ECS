##################################################
##################projeto 2
##################################################
##################################################
#Goncalo Freitas - 98012
##############################
#Tiago Alvim - 95584
##############################
#Vasco Costa - 97746
##############################
# fixar a seed
set.seed(20)

###########################################################################
#exercicio 1
####################
#b
###########################################################################
N = 10000000 
U = runif(N)
x = 10
I.all = x*exp(-((x*U)^2)/2) 
I.hat = mean(I.all)
print(paste("valor de I.hat -> ",I.hat))
#[1] 1.25389810116415

#fazer grafico
png(file="exer_1_hist.png",width=600,height=600)
hist(I.all, main="Histograma dos valores de Monte Carlo", 
	freq = TRUE, xlab = "g(u)" )

#variancia do metodo de monte carlo
print(paste("erro do metodo de monte carlo -> ",sqrt(var(I.all)/N) ))
#[1] 0.000854173476386978

#valor exato
exact.value = (pnorm(x)-0.5)*sqrt(2*pi)
print(paste("valor de I exato -> ",exact.value))
#[1] 1.2533141373155

# Bias: E(teta_hat-teta)
print(paste("Bias -> ",(I.hat-exact.value) ))
#[1] 0.000583963848648672

# Erro percentual
error = abs(exact.value - I.hat)/(exact.value) * 100
print(paste("erro percentual -> ", error))
#[1] 0.0465935738903797

#intervalo de confianca
desvio = sd(x*exp(-((x*U)^2)/2))

LI = I.hat - qnorm(0.975,0,1)*desvio/sqrt(N)
LS = I.hat + qnorm(0.975,0,1)*desvio/sqrt(N)
print(paste("limite inf -> ", LI))
#[1] 1.25222395191388
print(paste("limite sup -> ", LS))
# [1] 1.25557225041442

###########################################################################
#exercicio 2
####################
#a
###########################################################################
#i
##########################
set.seed(20)
a = 0.4
b = 0.1
N = 1000
obs = 100
t = numeric(N)
for (i in 1:N){
	x = numeric(obs)
	y = rnorm(obs)
	x[1] = y[1] #primeira observacao
	num = 0
	denom = 0
	for(n in 2:obs){
		x[n] = a*x[n-1]+b*x[n-1]*y[n-1]+y[n]
		num = num+x[n]*x[n-1]
		denom = denom+(x[n-1])**2	
	}
	t[i] = num/denom
}

a.hat = mean(t)
print(paste("valor de a teorico    -> ",a))
print(paste("valor de a estimado   -> ",a.hat))

print(paste("desvio padrao         -> ",sd(t) ))
print(paste("erro quadratico medio -> ",mean((t-a)^2 )))

print(paste("bias                  -> ",mean(t-a) ))
enviesamento = t-a

png(file="exer_2_boxplot.png",width=600,height=600)
boxplot(enviesamento,main="Boxplot dos enviesamentos",
				ylab="â-a"
)

##########################
#ii
##########################
#intervalo de confianca
#95%
alpha = 0.05
LI = mean(t) - qnorm(1-alpha/2,0,1)*sd(t)/sqrt(N)
LS = mean(t) + qnorm(1-alpha/2,0,1)*sd(t)/sqrt(N)
print(paste("limite inf -> ", LI))
print(paste("limite sup -> ", LS))

#bootstrap empirico e percentil
#N #usar o mesmo numero de a que se tem
mu = mean(t)
a.bootstrap_empirico  = numeric(N)
a.bootstrap_percentil = numeric(N)

for (b in 1:N){
	i = sample(1:N, size = N, replace = TRUE)
	x = t[i]
	a.bootstrap_empirico[b] = mean(x) - mu
	a.bootstrap_percentil[b] = mean(x)
}
png(file="exer_2_ii_emp.png",width=600,height=600)
hist(a.bootstrap_empirico,main="Bootstrap Empírico")

print(paste("media emp -> ", mean(a.bootstrap_empirico) ))
summary(a.bootstrap_empirico)
var(a.bootstrap_empirico)
q1 = quantile(a.bootstrap_empirico,prob=alpha/2) #sample quantile of order 0.025
q2 = quantile(a.bootstrap_empirico,prob=1-alpha/2)
print(paste("limite inf emp -> ", mu-q2 ))
print(paste("limite sup emp -> ", mu-q1 ))

#bootstrap percentil
png(file="exer_2_ii_percent.png",width=600,height=600)
hist(a.bootstrap_percentil,main="Bootstrap Percentil")
print(paste("media percentil -> ", mean(a.bootstrap_percentil) ))
summary(a.bootstrap_percentil)
var(a.bootstrap_percentil)
q1 = quantile(a.bootstrap_percentil,prob=alpha/2) #sample quantile of order 0.025
q2 = quantile(a.bootstrap_percentil,prob=1-alpha/2)
print(paste("limite inf percentil -> ", q1 ))
print(paste("limite sup percentil -> ", q2 ))

#amplitudes
print(paste("amplitude empirico -> ", mu-q1 - (mu-q2)))
print(paste("amplitude percentil -> ",q2-q1))
#mesmo valor


##########################
#b
##########################
#i
##########################
a = 0.4
b = 0.1
n = 100
sigma = 4/9 #assim fica na forma infinita
alpha = 0.05

vcrit = qchisq(alpha, n-1) #lower one tailed test
estat = function(x,sigma){(n-1)*var(x)/sigma}

N = 1000 # num de replicas
testrule = numeric(N)
for (j in 1:N) {
	x = rnorm(n, 0, sqrt(a/(1-b)) ) #sob H0 
	#decisão -- 1 (rejection) or 0 (not rejection)
	testrule[j] = as.integer(estat(x,sigma) < vcrit )
}
testrule
p.reject = mean(testrule)
p.reject
png(file="exer_2_b_i.png",width=600,height=600)
hist(testrule,main="Aceitações(0) e rejeições(1) para a = 0.4 e b = 0.1")

##########################
#ii
##########################
a = 0.1
b = 0.1
sigma = 4/9 #assim fica na forma infinita
alpha = 0.05

vcrit = qchisq(alpha, n-1)
estat = function(x,sigma){(n-1)*var(x)/sigma}

N = 1000 # num de replicas
testrule = numeric(N)
for (j in 1:N) {
	x = rnorm(n,0, sqrt(a/(1-b))  ) #sob H0 
	#decisão -- 1 (rejection) or 0 (not rejection)
	testrule[j] = as.integer(estat(x,sigma) < vcrit )
}
testrule
p.reject = mean(testrule)
p.reject
png(file="exer_2_b_ii.png",width=600,height=600)
hist(testrule,xlim=c(0,1),breaks=c(0,0.2,0.8,1),main="Aceitações(0) e rejeições(1) para a = 0.1 e b = 0.1")
