##################################################
##################projeto 1
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
####################
#a
###########################################################################
c = 2/sqrt(2*pi)*exp(1/2) #valor do c otimo, provado como sendo o limite superior

N = 1000 # numro de geracoes

k = 0 #observacoes aceites

j = 0 #iteracoes
x = numeric(N)

f = function(x){exp(-(x**2)/2 + x)*2/sqrt(2*pi)}
while(k < N){
	y = rexp(1) #Y, neste caso e a exponencial
	u = runif(1) #uniforme
	j = j+1
	if ( u < f(y) ){
		k = k+1
		x[k] = y
	}
}
x #valores gerados pelo metodo da rejicao
j #numero de iteracoes

print(paste(round(c*N,digits=0), "numero teorico de iteracoes necessarias"))
print(paste(j, "numero de iteracoes que foram mesmo precisas"))

png(file="exer_1.png",width=600,height=600)
hist(x,prob = TRUE, main = "1000 NPA's")
u = seq(0.0,4.0,0.01)
lines(u,(2/sqrt(2*pi))*exp(-(u**2)/2),col="red")
dev.off();

############################################################################
########b
############################################################################

N = 1000
k = 0 # valores aceites
j = 0 #iteracoes
x = numeric(N)
while(k<N){
	y = rnorm(1)
	j = j+1
	if( y>0){
		k=k+1
		x[k]=y
	}
}
print(paste(j, "numero de iteracoes que foram mesmo precisas"))
png(file="exer_1_b.png",width=600,height=600)
hist(x,prob = TRUE, main = "1000 NPA's")
u = seq(0.0,4.0,0.01)
lines(u,(2/sqrt(2*pi))*exp(-(u**2)/2),col="red")
dev.off();

##########################################################################
# exercicio 2
##########################################################################
n = 10000
U1 = runif(n)
U2 = runif(n)
Y = sqrt(U1)
X = U2*Y
XY = cbind(X,Y)

png(file="exer_2.png",width=600,height=600)
plot(XY, pch=".", xlim=c(0,1),ylim=c(0,1), main="Distribution f(x,y)=2, 0<x<y<1 (n=10000)")
dev.off();

##########################################################################
# 	exercicio 3
##########################################################################
### 	b
#################### gerais a toda a alinea b)
alpha = 3 #so por ser o numero do exercicio
n_total = c(100,1000,10000)

for (n in n_total){
	print(n)
	####################
	### 	i
	####################
	m = numeric(n)
	for (i in 0:n){
		tau = rgamma(n,shape = alpha, rate = alpha)
		x_norm = rnorm(n,0,1/tau) #ja espera o valor ao quadrado por isso 1/tau, nao a raiz
		m[i] = mean(x_norm)
	}
	
	####################
	### 	ii
	####################
	t_student = rt(n,2*alpha)

	###################
	###	iii
	###################
	print("Normal-Gama")
	print(summary(m))
	print(var(m))
	print("R")
	print(summary(t_student))
	print(var(t_student))
}

##########################################################################
# 	exercicio 4
##########################################################################
###	a
#######################
n = 5000

b_1 = rbeta(n,3,3)
b_2 = rbeta(n,2,3)
b_3 = rbeta(n,3,2)

png(file="exer_3.png",width=800,height=600)

par(mfrow=c(1,3)) 
hist(b_1,prob=TRUE,main="Beta(3,3)",xlab="x")
z = seq(-4,4,by=0.01)
lines(z,dbeta(z,3,3),col="red")


hist(b_2,prob=TRUE,main="Beta(2,3)",xlab="x")
z = seq(-4,4,by=0.01)
lines(z,dbeta(z,2,3),col="red")

 
hist(b_3,prob=TRUE,main="Beta(3,2)",xlab="x")
z = seq(-4,4,by=0.01)
lines(z,dbeta(z,3,2),col="red")

dev.off();

#######################
###	b- com a correcao da aula devido ao quadrado
#######################

num = 50
n = 5000

b_1 = numeric(num)
b_2 = numeric(num)
b_3 = numeric(num)

for (i in 0:num){
	beta_1 = rbeta(n,3,3)
	beta_2 = rbeta(n,2,3)
	beta_3 = rbeta(n,3,2)
	b_1[i] = ( 1/n*sum((beta_1-mean(beta_1))^3)/sd(beta_1)^3   )
	b_2[i] = ( 1/n*sum((beta_2-mean(beta_2))^3)/sd(beta_2)^3   )
	b_3[i] = ( 1/n*sum((beta_3-mean(beta_3))^3)/sd(beta_3)^3   )
}

summary(b_1)
var(b_1)

summary(b_2)
var(b_2)

summary(b_3)
var(b_3)
