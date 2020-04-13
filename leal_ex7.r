Lagrangian <- function (x,y){
  #sort data with respect to x
  o = order(x)
  x = x[o]
  y = y[o]
  
  n = length(x)
  
  #setup equation
  f = "f<-function(x)"
  for (i in 1:n){
    term = y[i]
    for (j in 1:(n)){
      if (i == j){
        next
      }
      num = paste("(x-",x[j],")",sep="")
      denom = paste("(",x[i],"-",x[j],")",sep="")
      mult = paste(num,denom,sep="/")
      term = paste(term,mult,sep=" * ")
    }
    if(i==n){
      f = paste(f,term)
    }
    else{
    f = paste(f,term,"+")
    }
  }
  poly = f
  f = eval(parse(text = f))
  return (list(func=f,poly=poly))
}

Neville <- function (x,y,given_x){
  n = length(x)
  rownames = c(1:n)
  colnames = c("i","xi","|x-xi|","Pi,0=f(xi)",paste("Pi,",1:(n-1),sep=""))

  #setup matrix
  m = matrix(nrow=length(rownames),ncol=length(colnames),dimnames = list(NULL,colnames))
  m[,"xi"]=x
  m[,"|x-xi|"]=c(abs(given_x - x[1:n]))
  m[,"Pi,0=f(xi)"]=y
  m = m[order(m[, "|x-xi|"]), ]
  m[,'i']=c(1:n)
  
  for (i in 1:(n-1)){
    for (j in 1:(n-i)){
      num = ((given_x-m[j,"xi"])*m[j+1,i+3])+((m[j+i,"xi"]-given_x)*m[j,i+3])
      denom = m[i+j,"xi"]-m[j,"xi"]
      m[j,(i+4)] = num/denom
    }
  }
  return(list(mat=m,value=m[1,ncol(m)]))
}

x = c(4,1,3,5)
y = c(1.3863,0,1.0986,1.6094)
x = c(1995,2000,2004,2005,2010,2015)
y = c(68349452,75505061,80824322,82079348,87940171,93440274)

lagrangian = Lagrangian(x,y)
print(lagrangian$f(2004))
print(lagrangian$poly)
n = Neville(x,y,2004)
print(n$mat)
print(n$value)

plot(1990:2020, lagrangian$f(1990:2020),pch=20,main="Philippine
population from 1995 to 2015", col = "red",xlab="Year",ylab="Population Count")
lagrangeModel = lm(y~poly(x,5,raw=TRUE))
lines(x,predict(lagrangeModel),col="green")