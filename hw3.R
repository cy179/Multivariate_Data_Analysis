multi.mean_and_cov = function(data, na.action=F){
  setClass('mean_and_cov', representation (mean = 'numeric', cov = 'matrix'))
  returnValue(new('mean_and_cov', mean = apply(data, 2, mean), cov = cov(data)))
}

test.TSquare = function(mu=rep(0,ncol(data)), data){
  n = nrow(data)
  p = ncol(data) 
  mean = apply(data, 2, mean)
  diff = mean - mu
  T.sq = n * t(diff) %*% solve(cov(data)) %*% (diff)
  F = (n - p)/((n - 1) * p) * T.sq;
  cat('T-square Test\n---------\n','T-square Value: ', T.sq, '\nP-value: ', 1 - pf(F, p, n - p),'\n', sep = '')
}

TF.convert = function(T.sq = 0, F=0, n, p, T_to_F = TRUE){
  if(T_to_F){
    F = (n - p)/((n - 1) * p) * T.sq
    df1 = p
    df2 = n-p
    r = data.frame(F,df1,df2)
    returnValue(r)
  }else{
    T.sq = F*(n-1)*p/(n-p)
    df1 = p
    df2 = n-1
    returnValue(data.frame(T.sq,df1,df2))
  }
}


T.ellipse = function(data, confident.level=0.95){
  n = nrow(data)
  p = ncol(data)
  cat('\n\nT.sq confidence ellipse')
  center = colMeans(data)
  cat('\ncenters\n')
  print(center)
  T.critical = TF.convert(F=qf(confident.level,p,(n-p)), n=n, p=p,T_to_F = FALSE)
  cat('\nThe Critical T Value\n')
  print(T.critical)
  S = cov(data)
  cat('\ncovariance matrix\n')
  print(S)
  ei = eigen(S)
  cat('\nwith eigenvalue and vectors\n')
  print(ei)
  length = 2*sqrt(ei$values) * sqrt(1/nrow(data)*T.critical$T.sq)
  cat('\nlength\n')
  print(length)
  axes = NULL
  for(i in c(1:length(length))){
    axes = cbind (axes, length[i]/2 * ei$vectors[,i])
  }
  cat('\naxes\n')
  print(axes)
  distance.simontaneous = sqrt(T.critical$T.sq/nrow(data))*sqrt(diag(S))
  CI.simontaneous = rbind (center - distance.simontaneous, center + distance.simontaneous)
  rownames(CI.simontaneous) = c('low','high')
  cat('\nsimontaneous confident interval:\n')
  print(t(CI.simontaneous))
  t.Bonferroni = (qt(df=nrow(data)-1,p=(1-confident.level)/(2*ncol(data)),lower.tail=F))
  cat('\nt value for Bonferroni: \n')
  print(t.Bonferroni)
  distance.Bonferroni = t.Bonferroni*sqrt(diag(S)/nrow(data))
  CI.Bonferroni = rbind(center-distance.Bonferroni, center+distance.Bonferroni)
  rownames(CI.Bonferroni) = c('low','high')
  cat('\nBonferroni confident interval:\n')
  print(t(CI.Bonferroni))
}

# question 5.1
mu =  matrix(c(7,11))
X = matrix(c(2,8,6,8,12,9,9,10),ncol = 2)
summary = multi.mean_and_cov(X)
print(summary)
test.TSquare(mu=c(7,11),X)

# question 5.3
test.TSquare(mu=c(7,11),X)
cov.xbar = (nrow(X)-1)/nrow(X) * summary@cov
cov.mu = cov.xbar + (summary@mean - mu) %*% t(summary@mean - mu)
T.sq = (nrow(X)-1) * (det(cov.mu)/det(cov.xbar)-1)

# question 5.4
x1 = c(3.7,5.7,3.8,3.2,3.1,4.6,2.4,7.2,6.7,5.4,3.9,4.5,3.5,4.5,1.5,8.5,4.5,6.5,4.1,5.5)
x2 = c(48.5,65.1,47.2,53.2,55.5,36.1,24.8,33.1,47.4,54.1,36.9,58.8,27.8,40.2,13.5,56.4,71.6,52.8,44.1,40.9)
x3 = c(9.3,8.0,10.9,12.0,9.7,7.9,14,7.6,8.5,11.3,12.7,12.3,9.8,8.4,10.1,7.1,8.2,10.9,11.2,9.4)
x = cbind(x1,x2,x3)
T.ellipse(x,0.9)

# question 5.5
T.ellipse(x)

# question 5.18
x1 = c(468,428,514,547,614,501,421,527,527,620,587,541,561,468,614,527,507,580,507,521,574,587,488,488,587,421,481,428,640,
       574,547,580,494,554,647,507,454,427,521,468,587,507,574,507,494,541,362,408,594,501,687,633,647,647,614,633,448,408,
       441,435,501,507,620,415,554,348,468,507,527,527,435,660,733,507,527,428,481,507,527,488,607,561,614,527,474,441,607)
x2 = c(41,39,53,67,61,67,46,50,55,72,63,59,53,62,65,48,32,64,59,54,52,64,51,62,56,38,52,40,65,61,64,64,53,51,58,65,52,57,
       66,57,55,61,54,53,41,47,36,28,68,25,75,52,67,65,59,65,55,51,36,60,54,42,71,52,69,28,49,54,47,47,50,70,73,45,62,37,
       48,61,66,41,69,59,70,49,41,47,67)
x3 = c(26,26,21,33,27,29,22,23,19,32,31,19,26,20,28,21,27,21,21,23,25,31,27,18,26,16,26,19,25,28,27,28,26,21,23,23,28,21,
       26,14,30,31,31,23,24,25,17,17,23,26,33,31,29,34,25,28,24,19,22,20,21,24,36,20,30,18,25,26,31,26,28,25,33,28,29,19,
       23,19,23,28,28,34,23,30,16,26,32)
X = cbind(x1,x2,x3)
mu = c(500,50,30)
multi.mean_and_cov(X)
test.TSquare(mu=mu,X)
T.ellipse(X)

# Problem 5.21
x1 = c(1.103,0.842,0.925,0.857,0.795,0.787,0.933,0.799,0.945,0.921,0.792,0.815,0.755,0.88,0.9,0.764,0.733,0.932,0.856,0.89,0.688,0.94,0.493,0.835,0.915)
x2 = c(1.052,0.859,0.873,0.744,0.809,0.779,0.88,0.851,0.876,0.906,0.825,0.751,0.724,0.866,0.838,0.757,0.748,0.898,0.786,0.95,0.532,0.85,0.616,0.752,0.936)
x3 = c(2.139,1.873,1.887,1.739,1.734,1.509,1.695,1.74,1.811,1.954,1.624,2.204,1.508,1.784,1.902,1.743,1.863,2.028,1.39,2.187,1.65,2.334,1.037,1.509,1.971)
x4 = c(2.238,1.741,1.809,1.547,1.715,1.474,1.656,1.777,1.759,2.009,1.657,1.846,1.458,1.811,1.606,1.794,1.869,2.032,1.324,2.087,1.378,2.225,1.268,1.422,1.869)
x5 = c(0.873,0.59,0.767,0.706,0.549,0.782,0.737,0.618,0.853,0.823,0.686,0.678,0.662,0.81,0.723,0.586,0.672,0.836,0.578,0.758,0.533,0.757,0.546,0.618,0.869)
x6 = c(0.872,0.744,0.713,0.674,0.654,0.571,0.803,0.682,0.777,0.765,0.668,0.546,0.595,0.819,0.677,0.541,0.752,0.805,0.61,0.718,0.482,0.731,0.615,0.664,0.868)
X = cbind(x1,x2,x3,x4,x5,x6)
T.ellipse(X)
