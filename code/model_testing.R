
start_time <- Sys.time()

source('function.R')
# install.packages('doParallel')
library(doParallel)

folder_name <- "../results/"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

folder_name <- "../results/Testing_result/"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}


folder_name <- "../results/Rolling_window_result/"
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}

data = read.csv('../data/factors_data.csv')


dim(data)
colnames(data)
T0 = nrow(data)

Factors = data[,-c(1,2)]

nmodel= 9

models1 = matrix(,nmodel,nmodel)
CAPM =  Factors[c('MKTRF')]
FF3 = Factors[c('MKTRF','SMB','HML')]
DHS3 = Factors[c('MKTRF','PEAD','FIN')]
HXZ4 = Factors[c('MKTRF','ME','IA','ROE')]
Q5 = Factors[c('MKTRF','ME','IA','ROE','REG')]
FF5 = Factors[c('MKTRF','SMB','HML','CMA','RMW')]
FF6 = Factors[c('MKTRF','SMB','HML','CMA','RMW','UMD')]
BS6 = Factors[c('MKTRF','SMB','IA','ROE','UMD','HMLM')]
FLWZ8 = Factors[c('MKTRF','SMB','REG','PEAD','HMLM','EPRD', 'STR', 'ILR')]

models1[1,c(1:1)] = c('MKTRF')
models1[2,c(1:3)] = c('MKTRF','SMB','HML')
models1[3,c(1:3)] = c('MKTRF','PEAD','FIN')
models1[4,c(1:4)] = c('MKTRF','ME','IA','ROE')
models1[5,c(1:5)] = c('MKTRF','ME','IA','ROE','REG')
models1[6,c(1:5)] = c('MKTRF','SMB','HML','CMA','RMW')
models1[7,c(1:6)] = c('MKTRF','SMB','HML','CMA','RMW','UMD')
models1[8,c(1:6)] = c('MKTRF','SMB','IA','ROE','UMD','HMLM')
models1[9,c(1:8)] = c('MKTRF','SMB','REG','PEAD','HMLM','EPRD', 'STR', 'ILR')


CAPM = as.matrix(CAPM)
FF3 = as.matrix(FF3)
DHS3 = as.matrix(DHS3)
HXZ4 = as.matrix(HXZ4)
Q5 = as.matrix(Q5)
FF5 = as.matrix(FF5)
FF6 = as.matrix(FF6)
BS6 = as.matrix(BS6)
FLWZ8 = as.matrix(FLWZ8)


RF1 = read.csv('../data/RF.csv')
Fama100 = read.csv('../data/Fama100.csv')

Fama100 = Fama100[,-c(1:2)]
Fama100 = Fama100 - (matrix(rep(RF1[,3], ncol(Fama100)),nrow(Fama100) , ncol(Fama100)))


Fama185 = read.csv('../data/Fama185.csv')
Fama185 = Fama185[,-c(1:2)]
Fama185 = Fama185 - (matrix(rep(RF1[,3], ncol(Fama185)),nrow(Fama185) , ncol(Fama185)))
Fama285 = cbind(Fama100, Fama185)


test_asset_name = c("Fama100",'Fama185',"Fama285") 
test_asset_name1 = list(Fama100,Fama185, Fama285) 
Model_name = c('CAPM','FF3','DHS3','HXZ4','Q5','FF5','FF6','BS6', 'FLWZ8')


T0 = nrow(CAPM)
L0 = T0 / 2

Model_data = list(CAPM, FF3, DHS3,HXZ4, Q5,FF5, FF6,BS6,FLWZ8)

for (t in 1:length(test_asset_name1)) {
  print(t)
  p_result = matrix(,length(Model_name)*3,5)
  test_asset1 = test_asset_name1[[t]]
  N = dim(test_asset1)[2]
  for(m in 1: length(Model_name)){
    
    Model1 = Model_data[[m]]
    
    F1 = Model1
    F2 = test_asset1
    F2 = as.matrix(t(F2))
    F1 = as.matrix(t(F1))
    N1 = dim(F2)[1]
    T1 = dim(F2)[2]
    result1 = Alpha_test(F1,F2,lambda=1)
    p_result[m,] = c(round(result1,3))
    
    for (z in 1:2) 
    {
print(z)
      F1 = Model1[(L0*(z-1)+1):(L0*z),]
      F2 = test_asset1[(L0*(z-1)+1):(L0*z),]
      F2 = as.matrix(t(F2))
      F1 = as.matrix(t(F1))
      N1 = dim(F2)[1]
      T1 = dim(F2)[2]
      result1 = Alpha_test(F1,F2,lambda)
      p_result[length(Model_name)*z + m,] = c(round(result1,3))

    }
  }
  
  files = '../results/Testing_result/'
  store_csv1 = paste(files,test_asset_name[t],'_testing','.csv',sep = "", collapse = "")
  
  
  rownames(p_result) = c(Model_name,Model_name,Model_name)
  colnames(p_result) = c('p(SSA)','p(MAX)','p(PY)','p(COM)','p(CGRS)')
  write.csv(p_result,store_csv1)
  
}



############### Rolling window tests

for (t in 1:length(test_asset_name1)) {
  print(t)
  test_asset1 = test_asset_name1[[t]]

  N = dim(test_asset1)[2]
  
  W0 = 240
  
  n_models <- length(Model_name)
  n_windows <- T0 - W0 + 1
  p_values <- list(
    SSA = matrix(nrow = n_windows, ncol = n_models),
    CGRS = matrix(nrow = n_windows, ncol = n_models),
    MAX = matrix(nrow = n_windows, ncol = n_models),
    PY = matrix(nrow = n_windows, ncol = n_models),
    COM = matrix(nrow = n_windows, ncol = n_models)
  )
  

  for (j in 1:length(Model_name)) {
    
    
    Model0 = Model_data[[j]]
    
    test_asset1 = as.matrix(test_asset1)
    for (i in 1:(T0 - W0+1))
    {
      print(c(j,i))
      F1 = Model0[(i):(i+ W0-1),]
      F2 = test_asset1[(i):(i+ W0-1),]
      F2 = as.matrix(t(F2))
      F1 = as.matrix(t(F1))
      result1 = Alpha_test(F1,F2,lambda)
      
      
      p_values$SSA[i, j] <- result1[1]
      p_values$CGRS[i, j] <- result1[2]
      p_values$MAX[i, j] <- result1[3]
      p_values$PY[i, j] <- result1[4]
      p_values$COM[i, j] <- result1[5]

    }

    
  }    

  
    reject_rate = matrix(,length(Model_name)*3,2*5)  
  thresholds = c(0.01, 0.05)
  for (z in 1:length(p_values)) {
    
    for (zz in 1:length(thresholds)) {
      p_vals1 = p_values[[z]]
      
      for (zzz in 1:ncol(p_vals1)) {
        reject_rate[zzz,((zz-1)*5 + z)] = length(which(p_vals1[-1,zzz] < thresholds[zz]))/length(p_vals1[-1,zzz])
        
        reject_rate[zzz + ncol(p_vals1),((zz-1)*5 + z)] = length(which(p_vals1[c(2:121),zzz] < thresholds[zz]))/120
        
        reject_rate[zzz + ncol(p_vals1)*2,((zz-1)*5 + z)] = length(which(p_vals1[c(122:241),zzz] < thresholds[zz]))/120
        
      }
    }  
  }
  
  

  store8 = paste('../results/Rolling_window_result','/','reject_rates','_','W',W0,'_N',N,'_',T0,".csv",sep="") 
  
  
  rownames(reject_rate) = c(Model_name,Model_name,Model_name)
  colnames(reject_rate) = c('SSA','MAX','PY','COM','CGRS',
                            'SSA','MAX','PY','COM','CGRS')
  
  write.csv(round(reject_rate,3),store8)
  
}


end_time <- Sys.time()
run_time <- end_time - start_time
print(paste("Runtime:", run_time))


