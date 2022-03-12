# beta一起估，轉換Y
rm(list=ls())         # 清除資料
library("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path)) # Rstudio改路徑
source("GEE41.r")       # 主程式需要與函數路徑一致
seed = 1234
set.seed(seed)
start = Sys.time()    # 時間測試
times = 200           # 模擬次數
n = 400               # 樣本數
period = 1            # 期數
D.break = TRUE        # 參數是否需要容許範圍 
D.R1.spatial = TRUE   # R1是否加入相關性矩陣估計
breaks = 10           # 參數容許範圍(倍數)
kk = 200              # 將距離分成幾份
tol.initial = 0.05    # beta起始值收斂範圍
tol.beta = 0.05       # beta收斂最大差值
tol.spatial = 0.05    # spatial收斂最大差值
tol.max.spatial = 500 # spatial容許最大數字
m.max = 20            # 迭代容許最大次數
D.jackknife = T       # 是否需要使用jackknife方法
K = 40                # jackknife切K組
tol.K = 0.8           # jackknife每K組須成功的比例
tol.D.times = 0       # 每次模擬，執行jackknife的次數
# (tol.D.times = 1，如果失敗，重新一次隨機分組，再做jackknife法)
# (如果參數估計小於想要組數，給系統重來幾次的設定)
D.print.table = T     # 模擬時輸出表格

## parameter for beta1
par1 = c(c1 = 0.5,tau1 = 2,kappa1 = 1,ak1 = 0.4)
par2 = c(c2 = 0.5,tau2 = 2,kappa2 = 1,ak2 = 0.8)
q1 = 1
q2 = 1
# c:時間相關性,tau:nugget,kappa:sill,ak:range的係數,q:power

## parameter for 
B1 = matrix(data = c(b10 = -0.3,b11 = -0.3,b12 = 0.5,b13 = 0.5),ncol = 1)
B2 = matrix(data = c(b20 =  0.3,b21 =  0.3,b22 = 0.3,b23 = 0.3),ncol = 1)
BB = rbind(B1,B2)
#---------------------------------------------------------------------------------
M = MM = matrix(data = NA,nr = length(BB),nc = times)
Q = QQ = matrix(data = NA,nr = length(BB),nc = times)
Zero = c()
Y.simulation = b.simulation = s.simulation = list()
all.s = all.b = list()
jackknife.indep = jackknife.spatial = list()
count.error.beta = count.error.spatial = 0
count.error.divergent = count.error.other = 0
a = 1
while(a <= times){
  
  Indicator = 1
  tryCatch(
    expr = {
      
      #======================================================================
      # 相關係數矩陣
      repeat{
        
        # 生成距離
        value = n*10
        S = cbind(sample(x = c(1:value),size = n,replace = F)/value,
                  sample(x = c(1:value),size = n,replace = F)/value)
        SSS = c()
        for( i in c(1:period)) SSS = rbind(SSS,S)
        distance = distance(x = SSS,method = "euclidean",mute.message = TRUE)
        Hmax = max(distance) # 最大距離
        par1[length(par1)] = par1[length(par1)] * Hmax
        par2[length(par2)] = par2[length(par2)] * Hmax 
        
        # 來自結構零之隨機變數的相關係數矩陣
        SIGMA1 = sigma.spatial(period = period,D = distance,par = par1,q = q1)
        
        # 一般卜瓦松之隨機變數的相關係數矩陣
        SIGMA2 = sigma.spatial(period = period,D = distance,par = par2,q = q2)
        det(SIGMA1);det(SIGMA2)
        
        break
        #if(is.singular.matrix(SIGMA1,tol = 1e-100)*
        #   is.singular.matrix(SIGMA2,tol = 1e-100) == 0)break
      }
      
      # 生成零膨脹卜瓦松之隨機變數
      zip = zeropois(n,B1,B2,SIGMA1,SIGMA2)
      YY   = zip$Y
      XX1  = zip$X1
      XX2  = zip$X2
      #=====================================================================
      Y  = YY
      X1 = XX1
      X2 = XX2
      SS = SSS

      # 參數估計  
      m = 1
      beta00 = matrix(data = c(BB/2),ncol = 1)
      if( period == 1 ){
        s10 = matrix(data = c(1,0.5,0.5*Hmax),ncol = 1)
        s20 = matrix(data = c(1,0.5,0.5*Hmax),ncol = 1)
      }else{
        s10 = matrix(data = c(0.5,1,0.5,0.5*Hmax),ncol = 1)
        s20 = matrix(data = c(0.5,1,0.5,0.5*Hmax),ncol = 1)
      }
      final.beta = beta00
      final.spatial = rbind(s10,s20)
      repeat{
        
        #==================================================
        if(max(final.spatial[,ncol(final.spatial)]) > tol.max.spatial){
          count.error.spatial = count.error.spatial + 1
          stop("\n","spatial is not under tolerance level")
        }
        
        if( D.print.table == TRUE ){
          
          cat("\014")
          cat("\n")
          cat("==========================================================","\n")
          information = c(times = a-1,
                          beta = count.error.beta,
                          spatial = count.error.spatial,
                          diverge = count.error.divergent,
                          all = count.error.other)
          if( a == 1 ){
            table = cbind(true = c(B1,B2),
                          indep = NA,spatial = NA,
                          sample = NA,sample = NA,
                          jackknife = NA,jackkife = NA)
          }else{
            table = cbind(true = c(B1,B2),
                          indep   = round(apply(M[ ,c(1:a)],1,mean,na.rm = T),4),
                          spatial = round(apply(MM[ ,c(1:a)],1,mean,na.rm = T),4),
                          sample = round(apply(M[ ,c(1:a)],1,var,na.rm = T),4),
                          sample = round(apply(MM[ ,c(1:a)],1,var,na.rm = T),4),
                          jackknife = round(apply(Q[ ,c(1:a)],1,mean,na.rm = T),4),
                          jackknife = round(apply(QQ[ ,c(1:a)],1,mean,na.rm = T),4))
          }
          rownames(table) = c(paste("B",rep(x = c(10:13),times = 1),sep = ""),
                              paste("B",rep(x = c(20:23),times = 1),sep = ""))
          show(information)
          show(table)
          cat("==========================================================","\n")
          cat(paste0("\n",
                     " Number:"," ",a,"-",m,
                     "\n",
                     "Process: ",a,"/",times," (",round(a/times*100,2),"% completed)",
                     "\n"))
          show(Sys.time()-start)
          
        }else{
          
          cat("\014")
          cat("\n")
          cat(paste0("\n",
                     " Number:"," ",a,"-",m,
                     "\n",
                     "Process: ",a,"/",times," (",round(a/times*100,2),"% completed)",
                     "\n"))
          show(Sys.time()-start)
          
        }
        
        #==================================================
        if(m == 1){ # 用working.independent估計beta1、beta2之起始值
          
          R1 = diag(x = 1,nrow = length(Y),ncol = length(Y))
          R2 = diag(x = 1,nrow = length(Y),ncol = length(Y))
          tol = tol.initial
          
        }else{ # 更新空間相關性矩陣估計beta1、beta2
          
          if( D.R1.spatial == TRUE ){
            
            SIGMA.beta1 = sigma.spatial(period = period,D = distance,
                                        par = s10,q = q1)
            R1 = SIGMA.beta1/s10[length(s10)-1]
            
          } 
          SIGMA.beta2 = sigma.spatial(period = period,D = distance,
                                      par = s20,q = q2)
          R2 = SIGMA.beta2/s20[length(s20)-1]
          tol = 100
          
        }
        
        result = newton(b = c(beta00),
                        Y = Y,X1 = X1,X2 = X2,
                        R1 = R1,R2 = R2,
                        stop.num = tol,D.break = D.break)
        (beta10 = result$beta)
        (z = result$z)
        
        # 估計空間相關性參數
        (s = spatial(par1 = s10,par2 = s20,result = result))
        (s11 = matrix(data = c(s$s11),ncol = 1))
        (s21 = matrix(data = c(s$s21),ncol = 1))
        final.beta = cbind(final.beta,beta10)
        final.spatial = cbind(final.spatial,rbind(s11,s21))
        dif.beta = max(abs(final.beta[ ,ncol(final.beta)] - 
                             final.beta[ ,(ncol(final.beta)-1)]))
        dif.spatial = max(abs(final.spatial[ ,ncol(final.spatial)] - 
                                final.spatial[ ,(ncol(final.spatial)-1)]))
        if( dif.beta < tol.beta & dif.spatial < tol.spatial ) break
        
        if( m == m.max ){
          
          count.error.divergent = count.error.divergent + 1
          stop("\n","beta is divergent")
          
        }else{
          
          m = m + 1
          beta00 = beta10
          s20 = s21
          if( D.R1.spatial == TRUE ) s10 = s11
          
        }
      } # 參數估計
      final.beta = final.beta[ ,-1] # 扣除起始值
      final.spatial = final.spatial[ ,-1] # 扣除起始值
      if( D.jackknife == TRUE){
        V2 = result$V2
        #LL = t(chol(V2))
        #Y.trans = solve(LL) %*% YY
      }
    }, # expr
    warning = function(msg){
      message("\n","Original warning message:")
      message(paste0(msg))
    }, # warning
    error = function(msg){
      message("\n","Original error message:")
      message(paste0(msg))
      assign(x = "Indicator",value = 0,envir = .GlobalEnv)
      assign(x = "count.error.other",value = {count.error.other + 1},envir = .GlobalEnv)
    } # error
  ) # trycatch
  
  if(Indicator == 1){
    
    Zero[a] = sum(Y == 0)/length(Y)
    Y.simulation[[a]] = Y
    b.simulation[[a]] = final.beta
    s.simulation[[a]] = final.spatial
    M[,a] = final.beta[ ,1]
    MM[,a] = final.beta[ ,ncol(final.beta)]
    if( D.jackknife == TRUE ){
      assign(x = "Indicator",value = 2,envir = .GlobalEnv) # 繼續往下做
    }else{
      assign(x = "a",value = {a + 1},envir = .GlobalEnv) # 此次模擬結束
    }
    
  } # 1
  
  # jackknife
  if(Indicator == 2){
    
    D.times = 0 # 第一次creatfolds
    repeat{
      
      # 隨機選擇刪除資料
      splits = createFolds(y = YY[1:n],k = K,list = TRUE)
      jk.splits = list()
      for( i in c(1:length(splits))){
        jk.splits[[i]] = rep(x = splits[[i]],times = period) + 
          n*(rep(x = c(0:(period-1)),each = length(splits[[i]])))
      }
      jk.beta.indep = matrix(data = NA,nrow = length(BB),ncol = length(jk.splits))
      jk.beta.spatial = jk.beta.indep
      
      for(k in c(1:length(jk.splits))){
      
        if(Indicator == 0) break # jackknife失敗，結束迴圈
        
        tryCatch(
          expr = {
            
            m = 1
            # 更新變數
            X1 = XX1[-jk.splits[[k]],]
            X2 = XX2[-jk.splits[[k]],]
            Y  = YY [-jk.splits[[k]]]
            SS = SSS[-jk.splits[[k]],]
            distance = distance(x = SS,method = "euclidean",mute.message = TRUE)
            beta00 = matrix(data = c(BB/2),ncol = 1)
            if( period == 1 ){
              s10 = matrix(data = c(1,0.5,0.5*Hmax),ncol = 1)
              s20 = matrix(data = c(1,0.5,0.5*Hmax),ncol = 1)
            }else{
              s10 = matrix(data = c(0.5,1,0.5,0.5*Hmax),ncol = 1)
              s20 = matrix(data = c(0.5,1,0.5,0.5*Hmax),ncol = 1)
            }
            final.beta = beta00
            final.spatial = rbind(s10,s20)
            repeat{
              
              #==============================================================
              if(max(final.spatial[,ncol(final.spatial)]) > tol.max.spatial){
                count.error.spatial = count.error.spatial + 1
                stop("\n","spatial is not under tolerance level")
              }
              
              # 判斷
              (success = sum(is.na(jk.beta.spatial[1,c(1:k)]) == 0))
              (fail = sum(is.na(jk.beta.spatial[1,c(1:k)]) == 1) - 1)
              (last = length(jk.splits) - ( success + fail ))
              if( (success + last) < (K * tol.K) ) {
                D.times = D.times + 1
                assign(x = "Indicator",value = {0},envir = .GlobalEnv) # 改0，重新模擬
                stop("\n",
                     "Original warning message:",
                     "\n",
                     "jackknife is under tolerance level")
              }
              
              if( D.print.table == TRUE ){
                
                cat("\014")
                cat("\n")
                cat("==========================================================","\n")
                information = c(times = a-1,
                                beta = count.error.beta,
                                spatial = count.error.spatial,
                                diverge = count.error.divergent,
                                all = count.error.other)
                if( a == 1 ){
                  
                  table = cbind(true = c(B1,B2),
                                indep = NA,spatial = NA,
                                sample = NA,sample = NA,
                                jackknife = NA,jackkife = NA)
                  
                }else{
                  
                  table = cbind(true = c(B1,B2),
                                indep = round(apply(M[ ,c(1:a)],1,mean,na.rm = T),4),
                                spatial = round(apply(MM[ ,c(1:a)],1,mean,na.rm = T),4),
                                sample = round(apply(M[ ,c(1:a)],1,var,na.rm = T),4),
                                sample = round(apply(MM[ ,c(1:a)],1,var,na.rm = T),4),
                                jackknife = round(apply(Q[ ,c(1:a)],1,mean,na.rm = T),4),
                                jackknife = round(apply(QQ[ ,c(1:a)],1,mean,na.rm = T),4))
                  
                }
                rownames(table) = c(paste("B",rep(x = c(10:13),times = 1),sep = ""),
                                    paste("B",rep(x = c(20:23),times = 1),sep = ""))
                show(information)
                show(table)
                cat("==========================================================","\n")
                cat(paste0("\n",
                           " Number:"," ",a,"-",k,"-",m," (",success,"-",fail,"-",K-(success+fail),")",
                           "\n",
                           "Process: ",a,"/",times," (",round(a/times*100,2),"% completed)",
                           "\n"))
                show(Sys.time()-start)
                
              }else{
                
                cat("\014")
                cat("\n")
                cat(paste0("\n"," Number:"," ",a,"-",m,"\n",
                           "Process: ",a,"/",times," (",round(a/times*100,2),
                           "% completed)","\n"))
                show(Sys.time()-start)
                
              }
              
              #==============================================================
              
              if(m == 1){ # 用working.independent估計beta1、beta2之起始值
                
                R1 = diag(x = 1,nrow = nrow(X1),ncol = nrow(X1))
                R2 = diag(x = 1,nrow = nrow(X1),ncol = nrow(X1))
                tol = tol.initial
                
              }else{ # 更新空間相關性矩陣估計beta1、beta2(c,tau,sigma,phi)
                
                if( D.R1.spatial == TRUE ){
                  
                  SIGMA.beta1 = sigma.spatial(period = period,
                                              D = distance,
                                              par = s10,
                                              q = q1)
                  R1 = SIGMA.beta1/s10[length(s10)-1]
                  
                }
                SIGMA.beta2 = sigma.spatial(period = period,
                                            D = distance,
                                            par = s20,
                                            q = q2)
                R2 = SIGMA.beta2/s20[length(s20)-1]
                #Lk = t(chol(V2[-jk.splits[[k]],-jk.splits[[k]]]))
                #Y  = (Lk) %*% Y.trans [-jk.splits[[k]]]
                tol = 100
                
              }
              result = newton(b = c(beta00),Y = Y,
                              X1 = X1,X2 = X2,R1 = R1,R2 = R2,
                              stop.num = tol,D.break = D.break)
              (beta10 = result$beta)
              (z = result$z)
              
              # 估計空間相關性參數
              (s = spatial(par1 = s10,par2 = s20,result = result))
              (s11 = matrix(data = c(s$s11),ncol = 1))
              (s21 = matrix(data = c(s$s21),ncol = 1))
              final.beta = cbind(final.beta,beta10)
              final.spatial = cbind(final.spatial,rbind(s11,s21))
              dif.beta = max(abs(final.beta[ ,ncol(final.beta)] - 
                                   final.beta[ ,(ncol(final.beta)-1)]))
              dif.spatial = max(abs(final.spatial[ ,ncol(final.spatial)] - 
                                      final.spatial[ ,(ncol(final.spatial)-1)]))
              if( dif.beta < tol.beta & dif.spatial < tol.spatial ) break
              
              if( m == m.max ){
                
                count.error.divergent = count.error.divergent + 1
                stop("\n","beta is divergent")
                
              }else{
                
                m = m + 1
                beta00 = beta10
                s20 = s21
                if( D.R1.spatial == TRUE ) s10 = s11
                
              }
            } # 參數估計
            jk.beta.indep[ ,k] = final.beta[ ,2]
            jk.beta.spatial[ ,k] = final.beta[ ,ncol(final.beta)]
          },
          warning = function(msg){NULL},
          error   = function(msg){
            message("\n","Original warning message:")
            message(paste0(msg,"\n"))
            assign(x = "count.error.other",
                   value = {count.error.other + 1},
                   envir = .GlobalEnv)
            return(NULL)
          }
        ) # trycatch
      } # for
      
      # 判斷
      jk.ncol = sum(is.na(jk.beta.spatial[1,]) == 0)
      if( jk.ncol >= K * tol.K ){ # jackknife成功
        
        assign(x = "Indicator",value = {3},envir = .GlobalEnv) # 改3
        break 
        
      }else{
        
        if( D.times < tol.D.times ){ # jackknife失敗，重新生成一組變數
          assign(x = "Indicator",value = {2},envir = .GlobalEnv) # 改2
        }else{ 
          break
        }
      } # 判斷
    } # repeat
  } # indicater = 2
  
  if(Indicator == 3){
    
    jk.ncol = sum(is.na(jk.beta.spatial[1,]) == 0)
    if( jk.ncol != K){
      jackknife.indep[[a]]   = jk.beta.indep[ ,-c(which(is.na(jk.beta.indep[1,]) == 1))]
      jackknife.spatial[[a]] = jk.beta.spatial[ ,-c(which(is.na(jk.beta.spatial[1,]) == 1))]
    }else{
      jackknife.indep[[a]]   = jk.beta.indep
      jackknife.spatial[[a]] = jk.beta.spatial
    }
    jk.var.indep   = apply(X = jackknife.indep[[a]],MARGIN = 1,FUN = var)*(jk.ncol-1)**2/(jk.ncol)
    jk.var.spatial = apply(X = jackknife.spatial[[a]],MARGIN = 1,FUN = var)*(jk.ncol-1)**2/(jk.ncol)
    Q [ ,a] = jk.var.indep
    QQ[ ,a] = jk.var.spatial
    assign(x = "a",value = {a + 1},envir = .GlobalEnv)
    
  } # 3
  
} # while

end = Sys.time()
time = end - start
mean = cbind(true = c(BB),indep = apply(M,1,mean,na.rm = T),spatial = apply(MM,1,mean,na.rm = T))
var = cbind(sample = apply(M,1,var,na.rm = T),sample = apply(MM,1,var,na.rm = T),
            jackknife = apply(X = Q,MARGIN = 1,FUN = mean,na.rm = T),
            jackknife = apply(X = QQ,MARGIN = 1,FUN = mean,na.rm = T))
MSE = ((apply(M,1,mean,na.rm = T)-BB)^2 + apply(X = Q,MARGIN = 1,FUN = mean,na.rm = T))/
  ((apply(MM,1,mean,na.rm = T)-BB)^2 + apply(X = QQ,MARGIN = 1,FUN = mean,na.rm = T))
compare = cbind(var = log(var[,1]/var[,2]),mse = log(MSE)) # 大於0表示spatial參數比較好
rownames(mean) = c(paste("B",c(10:13),sep = ""),paste("B",c(20:23),sep = ""))
rownames(var) = rownames(mean)
rownames(compare) = rownames(mean)
colnames(compare) = c("log.var","log.mse")
ll = list(information = c(seed = seed,times = a-1,sample = n,period = period),
          limitation = c(D.break = D.break,breaks = breaks,kk = kk),
          tol = c(initial = tol.initial,beta = tol.beta,alpha = tol.spatial,m.max = m.max),
          alpha = c(par1),lambda = c(par2),
          mean = round(mean,4),var = round(var,4),compare = round(compare,4),
          prop = c(mean = round(mean(Zero),4),var = round(var(Zero),4)),
          error = c(beta = count.error.beta,spatial = count.error.spatial,
                    divergent = count.error.divergent,other = count.error.other),
          time = time)
ll

#median(Zero)
#apply(all.b,1,median)
#filename = paste("C:\\Users\\user\\Desktop\\JK\\result",Sys.Date(),".txt")
#for(i in c( 1:length(ll) )){
#  write.table(x = ll[i],file = filename,
#              sep = " ",append = TRUE,
#              row.names = FALSE,col.names = TRUE)
#}
#sink(file = "sink-examp.txt",append = FALSE)
#ll
#sink()

