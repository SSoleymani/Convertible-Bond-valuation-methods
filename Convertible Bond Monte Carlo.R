# AFP Modeling Convertible Bonds with Correlation of Interest Rates and Equity Prices
# Ruixiang Lu
# Qichen Gao
# Siavash Soleymani
# Attika Raj

# Part I Modelling stock prices movement in correlation with interest rate movement

# Input: Stock prices S0
# Goal: Modeling the stock prices movement with volatility following Garch (1,1) process.

# clear environment and console
rm(list = ls())
cat("\014")
library(LSMonteCarlo)
### library(RQuantLib) # ###:shows RQuantLib codes

########### Functions ###########

# price bond with constant interest rate
mybond_price = function (faceval =100, couprate, discrate, maturity,freq = 2, compounding = c("same","continuous","annual")){
  
  times <- seq(to = maturity, by = (1/freq), length.out = ceiling(maturity * freq))
  if (compounding == "continuous") {
    pvfactors = exp(-discrate * times)
  }
  else if (compounding == "annual") {
    pvfactors = 1/(1 + discrate)^times
  }
  else {
    pvfactors = 1/(1 + discrate/freq)^(freq * times)
  }
  coupon <- couprate * faceval/(freq)
  cashflows <- rep(coupon, ceiling(maturity * freq))
  cashflows[length(cashflows)] = cashflows[length(cashflows)] + faceval
  price <- sum(cashflows * pvfactors)
  price = ifelse(maturity==0,faceval*(1+couprate/2),price)
  return(price)
}

# price bond with interest rate term structure as a data-frame
mybond_price_TermStructure = function (faceval =100, couprate, TermStructureDF, maturity,freq = 2, compounding = c("same","continuous","annual")){
  times <- seq(to = maturity, by = (1/freq), length.out = ceiling(maturity * freq))
  ts = approx(x = TermStructureDF$time, y = TermStructureDF$risk.free.rate, xout = times, rule = 2)
  discrate = ts$y
  if (compounding == "continuous") {
    pvfactors = exp(-discrate * times)
  }
  else if (compounding == "annual") {
    pvfactors = 1/(1 + discrate)^times
  }
  else {
    pvfactors = 1/(1 + discrate/freq)^(freq * times)
  }
  coupon <- couprate * faceval/(freq)
  cashflows <- rep(coupon, ceiling(maturity * freq))
  cashflows[length(cashflows)] = cashflows[length(cashflows)] + faceval
  price <- sum(cashflows * pvfactors)
  price = ifelse(maturity==0,faceval*(1+couprate/2),price)
  return(price)
}

mycoupon_pv = function (faceval =100, couprate, discrate, maturity,freq = 2, compounding = c("same","continuous","annual")){
  
  times <- seq(to = maturity, by = (1/freq), length.out = ceiling(maturity * freq))
  if (compounding == "continuous") {
    pvfactors = exp(-discrate * times)
  }
  else if (compounding == "annual") {
    pvfactors = 1/(1 + discrate)^times
  }
  else {
    pvfactors = 1/(1 + discrate/freq)^(freq * times)
  }
  coupon <- couprate * faceval/(freq)
  cashflows <- rep(coupon, ceiling(maturity * freq))
  # cashflows[length(cashflows)] = cashflows[length(cashflows)] + faceval
  price <- sum(cashflows * pvfactors)
  price = ifelse(maturity==0,faceval*couprate/2,price)
  return(price)
}

# clearing NA columns and naming cols as 0,0.5,1,..., T
clr.na.namecol = function(MAT){
  MAT = MAT[,colSums(is.na(MAT))<nrow(MAT)] #remove full NA columns
  colnames(MAT) = seq(0,T,length.out = ncol(MAT)) # name columns as Time in years
  return(MAT)
}

########### firm specific inputs ############
for (number in 1:nrow(conv_table)){
S0 = conv_table$Stock.Price[number]     # intial stock price
sigma = conv_table$Vol[number]   # annual volatility specific for each stock.
                 # there are three specific volatilities that we can use:
                      # 1) historical vol: in which we will use a price history as long as the term of the conv bond: 
                            # advantage: easy and straight forward, best when following two are not applicale or available.
                            # disadvantage : backward looking instead of forward looking.
                      # 2) Garch(1,1) method: 
                            # advantage: forward looking + path dependent results + flexible because of parameters.
                            # disadvantage: still not very diffrent results from historical if we use historical as the long term average: VL
                      # 3) Implied vol based on traded options:
                            # advantage: best correct answer!
                            # disadvantage: almost always no traded options are available, or if available maturities are significantly shorter! cannot extrapolate implied vol because of smile effect.
Delta=conv_table$Dividend[number]/100      # Setting intial delta as dividend rate (continously compounding)

########### Bond specific inputs ########### 
c=conv_table$Coupon[number]               # Assume 5.5% Bond Annual Coupon 
nt = conv_table$Conv.Ratio[number]         # conversion ratio
FaceVal=1000          # face value
strike = FaceVal/nt;
Exer_put=ifelse(is.na(conv_table$Put[number]) , 1, conv_table$Put[number])            # STRIKE PRICE OF PUT OPTION, SETTING TO CONSTANT NOW, NEED TO MODIFY!!!
Exer_call=ifelse(is.na(conv_table$Call[number]), 99999, conv_table$Call[number])      # STRIKE PRICE OF CALL OPTION, SETTING TO CONSTANT NOW, NEED TO MODIFY!!!
T=conv_table$matur_years[number]                   # Time to maturity of Bond 

########### Simpulation inputs ###########
dt=1/252                  # Time interval between two stock prices, decrease if path dependent features not needed, to save time and memory to simulate more paths
                          # If Conversion Decision: Daily: dt = 1/252  ; Weekly:  dt = 1/52 ; Monthly: dt = 1/12; yearly: dt = 1
h = round(T/dt)           # Total number of time steps 
paths = 10000             # Number of paths to generate

########### Interest rates: can be constant or term structure ###########
rf = approx(x = x.date,y = y.rate,xout = conv_table$matur_years[number],rule = 2)$y/100
#rf = 0.0116                # The risk free interest rate for simulating stock path under risk neutral measure and discount equity side (either flat or term structure)
credit_spread = conv_table$Spread[number]/10000    # can be either flat or term structure
r0 =  credit_spread + rf  # The total company specific discount rate used to discount cash payments of bond (either flat or term structure)

r0t=seq(from = r0, to = r0, length.out = h+1) # SETTING r AS CONSTANT FOR NOW, NEED TO MODIFY!!!
rft=seq(from = rf, to = rf, length.out = h+1) # SETTING r AS CONSTANT FOR NOW, NEED TO MODIFY!!!


########### simulation ###########

S=matrix(NA,paths, h+1)                 # Matrix of simulated stock prices
S[,1]=S0# Setting inital value

# # another way: uncomment if you wnat it!
# onescol = matrix(data = 1,nrow = paths, ncol = 1)
# ranmat = matrix(data = rnorm(paths*h),nrow = paths ,ncol = h)
# expranmat= exp((rt-0.5*sigma^2)*dt+sigma*sqrt(dt)*ranmat)
# S = S0*cbind(onescol,t(apply(X = expranmat,MARGIN = 1, FUN = cumprod)))

Norm=matrix(rnorm(paths* h ,0,1),paths) # Random normal matrix to model brownian motion

# revision needed! r0 should be r0(T) only!

for (i in 1:h){
  S[,i+1]=S[,i]+S[,i]*(rft[i]-Delta)*dt+S[,i]*Norm[,i]*sigma*sqrt(dt)
}# Simulating stock prices matrix that follows geometric brownian motion process

# # Another way: 
# S = fastGBM(Spot = S0,sigma = sigma,n = paths,m = h,r = rf,dr = Delta,mT = T)
# S = cbind(rep(S0,paths),S)


# # volatility and return checker : uncomment to check
# 
# S_ret = matrix(nrow = nrow(S),ncol= NCOL(S)-1)
# for (i in 1:ncol(S)-1){
#   S_ret[,i]= log(S[,i+1]/S[,i])
# }
# sd(S_ret[1,])/sqrt (dt) ; #should be equal to sigma
# 
# log(mean(S[,ncol(S)])/S0)/T # should be equal to r

########### Finding the Payoff at every time point to determine which action to take ###########

Vt = matrix(data = NA,nrow = nrow(S),ncol = NCOL(S)) # continuation value matrix

# at each point continuation value is equal to the price of the bond as the discounted remainning cash-flows and the embedded option
B = mybond_price(faceval = FaceVal,couprate = c,discrate = r0, maturity = T,freq = 2,compounding = "continuous" )
### O = AmericanOption(type = "call",underlying = S0,strike = strike ,dividendYield = Delta,riskFreeRate = r0,maturity = T,volatility = sigma,engine ="CrankNicolson")
### O = EuropeanOption(type = "call",underlying = S0,strike = strike ,dividendYield = Delta,riskFreeRate = r0,maturity = T,volatility = sigma)
O = EuCallBS(Spot = S0,sigma = sigma,Strike = strike,r = rf,dr = Delta,mT = T)
### Vt[,1] = B + O$value*nt
Vt[,1] = B + O*nt

for (i in 1:(2*T)){
  B = mybond_price(faceval = FaceVal,couprate = c,discrate = r0, maturity = T-i/2,freq = 2,compounding = "continuous" )
  for (j in 1:paths){
    # O = AmericanOption(type = "call",underlying = S[j,1+(1/dt)*i/2],strike = strike ,dividendYield = Delta,riskFreeRate = r0,maturity = T-i/2+0.01,volatility = sigma,engine ="CrankNicolson")
    # O = EuropeanOption(type = "call",underlying = S[j,1+(1/dt)*i/2],strike = strike ,dividendYield = Delta,riskFreeRate = r0,maturity = T-i/2+0.01,volatility = sigma)
    O = EuCallBS(Spot = S[j,1+(1/dt)*i/2],sigma = sigma,Strike = strike,r = rf,dr = Delta,mT = T-i/2)
    ### Vt[j,1+(1/dt)*i/2] = B + O$value*nt
    Vt[j,1+(1/dt)*i/2] = B + O*nt
  }
}

Vt_nodes = clr.na.namecol(Vt) # clear NAs on non-conversion dates, name columns as fraction of year

# 2 options for dicounting payoffs: (Which one is correct?)
    # 1) discount only cash payoff at last with r0 company, and nt*St with rf, which does not give the same answer as option + bond
          # mean(c(Vt_nodes[Vt_nodes[,11]<1027.6,11]*exp(-r0*T),Vt_nodes[Vt_nodes[,11]>1027.6,11]*exp(-rf*T))) 
    # 2) discount upto the cash payoff at last with r0 company, and the rest (ntSt - (F+c/2)) with rf! this is it! you just need to add PV (coupons) as well: 
          # (sum((Vt_nodes[Vt_nodes[,11]>1027.6,11]-1027.5)*exp(-rf*T))+paths*exp(-r0*T)*1027.5)/paths + PV (coupons) this is what gives the 

CVt=nt*S  # conversion value at each point
CVt_nodes = clr.na.namecol(CVt * (Vt>-1)) # only keep corresponding values on decision nodes

Kt=matrix(data = Exer_call,nrow = paths,ncol = ncol(Vt_nodes))  # Callability feature payoff at each node 
Pt=matrix(data = Exer_put,nrow = paths,ncol = ncol(Vt_nodes))   # Put option payoff at each node
REDt=FaceVal # redemption value at each point (can be kappa*FaceVal on very rare occasions; kappa is 1.1 for example at maturity)

rft_nodes = clr.na.namecol(rft*(Vt>-1))[1,]
r0t_nodes = clr.na.namecol(r0t*(Vt>-1))[1,]
Discount_rf = exp(-rft_nodes* seq(0,T,length.out = ncol(Vt_nodes)))
Discount_r0 = exp(-r0t_nodes* seq(0,T,length.out = ncol(Vt_nodes)))


M1=(CVt_nodes > Vt_nodes)   # Compare ntSt (conversion value) and Vt (continuation value)
M2=(CVt_nodes >= Pt)        # Compare ntSt (conversion value) and Pt (put payoff)
M3=(Pt > Vt_nodes)          # Compare Pt (put payoff) and Vt (continuation value)
M4=(Vt_nodes > Kt)          # Compare Vt (continuation value) and Kt (call payoff)
M5=(Kt >= CVt_nodes)       # Compare Kt (call payoff) and ntSt (conversion value)
M6=(REDt > CVt_nodes[,ncol(CVt_nodes)]) ## Compare FaceVal (redemption value) and ntSt (conversion value) at final node

C1=M1*M2      # Condition 1: Voluntary Conversion
C2=M3*(!M2)   # Condition 2: Put
C3=M4*M5      # Condition 3: Call
C4=M4*(!M5)   # Condition 4: Forced Conversion
C5=M6 


Action=matrix(data = 7,nrow(C1),ncol(C1))   # Constructing a matrix to record which action to take


# Assumption: Forced Conversion>> Call>> Put>> Voluntary Conversion>> Continuation>> Redemption
for (i in 1: (paths*ncol(C1))){
  if     (C4[i]==1) { Action[i]=1} # Forced Conversion
  else if(C3[i]==1) { Action[i]=2} # Call
  else if(C2[i]==1) { Action[i]=3} # Put
  else if(C1[i]==1) { Action[i]=4} # Voluntary Conversion
}


for(i in 1:paths){
  if(C5[i]==TRUE){Action[i,ncol(Action)]=6} #Redemption  
  else if(C5[i]==FALSE) {Action[i,ncol(Action)]=5} #conversion at end
}
# For each paths, we store the action number and the time of action
# Action number is the smallest number in each row
# Time of action is the respective column number
Final_Action = matrix(data = NA,nrow = paths,ncol = 4)
colnames(Final_Action) = c("action node","action number","payoff","PV payoff")


Final_Action[,1] = apply(Action,1,which.min)   # which node contains the min
Final_Action[,2] = apply(Action, 1, FUN = min) # Find min for each row in Matrix "Action"

for (i in 1:paths){
  if     (Final_Action[i,2]==1) {Final_Action[i,3]=CVt_nodes[i,Final_Action[i,1]]}   # If forced conversion,payoff=ntSt
  else if(Final_Action[i,2]==2) {Final_Action[i,3]=Kt[i,Final_Action[i,1]]}          # If called, payoff=Kt}
  else if(Final_Action[i,2]==3) {Final_Action[i,3]=Pt[i,Final_Action[i,1]]}          # If put, payoff=Pt
  else if(Final_Action[i,2]==4) {Final_Action[i,3]=CVt_nodes[i,Final_Action[i,1]]}   # If voluntary conversion, payoff=ntSt
  else if(Final_Action[i,2]==5) {Final_Action[i,3]=CVt_nodes[i,Final_Action[i,1]]}   # If coverted at end, payoff = nt*St
  else if(Final_Action[i,2]==6) {Final_Action[i,3]=Vt_nodes[i,Final_Action[i,1]]}    # If Redeemed, payoff = F+coupon/2
}

# put the following line in each condition of the above loop and replace proper discount rate!

for (i in 1:paths){
  Final_Action[i,4] = ifelse(Final_Action[i,2]==6 | Final_Action[i,2]==3 ,Discount_r0[Final_Action[i,1]] * Final_Action[i,3],Discount_rf[Final_Action[i,1]] * Final_Action[i,3]) + mycoupon_pv(couprate = c,discrate = r0,maturity = (Final_Action[i,1]-1)/2,freq = 2,compounding = "continuous")
}


# Part III Discounting to calculate the value of the bond at the beginning

conv_table$Price.wout.CIR[number] = mean(Final_Action[,4])/10

# result_mean = mean(Final_Action[,4])/10
# result_sd = sd(Final_Action[,4])
# result_upperB = result_mean + 1.96*result_sd/sqrt(paths)
# result_lowerB = result_mean - 1.96*result_sd/sqrt(paths)
# std_error = 1.96*result_sd/sqrt(paths)
# cat("Fair Value:", result_mean/10," with a standard error of", std_error/10)

}



