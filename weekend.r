library(gdata)
library(lpSolveAPI)

setwd("R:/projects/1 weekend thing")


power_per_site<-read.csv("sites.wind.power.per.csv", header=TRUE)
wind_sites<-read.csv("wind.sites.csv", header=TRUE)
demand<-read.csv("demand.csv", header=TRUE)
gen<-read.csv("gen.info.csv", header=TRUE)
import<-read.csv("import.limits.csv", header=TRUE)
trans<-read.csv("trans.limits.csv", header=TRUE)

west_energy <- gen$ann.gen.coal.gwh[1] + gen$ann.gen.nuc.gwh[1] + gen$ann.gen.hydro.gwh[1]
west_power <- (gen$ann.gen.coal.gwh[1] + gen$ann.gen.nuc.gwh[1] + gen$ann.gen.hydro.gwh[1])*1000/8760

east_energy <- gen$ann.gen.coal.gwh[2] + gen$ann.gen.nuc.gwh[2] + gen$ann.gen.hydro.gwh[2]
east_power <- (gen$ann.gen.coal.gwh[2] + gen$ann.gen.nuc.gwh[2] + gen$ann.gen.hydro.gwh[2])*1000/8760

downstate_energy <- gen$ann.gen.coal.gwh[3] + gen$ann.gen.nuc.gwh[3] + gen$ann.gen.hydro.gwh[3]
downstate_power <- (gen$ann.gen.coal.gwh[3] + gen$ann.gen.nuc.gwh[3] + gen$ann.gen.hydro.gwh[3])*1000/8760

state_energy <- west_energy + east_energy + downstate_energy
state_power <- west_power + east_power + downstate_power

west_cap <- gen$cap.gen.nonbase.mw[1] + import$ieso.mw[1] + import$hq.mw[1] + import$pjm.mw[1] + import$ne.mw[1]
east_cap <- gen$cap.gen.nonbase.mw[2] + import$ieso.mw[2] + import$hq.mw[2] + import$pjm.mw[2] + import$ne.mw[2]
down_cap <- gen$cap.gen.nonbase.mw[3] + import$ieso.mw[3] + import$hq.mw[3] + import$pjm.mw[3] + import$ne.mw[3]

demand$state <- demand$west + demand$east + demand$downstate

state_cap <- west_cap + east_cap + down_cap

# CASE 1
ranking <- as.data.frame(matrix(0,66,1))
ranking$site <- apply(power_per_site[2:67], 2, sum) 
ranking$V1 <- seq(1:66)
ranking <- as.data.frame(ranking[order(ranking$site, decreasing = TRUE),])

ranking$power <- wind_sites$n.t.max[ranking$V1]*ranking$site

ranking$cap_fac <- ranking$site/(3*8760)
ranking <- as.data.frame(ranking[order(ranking$cap_fac, decreasing = TRUE),])

total.turbine <- c(1,1000,2000,3000,4000,5000)
turbine.placement <- as.data.frame(matrix(0,66,7))
turbine.placement[1] <- ranking$site
days <- paste0('Case',5:11)
ranking[,days] <- 0
days <- paste0('Case_power',12:16)
ranking[,days] <- 0
for(a in seq(from=6, to=11, by=1))
{
  total=total.turbine[a-5]
 
  for(i in 1:66)
  {
    
    if(total==0)
    {
      total=0
    }
    else 
    {
      if(total>wind_sites[ranking[i,1],4])
      {
        
        ranking[i,a-1]=min(wind_sites[ranking[i,1],4],total)
        ranking[i,a+5]= ranking[i,a-1]*ranking$site[i]
        total=total-ranking[i,a-1]
      }
      else
      {
        if(total>0)
        {
          
          ranking[i,a-1]=total
          ranking[i,a+5]= ranking[i,a-1]*ranking$site[i]
          total=0
          
        }
        
      }
    }
  }
}
turbine.cap <- c(3, 3000, 6000, 9000, 12000, 15000)
apply(ranking[11:16],2,sum)/(turbine.cap*8760)


#CASE 1B
a = 0
sum<- array(0,dim=c(8760,13))

for(m in 5:10)
{
  #print(paste("case",m))

  for(i in 1:8760)
  {
    counter =0
    #print(paste("hour",i))
    for(k in 1:66)
    {
      
      if(ranking[k,m]!=0)
      {
        #print("inner")
       counter=counter+ranking[k,m]*power_per_site[i,ranking[k,1]+1]
      
      }
      else
      {
        sum[i,m-4]=counter
        if(m+2<13)
        {
          if((demand[i,5]-state_power)<sum[i,m-4])
          {
            print(i)
            sum[i,m+2]<- demand[i,5]-state_power
          }
          else
          {
            sum[i,m+2]<- sum[i,m-4]
          }
        }
        
        break
      }
      
    }
   
  }
}


turbine.cap <- c(3, 3000, 6000, 9000, 12000, 15000)
apply(sum[,7:12],2,sum)/(turbine.cap*8760)

# CASE 2
sites_no <- 66
sites <- wind_sites$site
turb_cap <- c(3, 3000, 6000, 9000, 12000, 15000)
times <- 8760


for (k in 1:length(turb_cap))
  {
  wind <- turb_cap[k]
  turb_no <-wind/3
  cat(sprintf("\r\n\r\n number of turbines %s",turb_no))
  optimize<- make.lp(nrow = 1+times, ncol = sites_no +times)
  optimize
  lp.control(optimize, sense = 'min')
  set.objfn(optimize,obj=c(rep(0,sites_no),rep(1,times)))
  set.type(optimize, c(1:sites_no),type = "integer")
  row.add.mode(optimize, "on")
    for(i in 1:times){
    add.constraint(optimize, xt =c(power_per_site[i, 1 +sites],1), type= ">=", rhs =demand[i, 5]-state_power,indices=c(sites,(sites_no+i)))
  }
  add.constraint(optimize, xt= rep(1,sites_no), type="<=", rhs=turb_no, indices=c(1:sites_no))
  row.add.mode(optimize, "off")
  set.bounds(optimize, lower=rep(0,sites_no + times), upper = c(wind_sites$n.t.max, rep(state_cap, times)))
 
  optimize
  solve(optimize)
 
  turbine.distribution <- get.variables(optimize)[1:sites_no]
  solution <- get.objective(optimize)
  print(solution)
  print("Distribution of Turbines")
  print(turbine.distribution)
  
  print("Total Wind Cap")
  wind_cap <- (sum(demand$state)- solution - state_power*times)/(wind*8760)
  print(wind_cap)
  
  #print("Optimized Solution")
  #print(solution)
 
  write.lp(optimize,'model.lp',type='lp')
}


