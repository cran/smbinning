#' Optimal Binning for Scoring Modeling
#'
#' \strong{Optimal Binning} categorizes a numeric characteristic into bins for ulterior usage in scoring modeling.
#' This process, also known as \emph{supervised discretization}, 
#' utilizes \href{http://cran.r-project.org/package=partykit}{Recursive Partitioning} to categorize 
#' the numeric characteristic.\cr
#' The especific algorithm is Conditional Inference Trees 
#' which initially excludes missing values (\code{NA}) to compute the cutpoints, adding them back later in the 
#' process for the calculation of the \emph{Information Value}.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @param x Continuous characteristic. At least 10 different values. Value \code{Inf} is not allowed.
#' Name of \code{x} must not have a dot.
#' @param p Percentage of records per bin. Default 5\% (0.05). 
#' This parameter only accepts values greater that 0.00 (0\%) and lower than 0.50 (50\%).
#' @return The command \code{smbinning} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'  
#' # Package application
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' result$ivtable # Tabulation and Information Value
#' result$iv # Information value
#' result$bands # Bins or bands
#' result$ctree # Decision tree from partykit
smbinning = function(df,y,x,p=0.05){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target variable is numeric
    return("Column name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Column name with a dot [.]")
    } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (tolower(y)=="default"){
    return("Field name 'default' not allowed")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (p<=0 | p>0.5){
    return("p must be greater than 0 and lower than 0.5 (50%)")
  } else if (!is.numeric(df[,j])){
    return("Characteristic (x) not found or it is not a number")
  } else if (length(unique(df[,j]))<10){
    return("Uniques values of x < 10")  
  } else { 
    ctree=ctree(formula(paste(y,"~",x)),
                data=df, 
                na.action=na.exclude,
                control=ctree_control(minbucket=ceiling(round(p*nrow(df)))))
    bins=width(ctree)
    if (bins<2){return("No significant splits")}
    # Append cutpoinstop()ts in a table (Automated)
    cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(ctree) # Number of nodes
    for (i in 1:n) {
      cutvct=rbind(cutvct,ctree[i]$node$split$breaks)
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    cutvct=ifelse(cutvct<0,trunc(10000*cutvct)/10000,ceiling(10000*cutvct)/10000) # Round to 4 dec. to avoid borderline cases
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Empty table
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '<= $cutpoint' as Cutpoint,
                  NULL as CntRec,
                  NULL as CntGood,
                  NULL as CntBad,
                  sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                  sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                  sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    cutpoint=max(df[,j],na.rm=T) # Calculte Max without Missing
    cutpoint=ifelse(cutpoint<0,trunc(10000*cutpoint)/10000,ceiling(10000*cutpoint)/10000) # Round to 4 dec. to avoid borderline cases
    maxcutpoint=max(cutvct) # Calculte Max cut point
    mincutpoint=min(df[,j],na.rm=T) # Calculte Min without Missing for later usage
    mincutpoint=ifelse(mincutpoint<0,trunc(10000*mincutpoint)/10000,ceiling(10000*mincutpoint)/10000) # Round to 4 dec. to avoid borderline cases 
    ivt=rbind(ivt,
              fn$sqldf(
                "select '> $maxcutpoint' as Cutpoint,
                NULL as CntRec,
                NULL as CntGood,
                NULL as CntBad,
                sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                NULL as PctRec,
               NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $x is not NULL and $y is not NULL")
              )
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } 
    
     else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))}
    
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert to table numeric
    options(warn=-1)
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    options(warn=0)

    # Complete Table 
    ivt[1,2]=ivt[1,5] # Nbr Records
    ivt[1,3]=ivt[1,6] # Nbr Goods
    ivt[1,4]=ivt[1,7] # Nbr Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,2]=ivt[i,5]-ivt[i-1,5]
                   ivt[i,3]=ivt[i,6]-ivt[i-1,6]
                   ivt[i,4]=ivt[i,7]-ivt[i-1,7]}
    
    ivt[2,2]=ivt[2,5]-ivt[1,5]
    ivt[2,3]=ivt[2,6]-ivt[1,6]
    ivt[2,4]=ivt[2,7]-ivt[1,7]
    
    # Missing row.  Update: Added "if" statement
        ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]

    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  bands=append(mincutpoint,cutvct)
  bands=append(bands,cutpoint)
  list(ivtable=ivt,iv=iv,ctree=ctree,bands=bands,x=x,col_id=j,cuts=cutvct)
  }

# Begin Summary IV 20160602 ###############################################
#' Information value Summary
#'
#' It gives the user the ability to calculate, in one step, the IV for each characteristic of the dataset.
#' This function also shows a progress bar so the user can see the status of the process.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @return The command \code{smbinning.sumiv} generates a table that lists each characteristic 
#' with its corresponding IV for those where the calculation is possible, otherwise it will generate a 
#' missing value (\code{NA}).
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'  
#' # Summary IV application
#' sumivt=smbinning.sumiv(chileancredit.train,y="FlagGB")
#' sumivt # Display table with IV by characteristic

smbinning.sumiv = function(df,y){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data not a data.frame")}  
  ncol=ncol(df)
  sumivt=data.frame(matrix(ncol=0,nrow=0)) # Empty table
  options(warn=-1) # Turn off warnings
  cat("","\n")
  pb = txtProgressBar(min = 0, max = 1, initial = 0, style=3, char = "-",width=50)
  # t0=Sys.time()
  # t1=0
  for (i in 1:ncol){
    smbnum=smbinning(df,y,colnames(df[i]))
    smbfac=smbinning.factor(df,y,colnames(df[i]))
    if (colnames(df[i])!=y) {
      if (is.numeric(df[,i]) & is.list(smbnum)) {sumivt=rbind(sumivt,data.frame(Char=colnames(df[i]), IV=smbnum$iv, Process="Numeric binning OK"))}
      else if (is.numeric(df[,i]) & !is.list(smbnum)) {sumivt=rbind(sumivt,data.frame(Char=colnames(df[i]), IV=NA, Process=smbnum))}
      else if (is.factor(df[,i]) & is.list(smbfac)) {sumivt=rbind(sumivt,data.frame(Char=colnames(df[i]), IV=smbfac$iv, Process="Factor binning OK"))}
      else if (is.factor(df[,i]) & !is.list(smbfac)) {sumivt=rbind(sumivt,data.frame(Char=colnames(df[i]), IV=NA,Process=smbfac))}
      else {sumivt=rbind(sumivt,data.frame(Char=colnames(df[i]),IV=NA,Process="Not numeric nor factor"))}
      # t1=Sys.time() # Recall system time
      # t1=round(t1-t0,2) # Compares against starting time t0
    }
    setTxtProgressBar(pb, i/ncol)
    # if(i<ncol) 
    # {cat(" | Time:",t1,"| Binning",substring(colnames(df[i]),1,10),"...")} 
    # else {cat(" | Time:",t1,"| Binning Done           ")} 
  }
  close(pb)
  options(warn=0) # Turn back on warnings
  sumivt=sumivt[with(sumivt,order(-IV)),]
  cat("","\n")
  return(sumivt)
}

# End Summary IV 20160602 #################################################


# Begin Plot Summary IV 20160602 ##########################################
#' Plot Information Value Summary
#'
#' It gives the user the ability to plot the Information Value by characteristic.
#' The chart only shows characteristics with a valid IV.
#' @param sumivt A data frame saved after \code{smbinning.sumiv}.
#' @param cex Optional parameter for the user to control the font size of the characteristics
#' displayed on the chart. The default value is 0.9
#' @return The command \code{smbinning.sumiv.plot} returns a plot that shows the IV
#' for each numeric and factor characteristic in the dataset.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'  
#' # Plotting smbinning.sumiv
#' sumivt=smbinning.sumiv(chileancredit.train,y="FlagGB")
#' sumivt # Display table with IV by characteristic
#' smbinning.sumiv.plot(sumivt,cex=0.8) # Plot IV summary table

smbinning.sumiv.plot=function(sumivt, cex=0.9){
  if (!is.data.frame(sumivt)){ # Check if data.frame
    return("Data not a data.frame")
  } else if (names(sumivt[1])!="Char" | names(sumivt[2])!="IV") {
    return("Not from smbinning.sumiv")}
  sumivtplot=sumivt
  sumivtplot=sumivtplot[complete.cases(sumivtplot$IV),]
  sumivtplot=sumivtplot[order(sumivtplot$IV),]
  sumivtplot=cbind(sumivtplot,Desc=ifelse(sumivtplot$IV>=0.3,"1:Strong",ifelse(sumivtplot$IV>=0.1,"2:Medium","3:Weak")))
  unique(sumivtplot$Desc)
  sumivtplot$Desc=as.factor(sumivtplot$Desc)
  smbsumivplot=dotchart(sumivtplot$IV, 
                        main="Information Value",
                        labels=sumivtplot$Char,
                        pch=ifelse(sumivtplot$IV>=0.3,21,ifelse(sumivtplot$IV>=0.1,21,1)),
                        color=ifelse(sumivtplot$IV>=0.3,"black",ifelse(sumivtplot$IV>=0.1,"black","black")),
                        bg=ifelse(sumivtplot$IV>=0.3,"black",ifelse(sumivtplot$IV>=0.1,"gray75","white")),
                        groups=sumivtplot$Desc,
                        cex=cex)
}

# End Plot Summary IV 20160602 ############################################


# Begin Exploratory Data Analysis 20160602 ################################
#' Exploratory Data Analysis (EDA)
#'
#' It shows basic statistics for each numeric, integer, and factor characteristic in a data frame.
#' @param df A data frame.
#' @param rounding Optional parameter to define the decimal points shown in the output table. Default is 3.
#' @param pbar Optional parameter that turns on or off a progress bar. Default value is 1 (On).
#' @return The command \code{smbinning.eda} generates two data frames that list each characteristic 
#' with basic statistics such as extreme values and quartiles;
#' and also percentages of missing values and outliers, among others. 
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'  
#' # EDA application
#' smbinning.eda(chileancredit.train,rounding=3)$eda # Table with basic statistics.
#' smbinning.eda(chileancredit.train,rounding=3)$edapct # Table with basic percentages.

smbinning.eda = function(df, rounding=3, pbar=1){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data not a data.frame")}  
  ncol=ncol(df)
  nrow=nrow(df)
  r=rounding
  eda=data.frame(matrix(ncol=0,nrow=0)) # Empty table
  options(scipen=999) # No scientific notation
  options(warn=-1) # Turn off warnings
  if (pbar==1){
    cat("","\n")
    pb = txtProgressBar(min = 0, max = 1, initial = 0, style=3, char = "-",width=50)
  }
  for (i in 1:ncol){
    # t1=round(Sys.time()-t0,2)
    Miss=sum(is.na(df[,i]))
    if (is.numeric(df[,i]) | is.integer(df[,i])) {
      q=unname(quantile(df[,i], na.rm=T))
      iqr=q[4]-q[2]
      iqrlow=q[2]-1.5*iqr
      iqrupp=q[4]+1.5*iqr
      Avg=mean(df[,i],na.rm=T)
      eda=rbind(eda,data.frame(Field=colnames(df[i]),
                               Type="Num/Int",
                               Recs=nrow,
                               Miss=Miss,
                               Unique=length(unique((df[,i][!is.na(df[,i])]))),
                               Min=round(q[1],r),
                               Q25=round(q[2],r),
                               Q50=round(q[3],r),
                               Avg=round(Avg,r), 
                               Q75=round(q[4],r),
                               Max=round(q[5],r),
                               StDv=round(sd(df[,i],na.rm=T),r),
                               Neg=nrow(subset(df, df[,i]<0 & !is.na(df[,i]))),
                               Zero=nrow(subset(df, df[,i]==0 & !is.na(df[,i]))),
                               Pos=nrow(subset(df, df[,i]>0 & !is.na(df[,i]))),
                               OutLo=nrow(subset(df,df[,i]<iqrlow)),
                               OutHi=nrow(subset(df,df[,i]>iqrupp))
      ))}
    else if (is.factor(df[,i])) {
      eda=rbind(eda,data.frame(Field=colnames(df[i]),
                               Type="Factor",
                               Recs=nrow,
                               Miss=Miss,
                               Unique=length(unique((df[,i][!is.na(df[,i])]))),
                               Min=NA,
                               Q25=NA,
                               Q50=NA,
                               Avg=NA,
                               Q75=NA,
                               Max=NA,
                               StDv=NA,
                               Neg=NA,
                               Zero=NA,
                               Pos=NA,
                               OutLo=NA,
                               OutHi=NA
      ))}
    else {
      eda=rbind(eda,data.frame(Field=colnames(df[i]),
                               Type="Other",
                               Recs=nrow,
                               Miss=Miss,
                               Unique=length(unique((df[,i][!is.na(df[,i])]))),
                               Min=NA,
                               Q25=NA,
                               Q50=NA,
                               Avg=NA,
                               Q75=NA,
                               Max=NA,
                               StDv=NA,
                               Neg=NA,
                               Zero=NA,
                               Pos=NA,
                               OutLo=NA,
                               OutHi=NA
      ))}
    if (pbar==1){setTxtProgressBar(pb, i/ncol)}
  }
  
  if (pbar==1){close(pb)}
  
  # Table with percentages (edapct)
  edapct=cbind(eda[,1:4],eda[,13:17])
  edapct[,4]=round(edapct[,4]/edapct[,3],r)
  edapct[,5]=round(edapct[,5]/edapct[,3],r)
  edapct[,6]=round(edapct[,6]/edapct[,3],r)
  edapct[,7]=round(edapct[,7]/edapct[,3],r)
  edapct[,8]=round(edapct[,8]/edapct[,3],r)
  edapct[,9]=round(edapct[,9]/edapct[,3],r)
  
  options(warn=0) # Turn back on warnings
  list(eda=eda,edapct=edapct)
}

# End Exploratory Data Analysis 20160602 ##################################


# Begin Custom Cutpoints 20150307 #########################################
#' Customized Binning
#'
#' It gives the user the ability to create customized cutpoints. In Scoring Modeling, the analysis
#' of a characteristic usually begins with intervals with the same length to understand its distribution,
#' and then intervals with the same proportion of cases to explore bins with a reasonable sample size.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @param x Continuous characteristic. At least 10 different values. Value \code{Inf} is not allowed. 
#' Name of \code{x} must not have a dot.
#' @param cuts Vector with the cutpoints selected by the user. It does not have a default so user must define it. 
#' @return The command \code{smbinning.custom} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#'  
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#'
#' # Remove exclusions from chileancredit dataset
#' TOB.train=
#'   subset(chileancredit,(FlagSample==1 & (FlagGB==1 | FlagGB==0)), select=TOB)
#' TOB.test=
#'   subset(chileancredit,(FlagSample==0 & (FlagGB==1 | FlagGB==0)), select=TOB)
#'   
#' # Custom cutpoints using percentiles (20% each)
#' TOB.Pct20=quantile(TOB.train, probs=seq(0,1,0.2), na.rm=TRUE)
#' TOB.Pct20.Breaks=as.vector(quantile(TOB.train, probs=seq(0,1,0.2), na.rm=TRUE))
#' Cuts.TOB.Pct20=TOB.Pct20.Breaks[2:(length(TOB.Pct20.Breaks)-1)]
#' 
#' # Package application and results
#' result=
#'   smbinning.custom(df=chileancredit.train,
#'                    y="FlagGB",x="TOB",cuts=Cuts.TOB.Pct20) # Run and save
#' result$ivtable # Tabulation and Information Value
smbinning.custom = function(df,y,x,cuts){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target vable is numeric
    return("Column name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Column name with a dot [.]")
  } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (tolower(y)=="default"){
    return("Field name 'default' not allowed")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (!is.numeric(df[,j])){
    return("Characteristic (x) not found or it is not a number")
  } else if (length(unique(df[,j]))<10){
    return("Uniques values of x < 10")  
  } else { 
    # Append cutpoints in a table (Automated)
    cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cuts) # Number of cutpoints
    if (n<1){return("No Bins")} # At least 1 cutpoint
    for (i in 1:n) {
      cutvct=rbind(cutvct,cuts[i])
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    cutvct=ifelse(cutvct<0,trunc(10000*cutvct)/10000,ceiling(10000*cutvct)/10000) # Round to 4 dec. to avoid borderline cases
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '<= $cutpoint' as Cutpoint,
                  NULL as CntRec,
                  NULL as CntGood,
                  NULL as CntBad,
                  sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                  sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                  sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    cutpoint=max(df[,j],na.rm=T) # Calculte Max without Missing
    cutpoint=ifelse(cutpoint<0,trunc(10000*cutpoint)/10000,ceiling(10000*cutpoint)/10000) # Round to 4 dec. to avoid borderline cases
    maxcutpoint=max(cutvct) # Calculte Max cut point
    mincutpoint=min(df[,j],na.rm=T) # Calculte Min without Missing for later usage
    mincutpoint=ifelse(mincutpoint<0,trunc(10000*mincutpoint)/10000,ceiling(10000*mincutpoint)/10000) # Round to 4 dec. to avoid borderline cases 
    ivt=rbind(ivt,
              fn$sqldf(
                "select '> $maxcutpoint' as Cutpoint,
                NULL as CntRec,
                NULL as CntGood,
                NULL as CntBad,
                sum(case when $x <= $cutpoint and $y in (1,0) then 1 else 0 end) as CntCumRec,
                sum(case when $x <= $cutpoint and $y=1 then 1 else 0 end) as CntCumGood,
                sum(case when $x <= $cutpoint and $y=0 then 1 else 0 end) as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $x is not NULL and $y is not NULL")
              )
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    }
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert to table numeric
    options(warn=-1)
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    options(warn=0)
    
    # Complete Table 
    ivt[1,2]=ivt[1,5] # Nbr Records
    ivt[1,3]=ivt[1,6] # Nbr Goods
    ivt[1,4]=ivt[1,7] # Nbr Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,2]=ivt[i,5]-ivt[i-1,5]
                   ivt[i,3]=ivt[i,6]-ivt[i-1,6]
                   ivt[i,4]=ivt[i,7]-ivt[i-1,7]}
    
    ivt[2,2]=ivt[2,5]-ivt[1,5]
    ivt[2,3]=ivt[2,6]-ivt[1,6]
    ivt[2,4]=ivt[2,7]-ivt[1,7]
    
    # Missing row
    ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]
    
    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  bands=append(mincutpoint,cutvct)
  bands=append(bands,cutpoint)
  list(ivtable=ivt,iv=iv,bands=bands,x=x,col_id=j,cuts=cutvct)
  }
# End Custom Cutpoints 20150307 ###############################################


# Begin Binning Factors 20150407 #############################################
#' Binning on Factor Variables
#'
#' It generates the output table for the uniques values of a given factor variable.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x A factor variable with at least 2 different values. Value \code{Inf} is not allowed. 
#' @param maxcat Specifies the maximum number of categories.  Default value is 10.
#' Name of \code{x} must not have a dot.
#' @return The command \code{smbinning.factor} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen.factor}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' str(chileancredit) # Quick description of the data
#' table(chileancredit$FlagGB) # Tabulate target variable
#' 
#' # Training and testing samples (Just some basic formality for Modeling) 
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' 
#' # Package application and results
#' result.train=smbinning.factor(df=chileancredit.train,
#'                                y="FlagGB",x="IncomeLevel")
#' result.train$ivtable
#' result.test=smbinning.factor(df=chileancredit.test,
#'                                y="FlagGB",x="IncomeLevel")
#' result.test$ivtable
#' 
#' # Plots
#' par(mfrow=c(2,2))
#' smbinning.plot(result.train,option="dist",sub="Income Level (Tranining Sample)")
#' smbinning.plot(result.train,option="badrate",sub="Income Level (Tranining Sample)")
#' smbinning.plot(result.test,option="dist",sub="Income Level (Test Sample)")
#' smbinning.plot(result.test,option="badrate",sub="Income Level (Test Sample)")

smbinning.factor = function(df,y,x,maxcat=10){
  # Check data frame and formats
  if (!is.data.frame(df)){ # Check if data.frame
    return("Data not a data.frame")
  } else if (is.numeric(y) | is.numeric(x)){ # Check if target vable is numeric
    return("Column name not string")
  } else if (grepl("[.]",y) | grepl("[.]",x)){ # Check if there is a dot
    return("Column name with a dot [.]")
  } else 
    i=which(names(df)==y) # Find Column for dependant
  j=which(names(df)==x) # Find Column for independant
  if (!is.numeric(df[,i])){ 
    return("Target (y) not found or it is not numeric")
  } else if (max(df[,i],na.rm=T)!=1){
    return("Maximum not 1")
  } else if (tolower(y)=="default"){
    return("Field name 'default' not allowed")
  } else if (fn$sqldf("select count(*) from df where cast($x as text)='Inf' or cast($x as text)='-Inf'")>0){
    return("Characteristic (x) with an 'Inf' value (Divided by Zero). Replace by NA")  
  } else if (min(df[,i],na.rm=T)!=0){
    return("Minimum not 0")
  } else if (!is.factor(df[,j])){
    return("Characteristic (x) not found or it is not a factor")
  } else if (length(unique(df[,j]))<=1){
    return("Characteristic (x) requires at leats 2 uniques categories")  
  } else if (length(unique(df[,j]))>maxcat){
    return("Too many categories")    
    } else { 
    # Append cutpoints in a table (Automated)
    # cutvct=data.frame(matrix(ncol=0,nrow=0)) # Shell
    cutvct=c()
    cuts=fn$sqldf("select distinct $x from df where $x is not NULL")
    cuts=as.vector(as.matrix(cuts))
    n=length(cuts) # Number of cutpoints
    if (n<1){return("No Bins")} # At least 1 cutpoint
    for (i in 1:n) {
      cutvct=rbind(cutvct,cuts[i])
    }
    cutvct=cutvct[order(cutvct[,1]),] # Sort / converts to a ordered vector (asc)
    # Build Information Value Table #############################################
    # Counts per not missing cutpoint
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(cutvct) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=cutvct[i]
      ivt=rbind(ivt,
                fn$sqldf(
                  "select '= ''$cutpoint''' as Cutpoint,
                  sum(case when $x = '$cutpoint' and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x = '$cutpoint' and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x = '$cutpoint' and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $x is not NULL and $y is not NULL")
                )
    }
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    y.na=fn$sqldf("select count(*) from df where $y is null")
    if(x.na>0){
      ivt=rbind(ivt,
                fn$sqldf(
                  "select 'Missing' as Cutpoint,
                  sum(case when $x is NULL and $y in (1,0) then 1 else 0 end) as CntRec,
                  sum(case when $x is NULL and $y=1 then 1 else 0 end) as CntGood,
                  sum(case when $x is NULL and $y=0 then 1 else 0 end) as CntBad,
                  NULL as CntCumRec,
                  NULL as CntCumGood,
                  NULL as CntCumBad,
                  NULL as PctRec,
                  NULL as GoodRate,
                  NULL as BadRate,
                  NULL as Odds,
                  NULL as LnOdds,
                  NULL as WoE,
                  NULL as IV
                  from df where $y is not NULL")
                )
    } else {
      ivt=rbind(ivt,
                c("Missing",0,0,0,NA,NA,NA,NA,NA,NA,NA,NA,NA))
    }
    # Total
    ivt=rbind(ivt,
              fn$sqldf(
                "select 'Total' as Cutpoint,
                count(*) as CntRec,
                sum(case when $y=1 then 1 else 0 end) as CntGood,
                sum(case when $y=0 then 1 else 0 end) as CntBad,
                NULL as CntCumRec,
                NULL as CntCumGood,
                NULL as CntCumBad,
                NULL as PctRec,
                NULL as GoodRate,
                NULL as BadRate,
                NULL as Odds,
                NULL as LnOdds,
                NULL as WoE,
                NULL as IV
                from df where $y is not NULL")
              )
    
    # Covert table to numeric
    options(warn=-1)
    ncol=ncol(ivt)
    for (i in 2:ncol){
      ivt[,i]=as.numeric(ivt[,i])
    }
    options(warn=0)

    # Complete Table: 1st row 
    ivt[1,5]=ivt[1,2] # Nbr Cum. Records
    ivt[1,6]=ivt[1,3] # Nbr Cum. Goods
    ivt[1,7]=ivt[1,4] # Nbr Cum. Bads
    
    # From 2nd row
    n=nrow(ivt)-2
    for (i in 2:n){ivt[i,5]=ivt[i,2]+ivt[i-1,5]
                   ivt[i,6]=ivt[i,3]+ivt[i-1,6]
                   ivt[i,7]=ivt[i,4]+ivt[i-1,7]}
    
    ivt[2,5]=ivt[2,2]+ivt[1,5]
    ivt[2,6]=ivt[2,3]+ivt[1,6]
    ivt[2,7]=ivt[2,4]+ivt[1,7]
    
    # Missing row
    ivt[i+1,5]=ivt[i,5]+ivt[i+1,2]
    ivt[i+1,6]=ivt[i,6]+ivt[i+1,3]
    ivt[i+1,7]=ivt[i,7]+ivt[i+1,4]
    
    # Calculating metrics
    options(scipen=999) # Remove Scientific Notation
    ivt[,8]=round(ivt[,2]/ivt[i+2,2],4) # PctRec
    ivt[,9]=round(ivt[,3]/ivt[,2],4) # GoodRate
    ivt[,10]=round(ivt[,4]/ivt[,2],4) # BadRate
    ivt[,11]=round(ivt[,3]/ivt[,4],4) # Odds
    ivt[,12]=round(log(ivt[,3]/ivt[,4]),4) # LnOdds
    G=ivt[i+2,3]
    B=ivt[i+2,4]
    LnGB=log(G/B) # IV Part 1
    ivt[,13]=round(log(ivt[,3]/ivt[,4])-LnGB,4) # WoE
    ivt[,14]=round(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),4) # Mg IV
    # ivt[i+2,14]=round(sum(ivt[,13]*(ivt[,3]/G-ivt[,4]/B),na.rm=T),4) -- Old Calculation
    # Calculates Information Value even with undefined numbers
    ivt[i+2,14]=0.0000
    for (k in 1:(nrow(ivt)-1))
    {
      if(is.finite(ivt[k,14])) {mgiv=ivt[k,14]} else {mgiv=0.0000}
      ivt[i+2,14]=ivt[i+2,14]+mgiv
    }
    iv=ivt[i+2,14]
    # End Inf. Value Table ###################################################### 
    }
  list(ivtable=ivt,iv=iv,x=x,col_id=j,cuts=cutvct)
  }

# End Binning Factors 20150407 #################################################


# Begin Plotting ##############################################################
#' Plots after binning
#'
#' It generates plots for distribution, bad rate, and weight of evidence after running \code{smbinning} 
#' and saving its output.
#' @param ivout An object generated by binning.
#' @param option Distribution ("dist"), Good Rate ("goodrate"), Bad Rate ("badrate"), and Weight of Evidence ("WoE").
#' @param sub Subtitle for the chart (optional).
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Plots
#' par(mfrow=c(2,2))
#' boxplot(chileancredit.train$TOB~chileancredit.train$FlagGB,
#'         horizontal=TRUE, frame=FALSE, col="lightgray",main="Distribution")
#' mtext("Time on Books (Months)",3)
#' smbinning.plot(result,option="dist",sub="Time on Books (Months)")
#' smbinning.plot(result,option="badrate",sub="Time on Books (Months)")
#' smbinning.plot(result,option="WoE",sub="Time on Books (Months)")
smbinning.plot=function(ivout,option="dist",sub=""){
  r=ifelse(ivout$ivtable[nrow(ivout$ivtable)-1,2]==0,2,1)
  if (option=="dist"){
    # Distribution
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,8])*1.25
    ch_dist=barplot(ivout$ivtable[1:x_upper,8],
                    names.arg=ivout$ivtable[1:x_upper,1],
                    axes=F, 
                    main="Percentage of Cases", 
                    ylim=c(0,y_upper),
                    col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_dist,y=ivout$ivtable[1:x_upper,8], label=round(ivout$ivtable[1:x_upper,8]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="goodrate"){
    # Good Rate
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,9],na.rm=T)*1.25
    ch_goodrate=barplot(ivout$ivtable[1:x_upper,9],
                        names.arg=ivout$ivtable[1:x_upper,1],
                        axes=F, 
                        main="Good Rate (%)",
                        ylim=c(0,y_upper),
                        col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_goodrate,y=ivout$ivtable[1:x_upper,9], label=round(ivout$ivtable[1:x_upper,9]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="badrate"){
    # Bad Rate
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,10],na.rm=T)*1.25
    ch_badrate=barplot(ivout$ivtable[1:x_upper,10],
                       names.arg=ivout$ivtable[1:x_upper,1],
                       axes=F, 
                       main="Bad Rate (%)",
                       ylim=c(0,y_upper),
                       col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_badrate,y=ivout$ivtable[1:x_upper,10], label=round(ivout$ivtable[1:x_upper,10]*100,1), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else if (option=="WoE") {
    # WoE
    x_upper=nrow(ivout$ivtable)-r
    y_upper=max(ivout$ivtable[1:x_upper,13],na.rm=T)*1.25
    y_lower=min(ivout$ivtable[1:x_upper,13],na.rm=T)*1.25
    ch_woe=barplot(ivout$ivtable[1:x_upper,13],
                   names.arg=ivout$ivtable[1:x_upper,1],
                   axes=F, 
                   main="Weight of Evidence", 
                   ylim=c(y_lower,y_upper),
                   col=gray.colors(length(unique(ivout$ivtable[1:x_upper,1]))))
    text(x=ch_woe,y=ivout$ivtable[1:x_upper,13], label=round(ivout$ivtable[1:x_upper,13],2), pos=3,cex=1)
    abline(h=0)
    mtext(sub,3)
  } else {
    return("Options are dist, goodrate, badrate, or WoE")
  }
}
# End Plotting ################################################################

# Begin Gen Characteristic #####################################################
#' Utility to generate a new characteristic from a numeric variable
#'
#' It generates a data frame with a new predictive characteristic after the binning process.
#' @param df Dataset to be updated with the new characteristic.
#' @param ivout An object generated after \code{smbinning}.
#' @param chrname Name of the new characteristic.
#' @return A data frame with the binned version of the characteristic analyzed with \code{smbinning}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Generate new binned characteristic into a existing data frame
#' chileancredit.train=
#' smbinning.gen(chileancredit.train,result,"gTOB") # Update training sample
#' chileancredit=
#'   smbinning.gen(chileancredit,result,"gTOB") # Update population
#' sqldf("select gTOB,count(*) as Recs 
#'       from chileancredit group by gTOB") # Check new field counts 
smbinning.gen=function(df,ivout,chrname="NewChar"){
  df=cbind(df,tmpname=NA)
  ncol=ncol(df)
  
  # Update 2016-08-20 MP
  #col_id=ivout$col_id
  col_id = which(names(df)==ivout$x) 
  
  # Updated 20160130
  b=ivout$bands
  df[,ncol][is.na(df[,col_id])]=0 # Missing
  df[,ncol][df[,col_id]<=b[2]]=1 # First valid
  # Loop goes from 2 to length(b)-2 if more than 1 cutpoint
  if (length(b)>3) {
  for (i in 2:(length(b)-2)) {
    df[,ncol][df[,col_id]>b[i] & df[,col_id]<=b[i+1]]=i
  }
  }
  df[,ncol][df[,col_id]>b[length(b)-1]]=length(b)-1 # Last
  df[,ncol]=as.factor(df[,ncol]) # Convert to factor for modeling
  blab=c(paste("01 <=",b[2]))
  if (length(b)>3) {
  for (i in 3:(length(b)-1)) {
    blab=c(blab,paste(sprintf("%02d",i-1),"<=",b[i]))
  }
  } else {i=2}
  blab=c(blab,paste(sprintf("%02d",i),">",b[length(b)-1]))
  
  # Are there ANY missing values
  # any(is.na(df[,col_id]))
  
  if (any(is.na(df[,col_id]))){
    blab=c("00 Miss",blab)
  }
  df[,ncol]=factor(df[,ncol],labels=blab)
  
  names(df)[names(df)=="tmpname"]=chrname
  return(df)
}
# End Gen Characteristic #######################################################


# Begin Gen Characteristic for factor variables ################################
#' Utility to generate a new characteristic from a factor variable
#'
#' It generates a data frame with a new predictive characteristic from a factor variable after the binning process.
#' @param df Dataset to be updated with the new characteristic.
#' @param ivout An object generated after \code{smbinning.factor}.
#' @param chrname Name of the new characteristic.
#' @return A data frame with the binned version of the characteristic analyzed with \code{smbinning.factor}.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=
#' smbinning.factor(df=chileancredit.train,y="FlagGB",x="IncomeLevel")
#' result$ivtable
#' 
#' # Generate new binned characteristic into a existing data frame
#' chileancredit=
#'   smbinning.factor.gen(chileancredit,result,"gInc") # Update population

smbinning.factor.gen=function(df,ivout,chrname="NewChar"){
  df=cbind(df,tmpname=NA)
  ncol=ncol(df)
  col_id=ivout$col_id
  # Updated 20160523
  b=ivout$cuts
  df[,ncol][is.na(df[,col_id])]=0 # Missing
  # df[,ncol][df[,col_id]<=b[2]]=1 # First valid
  # Loop through all factor values
  for (i in 1:length(b)) {
    df[,ncol][df[,col_id]==b[i]]=i
  }
  df[,ncol]=as.factor(df[,ncol]) # Convert to factor for modeling
  
  #Labeling
  blab=c(paste("01 =  '",b[1],"'"))
  for (i in 2:length(b)) {
    blab=c(blab,paste(sprintf("%02d",i),"=  '",b[i],"'"))
  }
  
  # Are there ANY missing values
  # any(is.na(df[,col_id]))
  
  if (any(is.na(df[,col_id]))){
    blab=c("00 Miss",blab)
  }

  # Some Make Up
  blab=gsub(" '","'",blab)
  blab=gsub("' ","'",blab)
  
  df[,ncol]=factor(df[,ncol],labels=blab)
  
  names(df)[names(df)=="tmpname"]=chrname
  return(df)
}
# End Gen Characteristic for factor variables #################################


# Begin: SQL Code #############################################################
#' SQL Code
#'
#' It outputs a SQL code to facilitate the generation of new binned characetristic 
#' in a SQL environment.
#' @param ivout An object generated by \code{smbinning}.
#' @return A text with the SQL code for binning.
#' @examples
#' # Package loading and data exploration
#' library(smbinning) # Load package and its data
#' data(chileancredit) # Load smbinning sample dataset (Chilean Credit)
#' chileancredit.train=subset(chileancredit,FlagSample==1)
#' chileancredit.test=subset(chileancredit,FlagSample==0)
#' result=smbinning(df=chileancredit.train,y="FlagGB",x="TOB",p=0.05) # Run and save result
#' 
#' # Generate SQL code
#' smbinning.sql(result)
smbinning.sql=function(ivout){
  lines=nrow(ivout$ivtable)-2
  sqlcodetable=as.list(matrix(ncol=0,nrow=0))
  sqlcodetable=rbind("case")
  for (k in 1:lines){
    sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,ivout$ivtable[k,1],"then","'",sprintf("%02d",k),":",ivout$x,gsub("'","",ivout$ivtable[k,1]),"'"))
  }
  sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,"Is Null then","'",ivout$x,"Is Null'"))
  sqlcodetable=rbind(sqlcodetable,paste("else '99: Error' end"))
  sqlcodetable
  # Some Make Up
  sqlcodetable=gsub(" '","'",sqlcodetable)
  sqlcodetable=gsub("' ","'",sqlcodetable)
  sqlcodetable=gsub("then","then ",sqlcodetable)
  sqlcodetable=gsub("'then","' then ",sqlcodetable)
  sqlcodetable=gsub("then  ","then ",sqlcodetable)
  sqlcodetable=gsub("else","else ",sqlcodetable)
  sqlcodetable=gsub("'end","' end ",sqlcodetable)
  sqlcodetable=gsub(" :",":",sqlcodetable)
  sqlcodetable=gsub("='","= '",sqlcodetable)
  return(gsub(",","",toString(sqlcodetable)))
}
# End: SQL Code ###############################################################


# Begin: Chilean Credit Data ##################################################
#' Chilean Credit Data 
#'
#' A simulated dataset based on six months of information collected by a Chilean Bank
#' whose objective was to develop a credit scoring model to determine the probability
#' of default within the next 12 months. The target 
#' variable is FlagGB, which represents the binary status of default (0) and not default(1).
#'
#' \itemize{
#'   \item CustomerId. Customer Identifier.
#'   \item TOB. Time on books in months since first account was open.
#'   \item IncomeLevel. Income level from 0 (Low) to 5 (High).
#'   \item Bal. Outstanding balance.
#'   \item MaxDqBin. Max. delinquency bin. 0:No Dq., 1:1-29 ... 6:150-179.
#'   \item MtgBal. Mortgage outstanding balance at the Credit Bureau.
#'   \item NonBankTradesDq. Number of non-bank delinquent trades.
#'   \item FlagGB. 1: Good, 0: Bad.
#'   \item FlagSample. Training and testing sample indicator (1:75\%,0:25\%).
#'   }
#'
#'
#' @format Data frame with 7,702 rows and 19 columns.
#' @name chileancredit
NULL
# End: Chilean Credit Data ####################################################
