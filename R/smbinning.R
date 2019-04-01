#' Optimal Binning for Scoring Modeling
#'
#' \strong{Optimal Binning} categorizes a numeric characteristic into bins for ulterior usage in scoring modeling.
#' This process, also known as \emph{supervised discretization}, 
#' utilizes \href{https://cran.r-project.org/package=partykit}{Recursive Partitioning} to categorize 
#' the numeric characteristic.\cr
#' The especific algorithm is Conditional Inference Trees 
#' which initially excludes missing values (\code{NA}) to compute the cutpoints, adding them back later in the 
#' process for the calculation of the \emph{Information Value}.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @param x Continuous characteristic. At least 5 different values. Value \code{Inf} is not allowed.
#' Name of \code{x} must not have a dot.
#' @param p Percentage of records per bin. Default 5\% (0.05). 
#' This parameter only accepts values greater that 0.00 (0\%) and lower than 0.50 (50\%).
#' @return The command \code{smbinning} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example: Optimal binning
#' result=smbinning(df=smbsimdf1,y="fgood",x="cbs1") # Run and save result
#' result$ivtable # Tabulation and Information Value
#' result$iv # Information value
#' result$bands # Bins or bands
#' result$ctree # Decision tree

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
  } else if (length(unique(df[,j]))<5){
    return("Uniques values < 5")  
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

# End smbinning ###########################################################


# Begin Custom Cutpoints 20150307 #########################################
#' Customized Binning
#'
#' It gives the user the ability to create customized cutpoints.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @param x Continuous characteristic. At least 5 different values. Value \code{Inf} is not allowed. 
#' Name of \code{x} must not have a dot.
#' @param cuts Vector with the cutpoints selected by the user. It does not have a default so user must define it. 
#' @return The command \code{smbinning.custom} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Custom cutpoints using percentiles (20% each)
#' cbs1cuts=as.vector(quantile(smbsimdf1$cbs1, probs=seq(0,1,0.2), na.rm=TRUE)) # Quantiles
#' cbs1cuts=cbs1cuts[2:(length(cbs1cuts)-1)] # Remove first (min) and last (max) values
#' 
#' # Example: Customized binning
#' result=smbinning.custom(df=smbsimdf1,y="fgood",x="cbs1",cuts=cbs1cuts) # Run and save
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
  } else if (length(unique(df[,j]))<5){
    return("Uniques values < 5")  
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


# Begin Exploratory Data Analysis 20160602 ################################
#' Exploratory Data Analysis (EDA)
#'
#' It shows basic statistics for each characteristic in a data frame.
#' The report includes:
#' \itemize{
#'   \item Field: Field name.
#'   \item Type: Factor, numeric, integer, other.
#'   \item Recs: Number of records.
#'   \item Miss: Number of missing records.
#'   \item Min: Minimum value.
#'   \item Q25: First quartile. It splits off the lowest 25\% of data from the highest 75\%.
#'   \item Q50: Median or second quartile. It cuts data set in half.
#'   \item Avg: Average value.
#'   \item Q75: Third quartile. It splits off the lowest 75\% of data from the highest 25\%.
#'   \item Max: Maximum value.
#'   \item StDv: Standard deviation of a sample.
#'   \item Neg: Number of negative values.
#'   \item Pos: Number of positive values.
#'   \item OutLo: Number of outliers. Records below \code{Q25-1.5*IQR}, where \code{IQR=Q75-Q25}. 
#'   \item OutHi: Number of outliers. Records above \code{Q75+1.5*IQR}, where \code{IQR=Q75-Q25}.
#'   }
#' @param df A data frame.
#' @param rounding Optional parameter to define the decimal points shown in the output table. Default is 3.
#' @param pbar Optional parameter that turns on or off a progress bar. Default value is 1.
#' @return The command \code{smbinning.eda} generates two data frames that list each characteristic 
#' with basic statistics such as extreme values and quartiles;
#' and also percentages of missing values and outliers, among others. 
#' @examples 
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example: Exploratory data analysis of dataset
#' smbinning.eda(smbsimdf1,rounding=3)$eda # Table with basic statistics
#' smbinning.eda(smbsimdf1,rounding=3)$edapct # Table with basic percentages

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


# Begin Binning Factors 20150407 #############################################
#' Binning on Factor Variables
#'
#' It generates a table with relevant metrics for all the categories of a given factor variable.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x A factor variable with at least 2 different values. Labesl with commas are not allowed. 
#' @param maxcat Specifies the maximum number of categories.  Default value is 10.
#' Name of \code{x} must not have a dot.
#' @return The command \code{smbinning.factor} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen.factor}.
#' @examples
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#'
#' # Binning a factor variable
#' result=smbinning.factor(smbsimdf1,x="inc",y="fgood", maxcat=11)
#' result$ivtable

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
  } else if (any(grepl(",", df[,j]))){
    return("Values contain comma")
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
    # End Inf. Value Table ##################################################### 
    }
  list(ivtable=ivt,iv=iv,x=x,col_id=j,cuts=cutvct)
  }
# End Binning Factors 20150407 #################################################


# Begin Custom Binning Factors 20171117 ########################################
#' Customized Binning on Factor Variables
#'
#' It gives the user the ability to combine categories and create new attributes for a given characteristic.
#' Once these new attribues are created in a list (called \code{groups}), the funtion generates a table for 
#' the uniques values of a given factor variable. 
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot.
#' @param x A factor variable with at least 2 different values. Value \code{Inf} is not allowed. 
#' @param groups Specifies customized groups created by the user.
#' Name of \code{x} must not have a dot.
#' @return The command \code{smbinning.factor.custom} generates an object containing the necessary information 
#' and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen.factor}.
#' @examples 
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example: Customized binning for a factor variable
#' # Notation: Groups between double quotes
#' result=smbinning.factor.custom(
#'   smbsimdf1,x="inc",
#'   y="fgood",
#'   c("'W01','W02'",        # Group 1
#'     "'W03','W04','W05'",  # Group 2
#'     "'W06','W07'",        # Group 3
#'     "'W08','W09','W10'")) # Group 4
#' result$ivtable

smbinning.factor.custom = function(df,y,x,groups){
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
  } else if (any(grepl(",", df[,j]))){
    return("Values contain comma")
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
  } else {
    ivt=data.frame(matrix(ncol=0,nrow=0)) # Shell
    n=length(groups) # Number of cutpoits
    for (i in 1:n) {
      cutpoint=groups[i]
      statement=paste("select '",gsub(",","/",gsub("'","",cutpoint)),"' as Cutpoint,
                      sum(case when $x in ($cutpoint) and $y in (1,0) then 1 else 0 end) as CntRec,
                      sum(case when $x in ($cutpoint) and $y=1 then 1 else 0 end) as CntGood,
                      sum(case when $x in ($cutpoint) and $y=0 then 1 else 0 end) as CntBad,
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
                      from df where $x is not NULL and $y is not NULL", sep="")
      ivt=rbind(ivt,fn$sqldf(statement))
    }
    # Missing Data
    x.na=fn$sqldf("select count(*) from df where $x is null")  
    #  y.na=fn$sqldf("select count(*) from df where $y is null")
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
    # End Inf. Value Table #####################################################
  }
  list(ivtable=ivt,iv=iv,x=x,col_id=j,groups=groups)
  }
# End Custom Binning Factors 20171117 ##########################################


# Begin Gen Characteristic for factor variables ################################
#' Utility to generate a new characteristic from a factor variable
#'
#' It generates a data frame with a new predictive characteristic from a factor variable after applying
#' \code{smbinning.factor} or \code{smbinning.factor.custom}.
#' @param df Dataset to be updated with the new characteristic.
#' @param ivout An object generated after \code{smbinning.factor} or \code{smbinning.factor.custom}.
#' @param chrname Name of the new characteristic.
#' @return A data frame with the binned version of the original characteristic.
#' @examples 
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' pop=smbsimdf1 # Set population
#' train=subset(pop,rnd<=0.7) # Training sample
#' 
#' # Binning a factor variable on training data
#' result=smbinning.factor(train,x="home",y="fgood")
#' 
#' # Example: Append new binned characteristic to population
#' pop=smbinning.factor.gen(pop,result,"g1home")
#' 
#' # Split training
#' train=subset(pop,rnd<=0.7) # Training sample
#' 
#' # Check new field counts
#' table(train$g1home)
#' table(pop$g1home)

# Updated 20170910
smbinning.factor.gen=function(df,ivout,chrname="NewChar"){
  df=cbind(df,tmpname=NA)
  ncol=ncol(df)
  col_id=ivout$col_id
  df[,ncol][is.na(df[,col_id])]=0 # Missing
  
  # Loop through all factor values
  # If smbinning.factor
  if(is.null(ivout$groups)){
    b=ivout$cuts
    for (i in 1:length(b)) {
      df[,ncol][df[,col_id]==b[i]]=i
    }
  } 
  # If smbinning.factor.custom
  else {
    for (i in 1:length(ivout$groups)) {
      gelements=as.list(strsplit(as.character(gsub("'","",ivout$groups[i])), ",")[[1]])
      df[,ncol][df[,col_id] %in% gelements]=i
    }
  }
  
  # Convert to factor for modeling
  df[,ncol]=as.factor(df[,ncol])
  
  # Labeling: If smbinning.factor
  if(is.null(ivout$groups)){
    blab=c(paste("01 =  '",b[1],"'"))
    for (i in 2:length(b)) {
      blab=c(blab,paste(sprintf("%02d",i),"=  '",b[i],"'"))
    }
  } 
  # Labeling: If smbinning.factor.custom
  else {
    blab=c(paste("01 in (",ivout$groups[1],")"))
    for (i in 2:length(ivout$groups)) {
      blab=c(blab,paste(sprintf("%02d",i),"in (",ivout$groups[i],")"))
    }
  }
  
  # If for any reason #bins in test sample are different, error
  if(any(is.na(df[,col_id]))==F & length(blab)>length(unique(df[,ncol])))
  {stop("Number of bins in dataset different from original result.\n  Likely due to splitting population in training/testing sample.")}
  
  if(any(is.na(df[,col_id]))==T & length(blab)>=length(unique(df[,ncol])))
  {stop("Number of bins in dataset different from original result.\n  Likely due to splitting population in training/testing sample.")}
  
    # Are there ANY missing values
  # any(is.na(df[,col_id]))
  
  if (any(is.na(df[,col_id]))){
    blab=c("00 Miss",blab)
  }
  
  # Some Make Up
  blab=gsub(" '","'",blab)
  blab=gsub("' ","'",blab)
  
  df[,ncol]=factor(df[,ncol],labels=blab) # Here is the error
  
  names(df)[names(df)=="tmpname"]=chrname
  return(df)
}
# End Gen Characteristic for factor variables #################################


# Begin Gen Characteristic #####################################################
#' Utility to generate a new characteristic from a numeric variable
#'
#' It generates a data frame with a new predictive characteristic after applying
#' \code{smbinning} or \code{smbinning.custom}.
#' @param df Dataset to be updated with the new characteristic.
#' @param ivout An object generated after \code{smbinning} or \code{smbinning.custom}.
#' @param chrname Name of the new characteristic.
#' @return A data frame with the binned version of the original characteristic.
#' @examples
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' pop=smbsimdf1 # Set population
#' train=subset(pop,rnd<=0.7) # Training sample
#' 
#' # Binning application for a numeric variable
#' result=smbinning(df=train,y="fgood",x="dep") # Run and save result
#' 
#' # Generate a dataset with binned characteristic
#' pop=smbinning.gen(pop,result,"g1dep")
#' 
#' # Check new field counts
#' table(pop$g1dep)

smbinning.gen=function(df,ivout,chrname="NewChar"){
  df=cbind(df,tmpname=NA)
  ncol=ncol(df)
  col_id=ivout$col_id
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
  
  # Labels computed with training sample
  blab=c(blab,paste(sprintf("%02d",i),">",b[length(b)-1]))
  
  # If for any reason #bins in test sample are different, error
  if(any(is.na(df[,col_id]))==F & length(blab)>length(unique(df[,ncol])))
  {stop("Number of bins in dataset different from original result.\n  Likely due to splitting population in training/testing sample.")}
  
  if(any(is.na(df[,col_id]))==T & length(blab)>=length(unique(df[,ncol])))
  {stop("Number of bins in dataset different from original result.\n  Likely due to splitting population in training/testing sample.")}
  
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


# Ini Logit Rank 20190329 #####################################################
#' Logistic Regression Ranking
#'
#' It runs all the possible logistic models for a given set of characteristics (\code{chr}) and then rank them
#' from highest to lowest performance based on AIC.
#' Important Note: This function may take time depending on the datset size and number of variables used in it.
#' The user should run it at the end of the modeling process once variables have been pre-selected in previous steps.   
#' @param df Data frame.
#' @param y Binary dependent variable.
#' @param chr Vector with the characteristics (independent variables).
#' @return The command \code{smbinning.logitrank} returns a table with the combination of characteristics 
#' and their corresponding AIC and deviance. The table is ordered by AIC from lowest (best) to highest.
#' @examples 
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example: Best combination of characteristics
#' smbinning.logitrank(y="fgood",chr=c("chr1","chr2","chr3"),df=smbsimdf3)

smbinning.logitrank = function(y,chr,df){
  f=c() # Initialize empty list of formulas
  att=c() # Initialize empty list of characteristics in each formula
  for(k in 1:length(chr)) {
    v=t(combn(chr,k))
    nrow=nrow(v)
    ncol=ncol(v)
    fnext=c() # Empty list for 1 set of combinations
    attnext=c() # Empty list for 1 set of characteristics
    for (j in 1:nrow){
      ftmp=paste0(y," ~ ",v[j,1])
      atttmp=v[j,1]
      if (ncol>1){
        for (i in 2:ncol){
          ftmp=paste0(ftmp,paste0("+",c(v[j,])[i]))
          atttmp=paste0(atttmp,paste0("+",c(v[j,])[i]))
        } # End columns
      } # End if more than 1 column
      fnext=c(ftmp,fnext)
      attnext=c(atttmp,attnext)
    } # End rows
    f=c(f,fnext)
    att=c(att,attnext)
  } 
  # List attributes
  chrsum=data.frame(character(0),numeric(0),numeric(0))
  model=glm(paste0(y," ~ 1"),family=binomial(link='logit'),data=df) # Intercept Only
  chrsum=rbind(chrsum,cbind(c("Intercept Only"),c(model$aic),c(model$deviance)))
  for(i in 1:length(f)) {
    model=glm(f[i],family=binomial(link='logit'),data=df)
    chrsum=rbind(chrsum,cbind(c(att[i]),c(model$aic),c(model$deviance)))
  }
  colnames(chrsum)=c("Characteristics","AIC","Deviance")
  chrsum$AIC=as.numeric(as.character(chrsum$AIC))
  chrsum$Deviance=as.numeric(as.character(chrsum$Deviance))
  chrsum=chrsum[order(chrsum$AIC),] 
  
  return(chrsum)
}


# Ini Logit Rank 20190329 #####################################################


# Begin Metrics 20171009 ######################################################
#' Performance Metrics for a Classification Model
#'
#' It computes the classic performance metrics of a scoring model, including AUC, KS and all the relevant ones
#' from the classification matrix at a specific threshold or cutoff.
#' @param dataset Data frame.
#' @param prediction Classifier. A value generated by a classification model (Must be numeric).
#' @param actualclass Binary variable (0/1) that represents the actual class (Must be numeric).
#' @param cutoff Point at wich the classifier splits (predicts) the actual class (Must be numeric). 
#' If not specified, it will be estimated by using the maximum value of Youden J (Sensitivity+Specificity-1).
#' If not found in the data frame, it will take the closest lower value.
#' @param report Indicator defined by user. 1: Show report (Default), 0: Do not show report.
#' @param plot Specifies the plot to be shown for overall evaluation. It has three options: 'auc' shows the ROC curve, 
#' 'ks' shows the cumulative distribution of the actual class and its maximum difference (KS Statistic), and 'none' (Default).
#' @param returndf Option for the user to save the data frame behind the metrics. 1: Show data frame, 0: Do not show (Default).
#' @return The command \code{smbinning.metrics} returns a report with classic performance metrics of a classification model.
#' @examples 
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example: Metrics Credit Score 1
#' smbinning.metrics(dataset=smbsimdf1,prediction="cbs1",actualclass="fgood",
#'                   report=1) # Show report
#' smbinning.metrics(dataset=smbsimdf1,prediction="cbs1",actualclass="fgood",
#'                   cutoff=600, report=1) # User cutoff
#' smbinning.metrics(dataset=smbsimdf1,prediction="cbs1",actualclass="fgood",
#'                   report=0, plot="auc") # Plot AUC
#' smbinning.metrics(dataset=smbsimdf1,prediction="cbs1",actualclass="fgood",
#'                   report=0, plot="ks") # Plot KS
#'
#' # Save table with all details of metrics
#' cbs1metrics=smbinning.metrics(
#'   dataset=smbsimdf1,prediction="cbs1",actualclass="fgood",
#'   report=0, returndf=1) # Save metrics details

smbinning.metrics=function(dataset,prediction,actualclass,cutoff=NA,report=1,plot="none", returndf=0){
  i=which(names(dataset)==actualclass) # Find Column for actualclass
  j=which(names(dataset)==prediction) # Find Column for prediction
  if(!is.data.frame(dataset)){ # Check if data.frame
    return("Data not a data.frame.")
  } else if (!is.na(cutoff) & !is.numeric(cutoff)){
    return("'cutoff' must be numeric.")
  } else if(!is.numeric(dataset[,which(names(dataset)==prediction)])){
    return("'prediction' not found.")
  } else if(max(dataset[,i],na.rm=TRUE)!=1 | min(dataset[,i],na.rm=TRUE)!=0){
    return("'actualclass' must be binary (0/1).")
  } else if(length(unique(na.omit(dataset[,c(actualclass)])))!=2){
    return("'actualclass' must be binary (0/1).")
  } else if(report!=1 & report!=0){ 
    return("'report' must be 0 (Deactivated) or 1 (Activated).")
  } else if(returndf!=1 & returndf!=0){ 
    return("'df' must be 0 (Deactivated) or 1 (Activated).")
  } else if(plot!="auc" & plot!="ks" & plot!="none"){ 
    return("'plot' options are: 'auc', 'ks' or 'none'.")
  } else if(!is.na(cutoff) & (max(dataset[,j],na.rm=TRUE)<cutoff | min(dataset[,j],na.rm=TRUE)>cutoff)){
   return("'cutoff' out of range.")
    }
  else {
    
    # Create table and its basic structure
    dataset=dataset[,c(prediction,actualclass)] # Only Score and class
    nmiss=nrow(dataset)-nrow(na.omit(dataset)) # Missing rows
    dataset=na.omit(dataset) # Complete Cases Only
    df=table(dataset[,c(prediction)],dataset[,c(actualclass)]) # Group by prediction
    df=as.data.frame.matrix(df) # Prediction becomes rowname 
    df$Prediction=rownames(df) # Bring rowname as a column with values
    rownames(df)=NULL
    names(df)[names(df)=="0"]="CntBad"
    names(df)[names(df)=="1"]="CntGood"
    df=df[,c("Prediction","CntGood","CntBad")] #Reorder
    df$CntTotal=df$CntBad+df$CntGood
    
    # Cumulative Ascending
    df$CumAscGood=cumsum(df$CntGood)
    df$CumAscBad=cumsum(df$CntBad)
    df$CumAscTotal=cumsum(df$CntTotal)
    
    # Totals for upcoming calculations
    SumRecords=sum(df$CntTotal)
    SumGoods=sum(df$CntGood)
    SumBads=sum(df$CntBad)
    
    # Cumulative Descending
    df$CumDescGood=NA
    c=which(names(df)=="CumDescGood")
    df[1,c]=sum(df$CntGood)
    for (i in 2:nrow(df)) {
      df[i,c]=df[1,c]-df[i-1,5]
    }
    
    df$CumDescBad=NA
    c=which(names(df)=="CumDescBad")
    df[1,c]=sum(df$CntBad)
    for (i in 2:nrow(df)) {
      df[i,c]=df[1,c]-df[i-1,6]
    }
    
    df$CumDescTotal=NA
    c=which(names(df)=="CumDescTotal")
    df[1,c]=sum(df$CntTotal)
    for (i in 2:nrow(df)) {
      df[i,c]=df[1,c]-df[i-1,7]
    }
    
    # Cumulative Percent (Column/Vertical Calculation)
    df$PctCumDescTotal=df$CumDescTotal/SumRecords
    
    # Good/Bad Rates (Row/Horizontal Calculation)
    df$GoodRateRow=df$CntGood/df$CntTotal
    df$BadRateRow=df$CntBad/df$CntTotal
    df$GoodRateDesc=df$CumDescGood/df$CumDescTotal
    df$BadRateDesc=df$CumDescBad/df$CumDescTotal
    
    # Cumulative Percentage
    df$PctCumAscGood=df$CumAscGood/SumGoods
    df$PctCumAscBad=df$CumAscBad/SumBads
    
    # Remove ascending stats to avoid confusion
    df$CumAscGood=NULL
    df$CumAscBad=NULL
    df$CumAscTotal=NULL
    
    # TN and FN
    df$FN=SumGoods-df$CumDescGood
    df$TN=SumBads-df$CumDescBad
    
    # Accuracy
    df$Accuracy=(df$CumDescGood+df$TN)/(df$TN+df$FN+df$CumDescGood+df$CumDescBad)
    
    # Specificity, Sensitivity
    df$Sensitivity=df$CumDescGood/(df$CumDescGood+df$FN)
    df$Specificity=df$TN/(df$TN+df$CumDescBad)
    
    # False Positive Rate
    df$FalPosRate=df$CumDescBad/(df$CumDescBad+df$TN)
    
    # Precision
    df$Precision=df$CumDescGood/(df$CumDescGood+df$CumDescBad)
    
    # Inverse Precision
    df$InvPrecision=df$TN/(df$FN+df$TN)
    
    # Youden Index
    df$YoudenJ=df$Sensitivity+df$Specificity-1
    # max(df$YoudenJ)
    
    # Optimal Cutoff
    optcut=df[df$YoudenJ==max(df$YoudenJ), ]$Prediction
    df$YoudenJ=NULL
    optcutcomment=" (Optimal)"
    
    # If cutoff is specified
    if(!is.na(cutoff)){
      if ((cutoff %in% df$Prediction)==FALSE){
        optcut=df[which.min(abs(as.numeric(df$Prediction)-cutoff)),1]
        optcutcomment=paste0(" (",cutoff," Not Found)")
      } else {
        optcut=cutoff
        optcutcomment=" (User Defined)"
        
      }
    }
    
    # For AUC Calculation (Trapezoid Method)
    df$TPR=df$Sensitivity
    df$FPR=1-df$Specificity
    df$MgAUC=0
    a=which(names(df)=="MgAUC")
    f=which(names(df)=="FPR")
    t=which(names(df)=="TPR")
    for (i in 1:nrow(df)-1) {
      df[i,a]=0.5*(df[i,t]+df[i+1,t])*(df[i,f]-df[i+1,f])
    }
    
    # AUC
    auc=sum(df$MgAUC)
    df$TPR=NULL
    df$FPR=NULL
    df$MgAUC=NULL
    
    # AUC Evaluation
    auceval=ifelse(auc<0.6,"Unpredictive",
                   ifelse(auc<0.7,"Poor",
                          ifelse(auc<0.8,"Fair",
                                 ifelse(auc<0.9,"Good","Excellent")))) 
    
    # KS
    df$MgKS=abs(df$PctCumAscGood-df$PctCumAscBad)
    ks=as.numeric(max(df$MgKS))
    scoreks=as.numeric(df[df$MgKS==ks, ]$Prediction)
    cgks=as.numeric(df[df$MgKS==ks, ]$PctCumAscGood)
    cbks=as.numeric(df[df$MgKS==ks, ]$PctCumAscBad)
    df$MgKS=NULL
    
    # KS Evaluation
    kseval=ifelse(ks<0.3,"Unpredictive",
                  ifelse(ks<0.4,"Fair",
                         ifelse(ks<0.5,"Good",
                                ifelse(ks<0.6,"Excellent", 
                                       ifelse(ks<0.7,"Awesome","That Good. Really?")))))
    
    # If report is activated (report = 1)
    if(report==1) {
      
      # Confusion Matrix Components
      tp=df[df$Prediction==optcut, ]$CumDescGood
      fp=df[df$Prediction==optcut, ]$CumDescBad
      fn=df[df$Prediction==optcut, ]$FN
      tn=df[df$Prediction==optcut, ]$TN
      p=SumGoods
      n=SumBads
      recsabovecutoff=df[df$Prediction==optcut, ]$CumDescTotal/SumRecords
      goodrate=df[df$Prediction==optcut, ]$GoodRateDesc
      badrate=df[df$Prediction==optcut, ]$BadRateDesc
      
      # Report on Metrics
      admetrics=character()
      admetrics=paste0(admetrics, "\n")
      admetrics=paste0(admetrics, "  Overall Performance Metrics \n")
      admetrics=paste0(admetrics, "  -------------------------------------------------- \n")
      admetrics=paste0(admetrics, "                    KS : ",sprintf("%.4f",round(ks,4))," (",kseval,")\n")
      admetrics=paste0(admetrics, "                   AUC : ",sprintf("%.4f",round(auc,4))," (",auceval,")\n")
      admetrics=paste0(admetrics, "\n")
      admetrics=paste0(admetrics, "  Classification Matrix \n")
      admetrics=paste0(admetrics, "  -------------------------------------------------- \n")
      admetrics=paste0(admetrics, "           Cutoff (>=) : ",round(as.numeric(optcut),4),optcutcomment,"\n")
      admetrics=paste0(admetrics, "   True Positives (TP) : ",tp,"\n")
      admetrics=paste0(admetrics, "  False Positives (FP) : ",fp,"\n")
      admetrics=paste0(admetrics, "  False Negatives (FN) : ",fn,"\n")
      admetrics=paste0(admetrics, "   True Negatives (TN) : ",tn,"\n")
      admetrics=paste0(admetrics, "   Total Positives (P) : ",p,"\n")
      admetrics=paste0(admetrics, "   Total Negatives (N) : ",n,"\n")
      admetrics=paste0(admetrics, "\n")
      admetrics=paste0(admetrics, "  Business/Performance Metrics \n")
      admetrics=paste0(admetrics, "  -------------------------------------------------- \n")
      admetrics=paste0(admetrics, "      %Records>=Cutoff : ",sprintf("%.4f",round(recsabovecutoff,4)),"\n")
      admetrics=paste0(admetrics, "             Good Rate : ",sprintf("%.4f",round(goodrate,4))," (Vs ",sprintf("%.4f",round(SumGoods/SumRecords,4))," Overall)\n")
      admetrics=paste0(admetrics, "              Bad Rate : ",sprintf("%.4f",round(badrate,4))," (Vs ",sprintf("%.4f",round(SumBads/SumRecords,4))," Overall)\n")
      admetrics=paste0(admetrics, "        Accuracy (ACC) : ",sprintf("%.4f",round((tp+tn)/(tp+fp+tn+fn),4)),"\n")
      admetrics=paste0(admetrics, "     Sensitivity (TPR) : ",sprintf("%.4f",round(tp/p,4)),"\n")
      admetrics=paste0(admetrics, " False Neg. Rate (FNR) : ",sprintf("%.4f",round(fn/p,4)),"\n")
      admetrics=paste0(admetrics, " False Pos. Rate (FPR) : ",sprintf("%.4f",round(fp/n,4)),"\n")
      admetrics=paste0(admetrics, "     Specificity (TNR) : ",sprintf("%.4f",round(tn/n,4)),"\n")
      admetrics=paste0(admetrics, "       Precision (PPV) : ",sprintf("%.4f",round(tp/(tp+fp),4)),"\n")
      admetrics=paste0(admetrics, "  False Discovery Rate : ",sprintf("%.4f",round(fp/(tp+fp),4)),"\n")
      admetrics=paste0(admetrics, "    False Omision Rate : ",sprintf("%.4f",round(fn/(fn+tn),4)),"\n")
      admetrics=paste0(admetrics, "  Inv. Precision (NPV) : ",sprintf("%.4f",round(tn/(fn+tn),4)),"\n")
      admetrics=paste0(admetrics,"\n")
      admetrics=paste0(admetrics, "  Note: ",nmiss," rows deleted due to missing data.\n")
      admetrics=paste0(admetrics,"\n")
      admetrics=gsub(", ","",admetrics)
      
      # Metric report
      cat(admetrics)
    } # End if for report
    
    # AUC Plot
    if(plot=="auc"){
      plot(NULL, 
           xlim=c(0,1), 
           ylim=c(0,1), 
           main="ROC Curve",
           xlab="1-Specificity",
           ylab="Sensitivity")
      grid()
      points(1-df$Specificity, df$Sensitivity,
             pch=21,col="grey25",bg="grey25")
      abline(0,1,col="grey50", lty=2)
      text(0.8,0.1,paste0("AUC : ",round(100*auc,2),"%"))
    }
    
    # KS Plot
    if(plot=="ks"){
      plot(NULL, 
           xlim=c(as.numeric(min(df$Prediction)),as.numeric(max(df$Prediction))), 
           ylim=c(0,1), 
           main="Cumulative Distribution of Goods/Bads",
           xlab=prediction,
           ylab="Cumulative Distribution")
      grid()
      text(as.numeric(min(df$Prediction)),0.95,"* Bads",pos=4, col="firebrick")
      text(as.numeric(min(df$Prediction)),0.90,"* Goods",pos=4, col="dodgerblue4")
      segments(scoreks,cgks,scoreks,cbks, lty=2, lwd = 2, col="grey50")
      points(df$Prediction, df$PctCumAscBad,
             pch=21,col="firebrick",bg="firebrick")
      points(df$Prediction, df$PctCumAscGood,
             pch=21,col="dodgerblue4",bg="dodgerblue4")
      text(0.95*as.numeric(max(df$Prediction)),0.1,paste0("KS : ",round(100*ks,2),"%"))
    }
    
    if(returndf==1){return(df)}
  } # Close else
} # Close function

# Ends Metrics 20171009 #######################################################


# Ini: Metrics Plot 20171022 #############################################
#' Visualization of a Classification Matrix
#'
#' It generates four plots after running and saving the output report from \code{smbinning.metrics}.
#' @param df Data frame generated with \code{smbinning.metrics}.
#' @param cutoff Value of the classifier that splits the data between positive (>=) and negative (<).
#' @param plot Plot to be drawn. Options are: 'cmactual' (default),'cmactualrates','cmmodel','cmmodelrates'.
#' @examples 
#' # Load library and its dataset
#' library(smbinning)
#' smbmetricsdf=smbinning.metrics(dataset=smbsimdf1, prediction="cbs1",
#'                                actualclass="fgood", returndf=1)
#' 
#' # Example 1: Plots based on optimal cutoff
#' smbinning.metrics.plot(df=smbmetricsdf,plot='cmactual')
#' 
#' # Example 2: Plots using user defined cutoff
#' smbinning.metrics.plot(df=smbmetricsdf,cutoff=600,plot='cmactual')
#' smbinning.metrics.plot(df=smbmetricsdf,cutoff=600,plot='cmactualrates')
#' smbinning.metrics.plot(df=smbmetricsdf,cutoff=600,plot='cmmodel')
#' smbinning.metrics.plot(df=smbmetricsdf,cutoff=600,plot='cmmodelrates')

smbinning.metrics.plot=function(df,cutoff=NA, plot="cmactual"){
  
  if(names(df)[1]!="Prediction" | names(df)[2]!="CntGood"){
    return("Data not from smbinning.metrics.")
  } else if (!is.na(cutoff) & !is.numeric(cutoff)){ # Check if target variable is numeric
    return("'cutoff' must be numeric.")
  } else if(!is.na(cutoff) & (max(df$Prediction)<cutoff | min(df$Prediction)>cutoff)){ # Check if target variable is numeric
    return("'cutoff' out of range.")
  } else if(plot!="cmactual" & plot!="cmactualrates" & plot!="cmmodel" & plot!="cmmodelrates"){ 
    return("'plot' options are: 'auc', 'ks' or 'none'.")
    }
    else {
  
  df$YoudenJ=df$Sensitivity+df$Specificity-1
  optcut=df[df$YoudenJ==max(df$YoudenJ), ]$Prediction
  df$YoudenJ=NULL
  
  # If cutoff is specified
  if(!is.na(cutoff)){
    if ((cutoff %in% df$Prediction)==FALSE){
      optcut=df[which.min(abs(as.numeric(df$Prediction)-cutoff)),1]
    } else {
      optcut=cutoff
    }
  }
  
  
  # Confusion Matrix Components
  tp=df[df$Prediction==optcut, ]$CumDescGood
  fp=df[df$Prediction==optcut, ]$CumDescBad
  fn=df[df$Prediction==optcut, ]$FN
  tn=df[df$Prediction==optcut, ]$TN
  
  p=sum(df$CntGood)
  n=sum(df$CntBad)
  
  # CM Metrics
  accuracy=(tp+tn)/(tp+fp+tn+fn)
  sensitivity=tp/p
  FNR=fn/p
  specificity=tn/n
  FPR=fp/n
  precision=tp/(tp+fp)
  invprecision=tn/(fn+tn)
  FDR=fp/(tp+fp)
  FOR=fn/(fn+tn)
  # For AUC Calculation (Trapezoid Method)
  df$TPR=df$Sensitivity
  df$FPR=1-df$Specificity
  df$MgAUC=0
  a=which(names(df)=="MgAUC")
  f=which(names(df)=="FPR")
  t=which(names(df)=="TPR")
  for (i in 1:nrow(df)-1) {
    df[i,a]=0.5*(df[i,t]+df[i+1,t])*(df[i,f]-df[i+1,f])
  }
  # AUC
  auc=sum(df$MgAUC)
  df$TPR=NULL
  df$FPR=NULL
  df$MgAUC=NULL
  
  # KS
  df$MgKS=abs(df$PctCumAscGood-df$PctCumAscBad)
  ks=as.numeric(max(df$MgKS))
  df$MgKS=NULL
  
  # Classification Matrix
  cmnbr=matrix(c(tp,fn,fp,tn),nrow=2,ncol=2)
  cmpctactual=matrix(c(sensitivity,FNR,FPR,specificity),nrow=2,ncol=2)
  cmpctmodel=matrix(c(precision,FDR,FOR,invprecision),nrow=2,ncol=2)

  # CM, where X-axis is actual class 
  if(plot=="cmactual"){
    cmnbrplot=cmnbr
    colnames(cmnbrplot)=c("Actual+\n[ TP | FN ]","Actual-\n[ FP | TN ]") # Actuals
    rownames(cmnbrplot)=c("Model+","Model-") # Model
    bpactual=
      barplot(cmnbrplot,
              ylim=c(0,round(1.25*max(cmnbrplot),0)),
              col=c("grey50","grey85"),
              beside = T,
              axes = T,
              border=NA,
              legend.text = T,
              args.legend = list(x = "top",border=NA,bty="n",horiz=F))
    mtext(side=3,"Classification Matrix",line=2,cex=1.2, font=2)
    mtext(side=3,paste0("Records by Actual Class. Cutoff >=",optcut),line=0.75,cex=1)
    mtext(side=1,"Actual Class",line=3,cex=1, font=2)
    abline(h=0) # Horizontal line
    text(x=bpactual, y=cmnbrplot, label=cmnbrplot, pos=3, cex=1)
  }
  
  
  # CM Percentages, where X-axis is actual class 
  if(plot=="cmactualrates"){
    cmpctactualplot=cmpctactual
    colnames(cmpctactualplot)=c("Actual+\n[ Sensitivity | FNR ]","Actual-\n[ FPR | Specificity ]") # Actuals
    rownames(cmpctactualplot)=c("Model+","Model-") # Model
    bpactualpct=
      barplot(cmpctactualplot,
              ylim=c(0,1),
              col=c("grey50","grey85"),
              beside = T,
              axes = T, 
              border=NA,
              legend.text = T,
              args.legend = list(x = "top",border=NA,bty="n",horiz=F))
    mtext(side=3,"Classification Matrix",line=2,cex=1.2, font=2)
    mtext(side=3,paste0("Percentage by Actual Class. Cutoff >=",optcut),line=0.75,cex=1)
    mtext(side=1,"Actual Class",line=3,cex=1, font=2)
    abline(h=0) # Horizontal line
    text(x=bpactualpct, y=cmpctactualplot, label=round(cmpctactualplot,4), pos=3, cex=1)
  }
  
  # CM Numbers, where X-axis is model prediction 
  if(plot=="cmmodel"){
    cmnbrplot=cmnbr
    colnames(cmnbrplot)=c("Actual+","Actual-") # Actuals
    rownames(cmnbrplot)=c("Model+\n[ TP | FP ]","Model-\n [ FN | TN ]") # Model
    bpmodel=
      barplot(t(cmnbrplot),
              ylim=c(0,round(1.25*max(cmnbrplot),0)),
              col=c("grey50","grey85"),
              beside = T,
              axes = T, 
              border=NA,
              legend.text = T,
              args.legend = list(x = "top",border=NA,bty="n",horiz=F))
    mtext(side=3,"Classification Matrix",line=2,cex=1.2, font=2)
    mtext(side=3,paste0("Records by Model Prediction. Cutoff >=",optcut),line=0.75,cex=1)
    mtext(side=1,"Model Prediction",line=3,cex=1, font=2)
    abline(h=0) # Horizontal line
    text(x=bpmodel, y=t(cmnbrplot), label=t(cmnbrplot), pos=3, cex=1)
  }
  
  # CM Numbers, where X-axis is model prediction 
  if(plot=="cmmodelrates"){
    cmpctmodelplot=cmpctmodel
    colnames(cmpctmodelplot)=c("Model+\n[ PPV | FDR ]","Model-\n[ FOR | NPV ]") # Actuals
    rownames(cmpctmodelplot)=c("Actual+","Actual-") # Model
    bpactualpct=
      barplot(cmpctmodelplot,
              ylim=c(0,1),
              col=c("grey50","grey85"),
              beside = T,
              axes = T, 
              border=NA,
              legend.text = T,
              args.legend = list(x = "top",border=NA,bty="n",horiz=F))
    mtext(side=3,"Classification Matrix",line=2,cex=1.2, font=2)
    mtext(side=3,paste0("Percentage by Model Prediction. Cutoff >=",optcut),line=0.75,cex=1)
    mtext(side=1,"Model Prediction",line=3,cex=1, font=2)
    abline(h=0) # Horizontal line
    text(x=bpactualpct, y=cmpctmodelplot, label=round(cmpctmodelplot,4), pos=3, cex=1)
  }
  
    } # end else
}

# End: Metrics Plot 20171022 #############################################


# Ini Monotonic Binning 20181016 #########################################
#' Monotonic Binning
#'
#' It gives the user the ability to impose a monotonic trend for good/bad rates per bin.
#' @param df A data frame.
#' @param y Binary response variable (0,1). Integer (\code{int}) is required.
#' Name of \code{y} must not have a dot. Name "default" is not allowed.
#' @param x Continuous characteristic. At least 5 different values. Value \code{Inf} is not allowed. 
#' Name of \code{x} must not have a dot.
#' @param p Percentage of records per bin. Default 5\% (0.05). 
#' @return The command \code{smbinning.monotonic} generates and object containing the necessary info and utilities for binning.
#' The user should save the output result so it can be used 
#' with \code{smbinning.plot}, \code{smbinning.sql}, and \code{smbinning.gen}.
#' @examples
#' # Load library and its dataset
#' library(smbinning) # Load package and its data
#' 
#' # Example 1: Monotonic Binning (Increasing Good Rate per Bin)
#' smbinning(df=smbsimdf2,y="fgood2",x="chr2",p=0.05)$ivtable # Run regular binning
#' smbinning.monotonic(df=smbsimdf2,y="fgood2",x="chr2",p=0.05)$ivtable # Run monotonic binning
#'  
#' # Example 2: Monotonic Binning (Decreasing Good Rate per Bin)
#' smbinning(df=smbsimdf2,y="fgood3",x="chr3",p=0.05)$ivtable # Run regular binning
#' smbinning.monotonic(df=smbsimdf2,y="fgood3",x="chr3",p=0.05)$ivtable # Run monotonic binning

smbinning.monotonic=function(df,y,x,p=0.05) { # Ini function monotonic
  
  i=which(names(df)==y) # Column for y
  j=which(names(df)==x) # Column for x
  
  result=smbinning(df, y, x, p) # Save result from usual binning
  c=cor(df[,i],df[,j],use = "complete.obs", method = c("pearson"))
  col=if(c>0) {9} else {10} # Increasing (Column 9 Good Rate) or Decreasing (Column 10 Bad Rate)?
  
  if(result$iv<0.1) {return("Not Meaningful (IV<0.1)")} 
  
  else { # Ini condition
    
    # Get relevant data
    ivtable=result$ivtable
    ratevalues=ivtable[,col] # Column 9 for increasing (Good Rate), 10 for decreasing (Bad Rate)
    ratevalues=ratevalues[1:(length(ratevalues)-2)] # Excludes total and missing
    ratecuts=result$cuts
    
    # Start WHILE Loop
    count=0
    iter=0
    
    while(count<1) {
      
      # If last bin not follow trend, change it the previous bin
      i=length(ratevalues)-1
      ratecuts[i]=ifelse(ratevalues[i+1]<ratevalues[i], ratecuts[i-1],ratecuts[i])
      # If next bin (i+1) lower than previous (i) then merge 
      for (i in 1:length(ratevalues)-1) {
        ratecuts[i]=ifelse(ratevalues[i+1]<ratevalues[i], ratecuts[i+1],ratecuts[i])
      }
      
      # Make unique bin values
      ratecuts=unique(ratecuts)
      ratecuts=ratecuts[!is.na(ratecuts)]
      
      result=smbinning.custom(df, y, x, ratecuts)
      ivtable=result$ivtable
      ratevalues=ivtable[,col]
      ratevalues=ratevalues[1:(length(ratevalues)-2)] # Excludes total and missing
      
      iter=iter+1 # Number of iterations
      
      count=if(all(ratevalues==cummax(ratevalues))==FALSE) {0} else{1}
      
    }
  } # End condition
  
  return(result)
  
} # End function monotonic

# End Monotonic Binning 20181016 #########################################


# Begin Plotting ##############################################################
#' Plots after binning
#'
#' It generates plots for distribution, bad rate, and weight of evidence after running \code{smbinning} 
#' and saving its output.
#' @param ivout An object generated by binning.
#' @param option Distribution ("dist"), Good Rate ("goodrate"), Bad Rate ("badrate"), and Weight of Evidence ("WoE").
#' @param sub Subtitle for the chart (optional).
#' @examples 
#' # Load library and its dataset
#' library(smbinning)
#' 
#' # Example 1: Numeric variable (1 page, 4 plots)
#' result=smbinning(df=smbsimdf1,y="fgood",x="cbs1") # Run and save result
#' par(mfrow=c(2,2))
#' boxplot(smbsimdf1$cbs1~smbsimdf1$fgood,
#'         horizontal=TRUE, frame=FALSE, col="lightgray",main="Distribution")
#' mtext("Credit Score",3)
#' smbinning.plot(result,option="dist",sub="Credit Score")
#' smbinning.plot(result,option="badrate",sub="Credit Score")
#' smbinning.plot(result,option="WoE",sub="Credit Score")
#' par(mfrow=c(1,1))
#' 
#' # Example 2: Factor variable (1 plot per page)
#' result=smbinning.factor(df=smbsimdf1,y="fgood",x="inc",maxcat=11)
#' smbinning.plot(result,option="dist",sub="Income Level")
#' smbinning.plot(result,option="badrate",sub="Income Level")
#' smbinning.plot(result,option="WoE",sub="Income Level")

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

# Ini: PSI 20170821 ###########################################################
#' Population Stability Index
#' 
#' Often models are developed using multiple periods in time for a number of reasons.
#' For example, to avoid seasonality, to increase the size of the population, and some others.
#' With a metrics like the Population Stability Index (PSI), users can check if there is
#' a significant variation in the distribution of a certain feature by partition (usually time)
#' using the first one as the reference.
#' @param df Data frame.
#' @param y Column name the indicates the different partitions.
#' @param x Feature to be evaluated in terms of stability (It must be factor).
#' @return Three crosstabs by feature and period that show the frequency (psicnt), 
#' percentage (psipct) and PSI (psimg), and a plot for the analyzed characteristic.
#' @examples 
#' # Load library and its dataset
#' library(smbinning)
#'
#' # Check stability for income
#' smbinning.psi(df=smbsimdf1,y="period",x="inc") 

smbinning.psi = function(df,y,x){
  i=which(names(df)==x)
  j=which(names(df)==y)
  if(!is.data.frame(df)){return(stop("Not a data frame"))}
  else if(identical(i,integer(0))){return(stop(paste("Characteristic",x,"not found")))}
  else if(identical(j,integer(0))){return(stop(paste("Characteristic",y,"not found")))}
  else if(class(df[,i])!="factor"){return(stop("x must be formatted as factor"))}
  else {
    psicnt=table(df[,i],df[,j],useNA = "ifany") # Table with counts including NA
    options(scipen=999) # Remove Scientific Notation
    psipct=prop.table(psicnt, margin=2) # Table with column percentage
    psimg=psipct # Shell for PSI table
    n=ncol(psipct) # Number of columns (Periods)
    m=nrow(psipct) # Number of rowa (Periods)
    
    psimg[,1]=0 # Step 1: PSI=0 for first column
    
    # PSI Period VS First Period
    for (k in 2:n){
      for (l in 1:m){
        if(psipct[l,1]>0 & psipct[l,k]>0) 
        {psimg[l,k]=round((psipct[l,k]-psipct[l,1])*log(psipct[l,k]/psipct[l,1]),8)}
        else {psimg[l,k]=0}
      }
    }
    
    psimg=rbind(psimg, PSI=colSums(psimg))
    psimg=as.table(psimg) # Table with Mg PSI
    psitable=psimg[nrow(psimg),] # Extract total PSI only
    psitable=as.data.frame(psitable)
    # Plot
    psitable$Partition=rownames(psitable) # Create column "Partition"
    rownames(psitable)=NULL  # Remove rownames
    names(psitable)=c("PSI","Partition") # Rename columns
    psitable=psitable[,c("Partition","PSI")] # Reorder
    
    maxpsi=max(psitable$PSI)
    psiylim=ifelse(maxpsi<0.25,0.25,maxpsi)
    plot(psitable$PSI,
         main = paste("Stability:",x),
         xlab = paste("Partition:",y),
         ylab="PSI",
         ylim=c(0,psiylim), 
         pch=19,
         col="black",
         xaxt = "n") # Remove x axis values
    abline(h=0.1,lty=2)
    abline(h=0.25,lty=2)
    axis(1,at=1:nrow(psitable),psitable$Partition) # Add period
    
    list(psicnt=psicnt,psipct=psipct,psimg=psimg)  
    }
}

# End: PSI 20170821 ###########################################################

# Begin Model Scaling 20170821 ################################################
#' Scaling
#'
#' It transforms the coefficients of a logistic regression into scaled points
#' based on the following three parameters pre-selected by the analyst: PDO, Score, and Odds.
#' @param logitraw Logistic regression (glm) that must have specified \code{family=binomial} and 
#' whose variables have been generated with \code{smbinning.gen} or \code{smbinning.factor.gen}.
#' @param pdo Points to double the oods.
#' @param odds Desired \code{odds} at the selected \code{score}.
#' @param score Score at which the desire \code{odds} occur.
#' @return A scaled model from a logistic regression built with binned variables, the parameters
#' used in the scaling process, the expected minimum and maximum score, and the original logistic model.
#' @examples 
#' # Load library and its dataset
#' library(smbinning)
#' 
#' # Sampling
#' pop=smbsimdf1 # Population
#' train=subset(pop,rnd<=0.7) # Training sample
#' 
#' # Generate binning object to generate variables
#' smbcbs1=smbinning(train,x="cbs1",y="fgood")
#' smbcbinq=smbinning.factor(train,x="cbinq",y="fgood")
#' smbcblineut=smbinning.custom(train,x="cblineut",y="fgood",cuts=c(30,40,50))
#' smbpmt=smbinning.factor(train,x="pmt",y="fgood")
#' smbtob=smbinning.custom(train,x="tob",y="fgood",cuts=c(1,2,3))
#' smbdpd=smbinning.factor(train,x="dpd",y="fgood")
#' smbdep=smbinning.custom(train,x="dep",y="fgood",cuts=c(10000,12000,15000))
#' smbod=smbinning.factor(train,x="od",y="fgood")
#' smbhome=smbinning.factor(train,x="home",y="fgood")
#' smbinc=smbinning.factor.custom(
#'   train,x="inc",y="fgood",
#'   c("'W01','W02'","'W03','W04','W05'","'W06','W07'","'W08','W09','W10'"))
#' 
# Generate new characteristics and update population dataset
#' pop=smbinning.gen(pop,smbcbs1,"g1cbs1")
#' pop=smbinning.factor.gen(pop,smbcbinq,"g1cbinq")
#' pop=smbinning.gen(pop,smbcblineut,"g1cblineut")
#' pop=smbinning.factor.gen(pop,smbpmt,"g1pmt")
#' pop=smbinning.gen(pop,smbtob,"g1tob")
#' pop=smbinning.factor.gen(pop,smbdpd,"g1dpd")
#' pop=smbinning.gen(pop,smbdep,"g1dep")
#' pop=smbinning.factor.gen(pop,smbod,"g1od")
#' pop=smbinning.factor.gen(pop,smbhome,"g1home")
#' pop=smbinning.factor.gen(pop,smbinc,"g1inc")
#' 
#' # Resample
#' train=subset(pop,rnd<=0.7) # Training sample
#' test=subset(pop,rnd>0.7) # Testing sample
#' 
#' # Run logistic regression
#' f=fgood~g1cbs1+g1cbinq+g1cblineut+g1pmt+g1tob+g1dpd+g1dep+g1od+g1home+g1inc
#' modlogisticsmb=glm(f,data = train,family = binomial())
#' summary(modlogisticsmb)
#' 
#' # Example: Scaling from logistic parameters to points
#' smbscaled=smbinning.scaling(modlogisticsmb,pdo=20,score=720,odds=99)
#' smbscaled$logitscaled # Scaled model
#' smbscaled$minmaxscore # Expected minimum and maximum Score
#' smbscaled$parameters # Parameters used for scaling
#' summary(smbscaled$logitraw) # Extract of original logistic regression
#' 
#' # Example: Generate score from scaled model
#' pop1=smbinning.scoring.gen(smbscaled=smbscaled, dataset=pop)
#' 
#' # Example Generate SQL code from scaled model
#' smbinning.scoring.sql(smbscaled)

smbinning.scaling=function(logitraw, pdo=20, score=720, odds=99) {
  if(missing(logitraw)){
    return(stop("Logistic model missing"))
  }
  else if(as.character(summary(logitraw)[1]$call)[1] != "glm") { 
    return(stop("Logistic model must be glm"))
  }
  else if(as.character(logitraw$family)[1] != "binomial") { 
    return(stop("Logistic model must be a binomial family"))
  }
  else if (names(logitraw$coefficients[1]) != "(Intercept)"){ 
    return(stop("Logistic model must have a constant"))
  }
  else if (!is.numeric(pdo)) {
    return(stop("PDO must be numeric"))
  }
  else if (!is.numeric(score)) {
    return(stop("Score must be numeric"))
  }
  else if (!is.numeric(odds)) {
    return(stop("Odds must be numeric"))
  }
  else if (all.equal(pdo, as.integer(pdo))!=TRUE | pdo<0) {
    return(stop("PDO must be positive integer"))
  }
  else if (all.equal(score, as.integer(score))!=TRUE | score<0) {
    return(stop("score must be positive integer"))
  }
  else if (all.equal(odds, as.integer(odds))!=TRUE | odds<0) {
    return(stop("Odds must be positive integer"))
  }
  else{
    
    # Number of characteristics
    nchr=length(logitraw$xlevels) 
    
    # Characteristics and bin names #
    chrname=list() # Generate empty list
    chrbinname=list() # Generate empty list
    chrbinlist <- list() # Generate empty list
    for (i in 1:nchr) {
      chrname[i]=names(data.frame(logitraw$xlevels[i])) # Name of the characteristics
      chrbin=data.frame(names(table(logitraw$xlevels[i]))) # Characteristic's bins
      n=as.character(chrname[i])
      binappend=list()
      for (j in 1:dim(chrbin)[1]) {
        bincut=as.character(chrbin[j,1]) # For example "00 Miss"
        chrbinname=c(chrbinname, paste(chrname[i],bincut,sep = "")) # For example "binnedage 01 <= 25"
        binappend=c(binappend, bincut)
      }
      chrbinlist[n] <- list(binappend)
    }
    
    # Bins and coefficients
    bincoeffname=rownames(as.data.frame(logitraw$coefficients)) # Getting bin names
    bincoeff=as.list(matrix(logitraw$coefficients)) # Getting bin coefficients
    # Creating a list of bins and their coefs
    names(bincoeff)=bincoeffname
    # Updating bins, their coefficients and adding base bins and their coefficients
    for (i in 1:length(chrbinname)) {
      # if not a base bin then next
      if (exists(as.character(chrbinname[i]), where=bincoeff)) {
        next
      }
      # if  a base bin then add it to bincoeff
      else {
        a=as.character(chrbinname[i])
        bincoeff[a]=0.0000
      }
    }
    
    # Naming
    bincoeffname=names(bincoeff)
    
    # Scaling parameters
    factor=pdo/(log(2))
    offset=(score)-((factor)*(log(odds)))
    intercept=unname(logitraw$coefficients["(Intercept)"])
    interceptweight=(intercept)*(factor)
    interceptbychr=(interceptweight)/(nchr)
    offsetbychr=(offset)/(nchr)
    scorebychr=(interceptbychr)+(offsetbychr)
    
    # Updating the bincoeff into a dataframe
    bincoeff=data.frame(sapply(bincoeff,function(x) x[[1]] ))
    colnames(bincoeff)="Coefficient"
    bincoeff=cbind(rownames(bincoeff), bincoeff)
    bincoeff$Weight=(bincoeff$Coefficient)*(factor)
    bincoeff$WeightScaled=(bincoeff$Weight)+(scorebychr)
    bincoeff[1,4]=0.000 #Make scaled constant equal to zero
    bincoeff$Points=round(bincoeff$WeightScaled,0)
    colnames(bincoeff)=c("FullName","Coefficient","Weight","WeightScaled","Points")
    # Sorting bincoeff
    FullName=unlist(chrbinname)
    FullName=c("(Intercept)",FullName)
    bincoeff$FullName=factor(bincoeff$FullName, levels=FullName)
    bincoeff=bincoeff[order(bincoeff$FullName),]
    #bincoeff=within(bincoeff, WeightScaled[FullName=='(Intercept)']==0)
    rownames(bincoeff) <- 1:dim(bincoeff)[1]
    # Create attributes
    Attribute=unlist(chrbinlist)
    # Add value for intercept to attributes
    Attribute=c("", Attribute)
    Attribute=unname(data.frame(Attribute))
    rownames(Attribute)=NULL
    # Creating Characteristic
    Characteristic=list()
    for (i in 1:length(FullName)) {
      # Characteristic=c(Characteristic, gsub(gsub("[[:punct:]]","",as.character(Attribute[[i,1]])),"",gsub("[[:punct:]]","",FullName[i])))
      Characteristic=c(Characteristic, gsub( "[0-9][0-9] .*$", "",  FullName[i]))
    }
    Characteristic=unlist(Characteristic)
    Characteristic=unname(data.frame(Characteristic))
    # Removing column FullName from the bincoeff
    drops=c("FullName")
    bincoeff=bincoeff[ , !(names(bincoeff) %in% drops)]
    # Adding attributes to bincoeff
    bincoeff=cbind(Attribute, bincoeff)
    # Adding characteristic to bincoeff
    bincoeff=cbind(Characteristic, bincoeff)
    
    # Get Min/Max Score
    chrpts=bincoeff
    chrpts=chrpts[,c("Characteristic","Points")]
    chrpts=chrpts[-1,] # Remove (intercept)
    minmaxscore=c(
      sum(aggregate(Points ~ Characteristic, chrpts, min)$Points),
      sum(aggregate(Points ~ Characteristic, chrpts, max)$Points)
    )
    
    # Converting bincoeff into a list
    bincoeff=list(bincoeff)
    
    # Creating a list for all the parameters
    parameters=list()
    parameters$pdo=pdo
    parameters$odds=odds
    parameters$factor=factor
    parameters$score=score
    parameters$offset=offset
    parameters$intercept=intercept
    parameters$interceptweight=interceptweight
    parameters$nchr=nchr
    parameters$interceptbychr=interceptbychr
    parameters$offsetbychr=offsetbychr
    parameters$scorebychr=scorebychr
    
    # Final output
    modelscaled=list()
    modelscaled$logitscaled=bincoeff
    modelscaled$parameters=parameters
    modelscaled$logitraw=logitraw
    modelscaled$minmaxscore=minmaxscore
    
    # Return
    return(modelscaled)
  }
}

# End Model Scaling 20170821 ################################################


# Begin Add Points and Score 20170925 ##########################################
#' Generation of Score and Its Weights
#'
#' After applying \code{smbinning.scaling} to the model, \code{smbinning.scoring} generates a data frame 
#' with the final Score and additional fields with the points assigned to each characteristic so the user
#' can see how the final score is calculated. Example shown on \code{smbinning.scaling} section.
#' @param smbscaled Object generated using \code{smbinning.scaling}.
#' @param dataset A data frame.
#' @return The command \code{smbinning.scoring} generates a data frame with the final scaled Score and its
#' corresponding scaled weights per characteristic.

smbinning.scoring.gen=function(smbscaled, dataset) {
  if(missing(smbscaled)){
    # Check if logitscaled is missing or not
    return(stop("Missing argument: smbscaled"))
  } 
  else if(missing(dataset)){
    # Check if dataset is missing or not
    return(stop("Missing argument: dataset"))
  } 
  else if (!is.data.frame(dataset)){ 
    # Check if data.frame
    return(stop("Not a data frame"))
  }
  else if (names(smbscaled)[1] != "logitscaled"){ 
    return(stop("Not from 'smbinning.scaling'"))
  } else {
    df=dataset
    logitraw=smbscaled$logitraw
    logitscaled=smbscaled$logitscaled
    chrattpts=as.data.frame(logitscaled)
    chrattpts=chrattpts[,c("Characteristic","Attribute","Points")]
    for (i in 1:length(logitraw$xlevels)) {
      chrname=names(logitraw$xlevels[i])
      chrattptstmp=subset(chrattpts,chrattpts$Characteristic==chrname) # Tmp df for char [i] with attr and points
      df=cbind(df,chrtmp=NA)
      colidx=which(names(df)==chrname) # Where is original characteristic
      # df=cbind(df,chrtmporiginal=NA)
      # df$chrtmporiginal=df[,colidx] # Populate temporary original characteristic
      for (j in 1:nrow(chrattptstmp)) {
        df=within(df,chrtmp[df[,colidx]==logitraw$xlevels[[i]][j]]<-chrattptstmp[j,][3])
      }
      # df$chrtmporiginal=NULL
      df$chrtmp=as.numeric(df$chrtmp)
      
      if(paste0(names(logitraw$xlevels[i]),"Points") %in% colnames(df))
      {stop("Column '",paste0(names(logitraw$xlevels[i]),"Points' already exists. Drop it or rename it."))}
      
      names(df)[names(df)=="chrtmp"]=paste(chrname,"Points",sep = "") # Rename tmp column name to Points.
    }
    
    # Create final score
    if("Score" %in% colnames(df)) {stop("Column 'Score' already exists. Drop it or rename it.")} # Stop process if column Score already exists.
    df=cbind(df,Score=0) # Add score column
    scorecolid=ncol(df) # Score column id
    
    # All characteristics
    crhpts=list()
    for (i in 1:length(logitraw$xlevels)) {
      crhpts=cbind(crhpts,paste(names(logitraw$xlevels[i]),"Points",sep=""))
    }
    nbrchar=length(crhpts) # Number of new generated chars
    colini=which(names(df)==crhpts[1]) # Where is the new characteristic
    colend=colini+nbrchar-1
    df[,scorecolid]=rowSums(df[,colini:colend]) # Sum rows to get Score
  }
  return(df)
}

# End Add Points and Score 20170925 ###########################################


# Begin Convert model into a SQL Statement 20171117 ##########################
#' Generation of SQL Code After Scaled Model
#'
#' After applying \code{smbinning.scaling} to the model, \code{smbinning.scoring.sql} generates a SQL code 
#' that creates and updates all variables present in the scaled model. Example shown on \code{smbinning.scaling} section.
#' @param smbscaled Object generated using \code{smbinning.scaling}.
#' @return The command \code{smbinning.scoring.sql} generates a SQL code to implement the model the model in SQL.

smbinning.scoring.sql=function (smbscaled) {
  if(missing(smbscaled)){
    return(stop("Argument 'smbscaled' is missing"))
  } 
  else if (names(smbscaled)[1] != "logitscaled"){ 
    return(stop("'smbscaled' must come from 'smbinning.scaling'"))
  }
  else {
    logitscaled=smbscaled$logitscaled
    # SQL code 1: Create table
    logitscaleddf=data.frame(logitscaled)
    logitscaleddf=logitscaleddf[-1,] # Remove (intercept)
    uniquechctrs=unique(logitscaleddf$Characteristic) # Characteristics 
    codecreate=list()
    codecreate=c(codecreate,paste("-- Replace 'TableName' with your table in SQL \n"))
    codecreate=c(codecreate, "alter table TableName add \n") # First line of the sql code
    for (i in 1:length(uniquechctrs)) {
      pointname=paste(uniquechctrs[[i]], "Points", sep="")
      codecreate=c(codecreate, paste(pointname, "int not null, \n", sep=' '))
    }
    codecreate=c(codecreate, paste("Score", "int not null \n", sep=' '))
    codecreate=c(codecreate, paste("go \n", sep=' '))
    
    sqlcreate=character()
    for (i in 1:length(codecreate)) {
      sqlcreate=paste(sqlcreate,codecreate[[i]], sep="")
    }
    
    # SQL code 2: Update table
    codeupdate=list()
    codeupdate=c(codeupdate,paste("\n-- Replace 'TableName' with your table in SQL \n"))
    codeupdate=c(codeupdate,paste("update TableName set \n"))
    for (i in 1:length(uniquechctrs)) {
      subdf=subset(logitscaleddf, logitscaleddf$Characteristic==uniquechctrs[i])
      pointname=paste(uniquechctrs[[i]], "Points", sep="")
      codeupdate=c(codeupdate, pointname, "=", "( \n", sep=" ")
      codeupdate=c(codeupdate, "case \n")
      for (j in 1:dim(subdf)[1]) {
        codeupdate=c(codeupdate, paste(" when",
                                       subdf[["Characteristic"]][j],
                                       paste("= ","'",gsub("'","''",subdf[["Attribute"]][j]),"'",sep=""),
                                       "then",
                                       subdf[["Points"]][j],
                                       " \n",
                                       sep=" "))
      }
      codeupdate=c(codeupdate, " else Null end \n ), \n")
    }
    codeupdate=c(codeupdate, "Score=( \n")
    for (i in 1:length(uniquechctrs)) {
      pointname=paste(uniquechctrs[[i]], "Points", sep="")
      if (i==1) {
        codeupdate=c(codeupdate, paste(" ", pointname, sep=" "))
      } else {
        codeupdate=c(codeupdate, paste(" +", pointname, sep=" "))
      }
    }
    codeupdate=c(codeupdate, "\n) \n")
    
    sqlupdate=character()
    for (i in 1:length(codeupdate)) {
      sqlupdate=paste(sqlupdate,codeupdate[[i]], sep="")
    }
  }
  return(cat(sqlcreate,sqlupdate))
}
# End Convert model intto a SQL Statement 20171117 ############################


# Begin: SQL Code #############################################################
#' SQL Code
#'
#' It outputs a SQL code to facilitate the generation of new binned characetristic 
#' in a SQL environment. User must define table and new characteristic name.
#' @param ivout An object generated by \code{smbinning}.
#' @return A text with the SQL code for binning.
#' @examples 
#' # Load library and its dataset
#' library(smbinning)
#' 
#' # Example 1: Binning a numeric variable
#' result=smbinning(df=smbsimdf1,y="fgood",x="cbs1") # Run and save result
#' smbinning.sql(result)
#' 
#' # Example 2: Binning for a factor variable
#' result=smbinning.factor(df=smbsimdf1,x="inc",y="fgood",maxcat=11)
#' smbinning.sql(result)
#' 
#' # Example 3: Customized binning for a factor variable
#' result=smbinning.factor.custom(
#'   df=smbsimdf1,x="inc",y="fgood",
#'   c("'W01','W02'","'W03','W04','W05'",
#'     "'W06','W07'","'W08','W09','W10'"))
#' smbinning.sql(result)

smbinning.sql=function(ivout){
  if(is.null(ivout$groups)) {
    lines=nrow(ivout$ivtable)-2
    sqlcodetable=as.list(matrix(ncol=0,nrow=0))
    sqlcodetable=rbind(" alter table  'TableName'  add  'NewCharName'\n go\n update  'TableName'  set  'NewCharName'\n case\n")
    for (k in 1:lines){
      sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,ivout$ivtable[k,1],"then","'",sprintf("%02d",k),":",ivout$x,gsub("'","",ivout$ivtable[k,1]),"'\n"))
    }
    k=nrow(ivout$ivtable)-1
    sqlcodetable=rbind(sqlcodetable,paste("when",ivout$x,"Is Null then","'",sprintf("%02d",k),":",ivout$x,"Is Null' \n"))
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
    return(cat(gsub(",","",toString(sqlcodetable))))
    # For customized factor 20170910
  } else {
    sqlcodetable=as.list(matrix(ncol=0,nrow=0))
    sqlcodetable=rbind(" alter table  'TableName'  add  'NewCharName'\n go\n update  'TableName'  set  'NewCharName'\n case\n")
    for (i in 1:length(ivout$groups)){
      sqlcodetable=rbind(sqlcodetable,paste(" when",ivout$x,"in (",ivout$groups[i],") then","'",sprintf("%02d",i),":",ivout$x,gsub("'","",ivout$groups[i]),"'\n"))
    }
    i=length(ivout$groups)+1
    sqlcodetable=rbind(sqlcodetable,paste(" when",ivout$x,"Is Null then","'",sprintf("%02d",i),":",ivout$x,"Is Null'\n"))
    sqlcodetable=rbind(sqlcodetable,paste(" else '99: Error' end"))
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
    return(cat(gsub(", ","",toString(sqlcodetable))))
  }
}
# End: SQL Code ###############################################################


# Begin Summary IV 20160602 ###############################################
#' Information Value Summary
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
#' # Load library and its dataset
#' library(smbinning)
#' 
#' # Test sample
#' test=subset(smbsimdf1,rnd>0.9) # Training sample
#' test$rnd=NULL
#' 
#' # Example: Information Value Summary
#' testiv=smbinning.sumiv(test,y="fgood")
#' testiv
#'
#' # Example: Plot of Information Value Summary
#' smbinning.sumiv.plot(testiv)

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
#' Example shown on \code{smbinning.sumiv} section.
#' @param sumivt A data frame saved after \code{smbinning.sumiv}.
#' @param cex Optional parameter for the user to control the font size of the characteristics
#' displayed on the chart. The default value is 0.9
#' @return The command \code{smbinning.sumiv.plot} returns a plot that shows the IV
#' for each numeric and factor characteristic in the dataset.

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


# Ini: Simulated Credit Data ##############################################
#' Simulated Credit Data 
#'
#' A simulated dataset where the target variable is fgood, 
#' which represents the binary status of default (0) and not default (1).
#'
#' \itemize{
#'   \item fgood: Default (0), Not Default (1).
#'   \item cbs1: Credit quality index (1-100).
#'   \item cbs2: Profitability index (1-100).
#'   \item cbinq: Number of inquiries.
#'   \item cbline: Number of credit lines.
#'   \item cbterm: Number of term loans.
#'   \item cblineut: Line utilization (0-100).
#'   \item cbtob: Number of years on file.
#'   \item cbdpd: Indicator of days past due on bureau (Yes, No).
#'   \item cbnew: Number of new loans.
#'   \item pmt: Type of payment (M: Manual, A: Autopay, P: Payroll).
#'   \item tob: Time on books (Years).
#'   \item dpd: Level of delinquency (No, Low, High).
#'   \item dep: Amount of deposits.
#'   \item dc: Number of transactions.
#'   \item od: Number of overdrafts.
#'   \item home: Home ownership indicator (Yes, No).
#'   \item inc: Level of income.
#'   \item dd: Number of electronic transfers.
#'   \item online: Indicator of online activity (Yes, No).
#'   \item rnd: Random number to select testing and training samples.
#'   \item period: Factor that indicates the year/month of the data (Based on rnd).
#'   }
#'
#' @format Data frame with 2,500 rows and 22 columns with 500 defaults.
#' @name smbsimdf1
NULL
# End: Simulated Credit Data ##############################################

# Begin: Monotonic Sample Data ############################################
#' Monotonic Binning Sample Data 
#'
#' A simulated dataset used to illustrate the application of monotonic binning.
#'
#' \itemize{
#'   \item fgood1: Default (0), Not Default (1) for Numeric Variable 1.
#'   \item chr1: Numeric variable 1.
#'   \item fgood2: Default (0), Not Default (1) for Numeric Variable 2.
#'   \item chr2: Numeric variable 2.
#'   \item fgood3: Default (0), Not Default (1) for Numeric Variable 3.
#'   \item chr3: Numeric variable 3.
#'   }
#'
#' @format Data frame with 2,500 rows and 6 columns.
#' @name smbsimdf2
NULL
# End: Monotonic Sample Data ##############################################


# Begin: Model Ranking Sample Data ########################################
#' Monotonic Binning Sample Data 
#'
#' A simulated dataset used to illustrate the application of model ranking.
#'
#' \itemize{
#'   \item fgood1: Default (0), Not Default (1) for Numeric Variable 1.
#'   \item chr1: Numeric variable 1.
#'   \item chr2: Numeric variable 2.
#'   \item chr3: Numeric variable 3.
#'   }
#'
#' @format Data frame with 1,000 rows and 4 columns.
#' @name smbsimdf3
NULL
# End: Model Ranking Sample Data ##########################################


