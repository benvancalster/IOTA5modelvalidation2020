#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#
#   Interim Analysis IOTA5 - Paper 2    #
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%#


##########################################
#### 1. Import data and load packages ####
##########################################

## Paths ##
FindMap    = function(x, string) x[grepl(string, x, fixed = T)]
setwd("C:/Users/u0123496/Documents/IOTA/IOTA5")
CurWd      = getwd()
AllMaps    = list.dirs(gsub("/WorkingDirectoryIOTA5_Jolien", "", CurWd), recursive = F)
Scripts    = FindMap(AllMaps, "Scripts")
CSVfiles   = FindMap(AllMaps, "CSV")    
Cleaning   = FindMap(AllMaps, "Datacleaning")
RdataFiles = FindMap(AllMaps, "Rdata")
ValData    = FindMap(AllMaps, "Database validation")
Data       = FindMap(AllMaps, "Data")
AddInfo    = FindMap(AllMaps, "AdditionalInfo")
Results    = FindMap(AllMaps, "Results Paper 2")

## Source functions and load packages ##
source(paste0(Scripts, "/AdditionalFunctions.R"))
source(paste0(Scripts, "/functions iota5_desktop.R"))
source(paste0(Scripts, "/Functions IOTA5.R"))
LibraryM(c("readxl", "plyr", "REMA", "mice", "zoo", "Rcpp", "doParallel", "magrittr", "RcryptAdj", "grDevices"))

#load("IOTA5 - Validation 25032020.Rdata")
#load("IOTA5 - Datasets.RData")

## Import data ##
LoadRcrypt("IOTA5_AllCenters_260218BoxSync.Rdata.gpg", path = RdataFiles)
NrPtns(iota5.CA) # 9655
length(table(iota5.CA$center)) # 38 centers
IOTA5 <- iota5.CA[iota5.CA$center != "BSP" & iota5.CA$center != "BCN",]
length(table(IOTA5$center)) # 36 centers
NrPtns(IOTA5) # 8599


## Update the patient status ##
iota5.firstobs = get("iota5.firstobs", envir = .GlobalEnv)
BeginStudy = as.Date("2012-01-01")
df.tmp    = iota5.firstobs(IOTA5)
IDbefore2012  = df.tmp[which(df.tmp$`Exam date` < BeginStudy),"Patient ID"]
IOTA5$Before2012 = sapply(IOTA5$`Patient ID`,
                          function(x) if(x %in% IDbefore2012) T else F)
IOTA5 = ddply(IOTA5, .(Before2012, `Patient ID`),
              function(x) {
                if(all(x$Before2012)) {
                  N = nrow(x)
                  x$PatientStatusCorrect = rep("old patient", N)
                  x$`Followup time` = sapply(x$`Followup time`,
                                             function(y) {
                                               if(is.na(y))
                                                 y
                                               else if(y == "")
                                                 NA
                                               else
                                                 y
                                             })
                  
                  if(is.na(iota5.firstobs(x, dominant = T)$`Followup time`)) {
                    FirstTmp  = iota5.firstobs(x)
                    After2012 = x[x$`Exam date` > as.Date("2012-01-01"), ]
                    After2012 = After2012[with(After2012, order(`Exam date`)),]
                    After2012 = After2012[! duplicated(After2012$`Patient ID`),]
                    FUtime = (as.yearmon(strptime(After2012$`Exam date`, format="%Y-%m-%d"))-
                                as.yearmon(strptime(FirstTmp$`Exam date`, format="%Y-%m-%d"))) * 12
                    x$`Followup time` = rep(FUtime, N)
                  } else {
                    if(length(unique(na.omit(x$`Followup time`))) > 1)
                      stop("More than 1 unique FU time.")
                    x$`Followup time` = rep(unique(na.omit(x$`Followup time`)), N)
                  }
                } 
                return(x)
              }, .parallel = T, .paropts = list(.packages = "zoo", .export = "iota5.firstobs"))



#######################################
#### 2. In- and Exclusion criteria ####
#######################################

## Remove minors ##
MinorFirstScan = iota5.firstobs(IOTA5)
MinorID = unname(unlist(MinorFirstScan[which(MinorFirstScan$`Patient age` < 18), "Patient ID"]))
IOTA5 = IOTA5[!IOTA5$`Patient ID` %in% MinorID, ]
NrPtns(IOTA5) # n = 8519

## Remove patients that withdrew consent ##
RmPtns = unname(unlist(dlply(IOTA5, .(`Patient ID`),
                             function(x){
                               if(!allNA(x$`Study Outcome`) &&
                                  any(na.omit(x$`Study Outcome`)=="patient withdrew consent")){ 
                                 return(unique(x$`Patient ID`))
                               } else {
                                 return(NULL)
                               }
                             })))
IOTA5 = IOTA5[!IOTA5$`Patient ID` %in% RmPtns,]
NrPtns(IOTA5) # n = 8494

## Select included centers ##
SelectedCenters <- c("AGR", "CIT", "FLI", "GBE", "IUK", "LBE", "MIT", "MPO", "MSW", "NCI", "NUK", "OIT", "PSP", "RIT", "SIT", "SSW", "TIT")
ExtraCenters <- c("AGR", "CIT", "FLI", "GBE", "IUK", "LBE", "MIT", "MPO", "MSW", "NCI", "NUK", "OIT", "PSP", "RIT", "SIT", "SSW", "TIT", "BIT", "PCR", "LPO")
iota5Centers <- IOTA5[IOTA5$center %in% SelectedCenters,]
NrPtns(iota5Centers) # n = 5717

## Select new patients ##
iota5Val <- subset(iota5Centers, iota5Centers$PatientStatusCorrect != 'old patient')
NrPtns(iota5Val) # n = 4905

## Check whether recruitment is between 1 January 2012 and 1 March 2015 ##
BeginStudy = as.Date("2012-01-01")
EndRecruit = as.Date("2015-03-01")
df.tmp    = iota5.firstobs(iota5Val)
subset(df.tmp, as.Date(df.tmp$`Exam date`, "%Y-%m-%d") < BeginStudy | as.Date(df.tmp$`Exam date`, "%Y-%m-%d") > EndRecruit, select = "Patient ID")
df.tmp[as.Date(df.tmp$`Exam date`, "%Y-%m-%d") < BeginStudy | as.Date(df.tmp$`Exam date`, "%Y-%m-%d") > EndRecruit, "Patient ID"]
df.tmp[as.Date(df.tmp$`Exam date`, "%Y-%m-%d") < BeginStudy, "Patient ID"]
df.tmp[as.Date(df.tmp$`Exam date`, "%Y-%m-%d") > EndRecruit, "Patient ID"]

# Multiple masses: only the dominant masses
iota5Val <- iota5Val[iota5Val$multiple.masses == 1 & iota5Val$`Dominant mass` == 'yes' | iota5Val$multiple.masses == 0,]
NrPtns(iota5Val) # 4905

# Make dataset with only the first observation
iota5ValAnalyses <- iota5.firstobs(iota5Val)
NrPtns(iota5ValAnalyses) #4905

iota5ValImp <- iota5.firstobs(iota5Val)

#########################
#### 3. Descriptives ####
#########################

# Binary SA
iota5Val$SAbinary = sapply(iota5Val$Subjective,
                           function(x) {
                             if(is.na(x))
                               NA
                             else if(x == "benign")
                               "benign"
                             else
                               "malignant"
                           })

# Time until surgery
iota5Val = ddply(iota5Val, .(`Patient ID`),
                 function(x){
                   N = nrow(x)
                   
                   TimeSurgery = as.Date(iota5.lastobs(x)$`Date of Surgery`, format = "%Y-%m-%d") - as.Date(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d")
                   
                   x$TimeSurgery = rep(TimeSurgery, N)
                   return(x)
                 })
iota5Val[!is.na(iota5Val$TimeSurgery) & iota5Val$TimeSurgery > 1200, ]

# Follow-up time
iota5Val = ddply(iota5Val, .(`Patient ID`),
                 function(x){
                   N = nrow(x)
                   
                   FU = (as.yearmon(strptime(iota5.lastobs(x)$`Exam date`, format = "%Y-%m-%d")) - as.yearmon(strptime(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d"))) * 12
                   
                   x$FollowUp = rep(FU, N)
                   return(x)
                 })

iota5Val = iota5.fu(iota5Val, unit.fu = "months") # Werkt nu wel zonder missings

# Cumulative follow-up time
iota5Val = ddply(iota5Val, .(`Patient ID`),
                 function(x){
                   
                   FUcum = (as.yearmon(strptime(x$`Exam date`, format = "%Y-%m-%d")) - as.yearmon(strptime(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d"))) * 12
                   
                   x$cumulativeFollowUp = FUcum
                   return(x)
                 })
iota5Val[, c('Patient ID', 'cumulativeFollowUp')]
class(iota5Val$cumulativeFollowUp)

# Updated variable for Postmenopausal
PMstatus = ddply(iota5.firstobs(iota5Val, dominant = T),
                 .(Postmenopausal),
                 function(x) {
                   if (all(x$Postmenopausal != "uncertain")) {
                     x$Postmenopausal2 = x$Postmenopausal
                   } else {
                     x$Postmenopausal2 = sapply(x$`Patient age`,
                                                function(y) {
                                                  if(is.na(y))
                                                    NA
                                                  else if(y < 50)
                                                    "no"
                                                  else
                                                    "yes"
                                                })
                   }
                   return(x)
                 })

iota5Val$Postmenopausal2 = sapply(iota5Val$`Patient ID`,
                                  function(x) {
                                    PMstatus[PMstatus$`Patient ID` == x, "Postmenopausal2"]
                                  })


#------------------------------#
#### 3.1 Clinical diagnosis ####
#------------------------------#

# Determine outcome
iota5Val = ddply(iota5Val, .(`Patient ID`),
                 function(x){
                   N = nrow(x)
                   lastObs = iota5.lastobs(x)
                   
                   ClinicalDiagnosis =
                     if(!allNA(x$`Study Outcome`)){
                       # A
                       if(lastObs$`Study Outcome` == "surgery performed"){
                         if(!is.na(lastObs$`Mass Outcome`)){
                           # 1
                           if(lastObs$`Mass Outcome` == "benign"){
                             "Benign: surgery performed"
                           } # Einde if benign
                           else{
                             # 2
                             if(!is.na(lastObs$TimeSurgery)){
                               # 2a
                               if(lastObs$TimeSurgery <= 120){
                                 "Malignant: surgery performed <= 120 days inclusion"
                               } # Einde if time until surgery <= 120 days
                               else{
                                 if(anyNA(x$Subjective)){
                                   "Surgery (malignant, > 120 days), Missing SA"
                                 } # Einde if een missing in Subjective
                                 else{
                                   # 2b
                                   if(all(grepl("malignant",x$Certainty))){
                                     "Malignant: surgery performed after 120 days inclusion, all SA borderline/malignant"
                                   } # Einde if subjective is always malignant
                                   # 2c
                                   else{
                                     "Uncertain: surgery performed after 120 days inclusion, SA not always borderline/malignant"
                                   } # Einde else (not always malignant/borderline)
                                 } # Einde else (no missings in Subjective)
                               } # Einde else (not within 120 days)
                             } # Einde if not missing time until surgery
                             else{
                               "Surgery and malignant, but missing time until surgery"
                             } # Einde else (missing time until surgery)
                           } # Einde else (not benign)
                         } # Einde if not missing mass outcome
                         else{
                           "Surgery, but missing mass outcome"
                         } # Einde else (missing mass outcome)
                       } # Einde if surgery performed
                       else{
                         # B
                         if(lastObs$`Study Outcome` == "cyst spontaneously resolved"){
                           "Benign: spontaneous resolution"
                         } # Einde if spontaneously resolved
                         else{
                           # C
                           if(lastObs$cumulativeFollowUp >= 10){
                             if(!anyNA(x$Subjective)){
                               # 1
                               if(all(grepl("benign",x$Certainty[x$cumulativeFollowUp <= 14]))
                                  & grepl("benign",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                                 "Benign: no surgery, all SA during first 14 months are benign"
                               } # Einde if SA is benign
                               else{
                                 # 2
                                 if(all(grepl("malignant",x$Certainty[x$cumulativeFollowUp <= 14]))
                                    & grepl("malignant",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                                   "Malignant: no surgery, all SA during first 14 months are malignant"
                                 } # Einde if SA is malignant
                                 else{
                                   # 3
                                   "Uncertain: no surgery, SA during first 14 months is inconsistent"
                                 } # Einde else (SA is not malignant)
                               } # Einde else (SA is not benign)
                             } # Einde if no missings in subjective
                             else{
                               "No surgery, missing SA"
                             } # Einde else (missings in subjective)
                           } # Einde if follow-up >= 10 months
                           else{
                             if(lastObs$cumulativeFollowUp == 0){
                               paste0("Uncertain: No follow-up, study outcome ", lastObs$`Study Outcome`)
                             } # einde if no follow-up
                             else{
                               paste0("Uncertain: Follow-up < 10 months, study outcome ", lastObs$`Study Outcome`)
                             } # Einde else (follow-up < 10 months)
                           } # Einde else (follow-up < 10 months)
                         } # Einde else (cyst not spontaneously resolved)
                       } # Einde else (no surgery performed)
                     }# Einde if no missing study outcome
                   else{
                     if(lastObs$cumulativeFollowUp >= 10){
                       if(!anyNA(x$Subjective)){
                         # 1
                         if(all(grepl("benign",x$Certainty[x$cumulativeFollowUp <= 14]))
                            & grepl("benign",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                           "Benign: no surgery, all SA during first 14 months are benign"
                         } # Einde if SA is benign
                         else{
                           # 2
                           if(all(grepl("malignant",x$Certainty[x$cumulativeFollowUp <= 14]))
                              & grepl("malignant",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                             "Malignant: no surgery, all SA during first 14 months are malignant"
                           } # Einde if SA is malignant
                           else{
                             # 3
                             "Uncertain: no surgery, SA during first 14 months is inconsistent"
                           } # Einde else (SA is not malignant)
                         } # Einde else (SA is not benign)
                       } # Einde if no missings in subjective
                       else{
                         "No surgery, missing SA"
                       } # Einde else (missings in subjective)
                     } # Einde if follow-up >= 10 months
                     else{
                       if(lastObs$cumulativeFollowUp == 0){
                         paste0("Uncertain: No follow-up, study outcome ", if(is.na(lastObs$`Study Outcome`)) "Missing" else lastObs$`Study Outcome`)
                       } # Einde if no follow-up
                       else{
                         paste0("Uncertain: Follow-up < 10 months, study outcome ", if(is.na(lastObs$`Study Outcome`)) "Missing" else lastObs$`Study Outcome`)
                       } # Einde else (follow-up < 10 months)
                     } # Einde else (follow-up < 10 months)
                   } # Einde else (missing in study outcome)
                   
                   x$ClinicalDiagnosis = rep(ClinicalDiagnosis, N)
                   return(x)
                 }# Einde function
)

NrPtns(iota5Val)
iotaValFirst = iota5.firstobs(iota5Val)
table(iotaValFirst$ClinicalDiagnosis)

# Number of examiners in the dataset
table(iota5Val$Operator)
length(unique(iota5Val$Operator)) # 76 (one already excluded, because same examiner had 2 different names)
length(unique(iota5ValAnalyses$Operator)) # 58 (one already excluded, because same examiner had 2 different names)
table(iota5ValAnalyses$Operator)

# Make dataset with only the first observation
iota5ValAnalyses <- iota5.firstobs(iota5Val)
NrPtns(iota5ValAnalyses) #4905

## Variable with the categories in short
iota5ValAnalyses$CD_short <- ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Benign: surgery performed", "B1",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Malignant: surgery performed <= 120 days inclusion", "M1",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Malignant: surgery performed after 120 days inclusion, all SA borderline/malignant", "M2",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Uncertain: surgery performed after 120 days inclusion, SA not always borderline/malignant", "U1",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Benign: spontaneous resolution", "B0",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Benign: no surgery, all SA during first 14 months are benign", "B2",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Malignant: no surgery, all SA during first 14 months are malignant", "M3",
                                    ifelse(iota5ValAnalyses$ClinicalDiagnosis == "Uncertain: no surgery, SA during first 14 months is inconsistent", "U2",
                                    ifelse(grepl("Uncertain: Follow-up < 10 months", iota5ValAnalyses$ClinicalDiagnosis), "U3", "U4")))))))))
table(iota5ValAnalyses$CD_short)


## Save excel file ##
head(iota5Val)
iota5ValCD <-  sapply(iota5Val, iconv, to = "UTF-8")
iota5ValCD <- data.frame(iota5ValCD)
AllDiagn = dlply(iota5ValCD, .(ClinicalDiagnosis))
require(openxlsx)
list_of_datasets <- list("BenignSurgery" = AllDiagn[[3]], "MalignantSurgeryWithin120" = AllDiagn[[5]], "MalignantSurgeryMore120" = AllDiagn[[6]], "UncertainSurgery" = AllDiagn[[16]],
                         "BenignSpontaneousResolution" = AllDiagn[[2]],
                         "BenignNoSurgery" = AllDiagn[[1]], "MalignantNoSurgery" = AllDiagn[[4]], "UncertainNoSurgery" = AllDiagn[[15]],
                         "UncertainNoFUlost" = AllDiagn[[11]], "UncertainNoFUmissing" = AllDiagn[[12]], "UncertainNoFUdied" = AllDiagn[[13]], "UncertainNoFUstopped" = AllDiagn[[14]],
                         "UncertainLimitedFUlost" = AllDiagn[[7]], "UncertainLimitedFUmissing" = AllDiagn[[8]], "UncertainLimitedFUdied" = AllDiagn[[9]], "UncertainLimitedFUstopped" = AllDiagn[[10]])
write.xlsx(list_of_datasets, file = "Results Paper 2/With cysts/IOTA5ClinicalDiagnosis 08-07-2019.xlsx")


## Binary CD
iota5ValAnalyses$binaryCD <- sapply(iota5ValAnalyses$ClinicalDiagnosis,
                                    function(x){
                                      if(grepl("Malignant", x))
                                        1
                                      else if(grepl("Benign", x))
                                        0
                                      else
                                        NA
                                    })


#-------------------------------#
#### 3.2 List of the centers ####
#-------------------------------#

IOTA5Lists = IOTA5

## Select new patients ##
IOTA5Lists <- subset(IOTA5Lists, IOTA5Lists$PatientStatusCorrect != 'old patient')
NrPtns(IOTA5Lists) # n = 7329

## Multiple masses: only the dominant masses
IOTA5Lists <- IOTA5Lists[IOTA5Lists$multiple.masses == 1 & IOTA5Lists$`Dominant mass` == 'yes' | IOTA5Lists$multiple.masses == 0,]
IOTA5Lists <- IOTA5Lists[!is.na(IOTA5Lists$`Patient ID`),]
NrPtns(IOTA5Lists) # 7329

## Time until surgery
IOTA5Lists = ddply(IOTA5Lists, .(`Patient ID`),
                   function(x){
                     N = nrow(x)
                     
                     TimeSurgery = as.Date(iota5.lastobs(x)$`Date of Surgery`, format = "%Y-%m-%d") - as.Date(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d")
                     
                     x$TimeSurgery = rep(TimeSurgery, N)
                     return(x)
                   })
IOTA5Lists[!is.na(IOTA5Lists$TimeSurgery) & IOTA5Lists$TimeSurgery > 1200, ]

## Follow-up time
IOTA5Lists = ddply(IOTA5Lists, .(`Patient ID`),
                   function(x){
                     N = nrow(x)
                     
                     FU = (as.yearmon(strptime(iota5.lastobs(x)$`Exam date`, format = "%Y-%m-%d")) - as.yearmon(strptime(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d"))) * 12
                     
                     x$FollowUp = rep(FU, N)
                     return(x)
                   })

## Cumulative follow-up time
IOTA5Lists = ddply(IOTA5Lists, .(`Patient ID`),
                   function(x){
                     FUcum = (as.yearmon(strptime(x$`Exam date`, format = "%Y-%m-%d")) - as.yearmon(strptime(iota5.firstobs(x)$`Exam date`, format = "%Y-%m-%d"))) * 12
                     
                     x$cumulativeFollowUp = FUcum
                     return(x)
                   })

## Updated variable for Postmenopausal
PMstatus = ddply(iota5.firstobs(IOTA5Lists, dominant = T),
                 .(Postmenopausal),
                 function(x) {
                   if (all(x$Postmenopausal != "uncertain")) {
                     x$Postmenopausal2 = x$Postmenopausal
                   } else {
                     x$Postmenopausal2 = sapply(x$`Patient age`,
                                                function(y) {
                                                  if(is.na(y))
                                                    NA
                                                  else if(y < 50)
                                                    "no"
                                                  else
                                                    "yes"
                                                })
                   }
                   return(x)
                 })

IOTA5Lists$Postmenopausal2 = sapply(IOTA5Lists$`Patient ID`,
                                    function(x) {
                                      PMstatus[PMstatus$`Patient ID` == x, "Postmenopausal2"]
                                    })

## Binary SA
IOTA5Lists$SAbinary = sapply(IOTA5Lists$Subjective,
                             function(x) {
                               if(is.na(x))
                                 NA
                               else if(x == "benign")
                                 0
                               else
                                 1
                             })

## Clinical Diagnosis
IOTA5Lists <- IOTA5Lists[!is.na(IOTA5Lists$`Patient ID`),]

IOTA5Lists = ddply(IOTA5Lists, .(`Patient ID`),
                   function(x){
                     N = nrow(x)
                     lastObs = iota5.lastobs(x)
                     
                     ClinicalDiagnosis =
                       if(!allNA(x$`Study Outcome`)){
                         # A
                         if(lastObs$`Study Outcome` == "surgery performed"){
                           if(!is.na(lastObs$`Mass Outcome`)){
                             # 1
                             if(lastObs$`Mass Outcome` == "benign"){
                               "Benign: surgery performed"
                             } # Einde if benign
                             else{
                               # 2
                               if(!is.na(lastObs$TimeSurgery)){
                                 # 2a
                                 if(lastObs$TimeSurgery <= 120){
                                   "Malignant: surgery performed <= 120 days inclusion"
                                 } # Einde if time until surgery <= 120 days
                                 else{
                                   if(anyNA(x$Subjective)){
                                     "Surgery (malignant, > 120 days), Missing SA"
                                   } # Einde if een missing in Subjective
                                   else{
                                     # 2b
                                     if(all(grepl("malignant",x$Certainty))){
                                       "Malignant: surgery performed after 120 days inclusion, all SA borderline/malignant"
                                     } # Einde if subjective is always malignant
                                     # 2c
                                     else{
                                       "Uncertain: surgery performed after 120 days inclusion, SA not always borderline/malignant"
                                     } # Einde else (not always malignant/borderline)
                                   } # Einde else (no missings in Subjective)
                                 } # Einde else (not within 120 days)
                               } # Einde if not missing time until surgery
                               else{
                                 "Surgery and malignant, but missing time until surgery"
                               } # Einde else (missing time until surgery)
                             } # Einde else (not benign)
                           } # Einde if not missing mass outcome
                           else{
                             "Surgery, but missing mass outcome"
                           } # Einde else (missing mass outcome)
                         } # Einde if surgery performed
                         else{
                           # B
                           if(lastObs$`Study Outcome` == "cyst spontaneously resolved"){
                             "Benign: spontaneous resolution"
                           } # Einde if spontaneously resolved
                           else{
                             # C
                             if(lastObs$cumulativeFollowUp >= 10){
                               if(!anyNA(x$Subjective)){
                                 # 1
                                 if(all(grepl("benign",x$Certainty[x$cumulativeFollowUp <= 14]))
                                    & grepl("benign",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                                   "Benign: no surgery, all SA during first 14 months are benign"
                                 } # Einde if SA is benign
                                 else{
                                   # 2
                                   if(all(grepl("malignant",x$Certainty[x$cumulativeFollowUp <= 14]))
                                      & grepl("malignant",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                                     "Malignant: no surgery, all SA during first 14 months are malignant"
                                   } # Einde if SA is malignant
                                   else{
                                     # 3
                                     "Uncertain: no surgery, SA during first 14 months is inconsistent"
                                   } # Einde else (SA is not malignant)
                                 } # Einde else (SA is not benign)
                               } # Einde if no missings in subjective
                               else{
                                 "No surgery, missing SA"
                               } # Einde else (missings in subjective)
                             } # Einde if follow-up >= 10 months
                             else{
                               if(lastObs$cumulativeFollowUp == 0){
                                 paste0("Uncertain: No follow-up, study outcome ", lastObs$`Study Outcome`)
                               } # Einde if no follow-up
                               else{
                                 paste0("Uncertain: Follow-up < 10 months, study outcome ", lastObs$`Study Outcome`)
                               } # Einde else (follow-up < 10 months)
                             } # Einde else (follow-up < 10 months)
                           } # Einde else (cyst not spontaneously resolved)
                         } # Einde else (no surgery performed)
                       }# Einde if no missing study outcome
                     else{
                       if(lastObs$cumulativeFollowUp >= 10){
                         if(!anyNA(x$Subjective)){
                           # 1
                           if(all(grepl("benign",x$Certainty[x$cumulativeFollowUp <= 14]))
                              & grepl("benign",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                             "Benign: no surgery, all SA during first 14 months are benign"
                           } # Einde if SA is benign
                           else{
                             # 2
                             if(all(grepl("malignant",x$Certainty[x$cumulativeFollowUp <= 14]))
                                & grepl("malignant",x$Certainty[x$cumulativeFollowUp >= 10][1])){
                               "Malignant: no surgery, all SA during first 14 months are malignant"
                             } # Einde if SA is malignant
                             else{
                               # 3
                               "Uncertain: no surgery, SA during first 14 months is inconsistent"
                             } # Einde else (SA is not malignant)
                           } # Einde else (SA is not benign)
                         } # Einde if no missings in subjective
                         else{
                           "No surgery, missing SA"
                         } # Einde else (missings in subjective)
                       } # Einde if follow-up >= 10 months
                       else{
                         if(lastObs$cumulativeFollowUp == 0){
                           paste0("Uncertain: No follow-up, study outcome ", if(is.na(lastObs$`Study Outcome`)) "Missing" else lastObs$`Study Outcome`)
                         } # Einde if no follow-up
                         else{
                           paste0("Uncertain: Follow-up < 10 months, study outcome ", if(is.na(lastObs$`Study Outcome`)) "Missing" else lastObs$`Study Outcome`)
                         } # Einde else (follow-up < 10 months)
                       } # Einde else (follow-up < 10 months)
                     } # Einde else (missing in study outcome)
                     
                     x$ClinicalDiagnosis = rep(ClinicalDiagnosis, N)
                     return(x)
                   }# Einde function
)
NrPtns(IOTA5Lists) # n = 7329

## Binary CD
IOTA5Lists$binaryCD <- sapply(IOTA5Lists$ClinicalDiagnosis,
                              function(x){
                                if(grepl("Malignant", x))
                                  1
                                else if(grepl("Benign", x))
                                  0
                                else
                                  NA
                              })


IOTA5Descrip <- iota5.firstobs(IOTA5Lists)
table(IOTA5Descrip$ClinicalDiagnosis)

## Oncological vs non-oncological
onco    = c("LPO", "SSW", "PCR", "CIT", "BCH", "OIT", "UDI", "BIT", "RIT", "LBE", "PSP",
            "AGR", "MPO", "CAI", "NCI", "DEP", "BAI", "KPO", "LIP", "VAS")
nononco = c("SIT", "MIT", "GBE", "MSW", "CEG", "FIT", "FLI", "BSP", "TUS", "NUK", "CRI", "TIT",
            "MCA", "MFR", "PFR", "RZT", "IUK") 

IOTA5Descrip$oncocenter <- sapply(IOTA5Descrip$center,
                                  function(x){
                                    if(x %in% onco)
                                      1
                                    else if(x %in% nononco)
                                      0
                                    else
                                      NA
                                  })


## Description for each center
library(bazar)
library(openxlsx) # creation of Excel .xlsx files

IOTA5Descrip$CD_short <- ifelse(IOTA5Descrip$ClinicalDiagnosis == "Benign: surgery performed", "B1",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Malignant: surgery performed <= 120 days inclusion", "M1",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Malignant: surgery performed after 120 days inclusion, all SA borderline/malignant", "M2",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Uncertain: surgery performed after 120 days inclusion, SA not always borderline/malignant", "U1",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Benign: spontaneous resolution", "B0",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Benign: no surgery, all SA during first 14 months are benign", "B2",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Malignant: no surgery, all SA during first 14 months are malignant", "M3",
                                ifelse(IOTA5Descrip$ClinicalDiagnosis == "Uncertain: no surgery, SA during first 14 months is inconsistent", "U2",
                                ifelse(grepl("Uncertain: Follow-up < 10 months", IOTA5Descrip$ClinicalDiagnosis), "U3", "U4")))))))))
table(IOTA5Descrip$CD_short)

CentersDescrip = ddply(IOTA5Descrip, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_oper_mal = sum(x$immed.oper == 1 & grepl("Malignant", x$ClinicalDiagnosis)),
                         n_oper_unc = sum(x$immed.oper == 1 & grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_Benign = sum(x$multiple.visits == 1 & grepl("Benign", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                         n_Maligne = sum(x$multiple.visits == 1 & grepl("Malignant", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                         n_Uncertain = sum(x$multiple.visits == 1 & grepl("Uncertain", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip
write.xlsx(CentersDescrip, "Results Paper 2/With cysts/List Centers - No inconsistency.xlsx")

CentersDescrip = ddply(IOTA5Descrip, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip
write.xlsx(CentersDescrip, "Results Paper 2/With cysts/List Centers - Manuscript.xlsx")


#----------------------------------#
#### 3.3 Descriptive statistics ####
#----------------------------------#

N <- nrow(iota5ValAnalyses)

## Patient age
iota5ValAnalyses[is.na(iota5ValAnalyses$`Patient age`), c("Patient ID", "Patient DOB", "Exam date", "Patient age")] # Missing for 2 patients: FLI44 and MSW585
library(eeptools)
iota5ValAnalyses[is.na(iota5ValAnalyses$`Patient age`), c("Patient age")] = floor(age_calc(as.Date(iota5ValAnalyses$`Patient DOB`[is.na(iota5ValAnalyses$`Patient age`)], format="%Y-%m-%d"), as.Date(iota5ValAnalyses$`Exam date`[is.na(iota5ValAnalyses$`Patient age`)], format="%Y-%m-%d"), units = "years"))
iota5ValAnalyses[iota5ValAnalyses$`Patient ID` == "FLI44" | iota5ValAnalyses$`Patient ID` == "MSW585", c("Patient ID", "Patient DOB", "Exam date", "Patient age")]

median(iota5ValAnalyses$`Patient age`, na.rm = T) # 48
quantile(iota5ValAnalyses$`Patient age`, na.rm = T, probs = 0.25) # 36
quantile(iota5ValAnalyses$`Patient age`, na.rm = T, probs = 0.75) # 62
min(iota5ValAnalyses$`Patient age`, na.rm = T) # 18
max(iota5ValAnalyses$`Patient age`, na.rm = T) # 98

mean(iota5ValAnalyses$`Patient age`, na.rm = T)

## Postmenopausal
Postmeno <- count(iota5.firstobs(iota5Val), 'Postmenopausal2') # no: 2754 and yes: 2151

Postmeno$percent <- round(Postmeno$freq / N * 100, 1)
Postmeno

## Gynaecological symptoms during the year preceding inclusion
Gynae <- count(iota5ValAnalyses$`Symptoms last year`) # no: 2340 and yes: 2565

Gynae$percent <- round(Gynae$freq / N * 100, 1)
Gynae

## Tumour Type using IOTA terminology
TumourType <- count(iota5ValAnalyses$`Tumour type`) # multilocular: 1011, multilocular-solid: 649, solid: 689, unclassified: 20, unilocular: 2140 and unilocular-solid: 396

TumourType$percent <- round(TumourType$freq / N * 100, 1)
TumourType

## Presence of solid components
iota5ValAnalyses$SolidComp <- ifelse(iota5ValAnalyses$`Solid: measurement 1` == 0 & iota5ValAnalyses$`Solid: measurement 2` == 0 & iota5ValAnalyses$`Solid: measurement 3` == 0, "no", "yes")
Solid <- count(iota5ValAnalyses$SolidComp) # no: 3168 and yes: 1737 
Solid <- count(iota5ValAnalyses$solidbin) # 0: 3171 and 1: 1734

Solid$percent <- round(Solid$freq / N * 100, 1)
Solid

## Observed CA-125
median(iota5ValAnalyses$CA125, na.rm = T) # 25
quantile(iota5ValAnalyses$CA125, probs = 0.25, na.rm = T) # 12
quantile(iota5ValAnalyses$CA125, probs = 0.75, na.rm = T) # 107
min(iota5ValAnalyses$CA125, na.rm = T) # 1
max(iota5ValAnalyses$CA125, na.rm = T) # 57900

mean(iota5ValAnalyses$CA125, na.rm = T)

## Maximum diameter of lesion (mm)
median(iota5ValAnalyses$`Lesion largest diameter`) # 55
quantile(iota5ValAnalyses$`Lesion largest diameter`, probs = 0.25) # 38
quantile(iota5ValAnalyses$`Lesion largest diameter`, probs = 0.75) # 83
min(iota5ValAnalyses$`Lesion largest diameter`) # 7
max(iota5ValAnalyses$`Lesion largest diameter`) # 751

iota5ValAnalyses[iota5ValAnalyses$`Lesion largest diameter` > 300, c("Patient ID", "Lesion largest diameter")]

## Bilateral masses
Bilateral <- count(iota5ValAnalyses$Bilaterality) # bilateral: 829 and unilateral: 4076

Bilateral$percent <- round(Bilateral$freq / N * 100, 1)
Bilateral

## Ultrasound examiner's subjective impression
SA <- count(iota5ValAnalyses$Certainty) # Certainly benign: 2488, certainly malignant: 592, probably benign: 1066, probably malignant: 392 and uncertain: 367

SA$percent <- round(SA$freq / N * 100, 1)
SA

## Ultrasound examiner's presumed diagnosis
PreDiagnosis <- count(iota5ValAnalyses$`Presumed diagnosis`)

PreDiagnosis$percent <- round(PreDiagnosis$freq / N * 100, 1)
PreDiagnosis

## Colour score
ColorScore <- count(iota5ValAnalyses$`Colour score`) # 1: 2031, 2: 1336, 3: 1099, 4: 439

ColorScore$percent <- round(ColorScore$freq / N * 100, 1)
ColorScore

## Max diameter of solid area (mm), median
median(iota5ValAnalyses$soldmax)

## Proportion solid tissue, median
median(iota5ValAnalyses$propsol)

## Irregular internal cyst walls
table(iota5ValAnalyses$irregbin)

## Acoustic shadows
table(iota5ValAnalyses$shadowsbin)

## Ascites
table(iota5ValAnalyses$ascitesbin)

## Number of papillations, mean
mean(iota5ValAnalyses$papnr)

## Papillations with blood flow
table(iota5ValAnalyses$papflow)

## Multilocular cyst
table(iota5ValAnalyses$multibin)

## >10 locules
table(iota5ValAnalyses$loc10)

## Abdominal metastases
table(iota5ValAnalyses$metasbin)

## Solid areas, but smaller than 7mm (b2)
nrow(iota5ValAnalyses[iota5ValAnalyses$soldmax > 0 & iota5ValAnalyses$soldmax < 7,])
table(Iota5ValNI$SR_b2)

## Smooth multilocular cyst <100mm (b4)
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular" & iota5ValAnalyses$`Lesion largest diameter` < 100 & iota5ValAnalyses$irregbin==0,])
table(Iota5ValNI$SR_b4)

## Irregular solid tumor (m1)
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="solid" & iota5ValAnalyses$irregbin==1,])
table(Iota5ValNI$SR_m1)

## At least 4 papillations (m3)
table(iota5ValAnalyses$papnr)
table(Iota5ValNI$SR_m3)

## Irregular multilocular-solid cyst ???100mm (m4)
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular-solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$`Lesion largest diameter`>=100,])
table(Iota5ValNI$SR_m4)


#### According to missingness ####
iota5ValAnalyses$Bmissing <- ifelse(grepl("Uncertain", iota5ValAnalyses$ClinicalDiagnosis), "CD missing", "CD present")

## Outcome by missingness of CA125/outcome
table(iota5ValAnalyses$CD_short, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$CD_short, iota5ValAnalyses$Bmissing)

## Management by missingness of CA125/outcome
table(iota5ValAnalyses$immed.oper, iota5ValAnalyses$CA125missing)
table(iota5ValAnalyses$multiple.visits, iota5ValAnalyses$CA125missing)
table(iota5ValAnalyses[iota5ValAnalyses$immed.oper == 0 & iota5ValAnalyses$multiple.visits == 0, c("CA125missing")])

table(iota5ValAnalyses$immed.oper, iota5ValAnalyses$Bmissing)
table(iota5ValAnalyses$multiple.visits, iota5ValAnalyses$Bmissing)
table(iota5ValAnalyses[iota5ValAnalyses$immed.oper == 0 & iota5ValAnalyses$multiple.visits == 0, c("Bmissing")])

## Age by missingness of CA125/outcome
mean(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("Patient age")]) ## 52
mean(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("Patient age")]) ## 47

mean(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("Patient age")]) ## 49
mean(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("Patient age")]) ## 51

## Postmenopausal by missingness of CA125/outcome
table(iota5ValAnalyses$Postmenopausal2, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$Postmenopausal2, iota5ValAnalyses$Bmissing)

## CA125 by missingness of CA125/outcome
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("CA125")]) ## 25
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("CA125")]) ## -

table(iota5ValAnalyses$CA125missing[iota5ValAnalyses$CA125missing == "CA125 present"])
table(iota5ValAnalyses$CA125missing[iota5ValAnalyses$CA125missing == "CA125 missing"])

median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("CA125")], na.rm = TRUE) ## 25
median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("CA125")], na.rm = TRUE) ## 19

table(iota5ValAnalyses$CA125missing[iota5ValAnalyses$Bmissing == "CD present"])
table(iota5ValAnalyses$CA125missing[iota5ValAnalyses$Bmissing == "CD missing"])

## Maximum diameter of lesion by missingness of CA125/outcome
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("Lesion largest diameter")]) ## 67
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("Lesion largest diameter")]) ## 47

median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("Lesion largest diameter")]) ## 56
median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("Lesion largest diameter")]) ## 46

## Maximum diameter of solid area by missingness of CA125/outcome
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("soldmax")]) ## 3
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("soldmax")]) ## 0

median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("soldmax")]) ## 0
median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("soldmax")]) ## 0

## Proportion solid tissue by missingness of CA125/outcome
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("propsol")]) ## 0.05
median(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("propsol")]) ## 0

median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("propsol")]) ## 0
median(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("propsol")]) ## 0

## Presence of solid areas by missingness of CA125/outcome
table(iota5ValAnalyses$solidbin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$solidbin, iota5ValAnalyses$Bmissing)

## Irregular internal cyst walls by missingness of CA125/outcome
table(iota5ValAnalyses$irregbin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$irregbin, iota5ValAnalyses$Bmissing)

## Acoustic shadows by missingness of CA125/outcome
table(iota5ValAnalyses$shadowsbin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$shadowsbin, iota5ValAnalyses$Bmissing)

## Ascites by missingness of CA125/outcome
table(iota5ValAnalyses$ascitesbin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$ascitesbin, iota5ValAnalyses$Bmissing)

## Number of papillations by missingness of CA125/outcome
mean(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 present", c("papnr")]) ## 0.36
mean(iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("papnr")]) ## 0.15

mean(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD present", c("papnr")]) ## 0.25
mean(iota5ValAnalyses[iota5ValAnalyses$Bmissing == "CD missing", c("papnr")]) ## 0.19

## Papillations with blood flow by missingness of CA125/outcome
table(iota5ValAnalyses$papflow, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$papflow, iota5ValAnalyses$Bmissing)

## Bilateral by missingness of CA125/outcome
table(iota5ValAnalyses$Bilaterality, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$Bilaterality, iota5ValAnalyses$Bmissing)

## Multilocular cyst by missingness of CA125/outcome
table(iota5ValAnalyses$multibin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$multibin, iota5ValAnalyses$Bmissing)

## >10 locules by missingness of CA125/outcome
table(iota5ValAnalyses$loc10, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$loc10, iota5ValAnalyses$Bmissing)

## Abdominal metastases by missingness of CA125/outcome
table(iota5ValAnalyses$metasbin, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$metasbin, iota5ValAnalyses$Bmissing)

## Unilocular cyst (b1) by missingness of CA125/outcome
table(iota5ValAnalyses$`Tumour type`, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$`Tumour type`, iota5ValAnalyses$Bmissing)

## Solid areas, but smaller than 7mm (b2) by missingness of CA125/outcome
nrow(iota5ValAnalyses[iota5ValAnalyses$soldmax > 0 & iota5ValAnalyses$soldmax < 7 & iota5ValAnalyses$CA125missing == "CA125 present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$soldmax > 0 & iota5ValAnalyses$soldmax < 7 & iota5ValAnalyses$CA125missing == "CA125 missing",])

table(Iota5ValNI$SR_b2, Iota5ValNI$CA125missing)

nrow(iota5ValAnalyses[iota5ValAnalyses$soldmax > 0 & iota5ValAnalyses$soldmax < 7 & iota5ValAnalyses$Bmissing == "CD present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$soldmax > 0 & iota5ValAnalyses$soldmax < 7 & iota5ValAnalyses$Bmissing == "CD missing",])

table(Iota5ValNI$SR_b2, Iota5ValNI$Bmissing)

## Smooth multilocular cyst <100mm (b4) by missingness of CA125/outcome
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular" & iota5ValAnalyses$`Lesion largest diameter` < 100 & iota5ValAnalyses$irregbin==0 & iota5ValAnalyses$CA125missing == "CA125 present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular" & iota5ValAnalyses$`Lesion largest diameter` < 100 & iota5ValAnalyses$irregbin==0 & iota5ValAnalyses$CA125missing == "CA125 missing",])

table(Iota5ValNI$SR_b4, Iota5ValNI$CA125missing)

nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular" & iota5ValAnalyses$`Lesion largest diameter` < 100 & iota5ValAnalyses$irregbin==0 & iota5ValAnalyses$Bmissing == "CD present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular" & iota5ValAnalyses$`Lesion largest diameter` < 100 & iota5ValAnalyses$irregbin==0 & iota5ValAnalyses$Bmissing == "CD missing",])

table(Iota5ValNI$SR_b4, Iota5ValNI$Bmissing)

## Color score by missingness of CA125/outcome
table(iota5ValAnalyses$`Colour score`, iota5ValAnalyses$CA125missing)

table(iota5ValAnalyses$`Colour score`, iota5ValAnalyses$Bmissing)

## Irregular solid tumor (m1) by missingness of CA125/outcome
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$CA125missing == "CA125 present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$CA125missing == "CA125 missing",])

table(Iota5ValNI$SR_m1, Iota5ValNI$CA125missing)

nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$Bmissing == "CD present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$Bmissing == "CD missing",])

table(Iota5ValNI$SR_m1, Iota5ValNI$Bmissing)

## At least 4 papillations (m3) by missingness of CA125/outcome
table(iota5ValAnalyses$papnr, iota5ValAnalyses$CA125missing)
table(Iota5ValNI$SR_m3, Iota5ValNI$CA125missing)

table(iota5ValAnalyses$papnr, iota5ValAnalyses$Bmissing)
table(Iota5ValNI$SR_m3, Iota5ValNI$Bmissing)

## Irregular multilocular-solid cyst ???100mm (m4) by missingness of CA125/outcome
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular-solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$`Lesion largest diameter`>=100 & iota5ValAnalyses$CA125missing == "CA125 present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular-solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$`Lesion largest diameter`>=100 & iota5ValAnalyses$CA125missing == "CA125 missing",])

table(Iota5ValNI$SR_m4, Iota5ValNI$CA125missing)

nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular-solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$`Lesion largest diameter`>=100 & iota5ValAnalyses$Bmissing == "CD present",])
nrow(iota5ValAnalyses[iota5ValAnalyses$`Tumour type`=="multilocular-solid" & iota5ValAnalyses$irregbin==1 & iota5ValAnalyses$`Lesion largest diameter`>=100 & iota5ValAnalyses$Bmissing == "CD missing",])

table(Iota5ValNI$SR_m4, Iota5ValNI$Bmissing)


#-------------------------------------#
#### 3.4 Malignancy vs tumour type ####
#-------------------------------------#

TableMal <- ddply(iota5ValAnalyses, .(`Tumour type`),
                  function(x){
                    data.frame(N = nrow(x),
                               n_with_outcome = nrow(x) - sum(is.na(x$binaryCD)),
                               n_malignant = count(x$binaryCD[x$binaryCD == 1])[1, 2],
                               malignancy_rate = round(count(x$binaryCD[x$binaryCD == 1])[1, 2] / (nrow(x) - sum(is.na(x$binaryCD))), 4),
                               check.names = F)
                  })
TableMal

## Description for Tumor type
DescripTT = ddply(iota5ValAnalyses, .(`Tumour type`),
                  function(x){
                    data.frame(N = nrow(x),
                               n_surgery = sum(x$immed.oper),
                               n_oper_mal = sum(x$immed.oper == 1 & grepl("Malignant", x$ClinicalDiagnosis)),
                               n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                               n_Benign = sum(x$multiple.visits == 1 & grepl("Benign", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                               n_Maligne = sum(x$multiple.visits == 1 & grepl("Malignant", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                               n_Uncertain = sum(x$multiple.visits == 1 & grepl("Uncertain", x$ClinicalDiagnosis) & x$CD_short != "U4"),
                               n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                               check.names = F)
                  })
DescripTT


################################
#### 4. Multiple Imputation ####
################################

#---------------------#
#### 4.1 Variables ####
#---------------------#

#### 4.1.1 Additional variables ####

## Onco- vs non-oncological centers
onco    = c("LPO", "SSW", "PCR", "CIT", "BCH", "OIT", "UDI", "BIT", "RIT", "LBE", "PSP",
            "AGR", "MPO", "CAI", "NCI", "DEP")
nononco = c("SIT", "MIT", "GBE", "MSW", "CEG", "FIT", "FLI", "BSP", "TUS", "NUK", "CRI",
            "TIT", "IUK") 

iota5ValAnalyses$oncocenter <- sapply(iota5ValAnalyses$center,
                                      function(x){
                                        if(x %in% onco)
                                          1
                                        else if(x %in% nononco)
                                          0
                                        else
                                          NA
                                      })

## Presumed endometrioma
iota5ValAnalyses$PresumedEndometrioma = sapply(iota5ValAnalyses$`Presumed diagnosis`,
                                               function(x){
                                                 if(is.na(x))
                                                   NA
                                                 else if(x == "endometrioma")
                                                   1
                                                 else
                                                   0
                                               })

# Origin diagnosis
iota5ValAnalyses$ClinicalDiagnosis
iota5ValAnalyses = ddply(iota5ValAnalyses, .(`Patient ID`),
                         function(x){
                           x$OriginCD =
                             if(is.na(x$`Study Outcome`)){
                               "Clinical diagnosis"
                             }
                           else{
                             if(x$`Study Outcome` == "surgery performed"){
                               if(grepl("Uncertain", x$ClinicalDiagnosis)){
                                 "Clinical diagnosis"
                               } else{
                                 "Histological analysis"
                               }
                             } 
                             else{
                               "Clinical diagnosis"
                             }
                           }
                           return(x)
                         })

# Outcome with 5 groups
table(iota5ValAnalyses$`Mass Outcome`)
table(iota5ValAnalyses$OriginCD, iota5ValAnalyses$`Mass Outcome`)
table(iota5ValAnalyses$OriginCD)
iota5ValAnalyses[iota5ValAnalyses$OriginCD == "Histological analysis" & is.na(iota5ValAnalyses$`Mass Outcome`), c('OriginCD', 'Mass Outcome')]

iota5ValAnalyses = ddply(iota5ValAnalyses, .(`Patient ID`),
                         function(x){
                           x$CD_5groups = 
                             if(!is.na(x$OriginCD)){
                               if(x$OriginCD == "Histological analysis"){
                                 if(x$`Mass Outcome` == 'benign'){
                                   "Benign"
                                 } # End if benign mass outcome
                                 else{
                                   if(grepl("invasive", x$`Mass Outcome`)){
                                     if(is.na(x$`FIGO stage`)){
                                       "Invasive FIGO stage missing"
                                     } # End if FIGO is missing
                                     else{
                                       if(x$`FIGO stage` == "I"){
                                         "Stage I invasive"
                                       } # End if FIGO is I
                                       else{
                                         "Stage II-IV invasive"
                                       } # End else (FIGO not I)
                                     } # End else (FIGO not missing)
                                   } # End if invasive mass outcome
                                   else{
                                     if(grepl("borderline", x$`Mass Outcome`)){
                                       "Borderline"
                                     } # End if borderline mass outcome
                                     else{
                                       "Metastatic"
                                     } # End else (no borderline mass outcome)
                                   } # End else (no invasive mass outcome)
                                 } # End else (no benign mass outcome)
                               } # End if Histological analysis
                               else{
                                 if(is.na(x$binaryCD)){
                                   NA
                                 } # End if binaryCD missing
                                 else {
                                   if(x$binaryCD == 0){
                                     "Benign"
                                   } # End if binaryCD == 0
                                   else{
                                     NA
                                   } # End else (binaryCD is not 0)
                                 } # End else (binaryCD is not missing)
                               } # End else (no histological analysis)
                             } # End if no missings in OriginCD
                           else{
                             NA
                           }# End else (missings in OriginCD)
                           return(x)
                         })
table(iota5ValAnalyses$CD_5groups)
iota5ValAnalyses$CD_5groups = factor(iota5ValAnalyses$CD_5groups, levels = c("Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"))

## Binary variable for multilocular
iota5ValAnalyses$multibin = sapply(iota5ValAnalyses$`Tumour type`,
                                   function(x){
                                     if(is.na(x))
                                       NA
                                     else if(grepl("multilocular", x))
                                       1
                                     else
                                       0
                                   })

## Binary variable for solid
iota5ValAnalyses$solidbin = sapply(iota5ValAnalyses$`Tumour type`,
                                   function(x){
                                     if(is.na(x))
                                       NA
                                     else if(grepl("solid", x))
                                       1
                                     else
                                       0
                                   })

## Binary variable for bilateral
iota5ValAnalyses$bilatbin = sapply(iota5ValAnalyses$Bilaterality,
                                   function(x){
                                     if(is.na(x))
                                       NA
                                     else if (x == "bilateral")
                                       1
                                     else
                                       0
                                   })

## Binary variable for Metastases
iota5ValAnalyses$metasbin = sapply(iota5ValAnalyses$Metastases,
                                   function(x){
                                     if(is.na(x))
                                       NA
                                     else if (x == "yes")
                                       1
                                     else
                                       0
                                   })

## Ordinal variable for Subjective assessment
iota5ValAnalyses$CertaintyOrdinal <- factor(iota5ValAnalyses$Certainty,
                                            levels = c("certainly benign", "probably benign", "uncertain", "probably malignant", "certainly malignant"),
                                            ordered = T)


## Variable indicating whether CA125 is missing
iota5ValAnalyses$CA125missing = factor(sapply(iota5ValAnalyses$CA125, function(x) as.numeric(is.na(x))),
                                       levels = 0:1, labels = c("CA125 present", "CA125 missing"))
iota5ValAnalyses[iota5ValAnalyses$CA125missing == "CA125 missing", c("CA125", "CA125missing")]

table(iota5ValAnalyses$CA125missing, iota5ValAnalyses$immed.oper)

## Variable indicating whether CD_5groups is missing
iota5ValAnalyses$CDmissing = factor(sapply(iota5ValAnalyses$CD_5groups, function(x) as.numeric(is.na(x))),
                                    levels = 0:1, labels = c("CD present", "CD missing"))
iota5ValAnalyses[iota5ValAnalyses$CDmissing == "CD missing", c("CD_5groups", "CDmissing")]

## Variable for the centers in analyses
iota5ValAnalyses$centerRE = sapply(iota5ValAnalyses$center,
                                   function(x){
                                     if(x == "AGR")
                                       "Athens (Greece)"
                                     else if(x == "CIT")
                                       "Milan 1 (Italy)"
                                     else if(x == "GBE")
                                       "Genk (Belgium)"
                                     else if(x == "LBE")
                                       "Leuven (Belgium)"
                                     else if(x == "MPO")
                                       "Katowice (Poland)"
                                     else if(x == "MSW")
                                       "Malmö (Sweden)"
                                     else if(x == "NCI")
                                       "Milan 2 (Italy)"
                                     else if(x == "OIT")
                                       "Monza (Italy)"
                                     else if(x == "PSP")
                                       "Pamplona (Spain)"
                                     else if(x == "RIT")
                                       "Rome (Italy)"
                                     else if(x == "SIT")
                                       "Cagliari (Italy)"
                                     else if(x == "SSW")
                                       "Stockholm (Sweden)"
                                     else if(x == "TIT")
                                       "Trieste (Italy)"
                                     else
                                       "Other"
                                   })
table(iota5ValAnalyses$center)
table(iota5ValAnalyses$centerRE)
table(iota5ValAnalyses$center, iota5ValAnalyses$centerRE)

## Variable for locules
iota5ValAnalyses$Nlocules <- ifelse(iota5ValAnalyses$Locules == "> 10", "> 10",
                                    ifelse(as.numeric(iota5ValAnalyses$Locules) == 1, "1",
                                    ifelse(as.numeric(iota5ValAnalyses$Locules) >= 2 & as.numeric(iota5ValAnalyses$Locules) <= 10, "2 - 10", "Other")))

## Binary variable for Personal history of ovarian cancer
iota5ValAnalyses$Personal_OvCAbin <- ifelse(iota5ValAnalyses$Personal_OvCA == "no", "No", "Yes")

## Log transformation of 'maximum diameter of lesion'
iota5ValAnalyses$Llesdmax <- log2(iota5ValAnalyses$`Lesion largest diameter`)

## Quadratic term for 'Proportion of solid tissue'
iota5ValAnalyses$Qpropsol <- iota5ValAnalyses$propsol^2


#### 4.1.2 Variables needed ####
iota5ValAnalyses$outcome1Correct = as.factor(iota5ValAnalyses$outcome1)
iota5ValAnalyses$binaryCDcorrect = as.factor(iota5ValAnalyses$binaryCD)

# Variables needed in imputation model
VarsImp = c('PresumedEndometrioma', 'CertaintyOrdinal', 'Patient age', 'oncocenter', 'Llesdmax', 'propsol', 'Qpropsol', 'Nlocules', 'papnr',
            'shadowsbin', 'ascitesbin', 'metasbin', 'bilatbin', 'Pelvic_pain', 'Personal_OvCAbin', 'Papillary height', 'papflow', 
            'ColorOrdinal', 'Echogenicity', 'CD_5groups', 'irregbin')
# Additional variables for the dataframe
VarsAdd = c('Lesion largest diameter', 'soldmax', 'multibin', 'solidbin', 'Tumour type', 'Patient ID', 'Mass ID', 'center', 'Study Outcome', 'Mass Outcome', 
            'CA125missing', 'outcome1', 'Postmenopausal2', 'binaryCDcorrect', 'centerRE', 'outcome1Correct', 'binaryCD', 'CDmissing', 'ClinicalDiagnosis', 
            'loc10', 'immed.oper', 'TimeSurgery', 'multiple.visits', 'cumulativeFollowUp', 'initial.policy', 'CD_short', 'FollowUp')

# Check for missingness in the variables (enkel variabelen van VarsImp)
iota5ValAnalyses[is.na(iota5ValAnalyses$PresumedEndometrioma), c("Patient ID", "PresumedEndometrioma")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$Certainty), c("Patient ID", "CertaintyOrdinal")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$`Patient age`), c("Patient ID", "Patient age")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$oncocenter), c("Patient ID", "oncocenter")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$`Lesion largest diameter`), c("Patient ID", "Llesdmax")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$propsol), c("Patient ID", "propsol")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$propsol), c("Patient ID", "Qpropsol")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$Nlocules), c("Patient ID", "Nlocules")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$loc10), c("Patient ID", "loc10")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$papnr), c("Patient ID", "papnr")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$shadowsbin), c("Patient ID", "shadowsbin")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$ascitesbin), c("Patient ID", "ascitesbin")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$metasbin), c("Patient ID", "metasbin")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$bilatbin), c("Patient ID", "bilatbin")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$Pelvic_pain), c("Patient ID", "Pelvic_pain")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$Personal_OvCA), c("Patient ID", "Personal_OvCAbin")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$`Papillary height`), c("Patient ID", "Papillary height")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$papflow), c("Patient ID", "papflow")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$ColorOrdinal), c("Patient ID", "ColorOrdinal")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$Echogenicity), c("Patient ID", "Echogenicity")] # No missings
iota5ValAnalyses[is.na(iota5ValAnalyses$CD_5groups), c("Patient ID", "CD_5groups", "ClinicalDiagnosis", "OriginCD")] # 490 Missings (454 uncertains and 5 malignants)
iota5ValAnalyses[is.na(iota5ValAnalyses$irregbin), c("Patient ID", "irregbin")] # No missings


#----------------------------------#
#### 4.2 Preparation imputation ####
#----------------------------------#

iota5ValAnalyses$ll_CA125 <- log(log(iota5ValAnalyses$CA125 + 1))

AllVars <- c('ll_CA125', VarsImp, VarsAdd)

PredMatr = matrix(0, length(AllVars), length(AllVars))
colnames(PredMatr) <- AllVars
rownames(PredMatr) <- AllVars

# PredMatr: A value of 1 indicates that the column variable is used as a predictor to impute the target (row) variable, and a 0 means that it is not used

# Prediction Matrix for CA125: the first row needs to be 1, except for ll_CA125. The rest of the matrix is zero
PredMatr[1, -1] <- 1
# Variables present in the dataset that are not necessary for the imputation set to zero
PredMatr[, c('Lesion largest diameter', 'soldmax', 'multibin', 'solidbin', 'Tumour type', 'Patient ID', 'Mass ID', 'center', 'Study Outcome', 'Mass Outcome', 
             'CA125missing', 'outcome1', 'Postmenopausal2', 'binaryCDcorrect', 'centerRE', 'outcome1Correct', 'binaryCD', 'CDmissing', 'ClinicalDiagnosis', 
             'loc10', 'immed.oper', 'TimeSurgery', 'multiple.visits', 'cumulativeFollowUp', 'initial.policy', 'CD_short', 'FollowUp')] <- 0
PredMatr

# Prediction Matrix for Clinical Diagnosis
PredMatr['CD_5groups', ] <- 1
PredMatr['CD_5groups', c('Lesion largest diameter', 'soldmax', 'multibin', 'solidbin', 'Tumour type', 'Patient ID', 'Mass ID', 'center', 'Study Outcome', 'Mass Outcome', 
                         'CA125missing', 'outcome1', 'Postmenopausal2', 'binaryCDcorrect', 'centerRE', 'outcome1Correct', 'binaryCD', 'CDmissing', 'ClinicalDiagnosis', 
                         'loc10', 'immed.oper', 'TimeSurgery', 'multiple.visits', 'cumulativeFollowUp', 'initial.policy', 'CD_short', 'FollowUp', 'CD_5groups')] <- 0
PredMatr


#----------------------#
#### 4.3 Imputation ####
#----------------------#

iota5Imp <- iota5ValAnalyses[, AllVars]

colnames(iota5Imp) = gsub(" ","\\.",colnames(iota5Imp))
rownames(PredMatr)   = colnames(PredMatr) = gsub(" ","\\.",colnames(PredMatr))

Init = mice(iota5Imp, m = 1, maxit = 0, predictorMatrix = PredMatr)
Init$loggedEvents 
Init$method
for(i in Init$loggedEvents$out)
  iota5Imp[, i] %<>% as.factor

Method = sapply(colnames(iota5Imp),
                function(x){
                  if(x == "ll_CA125")
                    "pmm"
                  else if(x == "CD_5groups")
                    "polyreg"
                  else
                    ""
                })
Method

iota5MI <- mice(iota5Imp, m = 100, me = Method,
                predictorMatrix = PredMatr, maxit = 50,
                seed = 1213, vis = "monotone", ridge = 1e-3)

iota5MI$loggedEvents
iota5MI$imp$ll_CA125
iota5MI$imp$CD_5groups
View(iota5MI$chainMean)

#-----------------------------#
#### 4.4 Check convergence ####
#-----------------------------#

## Convergence
tiff("Results Paper 2/With cysts/Convergence.tiff")
plot(iota5MI)
dev.off()

## Density plots
densityplot(iota5MI, ~ ll_CA125| .imp)

tiff("Results Paper 2/With cysts/Density CA125.tiff")
densityplot(iota5MI, ~ ll_CA125)
dev.off()

densityplot(iota5MI, ~ CD_5groups| .imp)

tiff("Results Paper 2/With cysts/Density CD.tiff")
densityplot(iota5MI, ~ CD_5groups)
dev.off()

## Boxplot CA125
Iota5Imputed                = mice::complete(iota5MI, "long")
Iota5Ni                     = complete(iota5MI)
Iota5Imputed$CA125 = exp(exp(Iota5Imputed$ll_CA125)) - 1
Iota5Ni$CA125 = exp(exp(Iota5Ni$ll_CA125)) - 1

tiff("Results Paper 2/With cysts/Boxplot CA125 imputed.tiff", width = 1200, height = 650)
boxplot(ll_CA125 ~ .imp, Iota5Imputed[Iota5Imputed$CA125missing == "CA125 missing",])
dev.off()


## Determine most commonly imputed value as outcome for patients in category M3
ImputeCD = dlply(Iota5Imputed[Iota5Imputed$CDmissing == "CD missing" & Iota5Imputed$CD_short == "M3",], .(Patient.ID),
                 function(x){
                   Outcome = table(x$CD_5groups)
                   names(Outcome)[which.max(Outcome)]
                 })
ImputeCD

ImputeCD$MSW666 # Stage II-IV invasive
ImputeCD$OIT112 # Borderline
ImputeCD$OIT215 # Borderline
ImputeCD$OIT299 # Borderline

Iota5Imputed[Iota5Imputed$Patient.ID == "MSW666", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "OIT112", c("CD_5groups")] <- "Borderline"
Iota5Imputed[Iota5Imputed$Patient.ID == "OIT215", c("CD_5groups")] <- "Borderline"
Iota5Imputed[Iota5Imputed$Patient.ID == "OIT299", c("CD_5groups")] <- "Borderline"

## adapt binaryCDcorrect according to CD_5groups
Iota5Imputed$binaryCDcorrect <- sapply(Iota5Imputed$CD_5groups,
                                       function(x){
                                         if(x == "Benign")
                                           0
                                         else
                                           1
                                       })
table(Iota5Imputed$binaryCDcorrect)

## Observed and imputed outcome
Observed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD present"])
round(Observed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD present",]) * 100, 2)
Imputed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD missing"])
round(Imputed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD missing",]) * 100, 2)

## Save dataset
Iota5ValI  = Iota5Imputed
Iota5ValNI = Iota5Imputed[Iota5Imputed$.imp == 1, ]

## Missingnes in CA125
table(Iota5ValNI$CA125missing) # Overall
table(Iota5ValNI$CA125missing, Iota5ValNI$center) # Per center
table(Iota5ValNI$CA125missing, Iota5ValNI$initial.policy) # Per suggested management
table(IOTA5Surg.NI$CA125missing) # Surgery
table(IOTA5FU.NI$CA125missing) # FU
table(Iota5ValNI$CA125missing, Iota5ValNI$Postmenopausal2) # Per menopausal status
table(Iota5ValNI$CA125missing, Iota5ValNI$oncocenter)# Per type of centre
table(Iota5ValNI$CA125missing, Iota5ValNI$center, Iota5ValNI$initial.policy)

IOTA5.Exclu <- iota5.firstobs(IOTA5Descrip[!(IOTA5Descrip$center %in% SelectedCenters),]) # Dataset with only the excluded centra
IOTA5.Exclu$CA125missing = factor(sapply(IOTA5.Exclu$CA125, function(x) as.numeric(is.na(x))),
                                  levels = 0:1, labels = c("CA125 present", "CA125 missing"))
table(IOTA5.Exclu$CA125missing, IOTA5.Exclu$center) # Per center
table(IOTA5.Exclu$CA125missing, IOTA5.Exclu$center, IOTA5.Exclu$initial.policy)


table(Iota5ValNI$CA125missing, Iota5ValNI$immed.oper)
## percentage of missing for FU: 78.69
## percentage of missing for immediately operated: 31.69
table(Iota5ValNI$CA125missing, Iota5ValNI$multiple.visits)
## percentage of missing for FU: 79.79
## percentage of missing for immediately operated: 35.80
table(Iota5ValNI$CA125missing, Iota5ValNI$initial.policy)
table(Iota5ValNI$CA125missing, Iota5ValNI$CertaintyOrdinal)


#-------------------------------------#
#### 4.5 Malignancy vs tumour type ####
#-------------------------------------#

TableMal <- ddply(Iota5ValI, .(Tumour.type, .imp),
                  function(x){
                    #x = x[!is.na(x$SAbinary)]
                    data.frame(N = nrow(x),
                               n_with_outcome = nrow(x) - sum(is.na(x$binaryCDcorrect)),
                               n_malignant = count(x$binaryCDcorrect[x$binaryCDcorrect == 1])[1, 2],
                               malignancy_rate = round(count(x$binaryCDcorrect[x$binaryCDcorrect == 1])[1, 2] / nrow(x) * 100, 4),
                               check.names = F)
                  })
TableMal
library(doBy)
summaryBy(cbind(N, n_with_outcome, n_malignant, malignancy_rate) ~ Tumour.type, data = TableMal)


#######################
#### 5. Validation ####
#######################

#---------------------------#
#### 5.1 Imputed dataset ####
#---------------------------#

#### 5.1.1 Model predictions ####
source(paste0(Scripts, "/FunctionsAllPredModels2.R"))


#### RMI ####
Iota5ValI = RMI(postmeno = Postmenopausal2, multilocular = multibin, solid = solidbin, metastases = metasbin, bilateral = bilatbin, ascites = ascitesbin, CA125 = CA125, patientid = Patient.ID, data = Iota5ValI, imputed = T, .imp = .imp)
Iota5ValI$RMI

## Histogram: back-to-back
library(scales)
RMI.hist <- ggplot(Iota5ValI, aes(x=RMI)) +
            geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '0'),aes(y = ..count.., fill = "Benign"), bins = 100, alpha=.75) +
            geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '1'),aes(y = -..count.., fill = "Malignant"), bins = 100, alpha=.75) +
            theme_minimal() +
            labs(x = "Risk of Malignancy Index (RMI)",
                 y = "Frequency",
                 fill = "Outcome") +
            scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
            scale_x_continuous(trans = scales::pseudo_log_trans(base = 10), breaks = c(0, 10, 100, 1000, 10000, 100000, 1000000), labels = c(0, '10', '100', '1000', '10000', '100000', '1000000')) + 
            scale_y_continuous(limits = c(-50000, 200000), breaks = c(-50000, 0, 50000, 100000, 150000, 200000), labels = c(500, 0, 500, 1000, 1500, 2000))
RMI.hist


#### LR2 ####
Iota5ValNI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5ValNI)
Iota5ValNI$LR2
Iota5ValI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5ValI)
Iota5ValI$LR2

## Histogram: back-to-back
LR2.hist <- ggplot(Iota5ValI, aes(x=LR2)) +
            geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
            geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
            theme_minimal() +
            labs(x = "Estimated risk",
                 y = "Frequency",
                 fill = "Outcome") +
            scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
            scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
            scale_y_continuous(limits = c(-25000, 90000), breaks = c(-25000, 0, 25000, 50000, 75000, 100000), labels = c(250, 0, 250, 500, 750, 1000))
LR2.hist


#### Simple Rules and Simple Rules Risk ####
Iota5ValNI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5ValNI)
Iota5ValNI$SRrisks
Iota5ValI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5ValI)
Iota5ValI$SRrisks

## Histogram: back-to-back
SRrisk.hist <- ggplot(Iota5ValI, aes(x=SRrisks)) +
               geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
               geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
               theme_minimal() +
               labs(x = "Estimated risk",
                    y = "Frequency",
                    fill = "Outcome") +
               scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
               scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
               scale_y_continuous(limits = c(-50000, 160000), breaks = c(-50000, 0, 50000, 100000 , 150000, 200000), labels = c(500, 0, 500, 1000, 1500, 2000))
SRrisk.hist


#### ADNEX ####
## With CA125
Iota5ValI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = T, woca = F, data = Iota5ValI)
Iota5ValI$pmalw

## Histogram: back-to-back
ADNEXw.hist <- ggplot(Iota5ValI, aes(x=pmalw)) +
               geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
               geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
               theme_minimal() +
               labs(x = "Estimated risk",
                    y = "Frequency",
                    fill = "Outcome") +
               scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
               scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
               scale_y_continuous(limits = c(-25000, 90000), breaks = c(-25000, 0, 25000, 50000, 75000, 100000), labels = c(250, 0, 250, 500, 750, 1000))
ADNEXw.hist


## Without CA125
Iota5ValNI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5ValNI)
Iota5ValNI$pmalwo
Iota5ValI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5ValI)
Iota5ValI$pmalwo

## Histogram: back-to-back
ADNEXwo.hist <- ggplot(Iota5ValI, aes(x=pmalwo)) +
                geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '0'),aes(y = ..count.., fill = "Benign"), binwidth = 0.01, alpha=.75) +
                geom_histogram(data=subset(Iota5ValI,binaryCDcorrect == '1'),aes(y = -..count.., fill = "Malignant"), binwidth = 0.01, alpha=.75) +
                theme_minimal() +
                labs(x = "Estimated risk",
                     y = "Frequency",
                     fill = "Outcome") +
                scale_fill_manual(name = "Outcome", values = c(Benign = "#00BFC4", Malignant = "#F8766D")) +
                scale_x_continuous(breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), labels = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)) +
                scale_y_continuous(limits = c(-25000, 85000), breaks = c(-25000, 0, 25000, 50000, 75000, 100000), labels = c(250, 0, 250, 500, 750, 1000))
ADNEXwo.hist


#### 5.1.2 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5ValI, method.MA = "SJ")

RMI.AUC$dataPlot <- RMI.AUC$dataPlot[c(1:9, 11:16, 10, 17:20),]
RMI.AUC$Plot     <- RMI.AUC$Plot[c(1:9, 11:16, 10, 17:20),]

RMI.AUC$Plot[17, "RRauc"] <- "   "
RMI.AUC$Plot[18, "RRauc"] <- "   "
RMI.AUC$Plot[18, "RRcenter"] <- "Meta-analysis"
RMI.AUC$Plot[19, "RRcenter"] <- "AUC (95% CI)"
RMI.AUC$Plot[20, "RRauc"] <- "        (0.74 to 0.96)"
RMI.AUC$Plot[20, "RRprev"] <- "   "

fpDrawBarCI <- function (lower_limit, estimate, upper_limit, size, col, y.offset = 0.5, ...) 
{
  size <- ifelse(is.unit(size), convertUnit(size, unitTo = "npc", valueOnly = TRUE), size) * 0.9
  grid.polygon(x = unit(c(lower_limit, upper_limit, upper_limit, 
                          lower_limit), "native"), 
               y = c(0.4, 0.4, 0.6, 0.6), 
               gp = gpar(fill = col, col = col))
}


tiff("Results Paper 2/With cysts/Supplementary/Forest RMI - PI.tiff", width = 25, height = 17.75, units = "cm", res = 300)
forestplot(RMI.AUC$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC$dataPlot$AUC,
           lower = RMI.AUC$dataPlot$LL,
           upper = RMI.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC$IncludedCenters)),FALSE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(RMI.AUC$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = RMI.AUC$Performance$AUC)
dev.off()


## Sensitivity and specificity
RMI.SS25  <- SS.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, threshold = 25, imp = .imp, data = Iota5ValI, method.MA = "SJ")
RMI.SS100 <- SS.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, threshold = 100, imp = .imp, data = Iota5ValI, method.MA = "SJ")
RMI.SS200 <- SS.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, threshold = 200, imp = .imp, data = Iota5ValI, method.MA = "SJ")
RMI.SS250 <- SS.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, threshold = 250, imp = .imp, data = Iota5ValI, method.MA = "SJ")

RMI.SS <- rbind(RMI.SS25$OverallPer, RMI.SS100$OverallPer, RMI.SS200$OverallPer, RMI.SS250$OverallPer)
RMI.SS

# Fixed sensitivity
FixedSens.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, Sensitivity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer

# Fixed specificity
FixedSpec.imp(pred = RMI, outcome = binaryCDcorrect, center = centerRE, Specificity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer


#### LR2 ####

## AUC
LR2.AUC <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5ValI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

LR2.AUC$dataPlot <- LR2.AUC$dataPlot[c(1:9, 11:16, 10, 17:20),]
LR2.AUC$Plot     <- LR2.AUC$Plot[c(1:9, 11:16, 10, 17:20),]

LR2.AUC$Plot[17, "RRauc"] <- "   "
LR2.AUC$Plot[18, "RRauc"] <- "   "
LR2.AUC$Plot[18, "RRcenter"] <- "Meta-analysis"
LR2.AUC$Plot[19, "RRcenter"] <- "AUC (95% CI)"
LR2.AUC$Plot[20, "RRauc"] <- "        (0.82 to 0.96)"
LR2.AUC$Plot[20, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest LR2 - PI.tiff", width = 25, height = 17.75, units = "cm", res = 300)
forestplot(LR2.AUC$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC$dataPlot$AUC,
           lower = LR2.AUC$dataPlot$LL,
           upper = LR2.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC$IncludedCenters)), FALSE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(LR2.AUC$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = LR2.AUC$Performance$AUC)
dev.off()


## Sensitivity and specificity
LR2.SS01  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS03  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.03, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS05  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.05, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS10  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.10, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS15  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.15, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS20  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.20, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS25  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.25, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS30  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.30, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS40  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.40, imp = .imp, data = Iota5ValI, method.MA = "SJ")
LR2.SS50  <- SS.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, threshold = 0.50, imp = .imp, data = Iota5ValI, method.MA = "SJ")

LR2.SS <- rbind(LR2.SS01$OverallPer, LR2.SS03$OverallPer, LR2.SS05$OverallPer, LR2.SS10$OverallPer, LR2.SS15$OverallPer, LR2.SS20$OverallPer, LR2.SS25$OverallPer, LR2.SS30$OverallPer, LR2.SS40$OverallPer, LR2.SS50$OverallPer)
LR2.SS

# Fixed sensitivity
FixedSens.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, Sensitivity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer

# Fixed specificity
FixedSpec.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, Specificity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer


#### Simple Rules Risk ####

## AUC
SRrisks.AUC <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5ValI, method.MA = "SJ", titleGraph = "Simple Rules risk: AUC per center")

SRrisks.AUC$dataPlot <- SRrisks.AUC$dataPlot[c(1:9, 11:16, 10, 17:20),]
SRrisks.AUC$Plot     <- SRrisks.AUC$Plot[c(1:9, 11:16, 10, 17:20),]

SRrisks.AUC$Plot[17, "RRauc"] <- "   "
SRrisks.AUC$Plot[18, "RRauc"] <- "   "
SRrisks.AUC$Plot[18, "RRcenter"] <- "Meta-analysis"
SRrisks.AUC$Plot[19, "RRcenter"] <- "AUC (95% CI)"
SRrisks.AUC$Plot[20, "RRauc"] <- "        (0.83 to 0.98)"
SRrisks.AUC$Plot[20, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest SRrisks - PI.tiff", width = 25, height = 17.75, units = "cm", res = 300)
forestplot(SRrisks.AUC$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC$dataPlot$AUC,
           lower = SRrisks.AUC$dataPlot$LL,
           upper = SRrisks.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC$IncludedCenters)), FALSE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(SRrisks.AUC$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = SRrisks.AUC$Performance$AUC)
dev.off()


## Sensitivity and specificity
SRrisks.SS01  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS03  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.03, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS05  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.05, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS10  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.10, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS15  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.15, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS20  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.20, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS25  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.25, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS30  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.30, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS40  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.40, imp = .imp, data = Iota5ValI, method.MA = "SJ")
SRrisks.SS50  <- SS.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, threshold = 0.50, imp = .imp, data = Iota5ValI, method.MA = "SJ")

SRrisks.SS <- rbind(SRrisks.SS01$OverallPer, SRrisks.SS03$OverallPer, SRrisks.SS05$OverallPer, SRrisks.SS10$OverallPer, SRrisks.SS15$OverallPer, SRrisks.SS20$OverallPer, SRrisks.SS25$OverallPer, SRrisks.SS30$OverallPer, SRrisks.SS40$OverallPer, SRrisks.SS50$OverallPer)
SRrisks.SS

# Fixed sensitivity
FixedSens.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, Sensitivity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer

# Fixed specificity
FixedSpec.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, Specificity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5ValI, method.MA = "SJ")

ADNEX.AUC$dataPlot <- ADNEX.AUC$dataPlot[c(1:9, 11:16, 10, 17:20),]
ADNEX.AUC$Plot     <- ADNEX.AUC$Plot[c(1:9, 11:16, 10, 17:20),]

ADNEX.AUC$Plot[17, "RRauc"] <- "   "
ADNEX.AUC$Plot[18, "RRauc"] <- "   "
ADNEX.AUC$Plot[18, "RRcenter"] <- "Meta-analysis"
ADNEX.AUC$Plot[19, "RRcenter"] <- "AUC (95% CI)"
ADNEX.AUC$Plot[20, "RRauc"] <- "        (0.83 to 0.98)"
ADNEX.AUC$Plot[20, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXw - PI.tiff", width = 25, height = 17.75, units = "cm", res = 300)
forestplot(ADNEX.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC$dataPlot$AUC,
           lower = ADNEX.AUC$dataPlot$LL,
           upper = ADNEX.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC$IncludedCenters)), FALSE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEX.AUC$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEX.AUC$Performance$AUC)
dev.off()


## Sensitivity and specificity
ADNEX.SS01  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS03  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.03, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS05  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.05, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS10  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.10, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS15  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.15, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS20  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.20, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS25  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.25, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS30  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.30, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS40  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.40, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEX.SS50  <- SS.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, threshold = 0.50, imp = .imp, data = Iota5ValI, method.MA = "SJ")

ADNEX.SS <- rbind(ADNEX.SS01$OverallPer, ADNEX.SS03$OverallPer, ADNEX.SS05$OverallPer, ADNEX.SS10$OverallPer, ADNEX.SS15$OverallPer, ADNEX.SS20$OverallPer, ADNEX.SS25$OverallPer, ADNEX.SS30$OverallPer, ADNEX.SS40$OverallPer, ADNEX.SS50$OverallPer)
ADNEX.SS

# Fixed sensitivity
FixedSens.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, Sensitivity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer

# Fixed specificity
FixedSpec.imp(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, Specificity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5ValI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

ADNEXwo.AUC$dataPlot <- ADNEXwo.AUC$dataPlot[c(1:9, 11:16, 10, 17:20),]
ADNEXwo.AUC$Plot     <- ADNEXwo.AUC$Plot[c(1:9, 11:16, 10, 17:20),]

ADNEXwo.AUC$Plot[17, "RRauc"] <- "   "
ADNEXwo.AUC$Plot[18, "RRauc"] <- "   "
ADNEXwo.AUC$Plot[18, "RRcenter"] <- "Meta-analysis"
ADNEXwo.AUC$Plot[19, "RRcenter"] <- "AUC (95% CI)"
ADNEXwo.AUC$Plot[20, "RRauc"] <- "        (0.82 to 0.98)"
ADNEXwo.AUC$Plot[20, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXwo - PI.tiff", width = 25, height = 17.75, units = "cm", res = 300)
forestplot(ADNEXwo.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC$dataPlot$AUC,
           lower = ADNEXwo.AUC$dataPlot$LL,
           upper = ADNEXwo.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC$IncludedCenters)), FALSE, TRUE, TRUE, TRUE),
           title = "",
           xlab = "AUC",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwo.AUC$IncludedCenters) + 4)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwo.AUC$Performance$AUC)
dev.off()


## Sensitivity and specificity
ADNEXwo.SS01  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS03  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.03, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS05  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.05, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS10  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.10, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS15  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.15, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS20  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.20, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS25  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.25, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS30  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.30, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS40  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.40, imp = .imp, data = Iota5ValI, method.MA = "SJ")
ADNEXwo.SS50  <- SS.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, threshold = 0.50, imp = .imp, data = Iota5ValI, method.MA = "SJ")

ADNEXwo.SS <- rbind(ADNEXwo.SS01$OverallPer, ADNEXwo.SS03$OverallPer, ADNEXwo.SS05$OverallPer, ADNEXwo.SS10$OverallPer, ADNEXwo.SS15$OverallPer, ADNEXwo.SS20$OverallPer, ADNEXwo.SS25$OverallPer, ADNEXwo.SS30$OverallPer, ADNEXwo.SS40$OverallPer, ADNEXwo.SS50$OverallPer)
ADNEXwo.SS

# Fixed sensitivity
FixedSens.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, Sensitivity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer

# Fixed specificity
FixedSpec.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, Specificity = 0.9, imp = .imp, data = Iota5ValI, method.MA = "SJ")$OverallPer


#### Summary plot for AUC ####

NA.forest <- RMI.AUC$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, RMI.AUC$Performance[1,], LR2.AUC$Performance[1,], SRrisks.AUC$Performance[1,], ADNEXwo.AUC$Performance[1,], ADNEX.AUC$Performance[1,])
Summary.AUC$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.PI <- rbind(NA.forest, RMI.AUC$Performance[2,], LR2.AUC$Performance[2,], SRrisks.AUC$Performance[2,], ADNEXwo.AUC$Performance[2,], ADNEX.AUC$Performance[2,])
Summary.PI$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.AUC
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(RMI.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(LR2.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(SRrisks.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEX.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC$Performance$UL[2], 2), nsmall = 2), ")"))
)


setEPS()
postscript("Results Paper 2/With cysts/Paper/Summary AUC - PI.eps", width = 12, height = 5.5, horizontal = FALSE, onefile = FALSE, paper = "special")
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()


#### Subjective assessment ####

Iota5ValI$SA <- ifelse(grepl("benign", Iota5ValI$Certainty, fixed = T), 0, 1)
Iota5ValNI$SA <- ifelse(grepl("benign", Iota5ValNI$Certainty, fixed = T), 0, 1)

## Sensitivity and specificity
SS.imp(pred = SA, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")$`OverallPer`

table(Iota5ValNI$SA, ifelse(Iota5ValNI$SA < 0.01, 0, 1))

#### Simple Rules ####

Iota5ValI$SRbin <- ifelse(Iota5ValI$SRconcl == "Benign", 0, 1)
Iota5ValNI$SRbin <- ifelse(Iota5ValNI$SRconcl == "Benign", 0, 1)

perIncon <- round(nrow(Iota5ValNI[Iota5ValNI$SRconcl == "Inconclusive",]) / nrow(Iota5ValNI) * 100, 2)
perIncon # 16.23 %

table(Iota5ValNI$SRbin, Iota5ValNI$SRconcl)
table(Iota5ValI$SRbin[Iota5ValI$.imp == 55], Iota5ValI$SRconcl[Iota5ValI$.imp == 55])

## Sensitivity and specificity
SS.imp(pred = SRbin, outcome = binaryCDcorrect, center = centerRE, threshold = 0.01, imp = .imp, data = Iota5ValI, method.MA = "SJ")$`OverallPer`


#### 5.1.3 Calibration of the risk of malignancy ####

SampleSize <- c(567, 166, 406, 139, 501, 794, 367, 98, 267, 334, 111, 681, 363, 111)

ColorCen <- c("red", "limegreen", "orange", "blue", "maroon", "red", "limegreen", "orange", "blue", "maroon", "red", "limegreen", "orange", "blue")

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI.tiff", width = 14, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = Iota5ValI, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()

## Curve per center
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI centers.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = Iota5ValI, center = centerRE, imputation.id = .imp, patientid = Patient.ID, CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2)
dev.off()


#### LR2 ####
LR2 <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID)

## Curve per center
tiff("Results Paper 2/With cysts/Supplementary/Calibration LR2 centers.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID, CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title = "")
dev.off()


#### Simple Rules Risk ####
SRrisk <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")

## Curve per center
tiff("Results Paper 2/With cysts/Supplementary/Calibration SRrisks centers.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID, CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title = "")
dev.off()


#### ADNEX with CA125 ####
ADNEXw <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = Iota5ValI, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")

## Curve per center
tiff("Results Paper 2/With cysts/Supplementary/Calibration ADNEXw centers.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = Iota5ValI, center = centerRE, imputation.id = .imp, patientid = Patient.ID, CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title = "")
dev.off()


#### ADNEX without CA125 ####
ADNEXwo <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")

## Curve per center
tiff("Results Paper 2/With cysts/Supplementary/Calibration ADNEXwo centers.tiff", width = 18.75, height = 14, units = "cm", res = 300)
RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = Iota5ValI, imputation.id = .imp, patientid = Patient.ID, CalibrLines = "centers", lty.centers = c(rep(1, 5), rep(4, 5), rep(3, 4)), col.centers = ColorCen, lwd.centers = round(sqrt(SampleSize)/6, 1), ncol.leg = 1, cex.leg = 1.2, title = "")
dev.off()


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2$Plot$y
p.LR2 = LR2$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk$Plot$y
p.SR = SRrisk$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw$Plot$y
p.ADNEXw = ADNEXw$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo$Plot$y
p.ADNEXwo = ADNEXwo$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo$SumSlo$UL, 2), nsmall = 2), ")"))

library(plotrix)

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
setEPS()
postscript("Results Paper 2/With cysts/Paper/Summary calibration - black white.eps", width = 5.5, height = 5.5, horizontal = FALSE, onefile = FALSE, paper = "special") #`width = 8.5, height = 8.5`
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 1, lty = 1, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.LR2, OverallCal.LR2, lty = 3, lwd = 2, col = "black")
lines(p.SR, OverallCal.SRrisk, lty = 2, lwd = 2, col = "black")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "gray50")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "black")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "black", "black", "gray50", "black"), lty = c(1,3,2,1,1), lwd = c(1,2,2,2,2), cex = 0.7, bty = "n", ncol = 1) #cex = 1.2
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7) #cex = 1.2
dev.off()


#### 5.1.4 Clinical utility ####

#### RMI ####

Iota5ValI$RMIdichom <- ifelse(Iota5ValI$RMI <= 200, 0, 1)

RMI.Imp <- list()

for(i in 1:100){
  RMI.Imp[[i]] <- list()
}

for(i in 1:100){
  RMI.Imp[[i]] <- DataWinBugs(pred = RMIdichom, outcome = binaryCDcorrect, data = Iota5ValI[Iota5ValI$.imp == i,], center = centerRE, sequence = 0.5)$`Results`[[1]]
}

library(data.table)
RMI.long <- rbindlist(RMI.Imp, fill = TRUE)

library(doBy)
RMI.sum <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB) ~ cbind(Cutoff, Center), data = RMI.long, FUN = mean)
RMI.sum


#### LR2 ####

LR2.Dec <- DataWinBugs.imp(pred = LR2, outcome = binaryCDcorrect, center = centerRE, data = Iota5ValI, imp = .imp)$`Results`
LR2.Dec
ggplot(data = LR2.Dec, aes(x = CutOff, y = NB.mean, colour = Center, group = Center)) + geom_line()


#### Simple Rules Risk ####

SRrisks.Dec <- DataWinBugs.imp(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, data = Iota5ValI, imp = .imp)$`Results`
SRrisks.Dec
ggplot(data = SRrisks.Dec, aes(x = CutOff, y = NB.mean, colour = Center, group = Center)) + geom_line()


#### Simple Rules ####

SR.Imp <- list()

for(i in 1:100){
  SR.Imp[[i]] <- list()
}

for(i in 1:100){
  SR.Imp[[i]] <- DataWinBugs(pred = SRbin, outcome = binaryCDcorrect, data = Iota5ValI[Iota5ValI$.imp == i,], center = centerRE, sequence = 0.5)$`Results`[[1]]
}

library(data.table)
SR.long <- rbindlist(SR.Imp, fill = TRUE)

library(doBy)
SR.long <- summaryBy(cbind(TP, TN, cases, controls, FP, FN, NB) ~ cbind(Cutoff, Center), data = SR.long, FUN = mean)
SR.long


#### ADNEX with CA125 ####

ADNEXw.Dec <- DataWinBugs.imp(pred = pmalw, outcome = binaryCDcorrect, data = Iota5ValI, center = centerRE, imp = .imp)$`Results` 
ggplot(data = ADNEXw.Dec, aes(x = CutOff, y = NB.mean, colour = Center, group = Center)) + geom_line()


#### ADNEX without CA125 ####

ADNEXwo.Dec <- DataWinBugs.imp(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, data = Iota5ValI, imp = .imp)$`Results`
ggplot(data = ADNEXwo.Dec, aes(x = CutOff, y = NB.mean, colour = Center, group = Center)) + geom_line()


#### Summary plot for Clinical utility ####
library(readxl)
Overview <- read_excel("Results Paper 2/With cysts/Results WinBugs - rounded.xlsx", sheet = "Overview")

setEPS()
postscript("Results Paper 2/With cysts/Paper/Summary Clinical utility - black white.eps", width = 5.5, height = 5.5, horizontal = FALSE, onefile = FALSE, paper = "special") #width = 8.5, height = 8.5
plot(Overview$`Risk threshold`, Overview$`Treat all`, type = "l", col = "black", xlab = "Risk threshold", ylab = "Net benefit", xlim = c(0.055, 0.5), ylim = c(-0.03, 0.22), lwd = 1, cex.lab = 1, cex.axis = 1, las = 1) #cex.lab = 1.5, cex.axis = 1.5
points(Overview$`Risk threshold`, Overview$`NB RMI`, type='l', col = "black", lty = 4, lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB LR2`, type = 'l', col = "black", lty = 3, lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB SRrisks`, type = 'l', col = "black", lty = 2, lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB SR`, type = 'l', col = "gray50", lty = 2, lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEX wo`, type = 'l', col = "gray50", lwd = 2)
points(Overview$`Risk threshold`, Overview$`NB ADNEXw`, type = 'l', col = "black", lwd = 2)
abline(a = 0, b = 0, lty = 1, lwd = 1, col = "gray50")
legend(x = 0.04, y = 0.22, legend = c( "RMI (Cut-off = 200)", "LR2", "SRRisk", "SR", "ADNEX without CA125", "ADNEX with CA125", "Treat all", "Treat none"), ncol = 3, #cex = 1.2, #bty = "n", #cex = 0.7, pt.cex = 1,
       col = c("black", "black", "black", "gray50", "gray50", "black", "black", 'gray50'), lty = c(4,3,2,2,1,1,1,1), lwd = c(2,2,2,2,2,2,1,1), cex = 0.7, bty = "n")
dev.off()



#### 5.1.5 Additional analysis for ADNEX ####

#### 5.1.5.1 Discrimination ####

#### ADNEX with CA125 ####

# Pair 1: benign vs borderline
Iota5ValI$pp1 <- Iota5ValI$pbenw / (Iota5ValI$pbenw + Iota5ValI$pborw)
AUC.imp(pred = pp1, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Borderline",])$Performance
# AUC = 0.8959

# Pair 2: benign vs stage I
Iota5ValI$pp2 <- Iota5ValI$pbenw / (Iota5ValI$pbenw + Iota5ValI$pst1w)
AUC.imp(pred = pp2, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Stage I invasive",])$Performance
# AUC = 0.9547

# Pair 3: benign vs stage II-IV
Iota5ValI$pp3 <- Iota5ValI$pbenw / (Iota5ValI$pbenw + Iota5ValI$pst2_4w)
AUC.imp(pred = pp3, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.9838

# Pair 4: benign vs metastatic
Iota5ValI$pp4 <- Iota5ValI$pbenw / (Iota5ValI$pbenw + Iota5ValI$pmetaw)
AUC.imp(pred = pp4, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.9496

# Pair 5: borderline vs stage I
Iota5ValI$pp5 <- Iota5ValI$pborw / (Iota5ValI$pborw + Iota5ValI$pst1w)
Iota5ValI$binbo <- ifelse(Iota5ValI$CD_5groups == "Borderline", 0, 1)
AUC.imp(pred = pp5, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Stage I invasive",])$Performance
# AUC = 0.7709

# Pair 6: borderline vs stage II-IV
Iota5ValI$pp6 <- Iota5ValI$pborw / (Iota5ValI$pborw + Iota5ValI$pst2_4w)
AUC.imp(pred = pp6, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.9221

# Pair 7: borderline vs metastatic
Iota5ValI$pp7 <- Iota5ValI$pborw / (Iota5ValI$pborw + Iota5ValI$pmetaw)
AUC.imp(pred = pp7, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.8813

# Pair 8: stage I vs stage II-IV
Iota5ValI$pp8 <- Iota5ValI$pst1w / (Iota5ValI$pst1w + Iota5ValI$pst2_4w)
Iota5ValI$bini <- ifelse(Iota5ValI$CD_5groups == "Stage I invasive", 0, 1)
AUC.imp(pred = pp8, outcome = bini, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage I invasive" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.8136

# Pair 9: stage I vs metastatic
Iota5ValI$pp9 <- Iota5ValI$pst1w / (Iota5ValI$pst1w + Iota5ValI$pmetaw)
AUC.imp(pred = pp9, outcome = bini, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage I invasive" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.7543

# Pair 10: stage II-IV vs metastatic
Iota5ValI$pp10 <- Iota5ValI$pst2_4w / (Iota5ValI$pst2_4w + Iota5ValI$pmetaw)
Iota5ValI$bin2 <- ifelse(Iota5ValI$CD_5groups == "Stage II-IV invasive", 0, 1)
AUC.imp(pred = pp10, outcome = bin2, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage II-IV invasive" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.7784


#### ADNEX without CA125 ####

# Pair 1: benign vs borderline
Iota5ValI$pp1wo <- Iota5ValI$pbenwo / (Iota5ValI$pbenwo + Iota5ValI$pborwo)
AUC.imp(pred = pp1wo, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Borderline",])$Performance
# AUC = 0.8938

# Pair 2: benign vs stage I
Iota5ValI$pp2wo <- Iota5ValI$pbenwo / (Iota5ValI$pbenwo + Iota5ValI$pst1wo)
AUC.imp(pred = pp2wo, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Stage I invasive",])$Performance
# AUC = 0.9519

# Pair 3: benign vs stage II-IV
Iota5ValI$pp3wo <- Iota5ValI$pbenwo / (Iota5ValI$pbenwo + Iota5ValI$pst2_4wo)
AUC.imp(pred = pp3wo, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.9714

# Pair 4: benign vs metastatic
Iota5ValI$pp4wo <- Iota5ValI$pbenwo / (Iota5ValI$pbenwo + Iota5ValI$pmetawo)
AUC.imp(pred = pp4wo, outcome = binaryCDcorrect, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Benign" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.9421

# Pair 5: borderline vs stage I
Iota5ValI$pp5wo <- Iota5ValI$pborwo / (Iota5ValI$pborwo + Iota5ValI$pst1wo)
Iota5ValI$binbo <- ifelse(Iota5ValI$CD_5groups == "Borderline", 0, 1)
AUC.imp(pred = pp5wo, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Stage I invasive",])$Performance
# AUC = 0.7764

# Pair 6: borderline vs stage II-IV
Iota5ValI$pp6wo <- Iota5ValI$pborwo / (Iota5ValI$pborwo + Iota5ValI$pst2_4wo)
AUC.imp(pred = pp6wo, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.9041

# Pair 7: borderline vs metastatic
Iota5ValI$pp7wo <- Iota5ValI$pborwo / (Iota5ValI$pborwo + Iota5ValI$pmetawo)
AUC.imp(pred = pp7wo, outcome = binbo, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Borderline" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.8791

# Pair 8: stage I vs stage II-IV
Iota5ValI$pp8wo <- Iota5ValI$pst1wo / (Iota5ValI$pst1wo + Iota5ValI$pst2_4wo)
Iota5ValI$bini <- ifelse(Iota5ValI$CD_5groups == "Stage I invasive", 0, 1)
AUC.imp(pred = pp8wo, outcome = bini, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage I invasive" | Iota5ValI$CD_5groups == "Stage II-IV invasive",])$Performance
# AUC = 0.7218

# Pair 9: stage I vs metastatic
Iota5ValI$pp9wo <- Iota5ValI$pst1wo / (Iota5ValI$pst1wo + Iota5ValI$pmetawo)
AUC.imp(pred = pp9wo, outcome = bini, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage I invasive" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.7522

# Pair 10: stage II-IV vs metastatic
Iota5ValI$pp10wo <- Iota5ValI$pst2_4wo / (Iota5ValI$pst2_4wo + Iota5ValI$pmetawo)
Iota5ValI$bin2 <- ifelse(Iota5ValI$CD_5groups == "Stage II-IV invasive", 0, 1)
AUC.imp(pred = pp10wo, outcome = bin2, imp = .imp, data = Iota5ValI[Iota5ValI$CD_5groups == "Stage II-IV invasive" | Iota5ValI$CD_5groups == "Metastatic",])$Performance
# AUC = 0.6643


#### 5.1.5.2 Calibration ####

library(VGAM)

#### ADNEX with CA125 ####

Iota5ValI$Groups <- ifelse(Iota5ValI$CD_5groups == "Benign", '1',
                           ifelse(Iota5ValI$CD_5groups == "Borderline", '2',
                           ifelse(Iota5ValI$CD_5groups == "Stage I invasive", '3',
                           ifelse(Iota5ValI$CD_5groups == "Stage II-IV invasive", '4', '5'))))

polADNEXw <- PolCal(outcome = Groups, k = 5, r = 1, p = c('pbenw', 'pborw', 'pst1w', 'pst2_4w', 'pmetaw'), LP = c('lp2w', 'lp3w', 'lp4w', 'lp5w'), data = Iota5ValI, imputation.id = .imp)    

tiff("Results Paper 2/With cysts/Paper/Calibration ADNEXw multinomial.tiff", width = 14, height = 14, units = "cm", res = 300) #width = 21.5, height = 21.5
plotlines5adapt(preds = polADNEXw$FIT, obs = polADNEXw$LP.predict, mtxt = "")
legend("topleft", legend = c("Ideal", "Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"),
       col = c("gray50", "red", "limegreen", "darkblue", "cyan", "magenta"), lty = c(2,1,1,1,1,1), lwd = 2, cex = 0.7, bty = "n")
dev.off()


#### ADNEX without CA125 ####

polADNEXwo <- PolCal(outcome = Groups, k = 5, r = 1, p = c('pbenwo', 'pborwo', 'pst1wo', 'pst2_4wo', 'pmetawo'), LP = c('lp2wo', 'lp3wo', 'lp4wo', 'lp5wo'), data = Iota5ValI, imputation.id = .imp)    

tiff("Results Paper 2/With cysts/Paper/Calibration ADNEXwo multinomial.tiff", width = 14, height = 14, units = "cm", res = 300) #width = 21.5, height = 21.5
plotlines5adapt(preds = polADNEXwo$FIT, obs = polADNEXwo$LP.predict, mtxt = "")
legend("topleft", legend = c("Ideal", "Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"),
       col = c("gray50", "red", "limegreen", "darkblue", "cyan", "magenta"), lty = c(2,1,1,1,1,1), lwd = 2, cex = 0.7, bty = "n")
dev.off()


#-----------------------------#
#### 5.2 Subgroup analyses ####
#-----------------------------#


#### 5.2.1 Patients who were operated within 120 days after the first scan and without follow-up scan ####

IOTA5Surg.Imp <- Iota5ValI[Iota5ValI$immed.oper == "1" & Iota5ValI$TimeSurgery <= 120,]
IOTA5Surg.NI  <- Iota5ValNI[Iota5ValNI$immed.oper == "1" & Iota5ValNI$TimeSurgery <= 120,]

table(IOTA5Surg.NI$centerRE, IOTA5Surg.NI$binaryCDcorrect)

CentersDescrip = ddply(IOTA5Surg.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.1.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.surg <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Surg.Imp, method.MA = "SJ")

forestplot(RMI.AUC.surg$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.surg$dataPlot$AUC,
           lower = RMI.AUC.surg$dataPlot$LL,
           upper = RMI.AUC.surg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.surg$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.surg$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.surg <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Surg.Imp, method.MA = "SJ")

forestplot(LR2.AUC.surg$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.surg$dataPlot$AUC,
           lower = LR2.AUC.surg$dataPlot$LL,
           upper = LR2.AUC.surg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.surg$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.surg$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.surg <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Surg.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.surg$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.surg$dataPlot$AUC,
           lower = SRrisks.AUC.surg$dataPlot$LL,
           upper = SRrisks.AUC.surg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.surg$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.surg$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.surg <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Surg.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.surg$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.surg$dataPlot$AUC,
           lower = ADNEX.AUC.surg$dataPlot$LL,
           upper = ADNEX.AUC.surg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.surg$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.surg$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.surg <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Surg.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.surg$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.surg$dataPlot$AUC,
           lower = ADNEXwo.AUC.surg$dataPlot$LL,
           upper = ADNEXwo.AUC.surg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.surg$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.surg$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.surg$Performance[1,]
NA.forest <- NA

Summary.AUC.surg <- rbind(NA.forest, RMI.AUC.surg$Performance[1,], LR2.AUC.surg$Performance[1,], SRrisks.AUC.surg$Performance[1,], ADNEXwo.AUC.surg$Performance[1,], ADNEX.AUC.surg$Performance[1,])
Summary.AUC.surg$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.surg.PI <- rbind(NA.forest, RMI.AUC.surg$Performance[2,], LR2.AUC.surg$Performance[2,], SRrisks.AUC.surg$Performance[2,], ADNEXwo.AUC.surg$Performance[2,], ADNEX.AUC.surg$Performance[2,])
Summary.AUC.surg.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.surg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.surg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.surg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.surg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.surg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.surg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.surg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.surg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.surg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.surg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.surg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.surg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.surg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.surg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.surg$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.surg$AUC, 3),
           lower = round(Summary.AUC.surg$LL, 3),
           upper = round(Summary.AUC.surg$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.81, 1))


#### 5.2.1.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Surgery.tiff", width = 14, height = 14, units = "cm", res = 300) 
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5Surg.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()

#### LR2 ####
LR2.surg <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5Surg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve")

#### Simple Rules Risk ####
SRrisk.surg <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5Surg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")

#### ADNEX with CA125 ####
ADNEXw.surg <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = IOTA5Surg.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")

#### ADNEX without CA125 ####
ADNEXwo.surg <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5Surg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")

#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.surg$Plot$y
p.LR2 = LR2.surg$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.surg$Plot$y
p.SR = SRrisk.surg$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.surg$Plot$y
p.ADNEXw = ADNEXw.surg$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.surg$Plot$y
p.ADNEXwo = ADNEXwo.surg$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.surg$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.surg$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.surg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.surg$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.surg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.surg$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.surg$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.surg$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.surg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.surg$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.surg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.surg$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.surg$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.surg$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.surg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.surg$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.surg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.surg$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.surg$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.surg$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.surg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.surg$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.surg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.surg$SumSlo$UL, 2), nsmall = 2), ")"))


## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Surgery.tiff", width = 14, height = 14, units = "cm", res = 300) # width = 21.5, height = 21.5
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) #cex = 1.2
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7) #x = 0.21, y = 0.02
dev.off()


#### 5.2.2 Patients who received at least one follow-up scan ####

IOTA5FU.Imp <- Iota5ValI[Iota5ValI$multiple.visits == 1 & Iota5ValI$CD_short != "U4",]
IOTA5FU.NI  <- Iota5ValNI[Iota5ValNI$multiple.visits == 1 & Iota5ValNI$CD_short != "U4",]

table(IOTA5FU.NI$center, IOTA5FU.NI$binaryCDcorrect) # Use the pooled data!!!

CentersDescrip = ddply(IOTA5FU.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.2.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
AUC.RMI.FU <- AUC.imp(pred = RMI, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp)$`Performance`[5:7]


#### LR2 ####

## AUC
AUC.LR2.FU <- AUC.imp(pred = LR2, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp)$`Performance`[5:7]


#### Simple Rules Risk ####

## AUC
AUC.SRrisks.FU <- AUC.imp(pred = SRrisks, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp)$`Performance`[5:7]


#### ADNEX with CA125 ####

## AUC
AUC.ADNEX.FU <- AUC.imp(pred = pmalw, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp)$`Performance`[5:7]


#### ADNEX without CA125 ####

## AUC
AUC.ADNEXwo.FU <- AUC.imp(pred = pmalwo, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp)$`Performance`[5:7]


#### Summary plot for AUC ####

NA.forest <- AUC.RMI.FU
NA.forest <- NA

Summary.AUC.FU <- rbind(NA.forest, AUC.RMI.FU, AUC.LR2.FU, AUC.SRrisks.FU, AUC.ADNEXwo.FU, AUC.ADNEX.FU)
Summary.AUC.FU$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.FU
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.RMI.FU$AUC, 2), nsmall = 2), " (", format(round(AUC.RMI.FU$LL, 2), nsmall = 2), "; ", format(round(AUC.RMI.FU$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.LR2.FU$AUC, 2), nsmall = 2), " (", format(round(AUC.LR2.FU$LL, 2), nsmall = 2), "; ", format(round(AUC.LR2.FU$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.SRrisks.FU$AUC, 2), nsmall = 2), " (", format(round(AUC.SRrisks.FU$LL, 2), nsmall = 2), "; ", format(round(AUC.SRrisks.FU$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwo.FU$AUC, 2), nsmall = 2), " (", format(round(AUC.ADNEXwo.FU$LL, 2), nsmall = 2), "; ", format(round(AUC.ADNEXwo.FU$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEX.FU$AUC, 2), nsmall = 2), " (", format(round(AUC.ADNEX.FU$LL, 2), nsmall = 2), "; ", format(round(AUC.ADNEX.FU$UL, 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.FU$AUC, 3),
           lower = round(Summary.AUC.FU$LL, 3),
           upper = round(Summary.AUC.FU$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.61, 1))


#### 5.2.2.2 Calibration of the risk of malignancy ####

#### RMI ####
png("Results Paper 2/With cysts/Supplementary/Calibration RMI - FollowUp.png", width = 14, height = 14, units = "cm", res = 300) #width = 21.5, height = 21.5
CalPool.RMI(pred = RMI, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp, title = "")
dev.off()

#### LR2 ####
LR2.FU <- CalPool.Imp(pred = LR2, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp, title = "LR2: Calibration curve")

#### Simple Rules Risk ####
SRrisk.FU <- CalPool.Imp(pred = SRrisks, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp, title = "Simple Rules Risk: Calibration curve")

#### ADNEX with CA125 ####
ADNEXw.FU <- CalPool.Imp(pred = pmalw, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp, title = "ADNEX with CA125: Calibration curve")

#### ADNEX without CA125 ####
ADNEXwo.FU <- CalPool.Imp(pred = pmalwo, outcome = binaryCDcorrect, imp = .imp, data = IOTA5FU.Imp, title = "ADNEX without CA125: Calibration curve")

#### Summary plot for Calibration ####
table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')
table[1, 2:3] <- c(paste0(format(round(as.numeric(LR2.FU$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(LR2.FU$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(LR2.FU$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(LR2.FU$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(LR2.FU$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(LR2.FU$Performance[2,3]), 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(as.numeric(SRrisk.FU$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(SRrisk.FU$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(SRrisk.FU$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(SRrisk.FU$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(SRrisk.FU$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(SRrisk.FU$Performance[2,3]), 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(as.numeric(ADNEXw.FU$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXw.FU$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXw.FU$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(ADNEXw.FU$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXw.FU$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXw.FU$Performance[2,3]), 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(as.numeric(ADNEXwo.FU$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXwo.FU$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXwo.FU$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(ADNEXwo.FU$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXwo.FU$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXwo.FU$Performance[2,3]), 2), nsmall = 2), ")"))

x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - FollowUp.tiff", width = 14, height = 14, units = "cm", res = 300) # width = 21.5, height = 21.5
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(LR2.FU$Predictions$p_pred, LR2.FU$Predictions$p_obs, lwd = 2, col = "blue")
lines(SRrisk.FU$Predictions$p_pred, SRrisk.FU$Predictions$p_obs, lwd = 2, col = "red")
lines(ADNEXwo.FU$Predictions$p_pred, ADNEXwo.FU$Predictions$p_obs, lwd = 2, col = "darkgreen")
lines(ADNEXw.FU$Predictions$p_pred, ADNEXw.FU$Predictions$p_obs, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 0.65, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX w/o CA125", "ADNEX w/ CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) #cex = 1.2
addtable2plot(x = -0.035, y = 0.8, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


#### 5.2.3 Suggested management surgery ####

IOTA5MSurg.Imp <- Iota5ValI[Iota5ValI$initial.policy == "surgery",]
IOTA5MSurg.NI  <- Iota5ValNI[Iota5ValNI$initial.policy == "surgery",]

table(IOTA5MSurg.NI$centerRE, IOTA5MSurg.NI$binaryCDcorrect)

CentersDescrip = ddply(IOTA5MSurg.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.3.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.Msurg <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5MSurg.Imp, method.MA = "SJ")

forestplot(RMI.AUC.Msurg$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.Msurg$dataPlot$AUC,
           lower = RMI.AUC.Msurg$dataPlot$LL,
           upper = RMI.AUC.Msurg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.Msurg$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.Msurg$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.Msurg <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5MSurg.Imp, method.MA = "SJ")

forestplot(LR2.AUC.Msurg$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.Msurg$dataPlot$AUC,
           lower = LR2.AUC.Msurg$dataPlot$LL,
           upper = LR2.AUC.Msurg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.Msurg$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.Msurg$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.Msurg <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5MSurg.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.Msurg$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.Msurg$dataPlot$AUC,
           lower = SRrisks.AUC.Msurg$dataPlot$LL,
           upper = SRrisks.AUC.Msurg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.Msurg$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.Msurg$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.Msurg <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5MSurg.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.Msurg$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.Msurg$dataPlot$AUC,
           lower = ADNEX.AUC.Msurg$dataPlot$LL,
           upper = ADNEX.AUC.Msurg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.Msurg$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.Msurg$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.Msurg <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5MSurg.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.Msurg$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.Msurg$dataPlot$AUC,
           lower = ADNEXwo.AUC.Msurg$dataPlot$LL,
           upper = ADNEXwo.AUC.Msurg$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.Msurg$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.Msurg$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.Msurg$Performance[1,]
NA.forest <- NA

Summary.AUC.Msurg <- rbind(NA.forest, RMI.AUC.Msurg$Performance[1,], LR2.AUC.Msurg$Performance[1,], SRrisks.AUC.Msurg$Performance[1,], ADNEXwo.AUC.Msurg$Performance[1,], ADNEX.AUC.Msurg$Performance[1,])
Summary.AUC.Msurg$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.Msurg.PI <- rbind(NA.forest, RMI.AUC.Msurg$Performance[2,], LR2.AUC.Msurg$Performance[2,], SRrisks.AUC.Msurg$Performance[2,], ADNEXwo.AUC.Msurg$Performance[2,], ADNEX.AUC.Msurg$Performance[2,])
Summary.AUC.Msurg.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.Msurg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.Msurg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.Msurg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.Msurg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.Msurg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.Msurg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.Msurg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.Msurg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.Msurg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.Msurg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.Msurg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.Msurg$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.Msurg$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.Msurg$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.Msurg$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.Msurg$AUC, 3),
           lower = round(Summary.AUC.Msurg$LL, 3),
           upper = round(Summary.AUC.Msurg$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.81, 1))


#### 5.2.3.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Suggested Surgery.tiff", width = 14, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5MSurg.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()

#### LR2 ####
LR2.Msurg <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5MSurg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve")

#### Simple Rules Risk ####
SRrisk.Msurg <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5MSurg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")

#### ADNEX with CA125 ####
ADNEXw.Msurg <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = IOTA5MSurg.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")

#### ADNEX without CA125 ####
ADNEXwo.Msurg <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5MSurg.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")

#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.Msurg$Plot$y
p.LR2 = LR2.Msurg$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.Msurg$Plot$y
p.SR = SRrisk.Msurg$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.Msurg$Plot$y
p.ADNEXw = ADNEXw.Msurg$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.Msurg$Plot$y
p.ADNEXwo = ADNEXwo.Msurg$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.Msurg$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.Msurg$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.Msurg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.Msurg$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.Msurg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.Msurg$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.Msurg$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.Msurg$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.Msurg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.Msurg$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.Msurg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.Msurg$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.Msurg$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.Msurg$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.Msurg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.Msurg$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.Msurg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.Msurg$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.Msurg$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.Msurg$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.Msurg$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.Msurg$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.Msurg$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.Msurg$SumSlo$UL, 2), nsmall = 2), ")"))


## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Suggested Surgery.tiff", width = 14, height = 14, units = "cm", res = 300) # width = 21.5, height = 21.5
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) #cex = 1.2
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7) # x = 0.21, y = 0.02
dev.off()


#### 5.2.4 Suggested management conservative ####

IOTA5MCon.Imp <- Iota5ValI[Iota5ValI$initial.policy == "conservative",]
IOTA5MCon.NI  <- Iota5ValNI[Iota5ValNI$initial.policy == "conservative",]

table(IOTA5MCon.NI$centerRE, IOTA5MCon.NI$binaryCDcorrect) ## Poolen

CentersDescrip = ddply(IOTA5MCon.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.4.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
AUC.RMI.Mcon <- AUC.imp(pred = RMI, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp)$`Performance`[5:7]


#### LR2 ####

## AUC
AUC.LR2.Mcon <- AUC.imp(pred = LR2, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp)$`Performance`[5:7]


#### Simple Rules Risk ####

## AUC
AUC.SRrisks.Mcon <- AUC.imp(pred = SRrisks, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp)$`Performance`[5:7]


#### ADNEX with CA125 ####

## AUC
AUC.ADNEX.Mcon <- AUC.imp(pred = pmalw, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp)$`Performance`[5:7]


#### ADNEX without CA125 ####

## AUC
AUC.ADNEXwo.Mcon <- AUC.imp(pred = pmalwo, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp)$`Performance`[5:7]


#### Summary plot for AUC ####

NA.forest <- AUC.RMI.Mcon
NA.forest <- NA

Summary.AUC.Mcon <- rbind(NA.forest, AUC.RMI.Mcon, AUC.LR2.Mcon, AUC.SRrisks.Mcon, AUC.ADNEXwo.Mcon, AUC.ADNEX.Mcon)
Summary.AUC.Mcon$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.Mcon
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(AUC.RMI.Mcon$AUC, 2), nsmall = 2), " (", format(round(AUC.RMI.Mcon$LL, 2), nsmall = 2), "; ", format(round(AUC.RMI.Mcon$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.LR2.Mcon$AUC, 2), nsmall = 2), " (", format(round(AUC.LR2.Mcon$LL, 2), nsmall = 2), "; ", format(round(AUC.LR2.Mcon$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.SRrisks.Mcon$AUC, 2), nsmall = 2), " (", format(round(AUC.SRrisks.Mcon$LL, 2), nsmall = 2), "; ", format(round(AUC.SRrisks.Mcon$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEXwo.Mcon$AUC, 2), nsmall = 2), " (", format(round(AUC.ADNEXwo.Mcon$LL, 2), nsmall = 2), "; ", format(round(AUC.ADNEXwo.Mcon$UL, 2), nsmall = 2), ")", sep = ""),
    paste(format(round(AUC.ADNEX.Mcon$AUC, 2), nsmall = 2), " (", format(round(AUC.ADNEX.Mcon$LL, 2), nsmall = 2), "; ", format(round(AUC.ADNEX.Mcon$UL, 2), nsmall = 2), ")", sep = ""))
)

forestplot(labeltext = tabletext,
           title = "Overall AUC of the prediction models",
           mean = round(Summary.AUC.Mcon$AUC, 3),
           lower = round(Summary.AUC.Mcon$LL, 3),
           upper = round(Summary.AUC.Mcon$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.61, 1))


#### 5.2.4.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Conservative management.tiff", width = 14, height = 14, units = "cm", res = 300) 
CalPool.RMI(pred = RMI, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp, title = "")
dev.off()

#### LR2 ####
LR2.Mcon <- CalPool.Imp(pred = LR2, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp, title = "LR2: Calibration curve")

#### Simple Rules Risk ####
SRrisk.Mcon <- CalPool.Imp(pred = SRrisks, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp, title = "Simple Rules Risk: Calibration curve")

#### ADNEX with CA125 ####
ADNEXw.Mcon <- CalPool.Imp(pred = pmalw, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp, title = "ADNEX with CA125: Calibration curve")

#### ADNEX without CA125 ####
ADNEXwo.Mcon <- CalPool.Imp(pred = pmalwo, outcome = binaryCDcorrect, imp = .imp, data = IOTA5MCon.Imp, title = "ADNEX without CA125: Calibration curve")

#### Summary plot for Calibration ####
table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')
table[1, 2:3] <- c(paste0(format(round(as.numeric(LR2.Mcon$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(LR2.Mcon$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(LR2.Mcon$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(LR2.Mcon$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(LR2.Mcon$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(LR2.Mcon$Performance[2,3]), 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(as.numeric(SRrisk.Mcon$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(SRrisk.Mcon$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(SRrisk.Mcon$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(SRrisk.Mcon$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(SRrisk.Mcon$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(SRrisk.Mcon$Performance[2,3]), 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(as.numeric(ADNEXw.Mcon$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXw.Mcon$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXw.Mcon$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(ADNEXw.Mcon$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXw.Mcon$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXw.Mcon$Performance[2,3]), 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(as.numeric(ADNEXwo.Mcon$Performance[1,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXwo.Mcon$Performance[1,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXwo.Mcon$Performance[1,3]), 2), nsmall = 2), ")"), paste0(format(round(as.numeric(ADNEXwo.Mcon$Performance[2,1]), 2), nsmall = 2), " (", format(round(as.numeric(ADNEXwo.Mcon$Performance[2,2]), 2), nsmall = 2), " to ", format(round(as.numeric(ADNEXwo.Mcon$Performance[2,3]), 2), nsmall = 2), ")"))

x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Conservative management.tiff", width = 14, height = 14, units = "cm", res = 300) 
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(LR2.Mcon$Predictions$p_pred, LR2.Mcon$Predictions$p_obs, lwd = 2, col = "blue")
lines(SRrisk.Mcon$Predictions$p_pred, SRrisk.Mcon$Predictions$p_obs, lwd = 2, col = "red")
lines(ADNEXwo.Mcon$Predictions$p_pred, ADNEXwo.Mcon$Predictions$p_obs, lwd = 2, col = "darkgreen")
lines(ADNEXw.Mcon$Predictions$p_pred, ADNEXw.Mcon$Predictions$p_obs, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 0.65, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX w/o CA125", "ADNEX w/ CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = -0.035, y = 0.8, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


#### 5.2.5 Premenopausal patients ####

IOTA5Pre.Imp <- Iota5ValI[Iota5ValI$Postmenopausal2 == "no",]
IOTA5Pre.NI  <- Iota5ValNI[Iota5ValNI$Postmenopausal2 == "no",]

table(IOTA5Pre.NI$centerRE, IOTA5Pre.NI$binaryCDcorrect)

CentersDescrip = ddply(IOTA5Pre.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.5.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.pre <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Pre.Imp, method.MA = "SJ")

forestplot(RMI.AUC.pre$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.pre$dataPlot$AUC,
           lower = RMI.AUC.pre$dataPlot$LL,
           upper = RMI.AUC.pre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.pre$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.pre$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.pre <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Pre.Imp, method.MA = "SJ")

forestplot(LR2.AUC.pre$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.pre$dataPlot$AUC,
           lower = LR2.AUC.pre$dataPlot$LL,
           upper = LR2.AUC.pre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.pre$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.pre$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.pre <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Pre.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.pre$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.pre$dataPlot$AUC,
           lower = SRrisks.AUC.pre$dataPlot$LL,
           upper = SRrisks.AUC.pre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.pre$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.pre$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.pre <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Pre.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.pre$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.pre$dataPlot$AUC,
           lower = ADNEX.AUC.pre$dataPlot$LL,
           upper = ADNEX.AUC.pre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.pre$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.pre$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.pre <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Pre.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.pre$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.pre$dataPlot$AUC,
           lower = ADNEXwo.AUC.pre$dataPlot$LL,
           upper = ADNEXwo.AUC.pre$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.pre$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.pre$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.pre$Performance[1,]
NA.forest <- NA

Summary.AUC.pre <- rbind(NA.forest, RMI.AUC.pre$Performance[1,], LR2.AUC.pre$Performance[1,], SRrisks.AUC.pre$Performance[1,], ADNEXwo.AUC.pre$Performance[1,], ADNEX.AUC.pre$Performance[1,])
Summary.AUC.pre$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.pre.PI <- rbind(NA.forest, RMI.AUC.pre$Performance[2,], LR2.AUC.pre$Performance[2,], SRrisks.AUC.pre$Performance[2,], ADNEXwo.AUC.pre$Performance[2,], ADNEX.AUC.pre$Performance[2,])
Summary.AUC.pre.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.pre
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.pre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.pre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.pre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.pre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.pre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.pre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.pre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.pre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.pre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.pre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.pre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.pre$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.pre$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.pre$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.pre$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.pre$AUC, 3),
           lower = round(Summary.AUC.pre$LL, 3),
           upper = round(Summary.AUC.pre$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.81, 1))


#### 5.2.5.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Premenopausal.tiff", width = 14, height = 14, units = "cm", res = 300) # width = 21.5, height = 21.5
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5Pre.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2.pre <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5Pre.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve")


#### Simple Rules Risk ####
SRrisk.pre <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5Pre.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXw.pre <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = IOTA5Pre.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwo.pre <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5Pre.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.pre$Plot$y
p.LR2 = LR2.pre$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.pre$Plot$y
p.SR = SRrisk.pre$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.pre$Plot$y
p.ADNEXw = ADNEXw.pre$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.pre$Plot$y
p.ADNEXwo = ADNEXwo.pre$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.pre$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.pre$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.pre$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.pre$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.pre$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.pre$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.pre$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.pre$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.pre$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.pre$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.pre$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.pre$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.pre$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.pre$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.pre$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.pre$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.pre$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.pre$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.pre$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.pre$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.pre$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.pre$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.pre$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.pre$SumSlo$UL, 2), nsmall = 2), ")"))


## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Premenopausal.tiff", width = 14, height = 14, units = "cm", res = 300) # width = 21.5, height = 21.5
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) #cex = 1.2
addtable2plot(x = 0.16, y = -0.01, table = table, display.colnames= TRUE, cex = 0.7) # x = 0.21, y = 0.02
dev.off()


#### 5.2.6 Postmenopausal patients ####

IOTA5Post.Imp <- Iota5ValI[Iota5ValI$Postmenopausal2 == "yes",]
IOTA5Post.NI  <- Iota5ValNI[Iota5ValNI$Postmenopausal2 == "yes",]

table(IOTA5Post.NI$centerRE, IOTA5Post.NI$binaryCDcorrect)

CentersDescrip = ddply(IOTA5Post.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.6.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.post <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Post.Imp, method.MA = "SJ")

forestplot(RMI.AUC.post$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.post$dataPlot$AUC,
           lower = RMI.AUC.post$dataPlot$LL,
           upper = RMI.AUC.post$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.post$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.post$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.post <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Post.Imp, method.MA = "SJ")

forestplot(LR2.AUC.post$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.post$dataPlot$AUC,
           lower = LR2.AUC.post$dataPlot$LL,
           upper = LR2.AUC.post$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.post$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.post$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.post <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Post.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.post$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.post$dataPlot$AUC,
           lower = SRrisks.AUC.post$dataPlot$LL,
           upper = SRrisks.AUC.post$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.post$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.post$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.post <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Post.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.post$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.post$dataPlot$AUC,
           lower = ADNEX.AUC.post$dataPlot$LL,
           upper = ADNEX.AUC.post$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.post$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.post$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.post <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Post.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.post$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.post$dataPlot$AUC,
           lower = ADNEXwo.AUC.post$dataPlot$LL,
           upper = ADNEXwo.AUC.post$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.post$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.post$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.post$Performance[1,]
NA.forest <- NA

Summary.AUC.post <- rbind(NA.forest, RMI.AUC.post$Performance[1,], LR2.AUC.post$Performance[1,], SRrisks.AUC.post$Performance[1,], ADNEXwo.AUC.post$Performance[1,], ADNEX.AUC.post$Performance[1,])
Summary.AUC.post$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.post.PI <- rbind(NA.forest, RMI.AUC.post$Performance[2,], LR2.AUC.post$Performance[2,], SRrisks.AUC.post$Performance[2,], ADNEXwo.AUC.post$Performance[2,], ADNEX.AUC.post$Performance[2,])
Summary.AUC.post.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.post
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.post$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.post$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.post$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.post$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.post$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.post$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.post$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.post$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.post$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.post$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.post$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.post$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.post$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.post$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.post$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.post$AUC, 3),
           lower = round(Summary.AUC.post$LL, 3),
           upper = round(Summary.AUC.post$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))


#### 5.2.6.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Postmenopausal.tiff", width = 14, height = 14, units = "cm", res = 300) #width = 21.5, height = 21.5
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5Post.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2.post <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5Post.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve", dostats = T)


#### Simple Rules Risk ####
SRrisk.post <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5Post.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXw.post <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, center = centerRE, data = IOTA5Post.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwo.post <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5Post.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.post$Plot$y
p.LR2 = LR2.post$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.post$Plot$y
p.SR = SRrisk.post$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.post$Plot$y
p.ADNEXw = ADNEXw.post$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.post$Plot$y
p.ADNEXwo = ADNEXwo.post$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.post$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.post$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.post$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.post$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.post$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.post$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.post$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.post$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.post$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.post$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.post$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.post$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.post$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.post$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.post$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.post$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.post$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.post$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.post$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.post$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.post$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.post$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.post$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.post$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Postmenopausal.tiff", width = 14, height = 14, units = "cm", res = 300) 
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) # cex.lab = 1.5, cex.axis = 1.5
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) 
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7) 
dev.off()


#### 5.2.7 Patients examined in oncology centers ####

IOTA5Onco.Imp <- Iota5ValI[Iota5ValI$oncocenter == 1,]
IOTA5Onco.NI  <- Iota5ValNI[Iota5ValNI$oncocenter == 1,]

table(IOTA5Onco.NI$center, IOTA5Onco.NI$binaryCDcorrect)

CentersDescrip = ddply(IOTA5Onco.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.7.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.onco <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Onco.Imp, method.MA = "SJ")

forestplot(RMI.AUC.onco$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.onco$dataPlot$AUC,
           lower = RMI.AUC.onco$dataPlot$LL,
           upper = RMI.AUC.onco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.onco$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.onco$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.onco <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Onco.Imp, method.MA = "SJ")

forestplot(LR2.AUC.onco$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.onco$dataPlot$AUC,
           lower = LR2.AUC.onco$dataPlot$LL,
           upper = LR2.AUC.onco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.onco$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.onco$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.onco <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Onco.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.onco$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.onco$dataPlot$AUC,
           lower = SRrisks.AUC.onco$dataPlot$LL,
           upper = SRrisks.AUC.onco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.onco$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.onco$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.onco <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Onco.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.onco$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.onco$dataPlot$AUC,
           lower = ADNEX.AUC.onco$dataPlot$LL,
           upper = ADNEX.AUC.onco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.onco$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.onco$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.onco <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5Onco.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.onco$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.onco$dataPlot$AUC,
           lower = ADNEXwo.AUC.onco$dataPlot$LL,
           upper = ADNEXwo.AUC.onco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.onco$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.onco$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.onco$Performance[1,]
NA.forest <- NA

Summary.AUC.onco <- rbind(NA.forest, RMI.AUC.onco$Performance[1,], LR2.AUC.onco$Performance[1,], SRrisks.AUC.onco$Performance[1,], ADNEXwo.AUC.onco$Performance[1,], ADNEX.AUC.onco$Performance[1,])
Summary.AUC.onco$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.onco.PI <- rbind(NA.forest, RMI.AUC.onco$Performance[2,], LR2.AUC.onco$Performance[2,], SRrisks.AUC.onco$Performance[2,], ADNEXwo.AUC.onco$Performance[2,], ADNEX.AUC.onco$Performance[2,])
Summary.AUC.onco.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.onco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.onco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.onco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.onco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.onco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.onco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.onco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.onco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.onco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.onco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.onco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.onco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.onco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.onco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.onco$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.onco$AUC, 3),
           lower = round(Summary.AUC.onco$LL, 3),
           upper = round(Summary.AUC.onco$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.81, 1))


#### 5.2.7.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Oncocenter.tiff", width = 14, height = 14, units = "cm", res = 300) 
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5Onco.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2.onco <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5Onco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve")


#### Simple Rules Risk ####
SRrisk.onco <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5Onco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXw.onco <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = IOTA5Onco.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwo.onco <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5Onco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.onco$Plot$y
p.LR2 = LR2.onco$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.onco$Plot$y
p.SR = SRrisk.onco$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.onco$Plot$y
p.ADNEXw = ADNEXw.onco$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.onco$Plot$y
p.ADNEXwo = ADNEXwo.onco$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.onco$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.onco$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.onco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.onco$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.onco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.onco$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.onco$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.onco$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.onco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.onco$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.onco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.onco$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.onco$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.onco$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.onco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.onco$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.onco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.onco$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.onco$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.onco$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.onco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.onco$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.onco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.onco$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Oncocenter.tiff", width = 14, height = 14, units = "cm", res = 300) 
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) 
addtable2plot(x = 0.17, y = 0.02, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


#### 5.2.8 Patients examined in other centers ####

IOTA5NonOnco.Imp <- Iota5ValI[Iota5ValI$oncocenter == 0,]
IOTA5NonOnco.NI  <- Iota5ValNI[Iota5ValNI$oncocenter == 0,]

CentersDescrip = ddply(IOTA5NonOnco.NI, .(center),
                       function(x){
                         data.frame(Included = if(x$center %in% SelectedCenters){
                           "yes"
                         } else{
                           "no"
                         },
                         Type =
                           if(is.na(x$oncocenter)){
                             "-"
                           } else if(x$oncocenter == 1){
                             "Oncological center"
                           } else{
                             "Non-oncological center"
                           },
                         N = nrow(x),
                         n_Benign = sum(grepl("Benign", x$ClinicalDiagnosis)),
                         n_Maligne = sum(grepl("Malignant", x$ClinicalDiagnosis)),
                         n_Uncertain = sum(grepl("Uncertain", x$ClinicalDiagnosis)),
                         n_surgery = sum(x$immed.oper & x$CD_short != "U4"),
                         n_followup = sum(x$multiple.visits == 1 & x$CD_short != "U4"),
                         n_nofollowup = sum(x$immed.oper != 1 & x$multiple.visits != 1 | x$CD_short == "U4"),
                         check.names = F)
                       })
CentersDescrip

#### 5.2.8.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMI.AUC.NonOnco <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5NonOnco.Imp, method.MA = "SJ")

forestplot(RMI.AUC.NonOnco$Plot,
           align = c("l", "c", "c"),
           mean = RMI.AUC.NonOnco$dataPlot$AUC,
           lower = RMI.AUC.NonOnco$dataPlot$LL,
           upper = RMI.AUC.NonOnco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMI.AUC.NonOnco$IncludedCenters)), TRUE),
           title = "RMI: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = RMI.AUC.NonOnco$Performance$AUC)


#### LR2 ####

## AUC
LR2.AUC.NonOnco <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5NonOnco.Imp, method.MA = "SJ")

forestplot(LR2.AUC.NonOnco$Plot,
           align = c("l", "c", "c"),
           mean = LR2.AUC.NonOnco$dataPlot$AUC,
           lower = LR2.AUC.NonOnco$dataPlot$LL,
           upper = LR2.AUC.NonOnco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2.AUC.NonOnco$IncludedCenters)), TRUE),
           title = "LR2: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = LR2.AUC.NonOnco$Performance$AUC)


#### Simple Rules Risk ####

## AUC
SRrisks.AUC.NonOnco <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5NonOnco.Imp, method.MA = "SJ")

forestplot(SRrisks.AUC.NonOnco$Plot,
           align = c("l", "c", "c"),
           mean = SRrisks.AUC.NonOnco$dataPlot$AUC,
           lower = SRrisks.AUC.NonOnco$dataPlot$LL,
           upper = SRrisks.AUC.NonOnco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisks.AUC.NonOnco$IncludedCenters)), TRUE),
           title = "Simple Rules Risk: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = SRrisks.AUC.NonOnco$Performance$AUC)


#### ADNEX with CA125 ####

## AUC
ADNEX.AUC.NonOnco <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5NonOnco.Imp, method.MA = "SJ")

forestplot(ADNEX.AUC.NonOnco$Plot,
           align = c("l", "c", "c"),
           mean = ADNEX.AUC.NonOnco$dataPlot$AUC,
           lower = ADNEX.AUC.NonOnco$dataPlot$LL,
           upper = ADNEX.AUC.NonOnco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEX.AUC.NonOnco$IncludedCenters)), TRUE),
           title = "ADNEX with CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEX.AUC.NonOnco$Performance$AUC)


#### ADNEX without CA125 ####

## AUC
ADNEXwo.AUC.NonOnco <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = IOTA5NonOnco.Imp, method.MA = "SJ")

forestplot(ADNEXwo.AUC.NonOnco$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwo.AUC.NonOnco$dataPlot$AUC,
           lower = ADNEXwo.AUC.NonOnco$dataPlot$LL,
           upper = ADNEXwo.AUC.NonOnco$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwo.AUC.NonOnco$IncludedCenters)), TRUE),
           title = "ADNEX without CA125: AUC per center",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           col = fpColors(box = "gray", lines = "black"),
           zero = ADNEXwo.AUC.NonOnco$Performance$AUC)


#### Summary plot for AUC ####

NA.forest <- RMI.AUC.NonOnco$Performance[1,]
NA.forest <- NA

Summary.AUC.NonOnco <- rbind(NA.forest, RMI.AUC.NonOnco$Performance[1,], LR2.AUC.NonOnco$Performance[1,], SRrisks.AUC.NonOnco$Performance[1,], ADNEXwo.AUC.NonOnco$Performance[1,], ADNEX.AUC.NonOnco$Performance[1,])
Summary.AUC.NonOnco$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.NonOnco.PI <- rbind(NA.forest, RMI.AUC.NonOnco$Performance[2,], LR2.AUC.NonOnco$Performance[2,], SRrisks.AUC.NonOnco$Performance[2,], ADNEXwo.AUC.NonOnco$Performance[2,], ADNEX.AUC.NonOnco$Performance[2,])
Summary.AUC.NonOnco.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.NonOnco
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMI.AUC.NonOnco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMI.AUC.NonOnco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(RMI.AUC.NonOnco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2.AUC.NonOnco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2.AUC.NonOnco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(LR2.AUC.NonOnco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisks.AUC.NonOnco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisks.AUC.NonOnco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(SRrisks.AUC.NonOnco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwo.AUC.NonOnco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwo.AUC.NonOnco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEXwo.AUC.NonOnco$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEX.AUC.NonOnco$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEX.AUC.NonOnco$Performance$LL[1], 2), nsmall = 2), "; ", format(round(ADNEX.AUC.NonOnco$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMI.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"))
)

forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC.NonOnco$AUC, 3),
           lower = round(Summary.AUC.NonOnco$LL, 3),
           upper = round(Summary.AUC.NonOnco$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.75, 0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.70, 1))


#### 5.2.8.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Non oncocenter.tiff", width = 14, height = 14, units = "cm", res = 300) 
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = IOTA5NonOnco.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2.NonOnco <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = IOTA5NonOnco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "LR2: Calibration curve")


#### Simple Rules Risk ####
SRrisk.NonOnco <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = IOTA5NonOnco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXw.NonOnco <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = IOTA5NonOnco.Imp, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwo.NonOnco <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = IOTA5NonOnco.Imp, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2.NonOnco$Plot$y
p.LR2 = LR2.NonOnco$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRrisk.NonOnco$Plot$y
p.SR = SRrisk.NonOnco$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXw.NonOnco$Plot$y
p.ADNEXw = ADNEXw.NonOnco$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwo.NonOnco$Plot$y
p.ADNEXwo = ADNEXwo.NonOnco$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2.NonOnco$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2.NonOnco$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2.NonOnco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2.NonOnco$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2.NonOnco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2.NonOnco$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRrisk.NonOnco$SumInt$Est, 2), nsmall = 2), " (", format(round(SRrisk.NonOnco$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRrisk.NonOnco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRrisk.NonOnco$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRrisk.NonOnco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRrisk.NonOnco$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXw.NonOnco$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXw.NonOnco$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.NonOnco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXw.NonOnco$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXw.NonOnco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXw.NonOnco$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwo.NonOnco$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.NonOnco$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.NonOnco$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwo.NonOnco$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwo.NonOnco$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwo.NonOnco$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Non oncocenter.tiff", width = 14, height = 14, units = "cm", res = 300) 
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) 
addtable2plot(x = 0.16, y = -0.01, table = table, display.colnames= TRUE, cex = 0.7)
dev.off()


#### 5.2.9 Summary forest plots of the subgroups ####

Summary.AUC <- rbind(NA.forest, Summary.AUC.surg[4:7], NA.forest, Summary.AUC.FU[1:4], 
                     NA.forest, Summary.AUC.Msurg[4:7], NA.forest, Summary.AUC.Mcon[1:4], 
                     NA.forest, Summary.AUC.pre[4:7], NA.forest, Summary.AUC.post[4:7], 
                     NA.forest, Summary.AUC.onco[4:7], NA.forest, Summary.AUC.NonOnco[4:7])

Summary.AUC.pre
tabletext <- cbind(
  c('Surgery within 120 days after the first', 'scan without any follow-up scan (n=2489)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'At least one follow-up scan (n=1958)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Suggested management surgery (n=2579)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Suggested management conservative (n=2326)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Premenopausal patients (n=2754)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Postmenopausal patients (n=2151)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Oncology centre (n=3094)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125',
    '', 'Other centre (n=1811)', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', '',
    paste(format(round(Summary.AUC$AUC[3], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[3], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[3], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[4], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[4], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[4], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[5], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[5], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[5], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[6], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[6], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[6], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[7], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[7], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[7], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[10], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[10], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[10], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[11], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[11], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[11], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[12], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[12], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[12], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[13], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[13], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[13], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[14], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[14], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[14], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[17], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[17], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[17], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[18], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[18], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[18], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[19], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[19], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[19], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[20], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[20], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[20], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[21], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[21], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[21], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[24], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[24], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[24], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[25], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[25], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[25], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[26], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[26], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[26], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[27], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[27], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[27], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[28], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[28], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[28], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[31], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[31], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[31], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[32], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[32], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[32], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[33], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[33], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[33], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[34], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[34], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[34], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[35], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[35], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[35], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[38], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[38], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[38], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[39], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[39], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[39], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[40], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[40], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[40], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[41], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[41], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[41], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[42], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[42], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[42], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[45], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[45], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[45], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[46], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[46], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[46], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[47], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[47], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[47], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[48], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[48], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[48], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[49], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[49], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[49], 2), nsmall = 2), ")", sep = ""),
    '', '', 
    paste(format(round(Summary.AUC$AUC[52], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[52], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[52], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[53], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[53], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[53], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[54], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[54], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[54], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[55], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[55], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[55], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(Summary.AUC$AUC[56], 2), nsmall = 2), " (", format(round(Summary.AUC$LL[56], 2), nsmall = 2), " to ", format(round(Summary.AUC$UL[56], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', '', 
    paste0("(", format(round(RMI.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.surg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.surg$Performance$UL[2], 2), nsmall = 2), ")"),
    '', '', 
    '', '', '', '', '',
    '', '', 
    paste0("(", format(round(RMI.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.Msurg$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.Msurg$Performance$UL[2], 2), nsmall = 2), ")"),
    '', '', 
    '', '', '', '', '',
    '', '', 
    paste0("(", format(round(RMI.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.pre$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.pre$Performance$UL[2], 2), nsmall = 2), ")"),
    '', '', 
    paste0("(", format(round(RMI.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.post$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.post$Performance$UL[2], 2), nsmall = 2), ")"),
    '', '', 
    paste0("(", format(round(RMI.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.onco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.onco$Performance$UL[2], 2), nsmall = 2), ")"),
    '', '', 
    paste0("(", format(round(RMI.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMI.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisks.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisks.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwo.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwo.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEX.AUC.NonOnco$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEX.AUC.NonOnco$Performance$UL[2], 2), nsmall = 2), ")"))
)


setEPS()
postscript("Results Paper 2/With cysts/Paper/Summary AUC - All Subgroups PI.eps", width = 16, height = 20, horizontal = FALSE, onefile = FALSE, paper = "special")
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .7,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graph.pos = 3,
           graphwidth = unit(9, "cm"),
           xticks = c(0.60, 0.70, 0.80, 0.9, 1), xlog = TRUE, clip = c(0.61, 1))
dev.off()



#-----------------------------------#
#### 5.3 Omit uncertain outcomes ####
#-----------------------------------#

Iota5Ni_cert <- Iota5ValNI[!is.na(Iota5ValNI$binaryCD), ]
table(Iota5Ni_cert$CD_short)
Iota5Imputed_cert <- Iota5ValI[!is.na(Iota5ValI$binaryCD), ]

#### 5.3.1 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMIcert.AUC <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5Imputed_cert, method.MA = "SJ")

RMIcert.AUC$dataPlot <- RMIcert.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
RMIcert.AUC$Plot <- RMIcert.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

RMIcert.AUC$Plot[17, "RRauc"] <- "   "
RMIcert.AUC$Plot[19, "RRauc"] <- "        (0.73 to 0.96)"
RMIcert.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest RMI - Without uncertain PI.tiff", width = 1000, height = 650)
forestplot(RMIcert.AUC$Plot,
           align = c("l", "c", "c"),
           mean = RMIcert.AUC$dataPlot$AUC,
           lower = RMIcert.AUC$dataPlot$LL,
           upper = RMIcert.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMIcert.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(RMIcert.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = RMIcert.AUC$Performance$AUC)
dev.off()


#### LR2 ####

## AUC
LR2cert.AUC <- AUC.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, method.MA = "SJ", titleGraph = "LR2: AUC per center")

LR2cert.AUC$dataPlot <- LR2cert.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
LR2cert.AUC$Plot <- LR2cert.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

LR2cert.AUC$Plot[17, "RRauc"] <- "   "
LR2cert.AUC$Plot[19, "RRauc"] <- "        (0.84 to 0.96)"
LR2cert.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest LR2 - Without uncertain PI.tiff", width = 1000, height = 650)
forestplot(LR2cert.AUC$Plot,
           align = c("l", "c", "c"),
           mean = LR2cert.AUC$dataPlot$AUC,
           lower = LR2cert.AUC$dataPlot$LL,
           upper = LR2cert.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2cert.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(LR2cert.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = LR2cert.AUC$Performance$AUC)
dev.off()


#### Simple Rules Risk ####

## AUC
SRrisksCert.AUC <- AUC.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, method.MA = "SJ", titleGraph = "Simple Rules risk: AUC per center")

SRrisksCert.AUC$dataPlot <- SRrisksCert.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
SRrisksCert.AUC$Plot <- SRrisksCert.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

SRrisksCert.AUC$Plot[17, "RRauc"] <- "   "
SRrisksCert.AUC$Plot[19, "RRauc"] <- "        (0.83 to 0.98)"
SRrisksCert.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest SRrisks - Without uncertain PI.tiff", width = 1000, height = 650)
forestplot(SRrisksCert.AUC$Plot,
           align = c("l", "c", "c"),
           mean = SRrisksCert.AUC$dataPlot$AUC,
           lower = SRrisksCert.AUC$dataPlot$LL,
           upper = SRrisksCert.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisksCert.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(SRrisksCert.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = SRrisksCert.AUC$Performance$AUC)
dev.off()


#### ADNEX with CA125 ####

## AUC
ADNEXcert.AUC <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5Imputed_cert, method.MA = "SJ")

ADNEXcert.AUC$dataPlot <- ADNEXcert.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
ADNEXcert.AUC$Plot <- ADNEXcert.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

ADNEXcert.AUC$Plot[17, "RRauc"] <- "   "
ADNEXcert.AUC$Plot[19, "RRauc"] <- "        (0.84 to 0.98)"
ADNEXcert.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXw - Without uncertain PI.tiff", width = 1000, height = 650)
forestplot(ADNEXcert.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXcert.AUC$dataPlot$AUC,
           lower = ADNEXcert.AUC$dataPlot$LL,
           upper = ADNEXcert.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXcert.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXcert.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXcert.AUC$Performance$AUC)
dev.off()


#### ADNEX without CA125 ####

## AUC
ADNEXwoCert.AUC <- AUC.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, method.MA = "SJ", titleGraph = "LR2: AUC per center")

ADNEXwoCert.AUC$dataPlot <- ADNEXwoCert.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
ADNEXwoCert.AUC$Plot <- ADNEXwoCert.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

ADNEXwoCert.AUC$Plot[17, "RRauc"] <- "   "
ADNEXwoCert.AUC$Plot[19, "RRauc"] <- "        (0.83 to 0.98)"
ADNEXwoCert.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXwo - Without uncertain PI.tiff", width = 1000, height = 650)
forestplot(ADNEXwoCert.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwoCert.AUC$dataPlot$AUC,
           lower = ADNEXwoCert.AUC$dataPlot$LL,
           upper = ADNEXwoCert.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwoCert.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwoCert.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwoCert.AUC$Performance$AUC)
dev.off()


#### Summary plot for AUC ####

NA.forest <- RMIcert.AUC$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, RMIcert.AUC$Performance[1,], LR2cert.AUC$Performance[1,], SRrisksCert.AUC$Performance[1,], ADNEXwoCert.AUC$Performance[1,], ADNEXcert.AUC$Performance[1,])
Summary.AUC$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.AUC.PI <- rbind(NA.forest, RMIcert.AUC$Performance[2,], LR2cert.AUC$Performance[2,], SRrisksCert.AUC$Performance[2,], ADNEXwoCert.AUC$Performance[2,], ADNEXcert.AUC$Performance[2,])
Summary.AUC.PI$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.AUC
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMIcert.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMIcert.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(RMIcert.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2cert.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2cert.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(LR2cert.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisksCert.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisksCert.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(SRrisksCert.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwoCert.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwoCert.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwoCert.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXcert.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXcert.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXcert.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMIcert.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMIcert.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2cert.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2cert.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisksCert.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisksCert.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwoCert.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwoCert.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXcert.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXcert.AUC$Performance$UL[2], 2), nsmall = 2), ")"))
)

tiff("Results Paper 2/With cysts/Paper/Summary AUC - Without uncertain PI.tiff", width = 31, height = 13.75, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()



#### 5.3.2 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Without uncertain.tiff", width = 14, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = Iota5Imputed_cert, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2cert <- RE.ValProb2(p = LR2, y = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, title = "LR2: Calibration curve")


#### Simple Rules Risk ####
SRriskCert <- RE.ValProb2(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXwCert <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = Iota5Imputed_cert, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwoCert <- RE.ValProb2(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = Iota5Ni_cert, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
p.LR2 = seq(min(Iota5Ni_cert$LR2), max(Iota5Ni_cert$LR2), length = 500)
X = cbind(1, Logit(p.LR2))
FE = fixef(LR2cert$FitCalSlope)
OverallCal.LR2 = inv.logit(X[order(p.LR2), ] %*% FE)
p.LR2 = p.LR2[order(p.LR2)]

# Simple Rules Risk
p.SR = seq(min(Iota5Ni_cert$SRrisks), max(Iota5Ni_cert$SRrisks), length = 500)
X = cbind(1, Logit(p.SR))
FE = fixef(SRriskCert$FitCalSlope)
OverallCal.SRrisk = inv.logit(X[order(p.SR), ] %*% FE)
p.SR = p.SR[order(p.SR)]

# ADNEX with CA125
p.ADNEXw = seq(min(Iota5Imputed_cert$pmalw[Iota5Imputed_cert$.imp == 1]), max(Iota5Imputed_cert$pmalw[Iota5Imputed_cert$.imp == 1]), length = 500)
X = cbind(1, Logit(p.ADNEXw))
FE = ADNEXwCert$ResultsBRR
OverallCal.ADNEXw = inv.logit(X[order(p.ADNEXw), ] %*% FE)
p.ADNEXw = p.ADNEXw[order(p.ADNEXw)]

# ADNEX without CA125
p.ADNEXwo = seq(min(Iota5Ni_cert$pmalwo), max(Iota5Ni_cert$pmalwo), length = 500)
X = cbind(1, Logit(p.ADNEXwo))
FE = fixef(ADNEXwoCert$FitCalSlope)
OverallCal.ADNEXwo = inv.logit(X[order(p.ADNEXwo), ] %*% FE)
p.ADNEXwo = p.ADNEXwo[order(p.ADNEXwo)]

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2cert$Performance[1,1], 2), nsmall = 2), " (", format(round(LR2cert$Performance[1,2], 2), nsmall = 2), " to ", format(round(LR2cert$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(LR2cert$Performance[2,1], 2), nsmall = 2), " (", format(round(LR2cert$Performance[2,2], 2), nsmall = 2), " to ", format(round(LR2cert$Performance[2,3], 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRriskCert$Performance[1,1], 2), nsmall = 2), " (", format(round(SRriskCert$Performance[1,2], 2), nsmall = 2), " to ", format(round(SRriskCert$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(SRriskCert$Performance[2,1], 2), nsmall = 2), " (", format(round(SRriskCert$Performance[2,2], 2), nsmall = 2), " to ", format(round(SRriskCert$Performance[2,3], 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXwCert$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwCert$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwCert$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwCert$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwCert$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwCert$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwoCert$Performance[1,1], 2), nsmall = 2), " (", format(round(ADNEXwoCert$Performance[1,2], 2), nsmall = 2), " to ", format(round(ADNEXwoCert$Performance[1,3], 2), nsmall = 2), ")"), paste0(format(round(ADNEXwoCert$Performance[2,1], 2), nsmall = 2), " (", format(round(ADNEXwoCert$Performance[2,2], 2), nsmall = 2), " to ", format(round(ADNEXwoCert$Performance[2,3], 2), nsmall = 2), ")"))


## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Without uncertain.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.16, y = -0.01, table = table, display.colnames= TRUE, cex = 0.7) 
dev.off()



#---------------------------------------------------------#
#### 5.4 Supplementary analysis: including all centers ####
#---------------------------------------------------------#

IOTA5AllCentra <- IOTA5Descrip[IOTA5Descrip$immed.oper == "1" & IOTA5Descrip$TimeSurgery <= 120,] # n = 3369
IOTA5AllCentra <- IOTA5AllCentra[!is.na(IOTA5AllCentra$`Patient ID`),]

## Variable with the categories in short
IOTA5AllCentra$CD_short <- ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Benign: surgery performed", "B1",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Malignant: surgery performed <= 120 days inclusion", "M1",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Malignant: surgery performed after 120 days inclusion, all SA borderline/malignant", "M2",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Uncertain: surgery performed after 120 days inclusion, SA not always borderline/malignant", "U1",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Benign: spontaneous resolution", "B0",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Benign: no surgery, all SA during first 14 months are benign", "B2",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Malignant: no surgery, all SA during first 14 months are malignant", "M3",
                                  ifelse(IOTA5AllCentra$ClinicalDiagnosis == "Uncertain: no surgery, SA during first 14 months is inconsistent", "U2",
                                  ifelse(grepl("Uncertain: Follow-up < 10 months", IOTA5AllCentra$ClinicalDiagnosis), "U3", "U4")))))))))
table(IOTA5AllCentra$CD_short, IOTA5AllCentra$ClinicalDiagnosis)

## Binary CD
IOTA5AllCentra$binaryCD <- sapply(IOTA5AllCentra$ClinicalDiagnosis,
                                  function(x){
                                    if(grepl("Malignant", x))
                                      1
                                    else if(grepl("Benign", x))
                                      0
                                    else
                                      NA
                                  })

## Complete age for all included patients
IOTA5AllCentra[is.na(IOTA5AllCentra$`Patient age`), c("Patient ID", "Patient DOB", "Exam date", "Patient age")]
library(eeptools)
IOTA5AllCentra[IOTA5AllCentra$`Patient ID` == "FLI44", c("Patient age")] = floor(age_calc(as.Date(IOTA5AllCentra$`Patient DOB`[IOTA5AllCentra$`Patient ID` == "FLI44"], format="%Y-%m-%d"), as.Date(IOTA5AllCentra$`Exam date`[IOTA5AllCentra$`Patient ID` == "FLI44"], format="%Y-%m-%d"), units = "years"))
IOTA5AllCentra[IOTA5AllCentra$`Patient ID` == "FLI44", c("Patient ID", "Patient DOB", "Exam date", "Patient age")]

# Look at some characteristics
table(IOTA5AllCentra$ClinicalDiagnosis) # 8 persons with surgery, but missing mass outcome
table(IOTA5AllCentra$`Mass Outcome`) # 269 borderline & 817 invasive malignant & 175 metastatic = 1261
table(IOTA5AllCentra$`FIGO stage`) # 1094 persons with a FIGO stage
IOTA5AllCentra[IOTA5AllCentra$ClinicalDiagnosis == "Surgery, but missing mass outcome", c('Patient ID', 'ClinicalDiagnosis', 'Mass Outcome')]


#### 5.4.1 Multiple imputation ####

#### 5.4.1.1 Variables ####

## Onco- vs non-oncological centers
onco    = c("LPO", "SSW", "PCR", "CIT", "BCH", "OIT", "UDI", "BIT", "RIT", "LBE", "PSP",
            "AGR", "MPO", "CAI", "NCI", "DEP", "LIP", "VAS", "BAI", "KPO")
nononco = c("SIT", "MIT", "GBE", "MSW", "CEG", "FIT", "FLI", "BSP", "TUS", "NUK", "CRI", 
            "TIT", "MCA", "RZT", "MFR", "PFR", "IUK") 

IOTA5AllCentra$oncocenter <- sapply(IOTA5AllCentra$center,
                                    function(x){
                                      if(x %in% onco)
                                        1
                                      else if(x %in% nononco)
                                        0
                                      else
                                        NA
                                    })

## Presumed endometrioma
IOTA5AllCentra$PresumedEndometrioma = sapply(IOTA5AllCentra$`Presumed diagnosis`,
                                             function(x){
                                               if(is.na(x))
                                                 NA
                                               else if(x == "endometrioma")
                                                 1
                                               else
                                                 0
                                             })



## Origin diagnosis
IOTA5AllCentra = ddply(IOTA5AllCentra, .(`Patient ID`),
                       function(x){
                         x$OriginCD =
                           if(is.na(x$`Study Outcome`)){
                             "Clinical diagnosis"
                           }
                         else{
                           if(x$`Study Outcome` == "surgery performed"){
                             if(grepl("Uncertain", x$ClinicalDiagnosis)){
                               "Clinical diagnosis"
                             } else{
                               "Histological analysis"
                             }
                           } 
                           else{
                             "Clinical diagnosis"
                           }
                         }
                         return(x)
                       })
IOTA5AllCentra[IOTA5AllCentra$ClinicalDiagnosis == "Surgery, but missing mass outcome", c('Patient ID', 'ClinicalDiagnosis', 'OriginCD', 'Mass Outcome', 'Study Outcome')]


## Outcome with 5 groups
IOTA5AllCentra = ddply(IOTA5AllCentra, .(`Patient ID`),
                       function(x){
                         x$CD_5groups = 
                           if(!is.na(x$OriginCD)){
                             if(x$OriginCD == "Histological analysis"){
                               if(!is.na(x$`Mass Outcome`)){
                                 if(x$`Mass Outcome` == 'benign'){
                                   "Benign"
                                 } # End if benign mass outcome
                                 else{
                                   if(grepl("invasive", x$`Mass Outcome`)){
                                     if(is.na(x$`FIGO stage`)){
                                       "Invasive FIGO stage missing"
                                     } # End if FIGO is missing
                                     else{
                                       if(x$`FIGO stage` == "I"){
                                         "Stage I invasive"
                                       } # End if FIGO is I
                                       else{
                                         "Stage II-IV invasive"
                                       } # End else (FIGO not I)
                                     } # End else (FIGO not missing)
                                   } # End if invasive mass outcome
                                   else{
                                     if(grepl("borderline", x$`Mass Outcome`)){
                                       "Borderline"
                                     } # End if borderline mass outcome
                                     else{
                                       "Metastatic"
                                     } # End else (no borderline mass outcome)
                                   } # End else (no invasive mass outcome)
                                 } # End else (no benign mass outcome)
                               } # End if mass outcome missing
                               else{
                                 NA
                               } # End else (missing mas outcome)
                             } # End if Histological analysis
                             else{
                               if(is.na(x$binaryCD)){
                                 NA
                               } # End if binaryCD missing
                               else {
                                 if(x$binaryCD == 0){
                                   "Benign"
                                 } # End if binaryCD == 0
                                 else{
                                   NA
                                 } # End else (binaryCD is not 0)
                               } # End else (binaryCD is not missing)
                             } # End else (no histological analysis)
                           } # End if no missings in OriginCD
                         else{
                           NA
                         }# End else (missings in OriginCD)
                         return(x)
                       })
table(IOTA5AllCentra$CD_5groups) # Invasive FIGO stage missing: 15

FIGO <- IOTA5AllCentra[IOTA5AllCentra$CD_5groups == "Invasive FIGO stage missing" & !is.na(IOTA5AllCentra$CD_5groups), c('Patient ID')]
IOTA5AllCentra[IOTA5AllCentra$CD_5groups == "Invasive FIGO stage missing" & !is.na(IOTA5AllCentra$CD_5groups), c('CD_5groups')] = NA
IOTA5AllCentra$CD_5groups = factor(IOTA5AllCentra$CD_5groups, levels = c("Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"))
table(IOTA5AllCentra$ClinicalDiagnosis, IOTA5AllCentra$CD_5groups)

## Binary variable for multilocular
IOTA5AllCentra$multibin = sapply(IOTA5AllCentra$`Tumour type`,
                                 function(x){
                                   if(is.na(x))
                                     NA
                                   else if(grepl("multilocular", x))
                                     1
                                   else
                                     0
                                 })

## Binary variable for solid
IOTA5AllCentra$solidbin = sapply(IOTA5AllCentra$`Tumour type`,
                                 function(x){
                                   if(is.na(x))
                                     NA
                                   else if(grepl("solid", x))
                                     1
                                   else
                                     0
                                 })

## Binary variable for bilateral
IOTA5AllCentra$bilatbin = sapply(IOTA5AllCentra$Bilaterality,
                                 function(x){
                                   if(is.na(x))
                                     NA
                                   else if (x == "bilateral")
                                     1
                                   else
                                     0
                                 })

## Binary variable for Metastases
IOTA5AllCentra$metasbin = sapply(IOTA5AllCentra$Metastases,
                                 function(x){
                                   if(is.na(x))
                                     NA
                                   else if (x == "yes")
                                     1
                                   else
                                     0
                                 })

## Ordinal variable for Subjective assessment
IOTA5AllCentra$CertaintyOrdinal <- factor(IOTA5AllCentra$Certainty,
                                          levels = c("certainly benign", "probably benign", "uncertain", "probably malignant", "certainly malignant"),
                                          ordered = T)

## Variable indicating whether CA125 is missing
IOTA5AllCentra$CA125missing = factor(sapply(IOTA5AllCentra$CA125, function(x) as.numeric(is.na(x))),
                                     levels = 0:1, labels = c("CA125 present", "CA125 missing"))
IOTA5AllCentra[IOTA5AllCentra$CA125missing == "CA125 missing", c("CA125", "CA125missing")]

table(IOTA5AllCentra$CA125missing, IOTA5AllCentra$immed.oper)

## Variable indicating whether CD_5groups is missing
IOTA5AllCentra$CDmissing = factor(sapply(IOTA5AllCentra$CD_5groups, function(x) as.numeric(is.na(x))),
                                  levels = 0:1, labels = c("CD present", "CD missing"))
IOTA5AllCentra[IOTA5AllCentra$CDmissing == "CD missing", c("CD_5groups", "CDmissing")]

## Variable for the centers in analyses
table(IOTA5AllCentra$center, IOTA5AllCentra$binaryCD)
IOTA5AllCentra$centerRE = sapply(IOTA5AllCentra$center,
                                 function(x){
                                   if(x == "AGR")
                                     "Athens (Greece)"
                                   else if(x == "BAI")
                                     "Bari (Italy)"
                                   else if(x == "BCH")
                                     "Beijing (China)"
                                   else if(x == "BIT")
                                     "Bologna (Italy)"
                                   else if(x == "CIT")
                                     "Milan 1 (Italy)"
                                   else if(x == "GBE")
                                     "Genk (Belgium)"
                                   else if(x == "KPO")
                                     "Krakow (Poland)"
                                   else if(x == "LBE")
                                     "Leuven (Belgium)"
                                   else if(x == "LPO")
                                     "Lublin (Poland)"
                                   else if(x == "MPO")
                                     "Katowice (Poland)"
                                   else if(x == "MSW")
                                     "Malmö (Sweden)"
                                   else if(x == "NCI")
                                     "Milan 2 (Italy)"
                                   else if(x == "OIT")
                                     "Monza (Italy)"
                                   else if(x == "PCR")
                                     "Prague (Czech Republic)"
                                   else if(x == "PSP")
                                     "Pamplona (Spain)"
                                   else if(x == "RIT")
                                     "Rome (Italy)"
                                   else if(x == "SIT")
                                     "Cagliari (Italy)"
                                   else if(x == "SSW")
                                     "Stockholm (Sweden)"
                                   else if(x == "TIT")
                                     "Trieste (Italy)"
                                   else if(x == "UDI")
                                     "Udine (Italy)"
                                   else
                                     "Other"
                                 })
table(IOTA5AllCentra$center)
table(IOTA5AllCentra$centerRE)
table(IOTA5AllCentra$center, IOTA5AllCentra$centerRE)

## Variable for locules
IOTA5AllCentra$Nlocules <- ifelse(IOTA5AllCentra$Locules == "> 10", "> 10",
                                  ifelse(as.numeric(IOTA5AllCentra$Locules) == 1, "1",
                                         ifelse(as.numeric(IOTA5AllCentra$Locules) >= 2 & as.numeric(IOTA5AllCentra$Locules) <= 10, "2 - 10", "Other")))

## Binary variable for Personal history of ovarian cancer
IOTA5AllCentra$Personal_OvCAbin <- ifelse(IOTA5AllCentra$Personal_OvCA == "no", "No", "Yes")

## Log transformation of 'maximum diameter of lesion'
IOTA5AllCentra$Llesdmax <- log2(IOTA5AllCentra$`Lesion largest diameter`)

## Quadratic term for 'Proportion of solid tissue'
IOTA5AllCentra$Qpropsol <- IOTA5AllCentra$propsol^2

#### 5.4.1.2 Variables needed ####
IOTA5AllCentra$outcome1Correct = as.factor(IOTA5AllCentra$outcome1)
IOTA5AllCentra$binaryCDcorrect = as.factor(IOTA5AllCentra$binaryCD)

# Check for missingness in the variables
IOTA5AllCentra[is.na(IOTA5AllCentra$PresumedEndometrioma), c("Patient ID", "PresumedEndometrioma")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$Certainty), c("Patient ID", "Certainty")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$`Patient age`), c("Patient ID", "Patient age")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$oncocenter), c("Patient ID", "oncocenter")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$`Lesion largest diameter`), c("Patient ID", "Lesion largest diameter")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$propsol), c("Patient ID", "propsol")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$Nlocules), c("Patient ID", "Nlocules")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$loc10), c("Patient ID", "loc10")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$papnr), c("Patient ID", "papnr")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$shadowsbin), c("Patient ID", "shadowsbin")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$ascitesbin), c("Patient ID", "ascitesbin")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$metasbin), c("Patient ID", "metasbin")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$bilatbin), c("Patient ID", "bilatbin")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$Pelvic_pain), c("Patient ID", "Pelvic_pain")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$Personal_OvCA), c("Patient ID", "Personal_OvCA")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$`Papillary height`), c("Patient ID", "Papillary height")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$papflow), c("Patient ID", "papflow")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$ColorOrdinal), c("Patient ID", "ColorOrdinal")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$Echogenicity), c("Patient ID", "Echogenicity")] # No missings
IOTA5AllCentra[is.na(IOTA5AllCentra$CD_5groups), c("Patient ID", "CD_5groups", "ClinicalDiagnosis", "OriginCD")] # 8 Missings (Surgery, but missing mass outcome)


#### 5.4.1.3 Preparation imputation ####

IOTA5AllCentra$ll_CA125 <- log(log(IOTA5AllCentra$CA125 + 1))

PredMatr


#### 5.4.1.4 Imputation ####

iota5Imp <- IOTA5AllCentra[, AllVars]

colnames(iota5Imp) = gsub(" ","\\.",colnames(iota5Imp))
rownames(PredMatr)   = colnames(PredMatr) = gsub(" ","\\.",colnames(PredMatr))

Init = mice(iota5Imp, m = 1, maxit = 0, predictorMatrix = PredMatr)
Init$loggedEvents 
Init$method
for(i in Init$loggedEvents$out)
  iota5Imp[, i] %<>% as.factor

Method = sapply(colnames(iota5Imp),
                function(x){
                  if(x == "ll_CA125")
                    "pmm"
                  else if(x == "CD_5groups")
                    "polyreg"
                  else
                    ""
                })
Method

iota5MI <- mice(iota5Imp, m = 100, me = Method,
                predictorMatrix = PredMatr, maxit = 50,
                seed = 1213, vis = "monotone", ridge = 1e-3)

iota5MI$loggedEvents


#### 5.4.1.5 Check convergence ####

## Convergence
tiff("Results Paper 2/With cysts/Convergence - All centers.tiff")
plot(iota5MI)
dev.off()

## Density plots
densityplot(iota5MI, ~ ll_CA125| .imp)

tiff("Results Paper 2/With cysts/Density CA125 - All centers.tiff")
densityplot(iota5MI, ~ ll_CA125)
dev.off()

densityplot(iota5MI, ~ CD_5groups| .imp)

tiff("Results Paper 2/With cysts/Density CD - All centers.tiff")
densityplot(iota5MI, ~ CD_5groups)
dev.off()

## Boxplot CA125
Iota5Imputed                = mice::complete(iota5MI, "long")
Iota5Ni                     = complete(iota5MI)
Iota5Imputed$CA125 = exp(exp(Iota5Imputed$ll_CA125)) - 1
Iota5Ni$CA125 = exp(exp(Iota5Ni$ll_CA125)) - 1

tiff("Results Paper 2/With cysts/Boxplot CA125 imputed - All centers.tiff", width = 1200, height = 650)
boxplot(ll_CA125 ~ .imp, Iota5Imputed[Iota5Imputed$CA125missing == "CA125 missing",])
dev.off()

## Determine most commonly imputed value as outcome for patients with FIGO stage missing
ImputeCD = dlply(Iota5Imputed[Iota5Imputed$CDmissing == "CD missing" & Iota5Imputed$Patient.ID %in% FIGO,], .(Patient.ID),
                 function(x){
                   Outcome = table(x$CD_5groups)
                   #names(Outcome)[which.max(Outcome)]
                 })
ImputeCD
table(IOTA5AllCentra$ClinicalDiagnosis)
IOTA5AllCentra[IOTA5AllCentra$ClinicalDiagnosis == "Surgery, but missing mass outcome", c('Patient ID')]
IOTA5AllCentra[IOTA5AllCentra$`Patient ID` %in% FIGO, c('Mass Outcome')]

ImputeCD$BAI64 # Stage II-IV invasive
ImputeCD$BCH52 # Stage II-IV invasive
ImputeCD$BCH89 # Stage II-IV invasive
ImputeCD$BIT710 # Borderline
ImputeCD$BIT971 # Borderline
ImputeCD$CAI128 # Borderline
ImputeCD$CAI140 # Stage I invasive
ImputeCD$KPO24 # Stage II-IV invasive
ImputeCD$PCR176 # Stage II-IV invasive
ImputeCD$PCR200 # Stage II-IV invasive
ImputeCD$PCR217 # Stage II-IV invasive
ImputeCD$PCR32 # Stage II-IV invasive
ImputeCD$PCR52 # Stage II-IV invasive
ImputeCD$PFR3 # Stage II-IV invasive
ImputeCD$PFR32 # Borderline

Iota5Imputed[Iota5Imputed$Patient.ID == "BAI64", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "BCH52", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "BCH89", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "BIT710", c("CD_5groups")] <- "Borderline"
Iota5Imputed[Iota5Imputed$Patient.ID == "BIT971", c("CD_5groups")] <- "Borderline"
Iota5Imputed[Iota5Imputed$Patient.ID == "CAI128", c("CD_5groups")] <- "Borderline"
Iota5Imputed[Iota5Imputed$Patient.ID == "CAI140", c("CD_5groups")] <- "Stage I invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "KPO24", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PCR176", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PCR200", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PCR217", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PCR32", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PCR52", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PFR3", c("CD_5groups")] <- "Stage II-IV invasive"
Iota5Imputed[Iota5Imputed$Patient.ID == "PFR32", c("CD_5groups")] <- "Borderline"


## adapt binaryCDcorrect according to CD_5groups
Iota5Imputed$binaryCDcorrect <- sapply(Iota5Imputed$CD_5groups,
                                       function(x){
                                         if(x == "Benign")
                                           0
                                         else
                                           1
                                       })
table(Iota5Imputed$binaryCDcorrect)

## Observed and imputed outcome
Observed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD present"])
round(Observed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD present",]) * 100, 2)
Imputed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD missing"])
round(Imputed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD missing",]) * 100, 2)


## Save dataset
Iota5AllI  = Iota5Imputed
Iota5AllNI   = Iota5Imputed[Iota5Imputed$.imp == 1, ]


#### 5.4.2 Model predictions ####
source(paste0(Scripts, "/FunctionsAllPredModels2.R"))

#### RMI ####
Iota5AllI = RMI(postmeno = Postmenopausal2, multilocular = multibin, solid = solidbin, metastases = metasbin, bilateral = bilatbin, ascites = ascitesbin, CA125 = CA125, patientid = Patient.ID, data = Iota5AllI, imputed = T, .imp = .imp)
Iota5AllI$RMI

#### LR2 ####
Iota5AllNI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5AllNI)
Iota5AllNI$LR2
Iota5AllI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5AllI)
Iota5AllI$LR2

#### Simple Rules and Simple Rules Risk ####
Iota5AllNI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5AllNI)
Iota5AllNI$SRrisks
Iota5AllI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5AllI)
Iota5AllI$SRrisks

#### ADNEX ####
## With CA125
Iota5AllI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = T, woca = F, data = Iota5AllI)
Iota5AllI$pmalw

## Without CA125
Iota5AllNI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5AllNI)
Iota5AllNI$pmalwo
Iota5AllI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5AllI)
Iota5AllI$pmalwo

#### 5.4.3 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMIall.AUC <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5AllI, method.MA = "SJ")

RMIall.AUC$dataPlot <- RMIall.AUC$dataPlot[c(1:8, 10:23, 9, 24:26),]
RMIall.AUC$Plot <- RMIall.AUC$Plot[c(1:8, 10:23, 9, 24:26),]

RMIall.AUC$Plot[24, "RRauc"] <- "   "
RMIall.AUC$Plot[26, "RRauc"] <- "        (0.73 to 0.94)"
RMIall.AUC$Plot[26, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest RMI - All centers PI.tiff", width = 1000, height = 650)
forestplot(RMIall.AUC$Plot,
           align = c("l", "c", "c"),
           mean = RMIall.AUC$dataPlot$AUC,
           lower = RMIall.AUC$dataPlot$LL,
           upper = RMIall.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMIall.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(RMIall.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = RMIall.AUC$Performance$AUC)
dev.off()


#### LR2 ####

## AUC
LR2all.AUC <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5AllI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

LR2all.AUC$dataPlot <- LR2all.AUC$dataPlot[c(1:8, 10:23, 9, 24:26),]
LR2all.AUC$Plot <- LR2all.AUC$Plot[c(1:8, 10:23, 9, 24:26),]

LR2all.AUC$Plot[24, "RRauc"] <- "   "
LR2all.AUC$Plot[26, "RRauc"] <- "        (0.77 to 0.95)"
LR2all.AUC$Plot[26, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest LR2 - All centers PI.tiff", width = 1000, height = 650)
forestplot(LR2all.AUC$Plot,
           align = c("l", "c", "c"),
           mean = LR2all.AUC$dataPlot$AUC,
           lower = LR2all.AUC$dataPlot$LL,
           upper = LR2all.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2all.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(LR2all.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = LR2all.AUC$Performance$AUC)
dev.off()


#### Simple Rules Risk ####

## AUC
SRrisksAll.AUC <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5AllI, method.MA = "SJ", titleGraph = "Simple Rules risk: AUC per center")

SRrisksAll.AUC$dataPlot <- SRrisksAll.AUC$dataPlot[c(1:8, 10:23, 9, 24:26),]
SRrisksAll.AUC$Plot <- SRrisksAll.AUC$Plot[c(1:8, 10:23, 9, 24:26),]

SRrisksAll.AUC$Plot[24, "RRauc"] <- "   "
SRrisksAll.AUC$Plot[26, "RRauc"] <- "        (0.79 to 0.96)"
SRrisksAll.AUC$Plot[26, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest SRrisks - All centers PI.tiff", width = 1000, height = 650)
forestplot(SRrisksAll.AUC$Plot,
           align = c("l", "c", "c"),
           mean = SRrisksAll.AUC$dataPlot$AUC,
           lower = SRrisksAll.AUC$dataPlot$LL,
           upper = SRrisksAll.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisksAll.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(SRrisksAll.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = SRrisksAll.AUC$Performance$AUC)
dev.off()


#### ADNEX with CA125 ####

## AUC
ADNEXwAll.AUC <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5AllI, method.MA = "SJ")

ADNEXwAll.AUC$dataPlot <- ADNEXwAll.AUC$dataPlot[c(1:8, 10:23, 9, 24:26),]
ADNEXwAll.AUC$Plot <- ADNEXwAll.AUC$Plot[c(1:8, 10:23, 9, 24:26),]

ADNEXwAll.AUC$Plot[24, "RRauc"] <- "   "
ADNEXwAll.AUC$Plot[26, "RRauc"] <- "        (0.79 to 0.97)"
ADNEXwAll.AUC$Plot[26, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXw - All centers PI.tiff", width = 1000, height = 650)
forestplot(ADNEXwAll.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwAll.AUC$dataPlot$AUC,
           lower = ADNEXwAll.AUC$dataPlot$LL,
           upper = ADNEXwAll.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwAll.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwAll.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwAll.AUC$Performance$AUC)
dev.off()


#### ADNEX without CA125 ####

## AUC
ADNEXwoAll.AUC <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5AllI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

ADNEXwoAll.AUC$dataPlot <- ADNEXwoAll.AUC$dataPlot[c(1:8, 10:23, 9, 24:26),]
ADNEXwoAll.AUC$Plot <- ADNEXwoAll.AUC$Plot[c(1:8, 10:23, 9, 24:26),]

ADNEXwoAll.AUC$Plot[24, "RRauc"] <- "   "
ADNEXwoAll.AUC$Plot[26, "RRauc"] <- "        (0.77 to 0.97)"
ADNEXwoAll.AUC$Plot[26, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXwo - All centers PI.tiff", width = 1000, height = 650)
forestplot(ADNEXwoAll.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwoAll.AUC$dataPlot$AUC,
           lower = ADNEXwoAll.AUC$dataPlot$LL,
           upper = ADNEXwoAll.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwoAll.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwoAll.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwoAll.AUC$Performance$AUC)
dev.off()


#### Summary plot for AUC ####

NA.forest <- RMIall.AUC$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, RMIall.AUC$Performance[1,], LR2all.AUC$Performance[1,], SRrisksAll.AUC$Performance[1,], ADNEXwoAll.AUC$Performance[1,], ADNEXwAll.AUC$Performance[1,])
Summary.AUC$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC.PI <- rbind(NA.forest, RMIall.AUC$Performance[2,], LR2all.AUC$Performance[2,], SRrisksAll.AUC$Performance[2,], ADNEXwoAll.AUC$Performance[2,], ADNEXwAll.AUC$Performance[2,])
Summary.AUC.PI$Model <- c('', 'RMI', 'LR2', 'Simple Rules risk', 'ADNEX w/o CA125', 'ADNEX w/ CA125')

Summary.AUC
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMIall.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMIall.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(RMIall.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2all.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2all.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(LR2all.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisksAll.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisksAll.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(SRrisksAll.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwoAll.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwoAll.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwoAll.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwAll.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwAll.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwAll.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMIall.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMIall.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2all.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2all.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisksAll.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisksAll.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwoAll.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwoAll.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwAll.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwAll.AUC$Performance$UL[2], 2), nsmall = 2), ")"))
)

tiff("Results Paper 2/With cysts/Paper/Summary AUC - All centers PI.tiff", width = 31, height = 13.75, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface = "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.80, 0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.80, 1))
dev.off()


#### 5.4.4 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - All centers.tiff", width = 14, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = Iota5AllI, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()

#### LR2 ####
LR2all <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = Iota5AllI, imputation.id = .imp, patientid = Patient.ID)


#### Simple Rules Risk ####
SRriskAll <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = Iota5AllI, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXwAll <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = Iota5AllI, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwoAll <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = Iota5AllI, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2all$Plot$y
p.LR2 = LR2all$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRriskAll$Plot$y
p.SR = SRriskAll$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXwAll$Plot$y
p.ADNEXw = ADNEXwAll$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwoAll$Plot$y
p.ADNEXwo = ADNEXwoAll$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2all$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2all$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2all$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2all$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2all$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2all$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRriskAll$SumInt$Est, 2), nsmall = 2), " (", format(round(SRriskAll$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRriskAll$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRriskAll$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRriskAll$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRriskAll$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXwAll$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwAll$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwAll$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwAll$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwAll$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwAll$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwoAll$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwoAll$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwoAll$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwoAll$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwoAll$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwoAll$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - All centers.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1) 
addtable2plot(x = 0.16, y = -0.01, table = table, display.colnames= TRUE, cex = 0.7) 
dev.off()



#----------------------------------------------------------------------#
#### 5.5 Sensitivity analysis: imputation uncertain and FU outcomes ####
#----------------------------------------------------------------------#

Iota5Sens <- iota5ValAnalyses

Iota5Sens[grepl("Uncertain", Iota5Sens$ClinicalDiagnosis), ("ClinicalDiagnosis")] = "Clinical Diagnosis"
Iota5Sens[Iota5Sens$CD_short == "M2", ("ClinicalDiagnosis")] = "Clinical Diagnosis"
Iota5Sens[Iota5Sens$CD_short == "B2", ("ClinicalDiagnosis")] = "Clinical Diagnosis"
Iota5Sens[Iota5Sens$CD_short == "M3", ("ClinicalDiagnosis")] = "Clinical Diagnosis"

Iota5Sens[Iota5Sens$ClinicalDiagnosis == "Clinical Diagnosis", ("CD_short")] = "CD"

table(Iota5Sens$ClinicalDiagnosis)
table(Iota5Sens$CD_short) # 465 B0, 2065 B1, 956 M1 and 1419 are based on FU data (CD)

## Binary CD
Iota5Sens$binaryCD <- sapply(Iota5Sens$ClinicalDiagnosis,
                             function(x){
                               if(grepl("Malignant", x))
                                 1
                               else if(grepl("Benign", x))
                                 0
                               else
                                 NA
                             })

table(Iota5Sens$binaryCD)

## Complete age for all included patients
Iota5Sens[is.na(Iota5Sens$`Patient age`), c("Patient ID", "Patient DOB", "Exam date", "Patient age")]

#### 5.5.1 Multiple imputation ####

#### 5.5.1.1 Variables ####

## Presumed endometrioma
Iota5Sens$PresumedEndometrioma = sapply(Iota5Sens$`Presumed diagnosis`,
                                        function(x){
                                          if(is.na(x))
                                            NA
                                          else if(x == "endometrioma")
                                            1
                                          else
                                            0
                                        })



## Origin diagnosis
Iota5Sens = ddply(Iota5Sens, .(`Patient ID`),
                  function(x){
                    x$OriginCD =
                      if(is.na(x$`Study Outcome`)){
                        "Clinical diagnosis"
                      }
                    else{
                      if(x$`Study Outcome` == "surgery performed"){
                        if(grepl("Clinical", x$ClinicalDiagnosis)){
                          "Clinical diagnosis"
                        } else{
                          "Histological analysis"
                        }
                      } 
                      else{
                        "Clinical diagnosis"
                      }
                    }
                    return(x)
                  })
table(Iota5Sens$OriginCD, Iota5Sens$CD_short)

## Outcome with 5 groups
Iota5Sens = ddply(Iota5Sens, .(`Patient ID`),
                  function(x){
                    x$CD_5groups = 
                      if(!is.na(x$OriginCD)){
                        if(x$OriginCD == "Histological analysis"){
                          if(!is.na(x$`Mass Outcome`)){
                            if(x$`Mass Outcome` == 'benign'){
                              "Benign"
                            } # End if benign mass outcome
                            else{
                              if(grepl("invasive", x$`Mass Outcome`)){
                                if(is.na(x$`FIGO stage`)){
                                  "Invasive FIGO stage missing"
                                } # End if FIGO is missing
                                else{
                                  if(x$`FIGO stage` == "I"){
                                    "Stage I invasive"
                                  } # End if FIGO is I
                                  else{
                                    "Stage II-IV invasive"
                                  } # End else (FIGO not I)
                                } # End else (FIGO not missing)
                              } # End if invasive mass outcome
                              else{
                                if(grepl("borderline", x$`Mass Outcome`)){
                                  "Borderline"
                                } # End if borderline mass outcome
                                else{
                                  "Metastatic"
                                } # End else (no borderline mass outcome)
                              } # End else (no invasive mass outcome)
                            } # End else (no benign mass outcome)
                          } # End if mass outcome missing
                          else{
                            NA
                          } # End else (missing mas outcome)
                        } # End if Histological analysis
                        else{
                          if(is.na(x$binaryCD)){
                            NA
                          } # End if binaryCD missing
                          else {
                            if(x$binaryCD == 0){
                              "Benign"
                            } # End if binaryCD == 0
                            else{
                              NA
                            } # End else (binaryCD is not 0)
                          } # End else (binaryCD is not missing)
                        } # End else (no histological analysis)
                      } # End if no missings in OriginCD
                    else{
                      NA
                    }# End else (missings in OriginCD)
                    return(x)
                  })
table(Iota5Sens$CD_5groups, Iota5Sens$CD_short)
Iota5Sens[Iota5Sens$CD_5groups == "Invasive FIGO stage missing" & !is.na(Iota5Sens$CD_5groups), c("Patient ID", "CD_5groups", "CD_short")]
Iota5Sens$CD_5groups = factor(Iota5Sens$CD_5groups, levels = c("Benign", "Borderline", "Stage I invasive", "Stage II-IV invasive", "Metastatic"))
table(Iota5Sens$ClinicalDiagnosis, Iota5Sens$CD_5groups)

## Binary variable for multilocular
Iota5Sens$multibin = sapply(Iota5Sens$`Tumour type`,
                            function(x){
                              if(is.na(x))
                                NA
                              else if(grepl("multilocular", x))
                                1
                              else
                                0
                            })

## Binary variable for solid
Iota5Sens$solidbin = sapply(Iota5Sens$`Tumour type`,
                            function(x){
                              if(is.na(x))
                                NA
                              else if(grepl("solid", x))
                                1
                              else
                                0
                            })

## Binary variable for bilateral
Iota5Sens$bilatbin = sapply(Iota5Sens$Bilaterality,
                            function(x){
                              if(is.na(x))
                                NA
                              else if (x == "bilateral")
                                1
                              else
                                0
                            })

## Binary variable for Metastases
Iota5Sens$metasbin = sapply(Iota5Sens$Metastases,
                            function(x){
                              if(is.na(x))
                                NA
                              else if (x == "yes")
                                1
                              else
                                0
                            })

## Ordinal variable for Subjective assessment
Iota5Sens$CertaintyOrdinal <- factor(Iota5Sens$Certainty,
                                     levels = c("certainly benign", "probably benign", "uncertain", "probably malignant", "certainly malignant"),
                                     ordered = T)

## Variable indicating whether CA125 is missing
Iota5Sens$CA125missing = factor(sapply(Iota5Sens$CA125, function(x) as.numeric(is.na(x))),
                                levels = 0:1, labels = c("CA125 present", "CA125 missing"))
Iota5Sens[Iota5Sens$CA125missing == "CA125 missing", c("CA125", "CA125missing")]

table(Iota5Sens$CA125missing, Iota5Sens$immed.oper)

## Variable indicating whether CD_5groups is missing
Iota5Sens$CDmissing = factor(sapply(Iota5Sens$CD_5groups, function(x) as.numeric(is.na(x))),
                             levels = 0:1, labels = c("CD present", "CD missing"))
Iota5Sens[Iota5Sens$CDmissing == "CD missing", c("CD_5groups", "CDmissing")]

## Variable for the centers in analyses
table(Iota5Sens$center, Iota5Sens$binaryCD)
Iota5Sens$centerRE = sapply(Iota5Sens$center,
                            function(x){
                              if(x == "AGR")
                                "Athens (Greece)"
                              else if(x == "BAI")
                                "Bari (Italy)"
                              else if(x == "BCH")
                                "Beijing (China)"
                              else if(x == "BIT")
                                "Bologna (Italy)"
                              else if(x == "CIT")
                                "Milan 1 (Italy)"
                              else if(x == "GBE")
                                "Genk (Belgium)"
                              else if(x == "KPO")
                                "Krakow (Poland)"
                              else if(x == "LBE")
                                "Leuven (Belgium)"
                              else if(x == "LPO")
                                "Lublin (Poland)"
                              else if(x == "MPO")
                                "Katowice (Poland)"
                              else if(x == "MSW")
                                "Malmö (Sweden)"
                              else if(x == "NCI")
                                "Milan 2 (Italy)"
                              else if(x == "OIT")
                                "Monza (Italy)"
                              else if(x == "PCR")
                                "Prague (Czech Republic)"
                              else if(x == "PSP")
                                "Pamplona (Spain)"
                              else if(x == "RIT")
                                "Rome (Italy)"
                              else if(x == "SIT")
                                "Cagliari (Italy)"
                              else if(x == "SSW")
                                "Stockholm (Sweden)"
                              else if(x == "TIT")
                                "Trieste (Italy)"
                              else if(x == "UDI")
                                "Udine (Italy)"
                              else
                                "Other"
                            })
table(Iota5Sens$center)
table(Iota5Sens$centerRE)
table(Iota5Sens$center, Iota5Sens$centerRE)

## Variable for locules
Iota5Sens$Nlocules <- ifelse(Iota5Sens$Locules == "> 10", "> 10",
                             ifelse(as.numeric(Iota5Sens$Locules) == 1, "1",
                                    ifelse(as.numeric(Iota5Sens$Locules) >= 2 & as.numeric(Iota5Sens$Locules) <= 10, "2 - 10", "Other")))

## Binary variable for Personal history of ovarian cancer
Iota5Sens$Personal_OvCAbin <- ifelse(Iota5Sens$Personal_OvCA == "no", "No", "Yes")

## Log transformation of 'maximum diameter of lesion'
Iota5Sens$Llesdmax <- log2(Iota5Sens$`Lesion largest diameter`)

## Quadratic term for 'Proportion of solid tissue'
Iota5Sens$Qpropsol <- Iota5Sens$propsol^2

#### 5.5.1.2 Variables needed ####
Iota5Sens$outcome1Correct = as.factor(Iota5Sens$outcome1)
Iota5Sens$binaryCDcorrect = as.factor(Iota5Sens$binaryCD)


# Check for missingness in the variables
Iota5Sens[is.na(Iota5Sens$PresumedEndometrioma), c("Patient ID", "PresumedEndometrioma")] # No missings
Iota5Sens[is.na(Iota5Sens$Certainty), c("Patient ID", "Certainty")] # No missings
Iota5Sens[is.na(Iota5Sens$`Patient age`), c("Patient ID", "Patient age")] # No missings
Iota5Sens[is.na(Iota5Sens$oncocenter), c("Patient ID", "oncocenter")] # No missings
Iota5Sens[is.na(Iota5Sens$`Lesion largest diameter`), c("Patient ID", "Lesion largest diameter")] # No missings
Iota5Sens[is.na(Iota5Sens$propsol), c("Patient ID", "propsol")] # No missings
Iota5Sens[is.na(Iota5Sens$Nlocules), c("Patient ID", "Nlocules")] # No missings
Iota5Sens[is.na(Iota5Sens$loc10), c("Patient ID", "loc10")] # No missings
Iota5Sens[is.na(Iota5Sens$papnr), c("Patient ID", "papnr")] # No missings
Iota5Sens[is.na(Iota5Sens$shadowsbin), c("Patient ID", "shadowsbin")] # No missings
Iota5Sens[is.na(Iota5Sens$ascitesbin), c("Patient ID", "ascitesbin")] # No missings
Iota5Sens[is.na(Iota5Sens$metasbin), c("Patient ID", "metasbin")] # No missings
Iota5Sens[is.na(Iota5Sens$bilatbin), c("Patient ID", "bilatbin")] # No missings
Iota5Sens[is.na(Iota5Sens$Pelvic_pain), c("Patient ID", "Pelvic_pain")] # No missings
Iota5Sens[is.na(Iota5Sens$Personal_OvCA), c("Patient ID", "Personal_OvCA")] # No missings
Iota5Sens[is.na(Iota5Sens$`Papillary height`), c("Patient ID", "Papillary height")] # No missings
Iota5Sens[is.na(Iota5Sens$papflow), c("Patient ID", "papflow")] # No missings
Iota5Sens[is.na(Iota5Sens$ColorOrdinal), c("Patient ID", "ColorOrdinal")] # No missings
Iota5Sens[is.na(Iota5Sens$Echogenicity), c("Patient ID", "Echogenicity")] # No missings
Iota5Sens[is.na(Iota5Sens$CD_5groups), c("Patient ID", "CD_5groups", "ClinicalDiagnosis", "OriginCD")] # 1415 Missings


#### 5.5.1.3 Preparation imputation ####

Iota5Sens$ll_CA125 <- log(log(Iota5Sens$CA125 + 1))

PredMatr


#### 5.5.1.4 Imputation ####

iota5Imp <- Iota5Sens[, AllVars]

colnames(iota5Imp) = gsub(" ","\\.",colnames(iota5Imp))
rownames(PredMatr)   = colnames(PredMatr) = gsub(" ","\\.",colnames(PredMatr))

Init = mice(iota5Imp, m = 1, maxit = 0, predictorMatrix = PredMatr)
Init$loggedEvents 
Init$method
for(i in Init$loggedEvents$out)
  iota5Imp[, i] %<>% as.factor

Method = sapply(colnames(iota5Imp),
                function(x){
                  if(x == "ll_CA125")
                    "pmm"
                  else if(x == "CD_5groups")
                    "polyreg"
                  else
                    ""
                })
Method

iota5MI <- mice(iota5Imp, m = 100, me = Method,
                predictorMatrix = PredMatr, maxit = 50,
                seed = 1213, vis = "monotone", ridge = 1e-3)

iota5MI$loggedEvents


#### 5.5.1.5 Check convergence ####

## Convergence
tiff("Results Paper 2/With cysts/Convergence - Sens CD.tiff")
plot(iota5MI)
dev.off()

## Density plots
densityplot(iota5MI, ~ ll_CA125| .imp)

tiff("Results Paper 2/With cysts/Density CA125 - Sens CD.tiff")
densityplot(iota5MI, ~ ll_CA125)
dev.off()

densityplot(iota5MI, ~ CD_5groups| .imp)

tiff("Results Paper 2/With cysts/Density CD - Sens CD.tiff")
densityplot(iota5MI, ~ CD_5groups)
dev.off()

## Boxplot CA125
Iota5Imputed                = mice::complete(iota5MI, "long")
Iota5Ni                     = complete(iota5MI)
Iota5Imputed$CA125 = exp(exp(Iota5Imputed$ll_CA125)) - 1
Iota5Ni$CA125 = exp(exp(Iota5Ni$ll_CA125)) - 1

tiff("Results Paper 2/With cysts/Boxplot CA125 imputed - Sens CD.tiff", width = 1200, height = 650)
boxplot(ll_CA125 ~ .imp, Iota5Imputed[Iota5Imputed$CA125missing == "CA125 missing",])
dev.off()


## adapt binaryCDcorrect according to CD_5groups
Iota5Imputed$binaryCDcorrect <- sapply(Iota5Imputed$CD_5groups,
                                       function(x){
                                         if(x == "Benign")
                                           0
                                         else
                                           1
                                       })
table(Iota5Imputed$binaryCDcorrect)

## Observed and imputed outcome
Observed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD present"])
round(Observed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD present",]) * 100, 2)
Imputed <- table(Iota5Imputed$CD_5groups[Iota5Imputed$CDmissing == "CD missing"])
round(Imputed/nrow(Iota5Imputed[Iota5Imputed$CDmissing == "CD missing",]) * 100, 2)

## Save dataset
Iota5SensI  = Iota5Imputed
Iota5SensNI   = Iota5Imputed[Iota5Imputed$.imp == 1, ]


#### 5.5.2 Model predictions ####

source(paste0(Scripts, "/FunctionsAllPredModels2.R"))

#### RMI ####
Iota5SensI = RMI(postmeno = Postmenopausal2, multilocular = multibin, solid = solidbin, metastases = metasbin, bilateral = bilatbin, ascites = ascitesbin, CA125 = CA125, patientid = Patient.ID, data = Iota5SensI, imputed = T, .imp = .imp)
Iota5SensI$RMI

#### LR2 ####
Iota5SensNI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5SensNI)
Iota5SensNI$LR2
Iota5SensI = LR2(age = Patient.age, ascites = ascitesbin, papflow = papflow, soldmax = soldmax, irregular = irregbin, shadows = shadowsbin, data = Iota5SensI)
Iota5SensI$LR2

#### Simple Rules and Simple Rules Risk ####
Iota5SensNI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5SensNI)
Iota5SensNI$SRrisks
Iota5SensI = SimpleRules(TumourType = Tumour.type, lesdmax = Lesion.largest.diameter, soldmax = soldmax, shadows = shadowsbin, colorscore = ColorOrdinal, irregular = irregbin, ascites = ascitesbin, papnr = papnr, oncocenter = oncocenter, data = Iota5SensI)
Iota5SensI$SRrisks

#### ADNEX ####
## With CA125
Iota5SensI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = T, woca = F, data = Iota5SensI)
Iota5SensI$pmalw

## Without CA125
Iota5SensNI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5SensNI)
Iota5SensNI$pmalwo
Iota5SensI = ADNEX(Age = Patient.age, Oncocenter = oncocenter, lesdmax = Lesion.largest.diameter, soldmaxorig = soldmax, loc10 = loc10, Papnr = papnr, Shadows = shadowsbin, Ascites = ascitesbin, ca125 = CA125, wica = F, woca = T, data = Iota5SensI)
Iota5SensI$pmalwo


#### 5.5.3 Discrimination between benign and malignant tumors ####

#### RMI ####

## AUC
RMIsens.AUC <- AUCimp.IOTA(pred = RMI, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5SensI, method.MA = "SJ")

RMIsens.AUC$dataPlot <- RMIsens.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
RMIsens.AUC$Plot <- RMIsens.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

RMIsens.AUC$Plot[17, "RRauc"] <- "   "
RMIsens.AUC$Plot[19, "RRauc"] <- "        (0.73 to 0.96)"
RMIsens.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest RMI - Sens FU PI.tiff", width = 1000, height = 650)
forestplot(RMIsens.AUC$Plot,
           align = c("l", "c", "c"),
           mean = RMIsens.AUC$dataPlot$AUC,
           lower = RMIsens.AUC$dataPlot$LL,
           upper = RMIsens.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(RMIsens.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(RMIsens.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = RMIsens.AUC$Performance$AUC)
dev.off()


#### LR2 ####

## AUC
LR2sens.AUC <- AUCimp.IOTA(pred = LR2, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5SensI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

LR2sens.AUC$dataPlot <- LR2sens.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
LR2sens.AUC$Plot <- LR2sens.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

LR2sens.AUC$Plot[17, "RRauc"] <- "   "
LR2sens.AUC$Plot[19, "RRauc"] <- "        (0.81 to 0.96)"
LR2sens.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest LR2 - Sens FU PI.tiff", width = 1000, height = 650)
forestplot(LR2sens.AUC$Plot,
           align = c("l", "c", "c"),
           mean = LR2sens.AUC$dataPlot$AUC,
           lower = LR2sens.AUC$dataPlot$LL,
           upper = LR2sens.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(LR2sens.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(LR2sens.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = LR2sens.AUC$Performance$AUC)
dev.off()


#### Simple Rules Risk ####

## AUC
SRrisksSens.AUC <- AUCimp.IOTA(pred = SRrisks, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5SensI, method.MA = "SJ", titleGraph = "Simple Rules risk: AUC per center")

SRrisksSens.AUC$dataPlot <- SRrisksSens.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
SRrisksSens.AUC$Plot <- SRrisksSens.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

SRrisksSens.AUC$Plot[17, "RRauc"] <- "   "
SRrisksSens.AUC$Plot[19, "RRauc"] <- "        (0.82 to 0.97)"
SRrisksSens.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest SRrisks - Sens FU PI.tiff", width = 1000, height = 650)
forestplot(SRrisksSens.AUC$Plot,
           align = c("l", "c", "c"),
           mean = SRrisksSens.AUC$dataPlot$AUC,
           lower = SRrisksSens.AUC$dataPlot$LL,
           upper = SRrisksSens.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(SRrisksSens.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(SRrisksSens.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = SRrisksSens.AUC$Performance$AUC)
dev.off()


#### ADNEX with CA125 ####

## AUC
ADNEXwSens.AUC <- AUCimp.IOTA(pred = pmalw, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5SensI, method.MA = "SJ")

ADNEXwSens.AUC$dataPlot <- ADNEXwSens.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
ADNEXwSens.AUC$Plot <- ADNEXwSens.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

ADNEXwSens.AUC$Plot[17, "RRauc"] <- "   "
ADNEXwSens.AUC$Plot[19, "RRauc"] <- "        (0.82 to 0.98)"
ADNEXwSens.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXw - Sens FU PI.tiff", width = 1000, height = 650)
forestplot(ADNEXwSens.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwSens.AUC$dataPlot$AUC,
           lower = ADNEXwSens.AUC$dataPlot$LL,
           upper = ADNEXwSens.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwSens.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwSens.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwSens.AUC$Performance$AUC)
dev.off()


#### ADNEX without CA125 ####

## AUC
ADNEXwoSens.AUC <- AUCimp.IOTA(pred = pmalwo, outcome = binaryCDcorrect, center = centerRE, imp = .imp, data = Iota5SensI, method.MA = "SJ", titleGraph = "LR2: AUC per center")

ADNEXwoSens.AUC$dataPlot <- ADNEXwoSens.AUC$dataPlot[c(1:9, 11:16, 10, 17:19),]
ADNEXwoSens.AUC$Plot <- ADNEXwoSens.AUC$Plot[c(1:9, 11:16, 10, 17:19),]

ADNEXwoSens.AUC$Plot[17, "RRauc"] <- "   "
ADNEXwoSens.AUC$Plot[19, "RRauc"] <- "        (0.81 to 0.98)"
ADNEXwoSens.AUC$Plot[19, "RRprev"] <- "   "

tiff("Results Paper 2/With cysts/Supplementary/Forest ADNEXwo - Sens FU PI.tiff", width = 1000, height = 650)
forestplot(ADNEXwoSens.AUC$Plot,
           align = c("l", "c", "c"),
           mean = ADNEXwoSens.AUC$dataPlot$AUC,
           lower = ADNEXwoSens.AUC$dataPlot$LL,
           upper = ADNEXwoSens.AUC$dataPlot$UL,
           is.summary = c(TRUE, FALSE, rep(FALSE, length(ADNEXwoSens.AUC$IncludedCenters)), FALSE, TRUE, TRUE),
           title = "",
           xlab = "AUC", # (95% CI)",
           xlog = TRUE,
           xticks = c(0.70, 0.80, 0.9, 1),
           clip = c(0.71, 1),
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55), xlab = gpar(cex = 1.5), label = gpar(cex = 1.5),
                            ticks = gpar(cex = 1.5), title = gpar(cex = 1.75)),
           graphwidth = unit(6, "cm"),
           col = fpColors(box = "gray", lines = "black", summary = "black"),
           fn.ci_sum = c(as.list(rep("fpDrawSummaryCI", length(ADNEXwoSens.AUC$IncludedCenters) + 3)),
                         fpDrawSummaryCI,
                         fpDrawBarCI),
           zero = ADNEXwoSens.AUC$Performance$AUC)
dev.off()


#### Summary plot for AUC ####

NA.forest <- RMIsens.AUC$Performance[1,]
NA.forest <- NA

Summary.AUC <- rbind(NA.forest, RMIsens.AUC$Performance[1,], LR2sens.AUC$Performance[1,], SRrisksSens.AUC$Performance[1,], ADNEXwoSens.AUC$Performance[1,], ADNEXwSens.AUC$Performance[1,])
Summary.AUC$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.AUC.PI <- rbind(NA.forest, RMIsens.AUC$Performance[2,], LR2sens.AUC$Performance[2,], SRrisksSens.AUC$Performance[2,], ADNEXwoSens.AUC$Performance[2,], ADNEXwSens.AUC$Performance[2,])
Summary.AUC.PI$Model <- c('', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')

Summary.AUC
tabletext <- cbind(
  c('Model', 'RMI', 'LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125'),
  c('AUC (95% CI)', 
    paste(format(round(RMIsens.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(RMIsens.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(RMIsens.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(LR2sens.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(LR2sens.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(LR2sens.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(SRrisksSens.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(SRrisksSens.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(SRrisksSens.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwoSens.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwoSens.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwoSens.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = ""),
    paste(format(round(ADNEXwSens.AUC$Performance$AUC[1], 2), nsmall = 2), " (", format(round(ADNEXwSens.AUC$Performance$LL[1], 2), nsmall = 2), " to ", format(round(ADNEXwSens.AUC$Performance$UL[1], 2), nsmall = 2), ")", sep = "")),
  c('95% PI', 
    paste0("(", format(round(RMIsens.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(RMIsens.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(LR2sens.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(LR2sens.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(SRrisksSens.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(SRrisksSens.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwoSens.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwoSens.AUC$Performance$UL[2], 2), nsmall = 2), ")"), 
    paste0("(", format(round(ADNEXwSens.AUC$Performance$LL[2], 2), nsmall = 2), " to ", format(round(ADNEXwSens.AUC$Performance$UL[2], 2), nsmall = 2), ")"))
)

tiff("Results Paper 2/With cysts/Paper/Summary AUC - Sens FU PI.tiff", width = 31, height = 13.75, units = "cm", res = 300)
forestplot(labeltext = tabletext,
           title = "",
           mean = round(Summary.AUC$AUC, 3),
           lower = round(Summary.AUC$LL, 3),
           upper = round(Summary.AUC$UL, 3),
           is.summary = c(FALSE, TRUE, TRUE, TRUE, TRUE,TRUE),
           xlab = "AUC (95% CI)",
           boxsize = .5,
           txt_gp = fpTxtGp(summary = gpar(cex = 1.55, fontface = "plain"), xlab = gpar(cex = 1.5, fontface= "plain"), label = gpar(cex = 1.5, fontface = "bold"),
                            ticks = gpar(cex = 1.5, fontface = "plain"), title = gpar(cex = 1.75)),
           graphwidth = unit(9, "cm"),
           graph.pos = 3,
           xticks = c(0.85, 0.9, 0.95, 1), xlog = TRUE, clip = c(0.85, 1))
dev.off()


#### 5.5.4 Calibration of the risk of malignancy ####

#### RMI ####
tiff("Results Paper 2/With cysts/Supplementary/Calibration RMI - Sens FU.tiff", width = 14, height = 14, units = "cm", res = 300)
RE.ValProbImp.RMI(p = RMI, y = binaryCDcorrect, data = Iota5SensI, center = centerRE, imputation.id = .imp, patientid = Patient.ID)
dev.off()


#### LR2 ####
LR2sens <- RE.ValProbImp(p = LR2, y = binaryCDcorrect, center = centerRE, data = Iota5SensI, imputation.id = .imp, patientid = Patient.ID)


#### Simple Rules Risk ####
SRriskSens <- RE.ValProbImp(p = SRrisks, y = binaryCDcorrect, center = centerRE, data = Iota5SensI, imputation.id = .imp, patientid = Patient.ID, title = "Simple Rules Risk: Calibration curve")


#### ADNEX with CA125 ####
ADNEXwSens <- RE.ValProbImp(p = pmalw, y = binaryCDcorrect, data = Iota5SensI, center = centerRE, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX with CA125: Calibration curve")


#### ADNEX without CA125 ####
ADNEXwoSens <- RE.ValProbImp(p = pmalwo, y = binaryCDcorrect, center = centerRE, data = Iota5SensI, imputation.id = .imp, patientid = Patient.ID, title = "ADNEX without CA125: Calibration curve")


#### Summary plot for Calibration ####

## Predict the outcome for the calibration curve
# LR2
OverallCal.LR2 = LR2sens$Plot$y
p.LR2 = LR2sens$Plot$x

# Simple Rules Risk
OverallCal.SRrisk = SRriskSens$Plot$y
p.SR = SRriskSens$Plot$x

# ADNEX with CA125
OverallCal.ADNEXw = ADNEXwSens$Plot$y
p.ADNEXw = ADNEXwSens$Plot$x

# ADNEX without CA125
OverallCal.ADNEXwo = ADNEXwoSens$Plot$y
p.ADNEXwo = ADNEXwoSens$Plot$x

table <- matrix(ncol = 3, nrow = 4)
colnames(table) <- c('Model', 'Intercept (95% CI)', 'Slope (95% CI)')
table[, 1] <- c('LR2', 'SRRisk', 'ADNEX without CA125', 'ADNEX with CA125')
table[1, 2:3] <- c(paste0(format(round(LR2sens$SumInt$Est, 2), nsmall = 2), " (", format(round(LR2sens$SumInt$LL, 2), nsmall = 2), " to ", format(round(LR2sens$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(LR2sens$SumSlo$Est, 2), nsmall = 2), " (", format(round(LR2sens$SumSlo$LL, 2), nsmall = 2), " to ", format(round(LR2sens$SumSlo$UL, 2), nsmall = 2), ")"))
table[2, 2:3] <- c(paste0(format(round(SRriskSens$SumInt$Est, 2), nsmall = 2), " (", format(round(SRriskSens$SumInt$LL, 2), nsmall = 2), " to ", format(round(SRriskSens$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(SRriskSens$SumSlo$Est, 2), nsmall = 2), " (", format(round(SRriskSens$SumSlo$LL, 2), nsmall = 2), " to ", format(round(SRriskSens$SumSlo$UL, 2), nsmall = 2), ")"))
table[4, 2:3] <- c(paste0(format(round(ADNEXwSens$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwSens$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwSens$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwSens$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwSens$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwSens$SumSlo$UL, 2), nsmall = 2), ")"))
table[3, 2:3] <- c(paste0(format(round(ADNEXwoSens$SumInt$Est, 2), nsmall = 2), " (", format(round(ADNEXwoSens$SumInt$LL, 2), nsmall = 2), " to ", format(round(ADNEXwoSens$SumInt$UL, 2), nsmall = 2), ")"), paste0(format(round(ADNEXwoSens$SumSlo$Est, 2), nsmall = 2), " (", format(round(ADNEXwoSens$SumSlo$LL, 2), nsmall = 2), " to ", format(round(ADNEXwoSens$SumSlo$UL, 2), nsmall = 2), ")"))

## Make graph
x = seq(0, 1, by = 0.05)
y = seq(0, 1, by = 0.05)
tiff("Results Paper 2/With cysts/Paper/Summary calibration - Sens FU.tiff", width = 14, height = 14, units = "cm", res = 300)
plot(x, y, xlim = c(0,1), ylim = c(0,1), type = "l", col = "gray50", lwd = 2, lty = 2, 
     xlab = "Estimated risk of malignancy", ylab = "Observed proportion of malignancy",
     main = "", cex.lab = 1, cex.axis = 1, las = 1) 
lines(p.LR2, OverallCal.LR2, lwd = 2, col = "blue")
lines(p.SR, OverallCal.SRrisk, lwd = 2, col = "red")
lines(p.ADNEXwo, OverallCal.ADNEXwo, lwd = 2, col = "darkgreen")
lines(p.ADNEXw, OverallCal.ADNEXw, lwd = 2, col = "darkorange")
legend(x = -0.035, y = 1, legend = c( "Ideal", "LR2", "SRRisk", "ADNEX without CA125", "ADNEX with CA125"),
       col = c("gray50", "blue", "red", "darkgreen", "darkorange"), lty = c(2,1,1,1,1), lwd = 2, cex = 0.7, bty = "n", ncol = 1)
addtable2plot(x = 0.16, y = -0.01, table = table, display.colnames= TRUE, cex = 0.7) 
dev.off()

save.image(file = "IOTA5 - Validation 25032020.RData") # Save the entire workspace

