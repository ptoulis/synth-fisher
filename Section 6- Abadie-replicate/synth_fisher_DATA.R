#' Alabama, 05/2004
#' Arkansas, 02/1993  (22->34, or 06/2003, 34->59)
#' California, 01/1989
#' Colorado, 01/2005
#' Connecticut, 04/1989
#' Delaware, 01/1991 (14->24, or larger >100% at 08/2003, 24->55)
#' Georgia, 07/2003
#' Idaho, 07/1994 (18-> 28, larger >100% at 06/2003, 28->57)
#' Illinois, 07/1989
#' Indiana, 07/2002
#' Iowa, 04/2007
#' Kansas, 07/2002
#' Kentucky, 06/2005
#' Louisiana, 08/2002
#' Maine, 11/1997 (37->74, could include 07/1991 bc of multiple increases, 16->37 in 10 years.)
#' Minnesota, 08/2005 (very large increase > 300%, or 06/91 but from 18->43 in several years of increase)
#' Mississippi, 05/2009
#' Missouri, *
#' Montana, 05/2003
#' Nebraska, 10/2002
#' Nevada, 07/1989
#' New Hampshire, 02/1990
#' New Mexico, 07/2003
#' North Carolina, 09/2005 (another increase in 08/1991 from $2 to $5.)
#' North Dakota, 05/1989
#' Ohio, 07/2002
#' Oklahoma, 01/2005
#' Pennsylvania, 08/1991
#' Rhode Island, 07/1997 (.->71, or 07/1993, 18->44 from 1982 to 1993)
#' South Carolina, 07/2010
#' South Dakota, 03/2003  (33->53, or 07/1995, 23->33)
#' Tennesee, 07/2002
#' Texas, 07/1990
#' Utah, 07/1991 (10->26.4 since 1980, or 07/1997 jump 26.5->51.5)
#' Vermont, 07/1995 (20->44, multiple increases since 1991)
#' Virginia, 09/2004
#' West Virginia, 05/2003
#' Wisconsin, 05/1992 (multiple 16->32 from 1980, starts 05/1982, then multiple increases)
#' Wyoming, 07/2003 (12->60, one more at 07/1989 but from 8->12)
#' 

AllStates = c("Alabama", "Arkansas", "California", "Colorado", "Connecticut", 
              "Delaware", "Georgia", "Idaho", "Illinois", "Indiana", 
              "Iowa", "Kansas", "Kentucky", "Louisiana",  "Maine", 
              "Minnesota", "Mississippi", "Missouri", "Montana", "Nebraska", 
              "Nevada", "New Hampshire", "New Mexico", "North Carolina", "North Dakota", 
              "Ohio", "Oklahoma", "Pennsylvania", "Rhode Island", "South Carolina", 
              "South Dakota", "Tennessee", "Texas", "Utah", "Vermont", 
              "Virginia", "West Virginia", "Wisconsin", "Wyoming")

adoption_Data = data.frame(state=AllStates, 
                           when=c("05/2004", "06/2003", "01/1989", "01/2005", "04/1989",
                                  "01/1991", "07/2003", "07/1994", "07/1989", "07/2002",
                                  "04/2007", "07/2002", "06/2005", "08/2002", "11/1997",
                                  "06/1991", "05/2009", "12/2014", "05/2003", "10/2002",
                                  "07/1989", "02/1990", "07/2003", "09/2005", "05/1989",
                                  "07/2002", "01/2005", "08/1991", "07/1997", "07/2010", 
                                  "03/2003", "07/2002", "07/1990", "07/1991", "07/1995", 
                                  "09/2004", "05/2003", "05/1992", "07/2003"))

adoption_Data_2 = data.frame(state=AllStates, 
                           when=c("05/2004", "02/1993", "01/1989", "01/2005", "04/1989",
                                  "08/2003", "07/2003", "06/2003", "07/1989", "07/2002",
                                  "04/2007", "07/2002", "06/2005", "08/2002", "07/1991",
                                  "08/2005", "05/2009", "12/2014", "05/2003", "10/2002",
                                  "07/1989", "02/1990", "07/2003", "09/2005", "05/1989",
                                  "07/2002", "01/2005", "08/1991", "07/1993", "07/2010", 
                                  "07/1995", "07/2002", "07/1990", "07/1997", "07/1995", 
                                  "09/2004", "05/2003", "05/1992", "07/1989"))

n = nrow(adoption_Data)
equal = which(as.character(adoption_Data[,2]) == as.character(adoption_Data_2[,2]))
rest = setdiff(1:n, equal)
M = matrix(0, nrow=1, ncol=n)
for(j in 1:length(rest)) {
  # all combinations of j
  A = t(combn(rest, j))
  for(i in 1:nrow(A)) {
    e0 = rep(0, n)
    e0[A[i,]] = 1
    M = rbind(M, e0)
  }
}
# 32 models in M
adoption_DataList = list()
for(i in 1:nrow(M)) {
  Mi = M[i,]
  adoption_DataList[[i]] = data.frame(state=AllStates, when=ifelse(Mi==0, as.character(adoption_Data[,2]), as.character(adoption_Data_2[,2])))
}

stopifnot(all(as.character(adoption_DataList[[1]][,2]) == as.character(adoption_Data[,2])))
stopifnot(all(as.character(adoption_DataList[[nrow(M)]][,2]) == as.character(adoption_Data_2[,2])))

print(sprintf("Found %d total specifications...", length(adoption_DataList)))
print("Different states")
print(adoption_Data[rest,])
print("Alternative specification --->")
print(adoption_Data_2[rest, ])

print("#####     DONE WITH DATA LOAD     #####")
print("#####    --------------------     #####")


load_smoke_data_plot = function() {
  load("smoking_with_dems.rdata")
  Smoking = smoking
  stopifnot(all(Smoking$cigsale==smoking$cigsale))
  Smoking$retprice = log(Smoking$retprice)
  Smoking$blue = smoking$dems
  rm(smoking)
  
  # load("smoking20.RData"
  kStates = unique(Smoking$state)
  kStatesId = unique(Smoking$unit.num)
  
  load("RTtobacco.RData")
  Legislation = data.frame(state=kStates)
  Legislation$time = sapply(kStates, function(st) {
    i = which(RTtobacco$state==st)
    RTtobacco[i, ]$time
  })
  rm(RTtobacco)
  # Nevada could be 89
  # Illinois also 89
  # Rhode Island 89
  # Wisco -> 92?
  # Arkansas -> 93
  # Pennsylvannia -> 91
  # South Dakota 95
  correct = list("New Hampshire"=1998, "Idaho"=1994, "Nevada"=2004, "Illinois"=1989, "Rhode Island"=1989, 
                 "Wisconsin"=1992, "Arkansas"=2003, "Pennsylvania"=1991, 
                 "South Dakota"=1995, "Utah"=1997)
  for(st in names(correct)) {
    i = which(Legislation$state==st)
    Legislation[i, 2] <- correct[[st]]
  }
  print("[INFO] -- Legislation times defined.")
  return(list(st=kStates, stId=kStatesId, smoke=Smoking, legal=Legislation))
}

Make_Fig1 = function() {
  D = load_smoke_data_plot()$smoke
  incl_states = c("California", "Illinois", "Kentucky", "Utah")
  D = subset(D, state %in% incl_states)
  
  g  = ggplot(data=D, aes(x=year, y=cigsale))
  g = g + geom_line(aes(lty=state), size=0.9)
  #g = g + geom_vline(xintercept = 1989) + geom_text(x=1987, y=215, label="T_CA", size=8)
  #g = g + geom_vline(xintercept = 1990) + geom_text(x=1991, y=190, label="T_IL", size=8)
  #g = g + geom_vline(xintercept = 2005) + geom_text(x=2010, y=190, label="T_KY", size=8)
  
  
  g = g+ theme(legend.text=element_text(size=18), 
               legend.title = element_text(size=20),
               axis.text.x = element_text(size=18),
               axis.title.x = element_text(size=20),
               axis.text.y = element_text(size=18), 
               axis.title.y = element_text(size=20),
               legend.key.size = unit(1.2, 'cm'), #change legend key size
               legend.key.width = unit(1.8, 'cm')
  ) 
  plot(g)
}

