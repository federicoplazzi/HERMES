#######################################################################################################################
#                                                                                                                     #
# HERMES is a straightforward index which tries to summarize the mitochondrial evolution pace using a single number.  #
# Several mitogenomic features are evaluated in a factor analysis framework; namely, in the current version, they are #
# the %URs, the Amount of Mitochondrial Identical Gene Arrangements (AMIGA) index, the absolute value of the SU skew, #
# the root-to-tip distance, the ML pairwise distance from a given outgroup, the %AT, the AT skew, the GC skew, the    #
# number of (annotated) genes, the length of the mtDNA, and the CAI.                                                  #
#                                                                                                                     #
# Copyright (C) 2020 Guglielmo Puccio, Federico Plazzi                                                                #
#                                                                                                                     #
# This program is free software: you can redistribute it and/or modify                                                #
# it under the terms of the GNU General Public License as published by                                                #
# the Free Software Foundation, either version 3 of the License, or	                                              #					      #
# (at your option) any later version.                                                                                 #
#                                                                                                                     #
# This program is distributed in the hope that it will be useful,                                                     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                       #
# GNU General Public License for more details.                                                                        #
#                                                                                                                     #
# You should have received a copy of the GNU General Public License                                                   #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                               #
#                                                                                                                     #
#######################################################################################################################

#HERMES version: 3.0

HERMES.version <- "3.0"

#Load library psych for factor analysis.

library(psych)

#Read the command line searching for the alpha argument.

HERMES.args <- commandArgs(trailingOnly=TRUE)

alpha <- 0.05
acceptable.KMO <- 0.6
acceptable.communality <- 0.5
order.option <- FALSE

#Read input file and parsing shades (note that the maximum allowed number of groups/shades is 657 - see the colors/colours function).

HERMES.file <- read.table(file="./Results/HERMES_variables.txt",header=TRUE,sep="\t",row.names=1)
if ("Shades" %in% colnames(HERMES.file)) {
	HERMES.shades <- HERMES.file[,"Shades"]
	HERMES.file <- HERMES.file[,colnames(HERMES.file)[colnames(HERMES.file) != "Shades"]]
	} else HERMES.shades <- rep(24,length(rownames(HERMES.file)))

#Read input file and remove variables with standard deviation equal to 0 (i.e., all identical values)

variable.names <- colnames(HERMES.file)
for (c in 1:length(variable.names)) {
	if (sd(HERMES.file[,variable.names[c]]) == 0) {
		HERMES.file <- HERMES.file[,colnames(HERMES.file)[colnames(HERMES.file) != variable.names[c]]]
		message(paste("The '",variable.names[c],"' column was dropped (standard deviation = 0)! Were all values identical?",sep=""))
		}
	}

#Compute the total number of possible variable combinations.

total.combinations <- 0
for (m in 2:dim(HERMES.file)[2]) total.combinations <- total.combinations+factorial(dim(HERMES.file)[2])/(factorial(m)*factorial(dim(HERMES.file)[2]-m))
landmarks <- round(seq(0.05,1,by=0.05)*total.combinations,digits=0)
explored.combinations <- 0
all.combinations <- list()

#Compute HERMES score with each possible combination. For each combination, the number of acceptable goodness-of-fit
#parameters is registered to the 'parameters' vector.

combination <- character()
KMO <- numeric()
TLI <- numeric()
SRMR <- numeric()
RMSEA <- numeric()
RMSEA.lower <- numeric()
RMSEA.upper <- numeric()
total.communality <- numeric()
parameters <- numeric()

for (m in 2:dim(HERMES.file)[2]) {
	combination.list <- combn(c(1:dim(HERMES.file)[2]),m,simplify=FALSE)
	combination.number <- length(combination.list)
	for (cc in 1:combination.number) {
		current.combination <- colnames(HERMES.file)[combination.list[[cc]]][1]
		for (p in 2:length(combination.list[[cc]])) current.combination <- paste(current.combination,colnames(HERMES.file)[combination.list[[cc]]][p],sep="_")
		explored.combinations <- explored.combinations+1

#Check whether factor analysis is possible with a given combination of variables.

		if (class(try(fa(HERMES.file[,combination.list[[cc]]],rotate="varimax",scores="tenBerge",fm="ml",cor="cor",normalize=TRUE,missing=TRUE,impute="median",alpha=alpha),silent=TRUE))[1] == "try-error") {
			message(paste("Factor analysis failed with combination '",current.combination,"' (",explored.combinations,"/",total.combinations,")!",sep=""))
			} else {

#If checkpoint is passed, perform factor analysis...

			good.parameters <- 0
			all.combinations[[length(all.combinations)+1]] <- combination.list[[cc]]
			HERMES.fa <- fa(HERMES.file[,combination.list[[cc]]],rotate="varimax",scores="tenBerge",fm="ml",cor="cor",normalize=TRUE,missing=TRUE,impute="median",alpha=alpha)
			combination <- c(combination,current.combination)
			KMO <- c(KMO,KMO(HERMES.file[,combination.list[[cc]]])$MSA)
			ifelse(KMO(HERMES.file[,combination.list[[cc]]])$MSA > acceptable.KMO,good.parameters <- good.parameters+1,NA)
			TLI <- c(TLI,HERMES.fa$TLI)
			ifelse(HERMES.fa$TLI > 0.95,good.parameters <- good.parameters+1,NA)
			SRMR <- c(SRMR,HERMES.fa$rms)
			ifelse(HERMES.fa$rms < 0.08,good.parameters <- good.parameters+1,NA)
			if (length(combination.list[[cc]]) > 3) {
				RMSEA <- c(RMSEA,as.numeric(HERMES.fa$RMSEA[1]))
				ifelse(as.numeric(HERMES.fa$RMSEA[1]) < 0.06,good.parameters <- good.parameters+1,NA)
				RMSEA.lower <- c(RMSEA.lower,as.numeric(HERMES.fa$RMSEA[2]))
				ifelse(as.numeric(HERMES.fa$RMSEA[2]) < 0.06,good.parameters <- good.parameters+1,NA)
				RMSEA.upper <- c(RMSEA.upper,as.numeric(HERMES.fa$RMSEA[3]))
				ifelse(as.numeric(HERMES.fa$RMSEA[3]) < 0.06,good.parameters <- good.parameters+1,NA)
				} else {
				RMSEA <- c(RMSEA,NA)
				RMSEA.lower <- c(RMSEA.lower,NA)
				RMSEA.upper <- c(RMSEA.upper,NA)
				}
			total.communality <- c(total.communality,mean(as.numeric(HERMES.fa$communality)))
			good.parameters <- good.parameters+length(as.numeric(HERMES.fa$communality)[as.numeric(HERMES.fa$communality) > acceptable.communality])
			parameters <- c(parameters,good.parameters)
			}

#Update the user on the work in progress.

		if (explored.combinations %in% landmarks) {
			message(paste("Explored ",explored.combinations," out of ",total.combinations," variable combinations (",round(explored.combinations/total.combinations*100,digits=0),"%).",sep=""))
			}
		}
	}

#Write the final result table, keeping only best-performing set of variables.

HERMES.data.frame <- data.frame(combination=combination,KMO=KMO,TLI=TLI,SRMR=SRMR,RMSEA=RMSEA,RMSEA.lower=RMSEA.lower,RMSEA.upper=RMSEA.upper,total.communality=total.communality,parameters=parameters)
best.HERMES.data.frame <- HERMES.data.frame[HERMES.data.frame$parameters == max(HERMES.data.frame$parameters),]

#In case of ties, select the best set of variables. Each goodness-of-fit test yields a value: the ratio is computed
#between the threshold and the value (for statistics that have a minimum threshold - KMO, TLI, and communalities) or
#between the value and the threshold (for statistics that have a maximum threshold - SRMR and RMSEA). The best set of
#variables is the one that minimizes the mean of such ratios.

if (dim(best.HERMES.data.frame)[1] == 1) {
	message(paste("Found one single best variable combination: it passed ",max(HERMES.data.frame$parameters)," tests!\nPlotting the HERMES index...",sep=""))
	best.HERMES.set <- as.numeric(rownames(best.HERMES.data.frame))
	write.table(best.HERMES.data.frame,file="./Results/best.HERMES.out",quote=FALSE,sep="\t",row.names=FALSE)
	} else {
	message(paste("Found ",dim(best.HERMES.data.frame)[1]," best-performing variable combinations: ",max(HERMES.data.frame$parameters)," tests were passed!\nSelecting the best set and computing HERMES scores...",sep=""))
	HERMES.performance <- numeric()
	for (bHr in 1:dim(best.HERMES.data.frame)[1]) {
		performances <- numeric()
		current.best.HERMES.combination <- all.combinations[[as.numeric(rownames(best.HERMES.data.frame)[bHr])]]
		HERMES.fa <- fa(HERMES.file[,current.best.HERMES.combination],rotate="varimax",scores="tenBerge",fm="ml",cor="cor",normalize=TRUE,missing=TRUE,impute="median",alpha=alpha)
		performances <- c(performances,acceptable.KMO/KMO(HERMES.file[,current.best.HERMES.combination])$MSA)
		performances <- c(performances,0.95/HERMES.fa$TLI)
		performances <- c(performances,HERMES.fa$rms/0.08)
		if (length(current.best.HERMES.combination) > 3) for (RMSEA in 1:3) performances <- c(performances,as.numeric(HERMES.fa$RMSEA[RMSEA])/0.06)
		for (c in 1:length(as.numeric(HERMES.fa$communality))) performances <- c(performances,acceptable.communality/as.numeric(HERMES.fa$communality)[c])
		HERMES.performance <- c(HERMES.performance,mean(performances,na.rm=TRUE))
		}
	best.HERMES.data.frame$HERMES.performance <- HERMES.performance
	combination.order <- order(best.HERMES.data.frame$HERMES.performance,decreasing=TRUE)
	best.HERMES.data.frame <- best.HERMES.data.frame[combination.order,]
	write.table(best.HERMES.data.frame,file="./Results/best.HERMES.out",quote=FALSE,sep="\t",row.names=FALSE)
	best.HERMES.set <- as.numeric(rownames(best.HERMES.data.frame[1,]))

#(Note that, in case of further ties, the first tested combination is arbitrarily selected between those minimizing this mean.)

	}

#Carry out the final HERMES analysis.

best.HERMES.combination <- all.combinations[[best.HERMES.set]]
HERMES.fa <- fa(HERMES.file[,best.HERMES.combination],rotate="varimax",scores="tenBerge",fm="ml",cor="cor",normalize=TRUE,missing=TRUE,impute="median",alpha=alpha)
HERMES.fa.scores <- as.numeric(HERMES.fa$scores)
HERMES.communalities <- paste(colnames(HERMES.file[,best.HERMES.combination]),"communality",sep=" ")
HERMES.KMO <- KMO(HERMES.file[,best.HERMES.combination])$MSA
HERMES.TLI <- HERMES.fa$TLI
HERMES.SRMR <- HERMES.fa$rms
HERMES.RMSEA <- as.numeric(HERMES.fa$RMSEA[1])
HERMES.RMSEA.lower <- as.numeric(HERMES.fa$RMSEA[2])
HERMES.RMSEA.upper <- as.numeric(HERMES.fa$RMSEA[3])
HERMES.total.communality <- mean(as.numeric(HERMES.fa$communality))
parameters <- c("KMO","TLI","SRMR","RMSEA",paste("RMSEA ",(1-alpha)*100,"% CI (lower bound)",sep=""),paste("RMSEA ",(1-alpha)*100,"% CI (upper bound)",sep=""),HERMES.communalities,"Total communality")
values <- c(HERMES.KMO,HERMES.TLI,HERMES.SRMR,HERMES.RMSEA,HERMES.RMSEA.lower,HERMES.RMSEA.upper,as.numeric(HERMES.fa$communality),HERMES.total.communality)
HERMES.parameters <- data.frame(parameters=parameters,values=values)
write.table(HERMES.parameters,file=paste("./Results/HERMES-v",HERMES.version,"_parameters.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
HERMES.fa.scores.data.frame <- data.frame(Species=rownames(HERMES.file[,best.HERMES.combination]),HERMES=HERMES.fa.scores,Shade=HERMES.shades)
HERMES.fa.scores.data.frame <- HERMES.fa.scores.data.frame[order(HERMES.fa.scores,decreasing=order.option),]
write.table(HERMES.fa.scores.data.frame,file=paste("./Results/HERMES-v",HERMES.version,"_scores.txt",sep=""),quote=FALSE,sep="\t",row.names=FALSE)

#Plot HERMES!

pdf(file=paste("./Results/HERMES-v",HERMES.version,"_plot.pdf",sep=""))
plot(HERMES.fa.scores.data.frame$HERMES,main="Hyper-Empirical Relative Mitochondrial Evolutionary Speed",xlab="",ylab="HERMES",col=colors()[HERMES.fa.scores.data.frame$Shade],bg=colors()[HERMES.fa.scores.data.frame$Shade],pch=21,xaxt="n")
axis(side=1,at=c(1:dim(HERMES.fa.scores.data.frame)[1]),labels=HERMES.fa.scores.data.frame$Species,las=2)
dev.off()
