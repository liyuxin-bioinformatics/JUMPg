#inputFile='all.tag.pks.acp.tab';
#inputDirectory='/home/yli4/development/JUMPg/JUMPg_v2.3.4/gnm_stage1_test1/intermediate/qc_accepted_PSMs';
#recoveryPct=99

setwd(inputDirectory);
library(MASS);
        ## Data loading
        tb=read.table(inputFile,head=T,sep="\t")
        tb=tb[tb$PPI==1,]
        tb$topTagEvalue[is.na(tb$topTagEvalue)]=0
        # specify training set
        tb$trainingSet=rep('notUsed',nrow(tb))
        tb$trainingSet[tb$accepted==1]='positive'
        tb$trainingSet[tb$matchedType=='unmatched']='negative'
        # Generation of a data matrix for LDA analysis
        subData=tb[tb$trainingSet!='notUsed',]
        X = as.data.frame(cbind(subData$topTagLength,subData$topTagEvalue, sqrt(subData$ms2PeaksCount), log10(subData$ms2Int)))
        group = subData$trainingSet
        # LAD analysis
        calclda <- function(variables,loadings)
        {
                as.data.frame(variables)
                numsamples <- nrow(variables)
                ld <- numeric(numsamples)
                numvariables <- length(variables)
                for (i in 1:numsamples)
                {
                        valuei <- 0
                        for (j in 1:numvariables)
                        {
                                valueij <- variables[i,j]
                                loadingj <- loadings[j]
                                valuei <- valuei + (valueij * loadingj)
                        }
                        ld[i] <- valuei
                }
                return(ld)
        }
        if (length(table(group))>=2 & table(group)[1]>1 & table(group)[2]>1)
        {
                # LDA modeling
                ldaModel = lda(group ~ ., X);
                scores = as.numeric(predict(ldaModel, X)$x);
                # apply same model to the whole dataset
                m=as.data.frame(cbind(tb$topTagLength,tb$topTagEvalue, sqrt(tb$ms2PeaksCount), log10(tb$ms2Int)))
                tb$ldaScore=calclda(m[,1:4], ldaModel$scaling[,1])
                # binning for ldaScore
                binN=100
                intv=(max(tb$ldaScore)-min(tb$ldaScore))/binN
                tb$binldaScore=floor((tb$ldaScore-min(tb$ldaScore))/intv)
                mp1=table(tb$binldaScore,tb$accepted)
                mp=matrix(nrow=nrow(mp1),ncol=6)
                colnames(mp)=c('fail','accept','successRate','identifiableSpectra','spectraUsed4search','IDlost')
                rownames(mp)=rownames(mp1)
                mp[,1:2]=mp1
                mp=as.data.frame(mp)
                mp$successRate=mp$accept/(mp$fail+mp$accept)
                #mp$identifiableSpectra=mp$fail*mp$successRate # based on successRate*fail
                mp$identifiableSpectra=mp$accept        # same as accept
                #spectraUsed4search
                #mp[,5]=mp[,1] # only fail
                mp[,5]=mp[,1]+mp[,2] # both fail and accept
                for (i in (nrow(mp)-1):1)
                {
                        #mp[i,5]=mp[i,1]+mp[i+1,5] # only fail
                        mp[i,5]=mp[i,1]+mp[i,2]+mp[i+1,5] # both fail and accept
                }
                #IDlost
                mp[,6]=mp$identifiableSpectra
                for (i in 2:nrow(mp))
                {
                        mp[i,6]=mp[i-1,6]+mp[i,4]
                }
                #mp$spectraUsed4searchPct=mp$spectraUsed4search/sum(mp[,1]) # only fail
                mp$spectraUsed4searchPct=mp$spectraUsed4search/(sum(mp[,1])+sum(mp[,2])) # both fail and accept
                mp$IDlostPct=mp$IDlost/sum(mp[,4])
                bin.lda=mp
                spectraSkip=100-tail(bin.lda$spectraUsed4searchPct[bin.lda$IDlostPct<(1-recoveryPct/100)]*100,1)
                # figure
                w=2
                par(las=1)
                pdf(file=paste(inputDirectory,'/PSM_recovery_vs_spctraUsed.pdf',sep=''),width=5,height=5)
                mp=bin.lda
                plot((1-mp$spectraUsed4searchPct)*100,(1-mp$IDlostPct)*100,xlab='% spectra skipped for database search',ylab='% expected identified PSMs',col='red',type='l',lwd=w)
                abline(h=recoveryPct,col='grey',lty=2)
                abline(v=spectraSkip,col='grey',lty=2)
                dev.off()
                write.table(tb,file=paste(inputDirectory,'/all.tag.pks.acp.qs.tab',sep=''),quote=F,sep="\t",row.names=F)
                write.table(bin.lda,file=paste(inputDirectory,'/threshold_table.txt',sep=''),quote=F,sep="\t")
        }
