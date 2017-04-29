#inputFile='feature_infor.txt'
#inputDirectory='/home/yli4/development/JUMPg/JUMPg_v2.3.3/gnm_stage1_test1/intermediate/sum_accepted_PSMs/misc'
#nTotLines=3724
#search_engine='jump'

setwd(inputDirectory);
        #suppressPackageStartupMessages(suppressWarnings(library(tcltk)))
        library(MASS);
        #library(qvalue);
        #search_engine=\$search_engine;
        ## Data loading
        #table5rows = read.table(inputFile, sep = "\t", header = T, nrows = 5);
        table5rows = read.table(inputFile, sep = "\t", header = T, nrows = 30);
        classes = sapply(table5rows, class);
        data = read.table(inputFile, sep = "\t", header = T, nrows = nTotLines, colClasses = classes);
        ## Initialization
        numSamples = dim(data)[1];
        pval = rep(1, numSamples);
        qval = rep(1, numSamples);
        fdr = rep(1, numSamples);
        ## Independent LDA model for each charge state
        chargeStates = c(1, 2, 3, 4);
        for (charge in chargeStates) {
                colName = paste("charge", charge, sep = "");
                subInd = which(data[colnames(data) == colName] == 1);
                if (length(subInd) > 30) {
                        subData = data[subInd, ];
                        ## Generation of a data matrix as in Du et al.
                        ## Three variables: log(XCorr)/log(peptideLength), sqrt(dCn) and ppm
                        if (!is.na(subData$ppm))
                        {
                        if ( search_engine == 'jump' )
                        {
                                if (sd(subData$dCn)>0)
                                {
                                        X = as.data.frame(cbind(subData$xcorr, sqrt(subData$dCn), abs(subData$ppm))); # maybe for Jscore?
                                }
                                else
                                {
                                        X = as.data.frame(cbind(subData$xcorr,abs(subData$ppm))); # dJn is constant
                                }
                        }
                        else
                        {
                                X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn), abs(subData$ppm)));
                                #X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn)));
                        }
                        } # end of if (sd(subData$ppm)>0)
                        else # ppm not considered
                        {
                        if ( search_engine == 'jump' )
                        {
                                if (sd(subData$dCn)>0)
                                {
                                        X = as.data.frame(cbind(subData$xcorr, sqrt(subData$dCn))); # for Jscore
                                }
                                else
                                {
                                        X = as.data.frame(cbind(subData$xcorr)); # dJn is constant
                                }
                        }
                        else
                        {
                                #X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn), abs(subData$ppm)));
                                X = as.data.frame(cbind((log(subData$xcorr)/log(subData$peptideLength)), sqrt(subData$dCn)));
                        }
                        } # end of else (sd(subData$ppm)>0)
                        group = subData$targetDecoy;
                        ## LDA modeling
                        if (length(table(group))>=2 & table(group)[1]>1 & table(group)[2]>1) {
                                ldaModel = lda(group ~ ., X);
                                ## Calculation of LDA scores of the samples
                                scores = as.numeric(predict(ldaModel, X)$x);
                                ## Let the scores of "target"s be always greater than "decoy"s
                                if (mean(scores[group == "target"]) < mean(scores[group == "decoy"])) {
                                        scores = -scores;
                                }
                                ## Null hypothesis distribution by fitting "decoy" scores to a normal distribution
H0 = fitdistr(scores[group == "decoy"], "normal");
                                p = 1 - pnorm(scores, mean = H0$estimate[1], sd = H0$estimate[2]);
                                pval[subInd] = p;
                                ## Calculate q-values and FDR (Benjamini-Hochberg)
                                #qval[subInd] = qvalue(p)$qvalue;
                                qval[subInd] = 1;
                                fdr[subInd] = p.adjust(p, method = "BH");
                        } else {
                                pval[subInd] = 1;
                                qval[subInd] = 1;
                                fdr[subInd] = 1;
                        }
                }
                else
                {
                        pval[subInd] = 1;
                        qval[subInd] = 1;
                        fdr[subInd] = 1;
                }
        }
        write.table(cbind(pval, qval, fdr), file = "LDA_Result.txt", sep = "\t");

