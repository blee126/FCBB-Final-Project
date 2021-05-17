library(stats)
library(logisticPCA)

df <- read.table(file = 'C:\\Users\\brend\\Downloads\\210414-010014_export_variant.tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df <- df[-1,]
df <- df[c('X.6','X.9','X.12', 'X.13', 'Gene_Ontology', 'X.16', 'GTEx_eQTLs', 'X.18', 'X.19')]
df$id <- paste(df$X.6, df$X.9)
df$X.12 <- as.numeric(df$X.12)

df2 <- read.table(file = 'C:\\Users\\brend\\Downloads\\210414-013040_export_variant.tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df2 <- df2[-1,]
df2 <- df2[c('X.6','X.9','X.12', 'X.13')]
df2$id <- paste(df2$X.6, df2$X.9)
df2$X.12 <- as.numeric(df2$X.12)

df$sc2 <- df2$X.12[match(df$id, df2$id)]
df$sc2[is.na(df$sc2)] <- 0
df$casesum <- df$X.12 + df$sc2

df$caselist <- df2$X.13[match(df$id, df2$id)]
df$caselist[is.na(df$caselist)] <- ""
df$caselist <- mapply(c, strsplit(df$X.13, ";"), strsplit(df$caselist, ";"))

df3 <- read.table(file = 'C:\\Users\\brend\\Downloads\\210414-030053_export_variant (1).tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df3 <- df3[-1,]
df3 <- df3[c('X.6','X.9','X.12', 'X.13')]
df3$id <- paste(df3$X.6, df3$X.9)
df3$X.12 <- as.numeric(df3$X.12)

df4 <- read.table(file = 'C:\\Users\\brend\\Downloads\\210414-015456_export_variant (1).tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df4 <- df4[-1,]
df4 <- df4[c('X.6','X.9','X.12', 'X.13')]
df4$id <- paste(df4$X.6, df4$X.9)
df4$X.12 <- as.numeric(df4$X.12)

df$sc3 <- df3$X.12[match(df$id, df3$id)]
df$sc3[is.na(df$sc3)] <- 0
df$sc4 <- df4$X.12[match(df$id, df4$id)]
df$sc4[is.na(df$sc4)] <- 0
df$controlsum <- df$sc3 + df$sc4

df$cl1 <- df3$X.13[match(df$id, df3$id)]
df$cl1[is.na(df$cl1)] <- ""
df$cl2 <- df4$X.13[match(df$id, df4$id)]
df$cl2[is.na(df$cl2)] <- ""
df$controllist <- mapply(c, strsplit(df$cl1, ";"), strsplit(df$cl2, ";"))

dnsaids <- df[grepl("aspirin", df$X.18, fixed=TRUE),]
#dnsaids <- df[grepl("aspirin|ibuprofen|naproxen|diclofenac|celecoxib|mefenamic|etoricoxib|indomethacin", df$X.18, fixed=TRUE),]

results <- apply(df, 1, 
      function(x) {
        ca <- as.numeric(x['casesum'])
        co <- as.numeric(x['controlsum'])
        tbl <- matrix(c(ca, 22-ca, co, 22-co), nrow=2, byrow=T)
        ft <- fisher.test(tbl, alternative="two.sided")
        return(c(ft$p.value, ft$estimate))
      })
d <- as.data.frame(matrix(unlist(results), nrow = 2))
df$pvalue <- unlist(d[1,])
df$oddsratio <- unlist(d[2,])
df <-  df[!(is.na(df$X.9) | df$X.9==""), ]
vardf <- df

df <- df[df$pvalue < 0.05 & df$casesum + df$controlsum > 8,]

gobp <- unlist(strsplit(df$Gene_Ontology, ';'))
commonbioprocess <- sort(table(gobp),decreasing=TRUE)
print(commonbioprocess[1:10])

gomf <- unlist(strsplit(df$X.16, ';'))
commonmolfunction <- sort(table(gomf),decreasing=TRUE)
print(commonmolfunction[1:10])

gtex <- unlist(lapply(strsplit(df$GTEx_eQTLs, '\\|'), unique))
commongtex <- sort(table(gtex),decreasing=TRUE)
print(commongtex[1:12])

commongenes <- sort(table(df$X.6),decreasing=TRUE)
print(commongenes[1:20])

#build table used to train/test ML models

#list of samples: index 1-22 are the cases, 23-43 are the controls
samples <- c(unique(unlist(vardf$caselist)), unique(unlist(vardf$controllist)))

#matrix where rows are the samples, columns are the variants to be included in ML models
mat <- matrix(0,43,30)
vars <- c("CYP2C9 c.430C>T", "CYP2C19 c.-806C>T", "CHRNE c.-888+1742A>G", "CYP4F2 c.1297G>A",
       "ALOX5AP c.241+2284A>T", "ALOX5 c.432-4727A>G", "ALOX5 c.349+779A>G", "ALOX5 c.982-5093T>C",
       "LTC4S c.-444A>C", "HLA-DPB1 c.*1227C>T", "HLA-DQB1 c.-2013A>G", "HLA-DQB1 c.121T>C",
       "ACOX2 c.1850+225G>T", "COX20 c.42+1860A>C", "COX10 c.500-4630G>T", "COX10 c.625-5106C>T",
       "TBXA2R c.795C>T", "TBXA2R c.-84+2263A>G", "CYSLTR1 c.-114-20950C>T", "CYSLTR1 c.-115+13590C>T",
       "CYSLTR2 c.-242-11829T>C", "PCSK9 c.1420G>A", "ZFYVE19 c.1126T>G", "GP6 c.709G>A",
       "TMEM8B c.1301G>A", "TMEM175 c.1178T>C", "CLDN8 c.451T>C", "CLDN7 c.590T>C",
       "IL10RB c.139A>G", "IL32 c.-89A>C")
#vars <- sample(vars,30)
colnames(mat) <- vars
#riskscores <- replicate(43,0)

for(col in vars) {
  ls <- unlist(vardf[grepl(col, vardf$id, fixed=TRUE), c('caselist', 'controllist')])
  indx <- grepl(paste(ls,collapse="|"), samples)
  mat[indx, col] <- 1
  #riskscores[indx] = riskscores[indx] + log(vardf$oddsratio[grepl(col, vardf$id, fixed=TRUE)][1])
}
#riskscores <- exp(riskscores/30)
#print(riskscores)
logsvd_model = logisticSVD(mat, k=6)
#write.table(file="pcs.csv", logsvd_model$A)

#create matrix for samples with asthma/nasal polyps but unknown sensitivity to aspirin
df0 <- read.table(file = 'C:\\Users\\brend\\Downloads\\210509-214351_export_variant.tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df0 <- df0[-1,]
df0$id <- paste(df0$X.6, df0$X.9)
df0$X.12 <- as.numeric(df0$X.12)

df02 <- read.table(file = 'C:\\Users\\brend\\Downloads\\210509-225523_export_variant.tsv', sep = '\t', header = TRUE, fill=TRUE, quote="")
df02 <- df02[-1,]
df02$id <- paste(df02$X.6, df02$X.9)
df02$X.12 <- as.numeric(df02$X.12)

df0$sc2 <- df02$X.12[match(df0$id, df02$id)]
df0$sc2[is.na(df0$sc2)] <- 0
df0$casesum <- df0$X.12 + df0$sc2

df0$caselist <- df02$X.13[match(df0$id, df02$id)]
df0$caselist[is.na(df0$caselist)] <- ""
df0$caselist <- mapply(c, strsplit(df0$X.13, ";"), strsplit(df0$caselist, ";"))

samples2 <- c(unique(unlist(df0$caselist)))
#matrix where rows are the samples, columns are the variants to be included in ML models
mat2 <- matrix(0,23,30)
colnames(mat2) <- vars

lst = c()
for(col2 in vars) {
  ls <- unlist(df0[grepl(col2, df0$id, fixed=TRUE), 'caselist'])
  indx <- grepl(paste(ls,collapse="|"), samples2)
  mat2[indx, col2] <- 1
}
mat2[,"HLA-DQB1 c.121T>C"] <- 0
mat2[,"CLDN7 c.590T>C"] <- 0
mat2[,"TMEM175 c.1178T>C"] <- 0
mat2[,"CHRNE c.-888+1742A>G"] <- 0
mat2[,"CYP2C19 c.-806C>T"] <- 0

#print(mat2)
write.table(file="pcsOUT.csv", predict(logsvd_model, mat2, type = "PCs"))
