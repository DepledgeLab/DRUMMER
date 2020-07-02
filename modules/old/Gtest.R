args<-commandArgs(TRUE)
input_path <-args[1]
output_path <- args[2]
#print(input_path)
#print(output_path)
#path = '/Users/mac/Desktop/Nanopore_prediction/adeno-m6A-script/[redux]_VA.m6A.complete (1).csv'
df <- read.csv(file = input_path,sep='\t')
library(RVAideMemoire)

for (row in 1:nrow(df)) {
  unmodifed_list <- c(df$A_unmod[row],df$C_unmod[row],df$G_unmod[row],df$T_unmod[row],df$N_unmod[row])
  modified_list <- c(df$A_mod[row],df$C_mod[row],df$G_mod[row],df$T_mod[row],df$N_mod[row])
  resulting_gtest <- G.test(as.matrix(cbind(unmodifed_list,modified_list)))
  #print(resulting_gtest)
  G_value <- resulting_gtest$statistic
  p_val <- resulting_gtest$p.value
  df$G[row] <- G_value
  df$p_value[row] <- p_val
  df$padj[row] <- p_val * nrow(df)
  }
#print(head(df))

write.csv(df,output_path,row.names = FALSE)
#out_path <-'/Users/mac/Desktop/Nanopore_prediction/adeno-m6A-script/[redux]_VA.m6A.complete_Gtest.csv'
#write.csv(df,output_path,row.names = FALSE)
