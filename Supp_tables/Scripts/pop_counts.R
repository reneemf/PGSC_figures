rm(list=ls())
source(file.path(Sys.getenv("PGSC_HOME", unset = "."), "config.R"))
library(gridExtra)
library(grid)
library(stringr)

setwd(paste0(fig_dir, "Supp_tables/pop_count/"))
base_cols <- c('#6699CC','#004488','#DDDDDD','#2b2d42')
col_names <- c("outer_row","header","inner_row","text")
names(base_cols) <- col_names
trans_cols <- trans_col(base_cols)
names(trans_cols) <- col_names

full_table <- matrix(data=NA,nrow=5,ncol=1)

for(context in context_list){
  table <- as.data.frame(read.table(paste0(home_dir, "pop_counts/pop_counts_", context, ".txt"), header = TRUE))
  colnames(table) <- str_to_sentence(gsub("\\.", " ", colnames(table)))
  full_table <- cbind(full_table, table[,c(2,3,1)])
}

full_table <- t(full_table[,-1])

write.csv(full_table, file = paste0("all_pop_counts.csv"))



