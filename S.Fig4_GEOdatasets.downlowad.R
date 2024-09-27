library(GEOquery)
library(tidyverse)

GEO.NO = 'GSE26511' #'GSE52903'
gset = getGEO(GEO.NO, destdir=".",getGPL = F)
gset=gset[[1]]
pdata=pData(gset)
exprSet=exprs(gset)
exprSet <- as.data.frame(exprSet)
identical(rownames(pdata),colnames(exprSet))
index = gset@annotation
GPL <-getGEO(index,destdir =".") 
GPL_anno <- Table(GPL)
probe2symbol_df <- GPL_anno %>% 
  select(1,11) 
colnames(probe2symbol_df) <- c("probe_id","symbol")
exprSet <- exprSet %>% 
  rownames_to_column("probe_id") %>% 
  inner_join(probe2symbol_df,by="probe_id") %>% 
  select(-probe_id) %>%  
  select(symbol,everything()) %>%  
  mutate(rowMean =rowMeans(.[,-1])) %>% 
  arrange(desc(rowMean)) %>% 
  distinct(symbol,.keep_all = T) %>% 
  select(-rowMean) %>% 
  filter(!is.na(symbol)) %>% 
  column_to_rownames("symbol") 
exprSet <- exprSet[rowSums(exprSet)>1,]
exprSet <- na.omit(exprSet)