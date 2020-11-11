
getReadsQuality = function(inputFile = "",
                           ...){
  if(!require("pacman")) install.packages("pacman");
  pacman::p_load(crayon, magrittr, tidyverse, rjson);
  
  if(!file.exists(inputFile)){
    stop(format(Sys.time(), usetz = TRUE), red(" inputFile ", inputFile, " is empty, please check it. quit now"), "\n")
  }
  
  sampleName = str_remove(basename(inputFile), "_report.json");
  
  jsonFile = fromJSON(file = inputFile);
  
  rbind.data.frame(jsonFile$read1_before_filtering$quality_curves %>% as.data.frame() %>% 
                     mutate(cycle = 1:nrow(.)) %>% 
                     tidyr::gather(key = "base", value = "score", A:mean ) %>% 
                     mutate("reads" = "read1",
                            "status" = "before_filter"),
                   jsonFile$read2_before_filtering$quality_curves %>% as.data.frame() %>% 
                     mutate(cycle = 1:nrow(.)) %>% 
                     tidyr::gather(key = "base", value = "score", A:mean ) %>% 
                     mutate("reads" = "read2",
                            "status" = "before_filter"),
                   
                   jsonFile$read1_after_filtering$quality_curves %>% as.data.frame() %>% 
                     mutate(cycle = 1:nrow(.)) %>% 
                     tidyr::gather(key = "base", value = "score", A:mean ) %>% 
                     mutate("reads" = "read1",
                            "status" = "after_filter"),
                   jsonFile$read2_after_filtering$quality_curves %>% as.data.frame() %>% 
                     mutate(cycle = 1:nrow(.)) %>% 
                     tidyr::gather(key = "base", value = "score", A:mean ) %>% 
                     mutate("reads" = "read2",
                            "status" = "after_filter")) %>% 
    mutate(base = factor(base, levels = c("A", "T", "C", "G", "mean")),
           status = factor(status, levels = c("before_filter", "after_filter"))) %>%
    ggplot(aes(x = cycle, y = score)) +
    geom_line(aes(color = base, linetype = status), size = 1) +
    scale_color_manual(values = c("green", "red", "blue", "orange", "black")) +
    scale_linetype_manual(values = c("before_filter" = "twodash", "after_filter" = "solid"))+
    facet_wrap(~reads) +
    labs(x = "Sequencing cycle",
         y = "Qaulity score",
         title = paste("Quality scores before and After filtering", sampleName, sep = " "));
}

getInsertSizeDistribution = function(inputFile = "",
                                     ...){
  if(!require("pacman")) install.packages("pacman");
  pacman::p_load(crayon, magrittr, tidyverse, rjson);
  
  if(!file.exists(inputFile)){
    stop(format(Sys.time(), usetz = TRUE), red(" inputFile ", inputFile, " is empty, please check it. quit now"), "\n")
  }
  
  sampleName = str_remove(basename(inputFile), "_report.json");
  
  jsonFile = fromJSON(file = inputFile);
  
  truelen = 0;
  if(length(jsonFile$insert_size$histogram) >= 2){
    hislen = length(jsonFile$insert_size$histogram);
    while(hislen >= 2){
      if(jsonFile$insert_size$histogram[hislen] != 0){
        truelen = hislen;
        break;
      }
      hislen = hislen -1;
    }
  } else {
    truelen = length(jsonFile$insert_size$histogram);
  }
  
  if(truelen == 0){
    stop(format(Sys.time(), usetz = TRUE), red(" histogram can not be made. quit now"), "\n")
  }
  
  
  jsonFile$insert_size %>% 
    as.data.frame() %>%
    dplyr::slice(1:truelen) %>% 
    mutate(base = 1:truelen,
           histogram2 = histogram *100 / sum(.$histogram)) %>% 
    ggplot(aes(x = base, y = histogram2)) +
    geom_col() +
    labs(x = "Insert size (bp)", y = "Read percentage (%)",
         title = paste("Insert size distribution", sampleName),
         caption = paste0("Insert size is estmiated on overlap analysis of paired-end reads. The insert size peak is ", jsonFile$insert_size$peak, " bp with ",
                          format(round(jsonFile$insert_size$unknown * 100 / sum(jsonFile$insert_size$histogram), 2), nsmall = 2),
                          "% reads of not overllapped"));
}

getDuplicationLevel = function(inputFile = "",
                                     ...){
  if(!require("pacman")) install.packages("pacman");
  pacman::p_load(crayon, magrittr, tidyverse, rjson);
  
  if(!file.exists(inputFile)){
    stop(format(Sys.time(), usetz = TRUE), red(" inputFile ", inputFile, " is empty, please check it. quit now"), "\n")
  }
  
  sampleName = str_remove(basename(inputFile), "_report.json");
  
  jsonFile = fromJSON(file = inputFile);
  
  jsonFile$duplication %>% 
    as.data.frame() %>% 
    mutate("Read_percentage" = histogram * 100 / sum(.$histogram),
           mean_gc = mean_gc * 100) %>% 
    mutate(base = 1:nrow(.)) %>% 
    dplyr::select(base, mean_gc, Read_percentage) %>% 
    tidyr::gather(key = "key", value = "value", mean_gc, Read_percentage) %>% 
    ggplot(aes(x = base, y = value, color = key)) +
    geom_line() +
    labs(x = "Duplication level", y = "Read percentage(%) / GC ratio",
         title = paste0("Duplication rate (", sample02_json$duplication$rate * 100, "%) ", sampleName),
         caption = "")
}
