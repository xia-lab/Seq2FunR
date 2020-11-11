seq2fun = function(suffixNameR1 = "_R1.fastq.gz",
                   suffixNameR2 = "",
                   outputDir = "outputSeq2fun",
                   genemap,
                   tfmi,
                   mode = "tGREEDY",
                   mismatch = 2,
                   minscore = 65,
                   minlength = 11,
                   maxtranslength = 65,
                   outputMappedCleanReads = FALSE,
                   profiling = FALSE,
                   nThreads = 8,
                   verbose = TRUE,
                   ...){
  
  if(!require("pacman")) install.packages("pacman");
  pacman::p_load(crayon, magrittr, stringr);
  
  fileName = c();
  fileNamesLeft = c();
  fileNamesRight = c();
  paired = TRUE;
  nSamples = 0;
  sampleTalbe = data.frame();
  
  if(suffixNameR2 == ""){
    cat(format(Sys.time(), usetz = TRUE), yellow(" suffixNameR2 is not specified, so is going to process your single end reads.\n
                                                 If your reads are paired end reads, please specify suffixNameR2 (eg. _R2.fastq.gz)"), "\n");
    paired = FALSE;
  } else {
    cat(format(Sys.time(), usetz = TRUE), yellow(" suffixNameR2 is specified, so is going to process your paired end reads"), "\n");
    paired = TRUE;
  }
  
  #check input dir;
    #get file infor
    fileNames = list.files();
    if(length(fileNames) == 0){
      stop(format(Sys.time(), usetz = TRUE), red(" current Dir is empty, please make sure there are reads files here, quit now"), "\n")
    }
    fileNamesLeft = fileNames[str_detect(fileNames, suffixNameR1)] %>% str_remove(suffixNameR1);
    
    if(length(fileNamesLeft) == 0){
      stop(format(Sys.time(), usetz = TRUE), red(" forward Reads in current Dir is not detected, please make sure there are reads files here or the suffixNameR1 ", suffixNameR1, " is correct. quit now"), "\n")
    }
    
    if(paired){
      fileNamesRight = fileNames[str_detect(fileNames, suffixNameR2)] %>% str_remove(suffixNameR2);
      if(length(fileNamesRight) == 0){
        stop(format(Sys.time(), usetz = TRUE), red(" reverse Reads in current Dir is not detected, please make sure there are reads files here or the suffixNameR2 ", suffixNameR2, " is correct. quit now"), "\n")
      } else {
        if(all(fileNamesLeft == fileNamesRight)){
          nSamples = length(fileNamesLeft);
          cat(format(Sys.time(), usetz = TRUE), yellow(" all read files (", nSamples, ") are paired\n"));
          fileLeftInput = paste0(fileNamesLeft, suffixNameR1);
          fileRightInput = paste0(fileNamesRight, suffixNameR2);
          sampleTalbe = data.frame("sampleNames" = file.path(outputDir, fileNamesLeft),
                                   "fileLeft" = fileLeftInput,
                                   "fileRight" = fileRightInput);
        } else {
          cat(format(Sys.time(), usetz = TRUE));
          stop(red(" there are ", length(fileNamesLeft), " forward Reads and ", length(fileNameRight), " reverse Reads; please check and ensure pairedend reads are paired"), "\n");
        }
      }
    } else {
      nSamples = length(fileNamesLeft);
      fileLeftInput = paste0(fileNamesLeft, suffixNameR1);
      cat(format(Sys.time(), usetz = TRUE), yellow(" detected ", nSamples, " single end samples"), "\n");
      sampleTalbe = data.frame("sampleNames" = file.path(outputDir, fileNamesLeft),
                               "fileLeft" = fileLeftInput);
    }
    
    write.table(sampleTalbe, 
                "sample.txt",
                col.names = F,
                row.names = F,
                quote = F,
                sep = "\t");
 
    
    
    
  # if(!dir.exists(inputDir)){
  #   stop(format(Sys.time(), usetz = TRUE), red(" inputDir ", inputDir, " does not exist, quit now"), "\n");
  # } else {
  #   #get file infor
  #   fileNames = list.files(inputDir);
  #   if(length(fileName) == 0){
  #     stop(format(Sys.time(), usetz = TRUE), red(" inputDir ", inputDir, " is empty, please make sure there are reads files here, quit now"), "\n")
  #   }
  #   fileNamesLeft = fileNames[str_detect(fileNames, suffixNameR1)] %>% str_remove(suffixNameR1);
  #   
  #   if(length(fileNameLeft) == 0){
  #     stop(format(Sys.time(), usetz = TRUE), red(" forward Reads R1 in inputDir ", inputDir, " is not detected, please make sure there are reads files here or the suffixNameR1 ", suffixNameR1, " is correct. quit now"), "\n")
  #   }
  #   
  #   if(paired){
  #     fileNamesRight = fileNames[str_detect(fileNames, suffixNameR2)] %>% str_remove(suffixNameR2);
  #     if(length(fileNamesRight) == 0){
  #       stop(format(Sys.time(), usetz = TRUE), red(" reverse Reads R2 in inputDir ", inputDir, " is not detected, please make sure there are reads files here or the suffixNameR2 ", suffixNameR2, " is correct. quit now"), "\n")
  #     } else {
  #       if(all(fileNamesLeft == fileNamesRight)){
  #         cat(format(Sys.time(), usetz = TRUE), yellow(" all read files (", nSamples, ") are paired\n"));
  #         fileLeftInput = paste0(fileNamesLeft, suffixNameR1);
  #         fileRightInput = paste0(fileNamesRight, suffixNameR2);
  #         
  #         write.table
  #       }
  #     }
  #   } else {
  #     nSamples = length(fileNamesLeft);
  #     fileLeftInput = paste0(fileNamesLeft, suffixNameR1);
  #     cat(format(Sys.time(), usetz = TRUE), yellow(" detected ", nSamples, " single end samples"), "\n");
  #     sampleTalbe = data.frame("sampleNames" = file.path(outputDir, fileNamesLeft),
  #                              "fileLeft" = fileLeftInput);
  #     write.table(sampleTalbe, 
  #                 "sample.txt",
  #                 col.names = F,
  #                 row.names = F,
  #                 quote = F,
  #                 sep = "\t");
  #   }
  # }
    
    if(!file.exists(genemap)){
      stop(red(" genemap ", genemap, " does not exit; please check and ensure the path and file are correct"), "\n");
    }
    
    if(!file.exists(tfmi)){
      stop(red(" tfmi ", tfmi, " does not exit; please check and ensure the path and file are correct"), "\n");
    }
    
    seq2funr(sampletable = "sample.txt",
             genemap = genemap,
             tfmi = tfmi,
             mode = mode,
             mismatch = mismatch,
             minscore = minscore,
             minlength = minlength,
             maxtranslength = maxtranslength,
             nThreads = nThreads,
             verbose = verbose,
             outputMappedCleanReads = outputMappedCleanReads,
             profiling = profiling);

}