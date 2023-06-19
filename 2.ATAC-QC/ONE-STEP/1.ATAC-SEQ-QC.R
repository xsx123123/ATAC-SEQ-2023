# AUTHOR:ZHANG JIAN
# DATE:2023.6.19
# VERSION:2.2.1V
# E-MAIL:1099651528@QQ.COM
# description: this is ATAC-seq QC scripts
# only for mus 
################################################################################
rm(list=ls())
################################################################################
# pre loading package
suppressMessages(library("crayon"))
suppressMessages(library("praise"))
# pre loading function
# Function print run condition style 1
print_color_note_type1 <- function(sologo){
  cat(underline(bold(cyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                               "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
  Sys.sleep(1)
}
# Function print run condition style 4
print_color_note_type3 <- function(sologo){
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                                 "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
  cat(underline(bold(cyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  Sys.sleep(1)
}
# Function check package install condition
install_package_check <- function(){
  list <- as.data.frame(installed.packages())
  # set need install package
  need_package_list <- c("TxDb.Mmusculus.UCSC.mm39.knownGene",
                         "BSgenome.Mmusculus.UCSC.mm39",
                         "crayon",
                         "praise",
                         "progress",
                         "ggplot2",
                         "ggalt",
                         "ATACseqQC",
                         "ChIPpeakAnno",
                         "base",
                         "parallel",
                         "Rsamtools",
                         "ggpubr")
  for (i in need_package_list){
    if(i %in% list$Package){
      cat("Analysis package",i,"install !!!!   ")
      cat(bold(green("PASS!!!")))
      cat("\n")
    }else{
      cat("Analysis package",i,"not install !!!!")
      cat(bold(red("FAIL!!!")))
      stop("you should install depend package!!! !!!!")
    }
    Sys.sleep(0.5)
  }
}
# check package install
print_color_note_type1("CHECK PACKAGE INSATLL DO !!! ")
install_package_check()
print_color_note_type3("CHECK PACKAGE INSATLL DONE !!! ")
################################################################################
# PATH 1 LOADING PACKAGE
suppressMessages(library("TxDb.Mmusculus.UCSC.mm39.knownGene"))
suppressMessages(library("BSgenome.Mmusculus.UCSC.mm39"))
suppressMessages(library("progress"))
suppressMessages(library("ggplot2"))
suppressMessages(library("ggpubr"))
suppressMessages(library("ggalt"))
suppressMessages(library("ATACseqQC"))
suppressMessages(library("ChIPpeakAnno"))
suppressMessages(library("base"))
suppressMessages(library("parallel"))
suppressMessages(library("Rsamtools"))
################################################################################
# PATH-2 define function
# Function1:create dir
create_dir <- function(list_dir){
  for (i in list_dir) {
    if(! dir.exists(i)){
      dir.create(i)
    }
  }
  cat("create dir finish!!!!")
}
# Function3:print run condition style 2
print_color_note_type2 <- function(sologo){
  cat("\n")
  cat("\n")
  cat(underline(bold(cyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                                "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
  cat(underline(bold(cyan("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  Sys.sleep(1)
}
# Function4:print run condition style 3
print_color_note_type2_warring <- function(sologo){
  cat("\n")
  cat("\n")
  cat(underline(bold(bgRed("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  factor <- paste0("(●´∀｀●)ﾉ ",as.character(Sys.time())," (●´∀｀●)ﾉ","\n",rep("                                "),sologo,"\n")
  cat(rep(" ",15),factor)
  cat("\n")
  cat(underline(bold(bgRed("#--------------------------------------------#####--------------------------------------------#"))))
  cat("\n")
  cat("\n")
  Sys.sleep(1)
}
# function6:GET FILE FULL NAME
get_bam_file_name <- function(bamFile){
  # bamFile <- "/data/jian/ATAC-seq-test/5.bam_file/B1C2C12siCtrl1_L3_Q0056W0150.sort.bam"
  name <- strsplit(bamFile,"/")[[1]][length(strsplit(bamFile,"/")[[1]])]
  name <- strsplit(name,".bam")[[1]][1]
  if(grepl(".sort", name, ignore.case = TRUE)){
    name <- strsplit(name,".sort")[[1]][1]
  }else{
    if (grepl(".rechrMt", name, ignore.case = TRUE)){
      name <- strsplit(name,".rechrMt")[[1]][1]
    }else{
      if (grepl(".shift", name, ignore.case = TRUE)){
        name <- strsplit(name,".shift")[[1]][1]
      }else{
        if (grepl(".fin", name, ignore.case = TRUE)){
          name <- strsplit(name,".fin")[[1]][1]
        }
      }
    }
  }
  return(name)
}
# Function7:build bam index by bamtools
build_index <- function(bam_dir,samtools_dir,num_threads){
    bam_index = paste0(samtools_dir," index -@ ",num_threads," ", bam_dir)
    cat("building index by samtools!!!!!","\n")
    system(bam_index)
}
# Function8:check BAM total QC condition
Origin_normal_QC_CHECK <- function(bamFile){
  # bamFile <- "/data/jian/ATAC-seq-test/5.bam_file/B1C2C12siCtrl1_L3_Q0056W0150.sort.bam"
  # check bai file is exit?
  if (file.exists(paste0(bamFile,".bai"))) {
    message("BAI file exists for the BAM file.")
  } else {
    message("BAI file does not exist for the BAM file.")
    build_index(bamFile,samtools_dir,num_threads)
  }
  # atac-seq QC
  ATACseqQC <- bamQC(bamFile, outPath = NULL,mitochondria="MT")
  name <- get_bam_file_name(bamFile)
  print_color_note_type2(paste0("Check ",name," ATAC-seq DO!!!"))
  # get normal
  mapping_quality <- ATACseqQC$MAPQ
  mapping_states <- ATACseqQC$idxstats
  # draw chr read rate
  # extert chr 1-MT chr READ rate
  extert_data <- mapping_states[c(1:as.numeric(rownames(mapping_states[grepl("MT",mapping_states$seqnames),]))),]
  extert_data$mapped <- (extert_data$mapped)/1000
  p <- ggplot(extert_data, aes(x=seqnames,y=mapped,fill=seqnames))+
    geom_bar(stat = 'identity') +
    scale_y_continuous(limits=c(0,max(extert_data$mapped)*1.2),n.breaks = 10,expand=c(0,0)) +
    xlab("chr ID") + ylab(expression("Read number 10"^"3")) + ggtitle(name) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size = 4),
          axis.title.x = element_text( size = 6),
          axis.title.y = element_text( size = 6),
          plot.title = element_text(hjust = 0.5,size = 6, face = "bold"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="#606c70",size=0.1,lineend = 0.04),
          axis.ticks.y=element_line(color="#606c70",size=0.1,lineend = 0.04),
          legend.position="none")
  # save plot
  ggsave(paste0("ATAC-seq"," chr read number ",name,".pdf"),plot = p,width = 6,height = 5,units="cm",device="pdf",path=origin_qc_result_dir)
  ggsave(paste0("ATAC-seq"," chr read number ",name,".png"),plot = p,width = 6,height = 5,units="cm",device="png",path=origin_qc_result_dir,dpi=2000)
  # extert other ATAC-seq QC
  other_ATAC_seq <- t(data.frame(totalQNAMEs=ATACseqQC$totalQNAMEs,
                                 duplicateRate=ATACseqQC$duplicateRate,
                                 mitochondriaRate=ATACseqQC$mitochondriaRate,
                                 properPairRate=ATACseqQC$properPairRate,
                                 unmappedRate=ATACseqQC$unmappedRate,
                                 hasUnmappedMateRate=ATACseqQC$hasUnmappedMateRate,
                                 notPassingQualityControlsRate=ATACseqQC$notPassingQualityControlsRate,
                                 nonRedundantFraction=ATACseqQC$nonRedundantFraction,
                                 PCRbottleneckCoefficient_1=ATACseqQC$PCRbottleneckCoefficient_1,
                                 PCRbottleneckCoefficient_2=ATACseqQC$PCRbottleneckCoefficient_2))
  colnames(other_ATAC_seq) <- name
  # save QC result
  write.csv(mapping_quality,file.path(origin_qc_result_dir,paste0(name,"-mapping_quality.csv")))
  write.csv(mapping_states,file.path(origin_qc_result_dir,paste0(name,"-mapping_states.csv")))
  write.csv(other_ATAC_seq,file.path(origin_qc_result_dir,paste0(name,"-other-ATAC-seq-QC.csv")))
  print_color_note_type2(paste0("Check ",name," ATAC-seq DONE!!!"))
}
# Function9:check BAM total QC condition
Clean_normal_QC_CHECK <- function(bamFile){
  # bamFile <- "/data/jian/ATAC-seq-test/5.bam_file/B1C2C12siCtrl1_L3_Q0056W0150.sort.bam"
  # check bai file is exit?
  if (file.exists(paste0(bamFile,".bai"))) {
    message("BAI file exists for the BAM file.")
  } else {
    message("BAI file does not exist for the BAM file.")
    build_index(bamFile,samtools_dir,num_threads)
  }
  # atac-seq QC
  ATACseqQC <- bamQC(bamFile, outPath = NULL,mitochondria="MT")
  name <- get_bam_file_name(bamFile)
  print_color_note_type2(paste0("Check ",name," ATAC-seq DO!!!"))
  # get normal
  mapping_quality <- ATACseqQC$MAPQ
  mapping_states <- ATACseqQC$idxstats
  # draw chr read rate
  # extert chr 1-MT chr READ rate
  extert_data <- mapping_states[c(1:as.numeric(rownames(mapping_states[grepl("Y",mapping_states$seqnames),]))),]
  extert_data$mapped <- (extert_data$mapped)/1000
  p <- ggplot(extert_data, aes(x=seqnames,y=mapped,fill=seqnames))+
    geom_bar(stat = 'identity') +
    scale_y_continuous(limits=c(0,max(extert_data$mapped)*1.2),n.breaks = 10,expand=c(0,0)) +
    xlab("chr ID") + ylab(expression("Read number 10"^"3")) + ggtitle(name) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 4),
          axis.text.y = element_text(size = 4),
          axis.title.x = element_text( size = 6),
          axis.title.y = element_text( size = 6),
          plot.title = element_text(hjust = 0.5,size = 6, face = "bold"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="#606c70",size=0.1,lineend = 0.04),
          axis.ticks.y=element_line(color="#606c70",size=0.1,lineend = 0.04),
          legend.position="none")
  # save plot
  ggsave(paste0("ATAC-seq"," chr read number ",name,".pdf"),plot = p,width = 6,height = 5,units="cm",device="pdf",path=clean_qc_result_dir)
  ggsave(paste0("ATAC-seq"," chr read number ",name,".png"),plot = p,width = 6,height = 5,units="cm",device="png",path=clean_qc_result_dir,dpi=2000)
  # extert other ATAC-seq QC
  other_ATAC_seq <- t(data.frame(totalQNAMEs=ATACseqQC$totalQNAMEs,
                                 duplicateRate=ATACseqQC$duplicateRate,
                                 mitochondriaRate=ATACseqQC$mitochondriaRate,
                                 properPairRate=ATACseqQC$properPairRate,
                                 unmappedRate=ATACseqQC$unmappedRate,
                                 hasUnmappedMateRate=ATACseqQC$hasUnmappedMateRate,
                                 notPassingQualityControlsRate=ATACseqQC$notPassingQualityControlsRate,
                                 nonRedundantFraction=ATACseqQC$nonRedundantFraction,
                                 PCRbottleneckCoefficient_1=ATACseqQC$PCRbottleneckCoefficient_1,
                                 PCRbottleneckCoefficient_2=ATACseqQC$PCRbottleneckCoefficient_2))
  colnames(other_ATAC_seq) <- name
  # save QC result
  write.csv(mapping_quality,file.path(clean_qc_result_dir,paste0(name,"-mapping_quality.csv")))
  write.csv(mapping_states,file.path(clean_qc_result_dir,paste0(name,"-mapping_states.csv")))
  write.csv(other_ATAC_seq,file.path(clean_qc_result_dir,paste0(name,"-other-ATAC-seq-QC.csv")))
  print_color_note_type2(paste0("Check ",name," ATAC-seq DONE!!!"))
}
# Function10:fragement size distribution
origin_ATAC_seq_fragsize <- function(bamFile){
  # bamFile <-  "/data/jian/ATAC-seq-test/6.rm_mt_read/B1C2C12siCtrl1_L3_Q0056W0150.rechrMt.bam"
  # check bai file is exit?
  if (file.exists(paste0(bamFile,".bai"))) {
    message("BAI file exists for the BAM file.")
  } else {
    message("BAI file does not exist for the BAM file.")
    build_index(bamFile,samtools_dir,num_threads)
  }
  name <- get_bam_file_name(bamFile)
  print_color_note_type2(paste0("Check ",name," ATAC-seq fragsize DO!!!"))
  # fragSize analysis & draw plot
  # data deal
  fragSize <- fragSizeDist(bamFile,"insertsize")
  fragSize <- as.data.frame(fragSize$insertsize)
  fragSize$Freq <- as.numeric(fragSize$Freq)
  fragSize$Var1 <- as.numeric(fragSize$Var1)
  colnames(fragSize) <- c("fragSize","Freq")
  # save fragment size
  write.csv(fragSize,file.path(origin_fragSize_save_dir,paste0(name,"-fragSize.csv")))
  fragSize_normal <- fragSize
  fragSize_normal$Freq <- (fragSize_normal$Freq)/1000
  # draw plot
  p1 <- ggplot(fragSize_normal,aes(x=fragSize,y=Freq)) + 
    geom_line(color="#29c6cd",linewidth=0.5,alpha=0.8) +
    geom_area(fill = "#29c6cd", alpha = 0.3) +
    labs(x = "Fragments length (bp)",
         y = expression("Normalized read density 10"^"3"),
         title=paste0(name," FragSize"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,(max(fragSize_normal$Freq)*1.1)),n.breaks = 8) +
    scale_x_continuous(limits = c(0,1000),n.breaks = 6) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5,size=8),
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.title.y = element_text(size=6),
          axis.line.x=element_line(linetype=1,color="black",size=0.4),       
          axis.line.y=element_line(linetype=1,color="black",size=0.4),
          panel.grid=element_blank()) +
    scale_color_brewer(palette="Paired")
  # 
  fragSize_log <- fragSize
  fragSize_log$Freq <- log10(fragSize_log$Freq)
  p2 <- ggplot(fragSize_log,aes(x=fragSize,y=Freq)) + 
    geom_line(color="#29c6cd",linewidth=0.5,alpha=0.8) +
    geom_area(fill = "#29c6cd", alpha = 0.3) +
    labs(x = "Fragments length (bp)",
         y = expression("log10(Normalized read density)"),
         title=paste0(name," FragSize"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,(max(fragSize_log$Freq)*1.1)),n.breaks = 8) +
    scale_x_continuous(limits = c(0,1000),n.breaks = 6) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5,size=8),
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.title.y = element_text(size=6),
          axis.line.x=element_line(linetype=1,color="black",size=0.4),       
          axis.line.y=element_line(linetype=1,color="black",size=0.4),
          panel.grid=element_blank()) +
    scale_color_brewer(palette="Paired")
  # merge plot
  p <- ggarrange(p1,p2,ncol=2,nrow=1,labels=c("A","B")) 
  ggsave(paste0("ATAC-seq","insert size",name,".pdf"),plot = p,width = 20,height = 8,units="cm",device="pdf",path=origin_fragSize_save_dir)
  ggsave(paste0("ATAC-seq","insert size",name,".png"),plot = p,width = 20,height = 8,units="cm",device="png",path=origin_fragSize_save_dir,dpi=2000)
  print_color_note_type2(paste0("Check ",name," ATAC-seq fragsize DONE!!!"))
}
# Function10:fragement size distribution
clean_ATAC_seq_fragsize <- function(bamFile){
  # bamFile <-  "/data/jian/ATAC-seq-test/6.rm_mt_read/B1C2C12siCtrl1_L3_Q0056W0150.rechrMt.bam"
  # check bai file is exit?
  if (file.exists(paste0(bamFile,".bai"))) {
    message("BAI file exists for the BAM file.")
  } else {
    message("BAI file does not exist for the BAM file.")
    build_index(bamFile,samtools_dir,num_threads)
  }
  name <- get_bam_file_name(bamFile)
  print_color_note_type2(paste0("Check ",name," ATAC-seq fragsize DO!!!"))
  # fragSize analysis & draw plot
  # data deal
  fragSize <- fragSizeDist(bamFile,"insertsize")
  fragSize <- as.data.frame(fragSize$insertsize)
  fragSize$Freq <- as.numeric(fragSize$Freq)
  fragSize$Var1 <- as.numeric(fragSize$Var1)
  colnames(fragSize) <- c("fragSize","Freq")
  # save fragment size
  write.csv(fragSize,file.path(origin_fragSize_save_dir,paste0(name,"-fragSize.csv")))
  fragSize_normal <- fragSize
  fragSize_normal$Freq <- (fragSize_normal$Freq)/1000
  # draw plot
  p1 <- ggplot(fragSize_normal,aes(x=fragSize,y=Freq)) + 
    geom_line(color="#29c6cd",linewidth=0.5,alpha=0.8) +
    geom_area(fill = "#29c6cd", alpha = 0.3) +
    labs(x = "Fragments length (bp)",
         y = expression("Normalized read density 10"^"3"),
         title=paste0(name," FragSize"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,(max(fragSize_normal$Freq)*1.1)),n.breaks = 8) +
    scale_x_continuous(limits = c(0,1000),n.breaks = 6) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5,size=8),
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.title.y = element_text(size=6),
          axis.line.x=element_line(linetype=1,color="black",size=0.4),       
          axis.line.y=element_line(linetype=1,color="black",size=0.4),
          panel.grid=element_blank()) +
    scale_color_brewer(palette="Paired")
  # 
  fragSize_log <- fragSize
  fragSize_log$Freq <- log10(fragSize_log$Freq)
  p2 <- ggplot(fragSize_log,aes(x=fragSize,y=Freq)) + 
    geom_line(color="#29c6cd",linewidth=0.5,alpha=0.8) +
    geom_area(fill = "#29c6cd", alpha = 0.3) +
    labs(x = "Fragments length (bp)",
         y = expression("log10(Normalized read density)"),
         title=paste0(name," FragSize"))+
    scale_y_continuous(expand = c(0,0),limits = c(0,(max(fragSize_log$Freq)*1.1)),n.breaks = 8) +
    scale_x_continuous(limits = c(0,1000),n.breaks = 6) +
    theme_classic() +
    theme(plot.title = element_text(hjust=0.5,size=8),
          axis.title.x = element_text(size=6),
          axis.text.x = element_text(size=6),
          axis.text.y = element_text(size=6),
          axis.title.y = element_text(size=6),
          axis.line.x=element_line(linetype=1,color="black",size=0.4),       
          axis.line.y=element_line(linetype=1,color="black",size=0.4),
          panel.grid=element_blank()) +
    scale_color_brewer(palette="Paired")
  # merge plot
  p <- ggarrange(p1,p2,ncol=2,nrow=1,labels=c("A","B")) 
  ggsave(paste0("ATAC-seq","insert size",name,".pdf"),plot = p,width = 20,height = 8,units="cm",device="pdf",path=clean_fragSize_save_dir)
  ggsave(paste0("ATAC-seq","insert size",name,".png"),plot = p,width = 20,height = 8,units="cm",device="png",path=clean_fragSize_save_dir,dpi=2000)
  print_color_note_type2(paste0("Check ",name," ATAC-seq fragsize DONE!!!"))
}
# Function11:BAM shift&split
shift_read_split_tss <- function(bamFile){
  # bamFile <- bam_file_list[1]
  name <- get_bam_file_name(bamFile)
  print_color_note_type2(paste0(name," ATAC-seq BAM shift & split & draw TSS PLOT DO!!!"))
  ## bamfile tags to be read in
  possibleTag <- combn(LETTERS, 2)
  possibleTag <- c(paste0(possibleTag[1, ], possibleTag[2, ]),
                   paste0(possibleTag[2, ], possibleTag[1, ]))
  bamTop100 <- scanBam(BamFile(bamFile, yieldSize = 100),
                       param = ScanBamParam(tag = possibleTag))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]
  # read bam file
  gal <- readBamFile(bamFile, tag=tags,
                     asMates=TRUE, bigFile=TRUE)
  # shift bam
  gal1 <- shiftGAlignmentsList(gal, outbam=file.path(shift_read_result_dir,paste0(name,".shift.bam")))
  # print run condition
  cat(bold(green(paste0("SHIFT BAM SAVE AT : ",file.path(shift_read_result_dir,paste0(name,".shift.bam"))))))
  # create split save dir
  split_save_dir <- file.path(split_bam_dir,name)
  create_dir(c(split_save_dir))
  # split shift bam file
  objs <- splitGAlignmentsByCut(gal1, txs=txs, genome=genome, outPath = split_save_dir)
  # print run condition
  cat(bold(green(paste0("SPLIT BAM SAVE AT : ",split_save_dir))))
  # bamFile_dir <- bam_file_list[1]
  # get split bam file path
  bamFiles <- file.path(split_save_dir,
                        c("NucleosomeFree.bam",
                          "mononucleosome.bam",
                          "dinucleosome.bam",
                          "trinucleosome.bam"))
  # get promoters range infor
  TSS <- promoters(txs, upstream=0, downstream=1)
  TSS <- unique(TSS)
  
  # 创建映射关系
  mapping <- c("chr1" = "1", "chr2" = "2", "chr3" = "3","chr4" = "4","chr5" = "5","chr6" = "6","chr7" = "7","chr8" = "8","chr9" = "9","chr10" = "10",
               "chr11" = "11", "chr12" = "12", "chr13" = "13","chr14" = "14","chr15" = "15","chr16" = "16","chr17" = "17","chr18" = "18","chr19" = "19",
               "chrX" = "X","chrY" = "Y")
  # 使用renameSeqlevels()函数重命名seqnames
  TSS <- renameSeqlevels(TSS, mapping)
  
  ## estimate the library size for normalization
  (librarySize <- estLibSize(bamFiles))
  librarySize_data <- as.data.frame(librarySize)
  # rename rownames
  for (i in 1:dim(librarySize_data)[1]){
    rownames(librarySize_data)[i] <- strsplit(rownames(librarySize_data)[i],"/")[[1]][length(strsplit(rownames(librarySize_data)[i],"/")[[1]])]
  }
  for (i in 1:dim(librarySize_data)[1]){
    rownames(librarySize_data)[i] <- strsplit(rownames(librarySize_data)[i],".bam")[[1]][1]
  }
  librarySize_data$nucleosome <- rownames(librarySize_data)
  # draw librarySize for NucleosomeFree & mononucleosome & dinucleosome & trinucleosome
  librarySize_data$nucleosome <- reorder(librarySize_data$nucleosome, -librarySize_data$librarySize)
  # draw plot
  p <- ggplot(librarySize_data,aes(x=nucleosome,y=librarySize,fill=nucleosome))+
    geom_bar(stat = 'identity',width=0.5) +
    scale_y_continuous(limits=c(0,max(librarySize_data$librarySize)*1.2),n.breaks = 10,expand=c(0,0)) +
    xlab("Nucleosome") + ylab(expression("Read number")) + ggtitle(name) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 4,hjust=0.9,vjust=0.9,angle=45),
          axis.text.y = element_text(size = 4),
          axis.title.x = element_text( size = 7),
          axis.title.y = element_text( size = 7),
          plot.title = element_text(hjust = 0.5,size = 6, face = "bold"),
          axis.line.x=element_line(linetype=1,color="black",size=0.2),       
          axis.line.y=element_line(linetype=1,color="black",size=0.2),
          axis.ticks.x=element_line(color="#606c70",size=0.1,lineend = 0.04),
          axis.ticks.y=element_line(color="#606c70",size=0.1,lineend = 0.04),
          legend.position="none")
  # save plot
  ggsave(paste0("ATAC-seq"," Nucleosome number ",name,".pdf"),plot = p,width = 6,height = 6,units="cm",device="pdf",path=TSS_dir)
  ggsave(paste0("ATAC-seq"," Nucleosome number ",name,".png"),plot = p,width = 6,height = 6,units="cm",device="png",path=TSS_dir,dpi=2000)
  ## calculate the signals around TSSs.
  NTILE <- 101
  dws <- 1010
  ups <- 1010
  seqlev <- c(1:19,"X","Y")
  sigs <- enrichedFragments(gal=objs[c("NucleosomeFree","mononucleosome","dinucleosome","trinucleosome")],
                            TSS=TSS,
                            seqlev=seqlev,
                            librarySize=librarySize,
                            TSS.filter=0.5,
                            n.tile = NTILE,
                            upstream = ups,
                            downstream = dws)
  ## log2 transformed signals
  sigs.log2 <- lapply(sigs, function(.ele) log2(.ele+1))
  #plot heatmap
  pdf(file.path(TSS_dir,paste0(name," metal TSS heatmap plot.pdf")),width=8, height=10)
  featureAlignedHeatmap(sigs.log2, reCenterPeaks(TSS, width=ups+dws),
                        zeroAt=.5, n.tile=NTILE)
  dev.off()
  # get signals normalized for nucleosome-free and nucleosome-bound regions.
  
  out <- featureAlignedDistribution(sigs,
                                    reCenterPeaks(TSS, width=ups+dws),
                                    zeroAt = .5, n.tile = NTILE, type = "l",
                                    ylab = "Averaged coverage")
  
  ## rescale the nucleosome-free and nucleosome signals to 0~1
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  out <- apply(out, 2, range01)
  pdf(file.path(TSS_dir,paste0(name," metal TSS plot.pdf")),width=6, height=4)
  matplot(out, type="l", xaxt="n",
          xlab="Position (bp)",
          ylab="Fraction of signal",
          main =paste0(name," TSS plot"))
  axis(1, at=seq(0, 100, by=10)+1,
       labels=c("-1K", seq(-800, 800, by=200), "1K"), las=2)
  abline(v=seq(0, 100, by=10)+1, lty=2, col="gray")
  dev.off()
  print_color_note_type2(paste0(name," ATAC-seq BAM shift & split & draw TSS PLOT DONE!!!"))
}
# Function12:print ProgressBar
ProgressBar <- function(){
  pb <- progress_bar$new(
    format = 'Waitting [:bar] :percent in :elapsed',
    total = 6, clear = FALSE, width = 100
  )
  for (i in 1:6) {
    pb$tick()
    Sys.sleep(0.5)
  }
}
# Function13:print run parameter
get_run_parameter <- function(){
  run_parameter <- c(paste0("ATAC-seq QC root_dir is :   ",root_dir),
                     paste0("ATAC-seq QC save dir is :   ",save_name),
                     paste0("ATAC-seq origion bam dir is :   ",origion_bam_dir),
                     paste0("ATAC-seq remove MT bam dir is :   ",rm_mt_bam),
                     paste0("ATAC-seq remove dup bam dir is :   ",rm_dup_bam),
                     paste0("ATAC-seq clean bam dir is :   ",clean_bam),
                     paste0("ATAC-seq QC use threads :   ",num_threads),
                     paste0("samtools dir is :   ",samtools_dir))
  run_parameter <- as.data.frame(run_parameter)
  print_color_note_type1("place check run parameter!!!")
  print(run_parameter,quote = FALSE, row.names = FALSE)
  print_color_note_type3("place check run parameter!!!")
  cat("\n")
  ProgressBar()
  cat("\n")
}
################################################################################
# PRINT NOTE INFOR
print_color_note_type2_warring(paste("THIS IS ",bold(green("!!! MOUSE !!!")) ," special ATAC-seq QC Rscripts"))
# PATH-3 set run parameter
args<-commandArgs(TRUE)
# PRINT RUN CONDITION
print_color_note_type1("set run parameter DO!!!")
# set root dir
root_dir <- args[1]
# set save dir 
save_name <- "13.ATAC-seq-QC"
# set origion bam dir
origion_bam_dir <- "5.bam_file"
# set remove MT chr bam dir
rm_mt_bam <- "6.rm_mt_read"
# set remove dup bam dir
rm_dup_bam <- "9.rm_dup_bam"
# set remove low quality bam dir
clean_bam <- "11.rm_low_quality_bam"
# set multithreads use core number
num_threads <- 10
# set samtools dir path
samtools_dir <- "/home/jian/biosoftware/samtools-1.16.1/samtools"
get_run_parameter()
# PRINT RUN CONDITION
print_color_note_type3("set run parameter DONE!!!")
################################################################################
# PATH-4 set save dir & create dir
# PRINT RUN CONDITION
print_color_note_type1("set save dir & create dir DO!!!")
# set save dir
save_dir <- file.path(root_dir,save_name)
# set Origin QC RESULT dir
origin_qc_result_dir <- file.path(save_dir,"Origin_QC_RESULT")
# set clean QC RESULT dir
clean_qc_result_dir <- file.path(save_dir,"Clean_QC_RESULT")
# set Origin QC fragSize RESULT dir
origin_fragSize_save_dir <- file.path(save_dir,"Origin_QC_fragSize")
# set clean QC fragSize RESULT dir
clean_fragSize_save_dir <- file.path(save_dir,"Clean_QC_fragSize")
# set shift read RESULT dir
shift_read_result_dir <- file.path(save_dir,"shift_read_RESULT")
# set split read  dir
split_bam_dir <- file.path(save_dir,"split_bam")
# set TSS  dir
TSS_dir <- file.path(save_dir,"TSS-plot")
# get save dir list
dir_list <- c(save_dir,origin_qc_result_dir,clean_qc_result_dir,
              origin_fragSize_save_dir,clean_fragSize_save_dir,
              shift_read_result_dir,split_bam_dir,TSS_dir)
# create dir
create_dir(dir_list)
# PRINT RUN CONDITION
print_color_note_type3("set save dir & create dir DONE!!!")
################################################################################
# PATH-5 ATAC-seq normal QC check
# PRINT RUN CONDITION
print_color_note_type1("ATAC-seq normal QC check DO!!!")
# get bam file full name
bam_file_list <- list.files(file.path(root_dir,origion_bam_dir),full.names=T,pattern="*.sort.bam$")
# use multthreads analysis ATAC-seq QC
mclapply(bam_file_list, Origin_normal_QC_CHECK, mc.cores = num_threads)
# get bam file full name
bam_file_list <- list.files(file.path(root_dir,clean_bam),full.names=T,pattern="*.fin.bam$")
# use multthreads analysis ATAC-seq QC
mclapply(bam_file_list, Clean_normal_QC_CHECK, mc.cores = num_threads)
# PRINT RUN CONDITION
print_color_note_type3("ATAC-seq normal QC check DONE!!!")
################################################################################
# PATH-6 ATAC-seq fragSize
# PRINT RUN CONDITION
print_color_note_type1("ATAC-seq fragSize DO!!!")
# fragSize for remove MT bam file
# get remove MT file full name
bam_file_list <- list.files(file.path(root_dir,rm_mt_bam),full.names=T,pattern="*.rechrMt.bam$")
# use mulitithreads
mclapply(bam_file_list, origin_ATAC_seq_fragsize, mc.cores = num_threads)
# get remove MT file full name
bam_file_list <- list.files(file.path(root_dir,clean_bam),full.names=T,pattern="*.fin.bam$")
# use mulitithreads
mclapply(bam_file_list, clean_ATAC_seq_fragsize, mc.cores = num_threads)
# PRINT RUN CONDITION
print_color_note_type3("ATAC-seq fragSize DONE!!!")
################################################################################
# PATH-6 ATAC-seq Splitting BAM files & draw TSS
# PRINT RUN CONDITION
print_color_note_type1("ATAC-seq Splitting BAM files & draw TSS DO!!!")
# Splitting BAM files for shift bam file
bam_file_list <- list.files(file.path(root_dir,clean_bam),full.names=T,pattern="*.fin.bam$")
# get transcription infor
txs <- transcripts(TxDb.Mmusculus.UCSC.mm39.knownGene)
# get genome infor
genome <- BSgenome.Mmusculus.UCSC.mm39
# use mulitithreads
mclapply(bam_file_list, shift_read_split_tss, mc.cores = num_threads)
# PRINT RUN CONDITION
print_color_note_type3("ATAC-seq Splitting BAM files & draw TSS DONE!!!")
################################################################################
