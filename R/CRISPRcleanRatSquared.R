#' @import ggplot2
#' @import ggrepel
#' @import Matrix
#' @import CRISPRcleanR
#' @import RSpectra
#' @import magrittr
#' @import dplyr
#' @import readr
#' @import tibble
#' @import stringr
#' @import CVXR
#' @import S4Vectors
#' @importFrom readxl read_xlsx
#' @importFrom GenomicRanges findOverlaps 
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ggpubr ggarrange
#' @importFrom igraph clusters
#' @importFrom igraph graph_from_edgelist
NULL

#' Load matched files for dual KO
#'
#' get_input_data() load files as indicated in param_file or param_list. 
#' Must include count1_file, library1_file and copy_number_file. 
#' Multiple count and library files can be passed (increasing index order). The output will keep separately the batches 
#' Additional used inputs are input_fold and out_fold (default = "./") and CL_name (default = "./")
#'
#' @param param_file .tsv file, the first column indicates the type of data, the
#'   second the file location. Can include multiple files that contain the name
#'   count or library BUT the length must coincide.
#' @param param_list list with input files' location to be passed. The names in the list must include count1_file, library1_file and copy_number_file
#'
#' @return a list 
#' - CNA: copy number, one row per gene 
#' - count: list of dual KO count (per batch), one row per guides combination 
#' - library: list of metadata for dual KO library (per batch)
#' - CL_name: cell line name
#' - out_fold: location to store results
#' @export
#' 
get_input_data <- function(
  param_file = NULL,
  param_list = NULL
  ) {
  
  if (is.null(param_file) & is.null(param_list)) {
    stop("one among param_file and param_list must be not NULL")
  }
  
  if (!is.null(param_file)) {
    input_info <- suppressWarnings(readr::read_table(
      param_file, 
      col_names = FALSE, 
      show_col_types = FALSE))
    
    for (i in 1:nrow(input_info)) {
      assign(input_info$X1[i], input_info$X2[i])    
    }
    count_files <- input_info$X1[grepl("count", input_info$X1)]
    library_files <- input_info$X1[grepl("library", input_info$X1)]
  }else{
    input_info <- param_list
    for (i in 1:length(input_info)) {
      assign(names(input_info)[i], input_info[[i]])    
    }
    count_files <- names(input_info)[grepl("count", names(input_info))]
    library_files <- names(input_info)[grepl("library", names(input_info))]
  }
  
  if (length(count_files) != length(library_files)) {
    stop("library and count files MUST match")
  }
  
  # check the correct input files were passed
  if (!exists("input_fold", inherits = F)) {input_fold <- "./"}
  if (!exists("out_fold", inherits = F)) {out_fold <- "./"}
  if (!exists("CL_name", inherits = F)) {CL_name <- ""}
  if (!exists("copy_number_file", inherits = F) | 
      !exists("count1_file", inherits = F) |
      !exists("library1_file", inherits = F)) {
    stop("Input MUST include copy_number_file, count1_file and library1_file entries")
  }
  
  n_files <- length(count_files)
  dual_library_list <- list()
  dual_count_list <- lapply(count_files, function(x)
    readr::read_table(sprintf('%s%s', input_fold, get(x)), 
                      col_types = readr::cols(.default = "?",
                                              sgRNA1_ID = "c", 
                                              sgRNA2_ID = "c"), 
                                show_col_types = FALSE))
  
  for (idx_file in seq_len(n_files)) {
    
    # load library and specify col types:
    if (grepl(".xls", library_files[idx_file])) {
      names_col_lib <- names(readxl::read_xlsx(
        path = sprintf('%s%s', input_fold, get(library_files[idx_file])),
        n_max = 0))

      col_types_lib <- ifelse(grepl("Chr", names_col_lib) | grepl("WGE_ID", names_col_lib),
                              "text", "guess")

      dual_library_list[[idx_file]] <- readxl::read_xlsx(
        path = sprintf('%s%s', input_fold, get(library_files[idx_file])),
        col_types = col_types_lib)
    }else{
      
      dual_library_list[[idx_file]] <- readr::read_tsv(
        sprintf('%s%s', input_fold, get(library_files[idx_file])),
        col_types = readr::cols(.default = "?", 
                                sgRNA1_Chr = "c", 
                                sgRNA2_Chr = "c", 
                                sgRNA1_WGE_ID = "c", 
                                sgRNA2_WGE_ID = "c"), 
        show_col_types = FALSE)
      }
    
    print(sprintf("############ load file n. %i ############", idx_file))
    
    dual_count_list[[idx_file]] <- dual_count_list[[idx_file]] %>%
      dplyr::mutate(ID = paste0(ID, "_lib", idx_file)) %>%
      dplyr::mutate(ID_lib = ID, .after = ID) %>%
      dplyr::mutate(lib = paste0("lib", idx_file), .after = ID_lib)
      
    dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>%
      dplyr::filter(!is.na(ID)) %>%
      dplyr::mutate(ID = paste0(ID, "_lib", idx_file))
    
    if (!identical(dual_count_list[[idx_file]]$ID, dual_library_list[[idx_file]]$ID)) {
      stop("count and library rows must match")
    }
    
    id_change <- grep("CTRL", dual_library_list[[idx_file]]$sgRNA1_WGE_ID)
    if (length(id_change) > 0) {
      print(sprintf("change CTRL to NONTARGET in position 1 for %i pairs", length(id_change)))
      dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_change] <- str_replace(
        string =  dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_change], 
        pattern = "CTRL", 
        replacement = "NONTARGET")
    }
    
    id_change <- grep("CTRL", dual_library_list[[idx_file]]$sgRNA2_WGE_ID)
    if (length(id_change) > 0) {
      print(sprintf("change CTRL to NONTARGET in position 2 for %i pairs", length(id_change)))
      dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_change] <- str_replace(
        string =  dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_change], 
        pattern = "CTRL", 
        replacement = "NONTARGET")
    }
    
    id_nontarget <- grep("NONTARGET", dual_count_list[[idx_file]]$Gene1)
    if (length(id_nontarget) > 0) {
      print(sprintf("assign missing sgRNA1_ID for %i NONTARGET pairs", length(id_nontarget)))
      dual_count_list[[idx_file]]$sgRNA1_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene1[id_nontarget]
      dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene1[id_nontarget]
    }
    
    id_nontarget <- grep("NONTARGET", dual_count_list[[idx_file]]$Gene2)
    if (length(id_nontarget) > 0) {
      print(sprintf("assign missing sgRNA2_ID for %i NONTARGET pairs", length(id_nontarget)))
      dual_count_list[[idx_file]]$sgRNA2_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene2[id_nontarget]
      dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene2[id_nontarget]
    }
    
    # if NA in count ID, use gene ID
    id_na <- which(is.na(dual_count_list[[idx_file]]$sgRNA1_ID))
    if (length(id_na) > 0) {
      print(sprintf("assign missing sgRNA1_ID for %i pairs", length(id_na)))
      dual_count_list[[idx_file]]$sgRNA1_ID[id_na] <- dual_count_list[[idx_file]]$Gene1[id_na]
      dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_na] <- dual_count_list[[idx_file]]$Gene1[id_na]
      dual_library_list[[idx_file]]$sgRNA1_Approved_Symbol[id_na] <- dual_count_list[[idx_file]]$Gene1[id_na]
    }
    
    id_na <- which(is.na(dual_count_list[[idx_file]]$sgRNA2_ID))
    if (length(id_na) > 0) {
      print(sprintf("assign missing sgRNA2_ID for %i pairs", length(id_na)))
      dual_count_list[[idx_file]]$sgRNA2_ID[id_na] <- dual_count_list[[idx_file]]$Gene2[id_na]
      dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_na] <- dual_count_list[[idx_file]]$Gene2[id_na]
      dual_library_list[[idx_file]]$sgRNA2_Approved_Symbol[id_na] <- dual_count_list[[idx_file]]$Gene2[id_na]
    }
    
    id_na_lib <- which(is.na(dual_library_list[[idx_file]]$sgRNA1_Approved_Symbol))
    if (length(id_na_lib) > 0) {
      print(sprintf("assign missing Gene1 name in library for %i pairs", length(id_na_lib)))
      dual_library_list[[idx_file]]$sgRNA1_Approved_Symbol[id_na_lib] <- dual_count_list[[idx_file]]$Gene1[id_na_lib]
    }
    
    id_na_lib <- which(is.na(dual_library_list[[idx_file]]$sgRNA2_Approved_Symbol))
    if (length(id_na_lib) > 0) {
      print(sprintf("assign missing Gene2 name in library for %i pairs", length(id_na_lib)))
      dual_library_list[[idx_file]]$sgRNA2_Approved_Symbol[id_na_lib] <- dual_count_list[[idx_file]]$Gene2[id_na_lib]
    }
    
    # same sequence with multiple ID, replace
    dual_library_seq <- data.frame(
      GI_ID = rep(dual_library_list[[idx_file]]$ID, 2),
      ID = c(dual_library_list[[idx_file]]$sgRNA1_WGE_ID, 
             dual_library_list[[idx_file]]$sgRNA2_WGE_ID), 
      GENES = c(dual_library_list[[idx_file]]$sgRNA1_Approved_Symbol, 
                dual_library_list[[idx_file]]$sgRNA2_Approved_Symbol), 
      CHR = c(dual_library_list[[idx_file]]$sgRNA1_Chr, 
                dual_library_list[[idx_file]]$sgRNA2_Chr), 
      START = c(dual_library_list[[idx_file]]$sgRNA1_Start, 
                dual_library_list[[idx_file]]$sgRNA2_Start), 
      END = c(dual_library_list[[idx_file]]$sgRNA1_End, 
                dual_library_list[[idx_file]]$sgRNA2_End), 
      LIBRARY = c(dual_library_list[[idx_file]]$sgRNA1_Library, 
                  dual_library_list[[idx_file]]$sgRNA2_Library), 
      SEQ = c(dual_library_list[[idx_file]]$sgRNA1_WGE_Sequence, 
              dual_library_list[[idx_file]]$sgRNA2_WGE_Sequence)
    ) %>%
      dplyr::mutate(CHR_START_END_SEQ = paste(CHR, START, END, SEQ, sep = "_"))
    
    multiple_ID <- dual_library_seq %>% 
      dplyr::group_by(CHR_START_END_SEQ) %>%
      dplyr::summarise(
        SEQ = unique(SEQ),
        n_ID = length(unique(ID)),
        ID_all = paste0(sort(unique(ID)), collapse = ","), 
        n_gene = length(unique(GENES)),
        gene_all = paste0(sort(unique(GENES)), collapse = ",")) %>%
      dplyr::filter(n_gene > 1)
    
    if (nrow(multiple_ID) > 0) {
      
      print(sprintf("reassign %i sgRNA1 or sgRNA2 with multiple genes/IDs but same CHR_START_END_SEQ", 
                    nrow(multiple_ID)))
      
      for (idx in seq_len(nrow(multiple_ID))) {
        
        # sgRNA1
        id_match <- which(dual_library_list[[idx_file]]$sgRNA1_WGE_Sequence == multiple_ID$SEQ[idx])
        if (length(id_match) > 0) {
          dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
          dual_library_list[[idx_file]]$sgRNA1_Approved_Symbol[id_match] <- multiple_ID$gene_all[idx] 
          dual_count_list[[idx_file]]$sgRNA1_ID[id_match] <- multiple_ID$ID_all[idx]   
          dual_count_list[[idx_file]]$Gene1[id_match] <- multiple_ID$gene_all[idx] 
        }
        # sgRNA2
        id_match <- which(dual_library_list[[idx_file]]$sgRNA2_WGE_Sequence == multiple_ID$SEQ[idx])
        if (length(id_match) > 0) {
          dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
          dual_library_list[[idx_file]]$sgRNA2_Approved_Symbol[id_match] <- multiple_ID$gene_all[idx] 
          dual_count_list[[idx_file]]$sgRNA2_ID[id_match] <- multiple_ID$ID_all[idx]   
          dual_count_list[[idx_file]]$Gene2[id_match] <- multiple_ID$gene_all[idx]   
        }
      }
    }
    
    
    dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>% 
      dplyr::mutate(COMB_ID =  paste(sgRNA1_WGE_Sequence, sgRNA2_WGE_Sequence, sep = "_"))
    
    if ("MyNotes" %in% colnames(dual_library_list[[idx_file]])) {
      dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>% 
        dplyr::select(-MyNotes)
    }
    
    # chr X -> 23, chr Y -> 24
    dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>% 
      dplyr::mutate(
        sgRNA1_Chr = dplyr::case_when(
          sgRNA1_Chr == "X" ~ "23",
          sgRNA1_Chr == "Y" ~ "24",
          # .default = as.character(sgRNA1_Chr)), 
          TRUE ~ as.character(sgRNA1_Chr)), 
        sgRNA2_Chr = dplyr::case_when(
          sgRNA2_Chr == "X" ~ "23",
          sgRNA2_Chr == "Y" ~ "24",
          TRUE ~ as.character(sgRNA2_Chr))
          #.default = as.character(sgRNA2_Chr))
        ) %>%
      dplyr::mutate(
        sgRNA1_Chr = as.numeric(sgRNA1_Chr), 
        sgRNA2_Chr = as.numeric(sgRNA2_Chr)
        )
  }
  
  CNA <- readr::read_table(sprintf('%s', copy_number_file), 
                         show_col_types = FALSE) %>%
    dplyr::mutate(
      CHROM = dplyr::case_when(
      CHROM == "chrX" ~ "chr23",
      CHROM == "chrY" ~ "chr24",
      TRUE ~ as.character(CHROM))
      # .default = as.character(CHROM))
    )
  
  if ("Sampleid" %in% colnames(CNA)) {
    CL_name_unif <- toupper(str_replace_all(CL_name, "[-|_]", ""))
    CNA <- CNA %>% 
      dplyr::filter(toupper(str_replace_all(Sampleid, "[-|_]", "")) %in% CL_name_unif)
    if (nrow(CNA) == 0) {
      stop(sprintf("no CN info available for %s", CL_name))
    }
  }
    
  return(list(CNA = CNA, 
              count = dual_count_list, 
              library =  dual_library_list, 
              CL_name = CL_name, 
              out_fold = out_fold))
  
}

#' Load matched files for dual KO
#'
#' get_input_data.v1() load files as indicated in param_file or param_list. 
#' Must include result_file, library_file and copy_number_file (tab separated). 
#' Only one file per type is passed (suitable for different batches previously combined in a unique file)
#' The result_file must have the first 13 column with the following structure:
#' - ID: guide pair ID string inside the batch (if any)
#' - lib: batch ID string
#' - ID_lib: combination from the columns ID and lib separated by an underscore
#' - Note1: class for guide pairs (e.g. PositiveControls, AnchorSingletons, etc)
#' - Note2: additional class for guide pairs (could be more granular like ESSENTIAL-NONTARGET)
#' - Gene_pair: HGNC symbol of the genes pair targeted by guide in position 1 and position 2 (separated by ~) 
#' - Gene1: HGNC symbol of the gene targeted by guide in position 1
#' - Gene2: HGNC symbol of the gene targeted by guide in position 2
#' - sgRNA1_WGE_ID: sgRNA ID for guide in position 1
#' - sgRNA1_WGE_Sequence: sgRNA sequence for guide in position 1
#' - sgRNA2_WGE_ID: sgRNA ID for guide in position 2
#' - sgRNA2_WGE_Sequence: sgRNA sequence for guide in position 2
#' - SEQ_pair: sgRNA sequence for guide in position 1 and 2 separated by ~
#' The remaining columns include log fold changes per cell line (<CL_name>_logFC)
#' *IMPORTANT*: the dual guides with the same SEQ_pair across batches (lib) will be merged taking the mean.
#' Additional used inputs are input_fold and out_fold (default = "./") and CL_name (default = "./")
#'
#' @param param_file .tsv file, the first column indicates the type of data, the
#' second the file location. 
#' @param param_list list with input files' location to be passed. The names in the list must include result_file, library_file and copy_number_file
#'
#'
#' @return a list 
#' - CNA: copy number, one row per gene 
#' - result: data.frame of dual KO logFC for the CL_name considered, one row per guides combination (unique identifier based on SEQ_pair)
#' - library: data.frame of metadata for dual KO library
#' - CL_name: cell line name
#' - out_fold: location to store results
#' 
#' @export
#'
#' @examples
get_input_data.v1 <- function(
  param_file = NULL,
  param_list = NULL
) {
  
  if (is.null(param_file) & is.null(param_list)) {
    stop("one among param_file and param_list must be not NULL")
  }
  
  if (!is.null(param_file)) {
    # read param_file which contains data location
    input_info <- suppressWarnings(readr::read_table(
      param_file, 
      col_names = FALSE, 
      show_col_types = FALSE))
    
    for (i in 1:nrow(input_info)) {
      assign(input_info$X1[i], input_info$X2[i])    
    }
    
  }else{
    # data location passed as input
    input_info <- param_list
    for (i in 1:length(input_info)) {
      assign(names(input_info)[i], input_info[[i]])    
    }
  }
  
  # check the correct input files were passed
  if (!exists("input_fold", inherits = F)) {input_fold <- "./"}
  if (!exists("out_fold", inherits = F)) {out_fold <- "./"}
  if (!exists("CL_name", inherits = F)) {CL_name <- ""}
  if (!exists("copy_number_file", inherits = F) | 
      !exists("result_file", inherits = F) |
      !exists("library_file", inherits = F)) {
    stop("Input MUST include copy_number_file, result_file and library_file entries")
  }
  
  dual_result <- readr::read_tsv(sprintf('%s%s', input_fold, result_file), 
                                 col_types = readr::cols(.default = "?",
                                                         sgRNA1_WGE_ID = "c", 
                                                         sgRNA2_WGE_ID = "c", 
                                                         sgRNA1_Chr = "c", 
                                                         sgRNA2_Chr = "c"), 
                                 show_col_types = FALSE)
  
  dual_library <- readr::read_tsv(sprintf('%s%s', input_fold, library_file), 
                                  col_types = readr::cols(.default = "?",
                                                          sgRNA1_WGE_ID = "c", 
                                                          sgRNA2_WGE_ID = "c", 
                                                          sgRNA1_Chr = "c", 
                                                          sgRNA2_Chr = "c"), 
                                  show_col_types = FALSE)
  
  if (!identical(dual_result$ID_lib, dual_library$ID_lib)) {
    warning("count and library rows must match, reorder library file")
    common_ID_lib <- intersect(dual_result$ID_lib, dual_library$ID_lib)
    dual_result <- dual_result[match(common_ID, dual_result$ID_lib),]
    dual_library <- dual_library[match(common_ID, dual_library$ID_lib),]
  }
  
  id_change <- grep("CTRL", dual_library$sgRNA1_WGE_ID)
  if (length(id_change) > 0) {
    print(sprintf("change CTRL to NONTARGET in position 1 for %i pairs", length(id_change)))
    dual_library$sgRNA1_WGE_ID[id_change] <- str_replace(
      string =  dual_library$sgRNA1_WGE_ID[id_change], 
      pattern = "CTRL", 
      replacement = "NONTARGET")
  }
  
  id_change <- grep("CTRL", dual_library$sgRNA2_WGE_ID)
  if (length(id_change) > 0) {
    print(sprintf("change CTRL to NONTARGET in position 2 for %i pairs", length(id_change)))
    dual_library$sgRNA2_WGE_ID[id_change] <- str_replace(
      string =  dual_library$sgRNA2_WGE_ID[id_change], 
      pattern = "CTRL", 
      replacement = "NONTARGET")
  }
  
  id_nontarget <- grep("NONTARGET", dual_result$Gene1)
  if (length(id_nontarget) > 0) {
    print(sprintf("assign missing sgRNA1_ID for %i NONTARGET pairs", length(id_nontarget)))
    dual_result$sgRNA1_WGE_ID[id_nontarget] <- dual_result$Gene1[id_nontarget]
    dual_library$sgRNA1_WGE_ID[id_nontarget] <- dual_result$Gene1[id_nontarget]
  }
  
  id_nontarget <- grep("NONTARGET", dual_result$Gene2)
  if (length(id_nontarget) > 0) {
    print(sprintf("assign missing sgRNA2_ID for %i NONTARGET pairs", length(id_nontarget)))
    dual_result$sgRNA2_WGE_ID[id_nontarget] <- dual_result$Gene2[id_nontarget]
    dual_library$sgRNA2_WGE_ID[id_nontarget] <- dual_result$Gene2[id_nontarget]
  }
  
  # if NA in count ID, use gene ID
  id_na <- which(is.na(dual_result$sgRNA1_WGE_ID))
  if (length(id_na) > 0) {
    print(sprintf("assign missing sgRNA1_ID for %i pairs", length(id_na)))
    dual_result$sgRNA1_WGE_ID[id_na] <- dual_result$Gene1[id_na]
    dual_library$sgRNA1_WGE_ID[id_na] <- dual_result$Gene1[id_na]
    dual_library$sgRNA1_Approved_Symbol[id_na] <- dual_result$Gene1[id_na]
  }
  
  id_na <- which(is.na(dual_result$sgRNA2_WGE_ID))
  if (length(id_na) > 0) {
    print(sprintf("assign missing sgRNA2_ID for %i pairs", length(id_na)))
    dual_result$sgRNA2_WGE_ID[id_na] <- dual_result$Gene2[id_na]
    dual_library$sgRNA2_WGE_ID[id_na] <- dual_result$Gene2[id_na]
    dual_library$sgRNA2_Approved_Symbol[id_na] <- dual_result$Gene2[id_na]
  }
  
  id_na_lib <- which(is.na(dual_library$sgRNA1_Approved_Symbol))
  if (length(id_na_lib) > 0) {
    print(sprintf("assign missing Gene1 name in library for %i pairs", length(id_na_lib)))
    dual_library$sgRNA1_Approved_Symbol[id_na_lib] <- dual_result$Gene1[id_na_lib]
  }
  
  id_na_lib <- which(is.na(dual_library$sgRNA2_Approved_Symbol))
  if (length(id_na_lib) > 0) {
    print(sprintf("assign missing Gene2 name in library for %i pairs", length(id_na_lib)))
    dual_library$sgRNA2_Approved_Symbol[id_na_lib] <- dual_result$Gene2[id_na_lib]
  }
  
  # same sequence with multiple ID, replace
  dual_library_seq <- data.frame(
    GI_ID = rep(dual_library$ID, 2),
    ID = c(dual_library$sgRNA1_WGE_ID, 
           dual_library$sgRNA2_WGE_ID), 
    GENES = c(dual_library$sgRNA1_Approved_Symbol, 
              dual_library$sgRNA2_Approved_Symbol), 
    CHR = c(dual_library$sgRNA1_Chr, 
            dual_library$sgRNA2_Chr), 
    START = c(dual_library$sgRNA1_Start, 
              dual_library$sgRNA2_Start), 
    END = c(dual_library$sgRNA1_End, 
            dual_library$sgRNA2_End), 
    LIBRARY = c(dual_library$sgRNA1_Library, 
                dual_library$sgRNA2_Library), 
    SEQ = c(dual_library$sgRNA1_WGE_Sequence, 
            dual_library$sgRNA2_WGE_Sequence)
  ) %>%
    dplyr::mutate(CHR_START_END_SEQ = paste(CHR, START, END, SEQ, sep = "_"))
  
  multiple_ID <- dual_library_seq %>% 
    dplyr::group_by(CHR_START_END_SEQ) %>%
    dplyr::summarise(
      SEQ = unique(SEQ),
      n_ID = length(unique(ID)),
      ID_all = paste0(sort(unique(ID)), collapse = ","), 
      n_gene = length(unique(GENES)),
      gene_all = paste0(sort(unique(GENES)), collapse = ",")) %>%
    dplyr::filter(n_gene > 1)
  
  if (nrow(multiple_ID) > 0) {
    
    print(sprintf("reassign %i sgRNA1 or sgRNA2 with multiple genes/IDs but same CHR_START_END_SEQ", 
                  nrow(multiple_ID)))
    
    for (idx in seq_len(nrow(multiple_ID))) {
      
      # sgRNA1
      id_match <- which(dual_library$sgRNA1_WGE_Sequence == multiple_ID$SEQ[idx])
      if (length(id_match) > 0) {
        dual_library$sgRNA1_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
        dual_library$sgRNA1_Approved_Symbol[id_match] <- multiple_ID$gene_all[idx] 
        dual_result$sgRNA1_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
        dual_result$Gene1[id_match] <- multiple_ID$gene_all[idx] 
      }
      # sgRNA2
      id_match <- which(dual_library$sgRNA2_WGE_Sequence == multiple_ID$SEQ[idx])
      if (length(id_match) > 0) {
        dual_library$sgRNA2_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
        dual_library$sgRNA2_Approved_Symbol[id_match] <- multiple_ID$gene_all[idx] 
        dual_result$sgRNA2_WGE_ID[id_match] <- multiple_ID$ID_all[idx]   
        dual_result$Gene2[id_match] <- multiple_ID$gene_all[idx]   
      }
    }
    dual_result$Gene_pair <- paste(dual_result$Gene1, dual_result$Gene2, sep = "~")
  }
  
  dual_library <- dual_library %>% 
    dplyr::mutate(COMB_ID =  paste(sgRNA1_WGE_Sequence, sgRNA2_WGE_Sequence, sep = "_"))
  
  # merge common_pairs across libraries via mean
  unique_fun <- function(x){
    collap <- paste0(unique(x), collapse = ",")
    collap[collap == "NA"] <- NA
    return(collap)
  }
  
  tmp_lib1 <- dual_library %>%
    dplyr::group_by(SEQ_pair) %>%
    dplyr::summarise_all(unique_fun)
  
  tmp_lib2 <- dual_library %>%
    dplyr::group_by(SEQ_pair) %>%
    dplyr::summarise(n_pairs = dplyr::n())
  
  dual_library <- dplyr::full_join(tmp_lib1, 
                       tmp_lib2, 
                       by = "SEQ_pair")
  
  # chr X -> 23, chr Y -> 24
  dual_library <- dual_library %>% 
    dplyr::mutate(
      sgRNA1_Chr = dplyr::case_when(
        sgRNA1_Chr == "X" ~ "23",
        sgRNA1_Chr == "Y" ~ "24",
        # .default = as.character(sgRNA1_Chr)), 
        TRUE ~ as.character(sgRNA1_Chr)), 
      sgRNA2_Chr = dplyr::case_when(
        sgRNA2_Chr == "X" ~ "23",
        sgRNA2_Chr == "Y" ~ "24",
        TRUE ~ as.character(sgRNA2_Chr))
      #.default = as.character(sgRNA2_Chr))
    ) %>%
    dplyr::mutate(
      sgRNA1_Chr = as.numeric(sgRNA1_Chr), 
      sgRNA2_Chr = as.numeric(sgRNA2_Chr), 
      sgRNA1_Start = as.numeric(sgRNA1_Start), 
      sgRNA2_Start = as.numeric(sgRNA2_Start), 
      sgRNA1_End = as.numeric(sgRNA1_End), 
      sgRNA2_End = as.numeric(sgRNA2_End), 
    )
  
  # if CL_name in columns dual_results, get only that
  new_name <- sprintf("%s_logFC", CL_name)
  if (CL_name %in% colnames(dual_result)) {
    dual_result_CL <- dual_result[, c(1:13, which(colnames(dual_result) == CL_name))] %>%
      dplyr::rename(Note = Note1, 
                    MyNote = Note2, 
                    !!new_name := CL_name)
    
  }else{
    stop("CL is not present in the result table!")
  }
  
  if ("Gene_pair" %in% colnames(dual_result_CL)) {
    dual_result_CL <- dual_result_CL 
      dplyr::rename(Gene_Pair = Gene_pair)
  }
  
  tmp_res1 <- dual_result_CL %>%
    dplyr::group_by(SEQ_pair) %>%
    dplyr::summarise_if(is.character, unique_fun)
  
  tmp_res2 <- dual_result_CL %>%
    dplyr::group_by(SEQ_pair) %>%
    dplyr::summarise(logFC =  mean(!!sym(sprintf("%s_logFC", CL_name))),  
                     n_pairs = dplyr::n())
  colnames(tmp_res2)[colnames(tmp_res2) == "logFC"] <- sprintf("%s_logFC", CL_name)
  dual_result_CL <- dplyr::full_join(tmp_res1, 
                       tmp_res2, 
                       by = "SEQ_pair")
  
  if (!identical(dual_result_CL$ID_lib, dual_library$ID_lib)) {
   stop("count and library rows SHOULD match, check what caused this...")
  }
  
  # load CNA
  CNA <- readr::read_table(sprintf('%s', copy_number_file), 
                           show_col_types = FALSE) %>%
    dplyr::mutate(
      CHROM = dplyr::case_when(
        CHROM == "chrX" ~ "chr23",
        CHROM == "chrY" ~ "chr24",
        TRUE ~ as.character(CHROM))
      # .default = as.character(CHROM))
    )
  
  if ("Sampleid" %in% colnames(CNA)) {
    CL_name_unif <- toupper(str_replace_all(CL_name, "[-|_]", ""))
    CNA <- CNA %>% 
      dplyr::filter(toupper(str_replace_all(Sampleid, "[-|_]", "")) %in% CL_name_unif)
    if (nrow(CNA) == 0) {
      stop(sprintf("no CN info available for %s", CL_name))
    }
  }
  
  return(list(CNA = CNA, 
              result = dual_result_CL, 
              library =  dual_library, 
              CL_name = CL_name, 
              out_fold = out_fold))
  
}

#' Run complete CRISPRCleanR analysis on genome-wide CRISPR-cas9 single KO
#'
#' This function performs a complete analysis using CRISPRcleanR on genome-wide CRISPR-cas9 KO data
#' It includes normalization, fold change computation and CRISPRcleanR correction.
#'
#' @param filename_single A string specifying the path of a tsv file containing the raw sgRNA counts. This must be a tab delimited file with one row per sgRNA and the following columns/headers: 
#' - sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
#' - gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration;
#' followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
#' @param min_reads_single This parameter defines a filter threshold value for sgRNAs, based on their average counts in the control sample. 
#' Specifically, it indicates the minimal number of counts that each individual sgRNA needs to have in the controls (on average) in order to be included in the output.
#' @param EXPname A string specifying the name of the experiment. 
#' @param libraryAnnotation_single A data frame containing the sgRNA annotations, with a named row for each sgRNA, and columns for targeted genes, genomic coordinates and possibly other information. 
#' @param display A logic value specifying whether figures containing boxplots with the count values pre/post normalisation and log fold-changes should be visualised (deafult is FALSE).
#' @param outdir Output directory for storing results (default is "./").
#'
#' @return A list containing the following components:
#' - FC: Data frame containing corrected log fold changes
#' - segment: Data frame with segmented results
#' - library: Updated library annotation
#'
#' @examples
#' @export
ccr.run_complete <- function(
  
  filename_single, 
  min_reads_single = 30,
  EXPname,
  libraryAnnotation_single,  
  display = FALSE, 
  outdir = "./"
  
  ) {
  
  normANDfcs <- ccr.NormfoldChanges(
    filename = filename_single,
    min_reads = min_reads_single,
    EXPname = EXPname,
    libraryAnnotation = libraryAnnotation_single, 
    display = display, 
    outdir = outdir
    )
  
  # order and add sequence
  gwSortedFCs <- ccr.logFCs2chromPos(
    normANDfcs$logFCs, 
    libraryAnnotation_single
    )
  
  gwSortedFCs <- gwSortedFCs %>%
    dplyr::mutate(SEQ = libraryAnnotation_single$seq[match(rownames(gwSortedFCs),
                                                           rownames(libraryAnnotation_single))]) %>%
    dplyr::rename(avgFC_notcentered = avgFC) %>%
    dplyr::mutate(avgFC = avgFC_notcentered - median(avgFC_notcentered)) %>%
    dplyr::relocate(avgFC, .after = genes) %>%
    dplyr::relocate(avgFC_notcentered, .after = last_col())
  
  # get CRISPRCleanR single correction
  single_correctedFCs <- ccr.GWclean(gwSortedFCs,
                                     display = display,
                                     label = sprintf("single_%s", EXPname)) 
  
  single_correctedFCs$corrected_logFCs <- single_correctedFCs$corrected_logFCs %>% 
    dplyr::mutate(correction = correctedFC - avgFC)

  # save segments
  single_ccr_segments <- single_correctedFCs$segments
  # add mean and sd
  guideIdx_s_e <- lapply(single_ccr_segments$guideIdx, function(x) 
    as.numeric(str_trim(str_split_1(x, pattern = ","))))
  
  # NOTE: cannot use avg.logFC (from segment), 
  # mean_logFC correspond to ccr correction
  # TODO: remove avg.logFC
  single_ccr_segments$mean_logFC <- sapply(guideIdx_s_e, function(x) 
    mean(single_correctedFCs$corrected_logFCs$avgFC[x[1]:x[2]]))
  
  single_ccr_segments$sd_logFC <- sapply(guideIdx_s_e, function(x) 
    sd(single_correctedFCs$corrected_logFCs$avgFC[x[1]:x[2]]))
  
  # update library with match
  single_correctedFCs <- single_correctedFCs$corrected_logFCs
    
  libraryAnnotation_single <- libraryAnnotation_single %>% 
    dplyr::filter(CODE %in% rownames(single_correctedFCs)) 
  libraryAnnotation_single <- libraryAnnotation_single[rownames(single_correctedFCs),]
  
  # X --> 23, Y --> 24 
  libraryAnnotation_single <- libraryAnnotation_single %>% 
    dplyr::mutate(
      CHRM = dplyr::case_when(
        CHRM == "X" ~ "23",
        CHRM == "Y" ~ "24",
        TRUE ~ as.character(CHRM))
#        .default = as.character(CHRM))
      ) %>%
    dplyr::mutate(CHRM = as.numeric(CHRM))
  
  return(list(
    FC = single_correctedFCs, 
    segment = single_ccr_segments, 
    library = libraryAnnotation_single
    ))
  
}

#' Get summary statistics for singletons
#'
#' This function calculates summary statistics for singletons in a dual KO screen.
#' It provides information such as median, mean, standard deviation, and interquartile range for the fold changes
#' of singletons at both positions (Position 1 and Position 2).
#'
#' @param dual_FC Data frame containing dual-KO fold changes.
#' @param corrected Logical, indicating whether to use ccr2 corrected fold changes (default is FALSE).
#'
#' @return A list containing two data frames with summary statistics for singletons at Position 1 and Position 2.
#' Each data frame includes columns such as target gene, target gene ID, median fold change, interquartile range, mean fold change,
#' standard deviation, sample size, related singletons' genes, related singletons' gene IDs, and additional information.
#'
#' @examples
#' @export
ccr2.get_summary_singletons <- function(dual_FC, 
                                       corrected = FALSE){
  
  name_var <- "avgFC"
  if (corrected) {
    name_var <- "correctedFC"
  }
  
  singletons_pos1 <- dual_FC %>%
    dplyr::filter(grepl("Singletons", info_subtype)) %>%
    dplyr::group_by(sgRNA1_WGE_ID) %>%
    dplyr::summarise(n_ID = dplyr::n()) %>%
    dplyr::arrange(desc(n_ID)) %>%
    dplyr::filter(n_ID > 30) %>% # how to set this?? This threshold is used to 
    # select sgRNA1_WGE_ID that are NON-ESSENTIAL genes / INTERGENIC regions 
    # used as negative controls in the singletons pairs. 
    # Only "high" counts will be associated to a lot of matches having those guides
    dplyr::pull(sgRNA1_WGE_ID) %>%
    sort()
  
  singletons_pos2 <- dual_FC %>%
    dplyr::filter(grepl("Singletons", info_subtype)) %>%
    dplyr::group_by(sgRNA2_WGE_ID) %>%
    dplyr::summarise(n_ID = dplyr::n()) %>%
    dplyr::arrange(desc(n_ID)) %>%
    dplyr::filter(n_ID > 30) %>% # how to set this?? This threshold is used to 
    # select sgRNA2_WGE_ID that are NON-ESSENTIAL genes / INTERGENIC regions 
    # used as negative controls in the singletons pairs. 
    # Only "high" counts will be associated to a lot of matches having those guides
    dplyr::pull(sgRNA2_WGE_ID) %>%
    sort()
  
  if (!identical(singletons_pos1, singletons_pos2)) {
    warning("Different singletons in position 1 and position 2: take the union")
  }
  
  singleton_id <- union(singletons_pos1, singletons_pos2)
  guide_pos1_id <- setdiff(unique(dual_FC$sgRNA1_WGE_ID), singleton_id)
  guide_pos2_id <- setdiff(unique(dual_FC$sgRNA2_WGE_ID), singleton_id)
  
  median_FC_pos1 <- dual_FC %>%
    dplyr::filter(sgRNA1_WGE_ID %in% guide_pos1_id & sgRNA2_WGE_ID %in% singleton_id) %>%
    dplyr::group_by(sgRNA1_WGE_ID) %>% 
    dplyr::summarise(
      target_gene = unique(Gene1), 
      target_gene_ID = unique(sgRNA1_WGE_ID), 
      median = median(!!sym(name_var), na.rm = T), 
      iqr = stats::IQR(!!sym(name_var), na.rm = T), 
      mean = mean(!!sym(name_var), na.rm = T), 
      sd = sd(!!sym(name_var), na.rm = T), 
      n = sum(!is.na(!!sym(name_var))), 
      singletons_gene = paste0(sort(unique(Gene2)), collapse = ","), 
      singletons_gene_ID = paste0(sort(sgRNA2_WGE_ID), collapse = ","), 
      info =  paste0(sort(unique(info)), collapse = ","), 
      info_subtype =  paste0(sort(unique(info_subtype)), collapse = ",")) %>%
    dplyr::mutate(position_target_gene = "Position 1") %>%
    dplyr::select(-sgRNA1_WGE_ID)
  
  median_FC_pos2 <- dual_FC %>%
    dplyr::filter(sgRNA2_WGE_ID %in% guide_pos2_id & sgRNA1_WGE_ID %in% singleton_id) %>%
    dplyr::group_by(sgRNA2_WGE_ID) %>% 
    dplyr::summarise(
      target_gene = unique(Gene2), 
      target_gene_ID = unique(sgRNA2_WGE_ID), 
      median = median(!!sym(name_var), na.rm = T), 
      iqr = stats::IQR(!!sym(name_var), na.rm = T), 
      mean = mean(!!sym(name_var), na.rm = T), 
      sd = sd(!!sym(name_var), na.rm = T), 
      n = sum(!is.na(!!sym(name_var))), 
      singletons_gene = paste0(sort(unique(Gene1)), collapse = ","), 
      singletons_gene_ID = paste0(sort(sgRNA1_WGE_ID), collapse = ","), 
      info =  paste0(sort(unique(info)), collapse = ","), 
      info_subtype =  paste0(sort(unique(info_subtype)), collapse = ",")) %>%
    dplyr::mutate(position_target_gene = "Position 2") %>%
    dplyr::select(-sgRNA2_WGE_ID)
  
  return(list(guide1 = median_FC_pos1, guide2 = median_FC_pos2))
  
}

#' Scale fold changes based on positive and negative controls in dual KO screens
#' This function scales fold changes in a dual KO screens based on the median values of positive and negative controls.
#' It calculates a scaled version of fold changes, making them such that negative controls are centered at 0 and positive controls at 1.
#'
#' @param dual_FC Data frame containing dual KO log fold changes.
#' @param corrected Logical, indicating whether to use corrected fold changes (default is FALSE).
#'
#' @return A data frame with scaled fold changes, where each fold change is centered around the median value of negative controls.
#' If `corrected` is TRUE, the resulting data frame includes a column named `correctedFC_scaled`; otherwise, it includes `avgFC_scaled`.
#'
#' @examples
#' @export
ccr2.scale_pos_neg <- function(dual_FC, corrected = FALSE){
  
  dual_FC_original <- dual_FC
  dual_FC <- dual_FC[!is.na(dual_FC$correction),]
  
  name_var <- "avgFC"
  if (corrected) {
    name_var <- "correctedFC"
  }
  
  median_neg <-  median(dual_FC[dual_FC$info == "NegativeControls", name_var, drop = T], na.rm = TRUE)
  median_pos <-  median(dual_FC[dual_FC$info == "PositiveControls", name_var, drop = T], na.rm = TRUE)
  dual_FC <- dual_FC %>%
    dplyr::mutate(tmp = get(name_var)) %>%
    dplyr::mutate(tmp_scaled = (tmp - median_neg)/(median_neg - median_pos)) %>%
    dplyr::select(-tmp)
  
  # add non-target / non-target
  dual_FC_missing <- dual_FC_original %>%
    dplyr::filter(!ID %in% dual_FC$ID) %>%
    dplyr::mutate(tmp_scaled = NA)
  dual_FC <- rbind(dual_FC, dual_FC_missing)
  
  if (corrected) {
    dual_FC <- dual_FC %>% 
      dplyr::rename(correctedFC_scaled = tmp_scaled)
  }  else {
    dual_FC <- dual_FC %>% 
      dplyr::rename(avgFC_scaled = tmp_scaled)
  }
  
  return(dual_FC)
  
}


#' Match single KO data with singletons and generate comparison plot
#'
#' This function matches single KO data with singletons median results from dual KO (output of ccr2.get_summary_singletons) 
#' and generates a comparison plot between single KO fold changes (FC) and median FC of singletons for the paired genes.
#'
#' @param match_dual_single_seq Data frame containing the matching information between sgRNA in dual KO and single KO libraries (output of ccr2.matchDualandSingleSeq).
#' @param single_FC Data frame containing single KO fold changes.
#' @param singletons_summary_FC List of 2 data frames with summary statistics for singletons from dual KO (output of ccr2.get_summary_singletons).
#' @param saveToFig Logical, indicating whether to save the comparison plot as a figure (default is FALSE).
#' @param display Logical, indicating whether to display the comparison plot (default is TRUE).
#' @param saveFormat If saveToFig is TRUE, the format in which to save the figure (e.g., "png", "pdf", "jpeg").
#' @param EXPname Name of the experiment.
#' @param outdir Output directory for storing results (default is "./").
#'
#' @return A data frame containing matched information between single KO data and singletons by gene name.
#'
#' @examples
#' @export
ccr2.match_singletons_singleFC <- function(match_dual_single_seq, 
                                           single_FC, 
                                           singletons_summary_FC, 
                                           saveToFig = FALSE, 
                                           display = TRUE, 
                                           saveFormat = NULL, 
                                           EXPname = "", 
                                           outdir = "./") {
  
  if (saveToFig) {
    display <- TRUE
    file_name <- sprintf("%s%s_singleFC_vs_singletonsFC.%s", outdir, EXPname, saveFormat)
  }
  
  singletons_summary_FC <- do.call(rbind, singletons_summary_FC)
  
  match_dual_single_seq <- match_dual_single_seq %>% 
    dplyr::filter(ID_single %in% rownames(single_FC))
  
  single_info <- single_FC[match_dual_single_seq$ID_single, ] %>%
    dplyr::mutate(target_gene_ID = match_dual_single_seq$ID)
  
  merged_FC <- dplyr::inner_join(singletons_summary_FC, 
                                 single_info, 
                                 by = "target_gene_ID") %>%
    dplyr::mutate(n = as.character(n))

  if ( display ) {
    
    pl <- ggplot(merged_FC, aes(x = avgFC, y = median, color = n)) + 
      facet_wrap(.~position_target_gene) +
      geom_hline(yintercept = 0, linetype = "dashed") + 
      geom_vline(xintercept = 0, linetype = "dashed") + 
      geom_point(alpha = 0.5) +
      theme_bw() + 
      theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
      xlab("single FC (ProjectScore)") + 
      ylab("singletons FC (median)") + 
      guides(color = guide_legend(title = "n. singletons guides")) + 
      ggtitle("Genes paired with singletons guides")
    print(pl)
    
    if (saveToFig) {
      ggsave(filename = file_name, plot = pl, width = 8, height = 5)
    }
  }
 
  return(merged_FC)
  
}

#' Compute Bliss z-score for the dual KO
#'
#' This function computes Bliss z-score for the dual KO screen, considering the expected fold changes
#' based on the medians of singletons associated with each guide pair.
#'
#' @param dual_KO Data frame containing dual KO fold changes.
#' @param corrected Logical, indicating whether to use CRISPRcleanR^2 corrected fold changes (default is FALSE).
#'
#' @return A data frame with Bliss z-score computed for each guide pair in the dual KO screen.
#' If `corrected` is TRUE, the resulting data frame includes columns named `expected_pairFC_corrected` and `bliss_zscore_corrected`;
#' otherwise, it includes `expected_pairFC` and `bliss_zscore`.
#'
#' @examples
#'
#' @export
ccr2.compute_bliss <- function(dual_FC,
                               corrected = FALSE) { 
  
  dual_FC <- dual_FC %>%
    dplyr::mutate(tmp_expected = NA, 
                  tmp_bliss = NA)
  
  name_var <- "avgFC"
  if (corrected) {
    name_var <- "correctedFC"
  }
  
  id_comp <- which(!grepl("Singletons", dual_FC$info))
  
  singletons_summary <- ccr2.get_summary_singletons(
    dual_FC = dual_FC, 
    corrected = corrected)
  
  for (id in id_comp) {
    
    id_pair <- c(dual_FC$sgRNA1_WGE_ID[id], dual_FC$sgRNA2_WGE_ID[id])
    pos1 <- singletons_summary$guide1$median[singletons_summary$guide1$target_gene_ID == id_pair[1]]
    pos2 <- singletons_summary$guide2$median[singletons_summary$guide2$target_gene_ID == id_pair[2]]
    
    if (length(pos1) > 0 & length(pos2) > 0) {
      dual_FC$tmp_expected[id] <- pos1 + pos2 
    }
  }
  
  loess_model <- loess(formula = as.formula(paste0(name_var, "~tmp_expected")), 
                       data = dual_FC)
  bliss_zscore <- residuals(loess_model)
  dual_FC$tmp_bliss <- NA
  dual_FC$tmp_bliss[as.integer(names(bliss_zscore))] <- bliss_zscore
  
  if (corrected) {
    dual_FC <- dual_FC %>% 
      dplyr::rename(expected_pairFC_corrected = tmp_expected, 
                    bliss_zscore_corrected = tmp_bliss)
  }  else {
    dual_FC <- dual_FC %>% 
      dplyr::rename(expected_pairFC = tmp_expected, 
                    bliss_zscore = tmp_bliss)
  }
  
  return(dual_FC)
  
}


#' Match dual library sequence to single library sequences
#'
#' This function matches a sequence from a dual KO library to sequences in a single KO library based on chromosome and sequence.
#' It matches one sequence in dual KO per time
#'
#' @param dual_seq Data frame with one row containing the sequence information from the dual library.
#' @param single_library_seq Data frame containing the sequence information from the single library for all sgRNA in the library.
#'
#' @return A named logical vector indicating matches between the single library sequences and the one dual library sequence.
#' The names of the vector represent the CHR_GENES_SEQ values from the single library, and the values are logical indicating matches.
#'
#' @examples
match_dual_to_single <- function(dual_seq, single_library_seq){
  
  single_nchar <- sort(unique(nchar(single_library_seq$SEQ)))
  if (length(single_nchar) > 1) {
    stop("seq length in single library MUST always be the same")
  }
  if (nrow(dual_seq) > 1) {
    stop("only one seq from dual library MUST be passed")
  }
  
  match_vect <- rep(0, nrow(single_library_seq))
  names(match_vect) <- single_library_seq$CHR_GENES_SEQ
  
  if (!is.na(dual_seq$CHR)) {
    single_library_seq_chr <- single_library_seq %>%
      dplyr::filter(CHR == dual_seq$CHR)
    
    if (nchar(dual_seq$SEQ) >= nchar(dual_seq$SEQ)) {
      match_id <- sapply(single_library_seq_chr$SEQ, function(x) 
        grepl(pattern = x, x = dual_seq$SEQ))
    }else{
      match_id <- grepl(pattern = dual_seq$SEQ, 
                        x = single_library_seq_chr$SEQ)
    }
    match_vect[single_library_seq_chr$CHR_GENES_SEQ] <- match_id 
  }
  
  return(match_vect)
  
}  

#' Match sequences between sgRNA in dual and single libraries
#'
#' This function matches sequences between a dual library and a single library in a CRISPRcleanR^2 analysis.
#' It creates a data frame with all sgRNAs from dual library (no repetitions, union of those in position 1 and position2) and 
#' includes the matched single library IDs and sequences, if a match exists.
#'
#' @param dual_library Data frame containing information from the dual library (one row per guide pair).
#' @param single_library Data frame containing information from the single library (one row per guide).
#'
#' @return A data frame with one row per sgRNA in dual library (union of position 1 and position 2), including matched single library IDs and sequences.
#'
#' @examples
#' @export
ccr2.matchDualandSingleSeq <- function(dual_library, single_library) {
  
  dual_library_seq <- data.frame(
    CHR = c(dual_library$sgRNA1_Chr, 
            dual_library$sgRNA2_Chr), 
    GENES = c(dual_library$sgRNA1_Approved_Symbol, 
              dual_library$sgRNA2_Approved_Symbol), 
    ID = c(dual_library$sgRNA1_WGE_ID, 
           dual_library$sgRNA2_WGE_ID), 
    LIBRARY = c(dual_library$sgRNA1_Library, 
                dual_library$sgRNA2_Library), 
    SEQ = c(dual_library$sgRNA1_WGE_Sequence, 
            dual_library$sgRNA2_WGE_Sequence)
    ) %>%
    dplyr::mutate(CHR_GENES_SEQ = paste(CHR, GENES, SEQ, sep = "_"))
  
  dual_library_seq <- dual_library_seq %>% 
    dplyr::filter(!duplicated(CHR_GENES_SEQ)) %>%
    dplyr::mutate(ID_single = NA, SEQ_single = NA)

  if (!"hgcn_additional_name" %in% colnames(single_library)) {
    single_library$hgcn_additional_name <- single_library$GENES
  }
    
  single_library_seq <- data.frame(
    CHR = single_library$CHRM,
    GENES = single_library$GENES, 
    ID = rownames(single_library),
    SEQ = single_library$seq, 
    ADDITIONAL_GENE_NAME = single_library$hgcn_additional_name)  %>%
    dplyr::mutate(CHR_GENES_SEQ = paste(CHR, GENES, SEQ, sep = "_"))
    
  # filter for genes and non-targeting
  single_library_seq <- single_library_seq %>%
      dplyr::filter(!duplicated(CHR_GENES_SEQ), 
                    GENES %in% dual_library_seq$GENES | 
                    ADDITIONAL_GENE_NAME %in% dual_library_seq$GENES | 
                    GENES == "NON-TARGETING")
 
  # Create matrix to match single x dual SEQ
  # ALTERNATIVE WITH LIST, TODO: RM
  # match_matrix_tmp <- list()
  # for (id in 1:nrow(dual_library_seq)) {
  #   match_matrix_tmp[[id]] <- match_dual_to_single(dual_library_seq[id,], single_library_seq)
  # }
  # match_matrix <- do.call(cbind, match_matrix_tmp)
  # colnames(match_matrix) <- dual_library_seq$CHR_GENES_SEQ
  
  match_matrix <- sapply(
   seq_len(nrow(dual_library_seq)), function(id)
     match_dual_to_single(dual_library_seq[id,], single_library_seq)
   )
  colnames(match_matrix) <- dual_library_seq$CHR_GENES_SEQ
  
  if (any(rowSums(match_matrix) > 1)) {
    id_mol <- which(rowSums(match_matrix) > 1)
    print(sprintf("Same single seq matches multiple double seq (n. SEQ = %i)", 
                  length(id_mol)))
    match_matrix[id_mol, ] <- 0
  }
  
  if (any(colSums(match_matrix) > 1)) {
    id_mol <- which(colSums(match_matrix) > 1)
    print(sprintf("Multiple single seq match same double seq (n. SEQ = %i)", 
                  length(id_mol)))
    match_matrix[,id_mol] <- 0
  }
  
  # add info on the final table:
  match_id_dual <- which(colSums(match_matrix) == 1)
  match_id_single <- unlist(apply(match_matrix, 2, function(x) which(x == 1)))
  dual_library_seq$ID_single[match_id_dual] <- single_library_seq$ID[match_id_single]
  dual_library_seq$SEQ_single[match_id_dual] <- single_library_seq$SEQ[match_id_single]

  return(dual_library_seq)
}


#' Normalize fold changes in dual KO screens
#'
#' This function normalizes the fold changes in a dual KO screen, removing pairs with low read counts and
#' calculating log2 fold changes relative to a control sample.
#'
#' @param Dframe Data frame containing the raw counts from dual KO screen, raw counts start from the 9th column. The 9th column include the control raw counts.
#' @param min_reads Minimum number of reads required in the control sample (e.g. plasmid) for a pair to be retained (default is 30).
#'
#' @return A list containing log2 fold changes, normalized counts, and normalization factors.
#'
#' @examples
#' @export
ccr2.NormfoldChanges <- function(
  Dframe, 
  min_reads = 30
) {
  
  # remove column already normalized and logFC
  counts <- Dframe %>% 
    dplyr::select(-dplyr::ends_with("_norm"), -dplyr::ends_with("_logFC"))
  # TODO: check that the correct columns are selected (possibly not by number)
  info <- counts[, 1:8]
  counts <- counts[, 9:ncol(counts)]
  # remove pairs with less than min_reads in plasmid
  id_keep <-  counts[, 1] >= min_reads
  info <- info[id_keep, ]
  info <- info %>% 
    dplyr::mutate(
      COMB_ID = dplyr::case_when(
        !is.na(sgRNA1_ID) & !is.na(sgRNA2_ID) ~ paste(sgRNA1_ID, sgRNA2_ID, sep = "_"),
        !is.na(sgRNA1_ID) & is.na(sgRNA2_ID) ~ paste(sgRNA1_ID, Gene2, sep = "_"), 
        is.na(sgRNA1_ID) & !is.na(sgRNA2_ID) ~ paste(Gene1, sgRNA2_ID, sep = "_"), 
        is.na(sgRNA1_ID) & is.na(sgRNA2_ID) ~ paste(Gene1, Gene2, sep = "_")))
  counts <- counts[id_keep, ]
  
  norm_fact_vector <- colSums(counts)
  norm_fact <- t(matrix(rep( norm_fact_vector, nrow(counts)), 
                        ncol(counts), nrow(counts)))
  norm_counts <- counts/norm_fact*10000000
  
  name_control <- colnames(norm_counts)[1]
  norm_controls <- norm_counts[, 1]
  norm_counts <- norm_counts[, -1]
  
  log_fc <- apply(norm_counts, 2, function(x) log2((x + 0.5)/(norm_controls + 0.5)))
  colnames(log_fc) <- paste(colnames(log_fc), "logFC", sep = "_")
  log_fc <- cbind(info, log_fc)
  norm_counts <- cbind(info, norm_controls , norm_counts)
  colnames(norm_counts)[ncol(info) + 1] <- name_control
  
  return(list(logFCs = log_fc, 
              norm_counts = norm_counts, 
              factor = norm_fact_vector))
}

#' Compute mean of logFCs across replicates, annotate sgRNA pairs with chromosomal positions
#'
#' This function computes the mean of logFCs across replicates (if any) and annotate each 
#' sgRNA pairs with chromosomal positions (BP) as the median point of the sgRNA 
#' in position 1/position 2 start and end. It provides additional information about the sgRNAs and their genomic locations.
#'
#' @param dual_FC data frame containing dual KO log fold changes,
#' @param dual_library data frame containing information about the sgRNAs and their genomic locations.
#'
#' @return A data frame with avg log fold changes and information about the sgRNAs' genomic locations.
#'
#' @examples
#' @export
ccr2.logFCs2chromPos <- function(
  dual_FC, 
  dual_library
) {
  
  dual_library <-  dual_library %>% 
    dplyr::select(ID, ID_lib, lib, sgRNA1_Library, 
                  sgRNA1_WGE_ID, sgRNA1_WGE_Sequence, 
                  sgRNA1_Chr, sgRNA1_Start, sgRNA1_End, 
                  sgRNA2_Library, 
                  sgRNA2_WGE_ID, sgRNA2_WGE_Sequence, 
                  sgRNA2_Chr, sgRNA2_Start, sgRNA2_End)
  
  # average FC
  # # remove scaled
  # dual_FC <- dual_FC %>% dplyr::select(-dplyr::ends_with("_Scaled_logFC"))
  avg_FC <- dual_FC %>% 
    # dplyr::filter(.[[9]] >= min_reads) %>%
    dplyr::select(ID, ID_lib, lib, MyNote, Note, Gene_Pair, Gene1, Gene2, 
                  # dplyr::ends_with("_Scaled_logFC")) %>%
                  dplyr::ends_with("_logFC")) %>%
    dplyr::rename(info = Note, info_subtype = MyNote) %>% 
    dplyr::mutate(avgFC = rowMeans(dplyr::select(.,ends_with("_logFC"))))
  #dplyr::mutate(avgFC = rowMeans(dplyr::select(., ends_with("_Scaled_logFC"))))
  avg_FC <-  avg_FC %>% dplyr::select(-dplyr::ends_with("_logFC"))
  
  # merge with position info
  combined <- dplyr::left_join(avg_FC, dual_library, by = c("ID_lib", "ID", "lib") ) %>%
    dplyr::mutate(sgRNA1_BP = sgRNA1_Start + (sgRNA1_End - sgRNA1_Start)/2, .after = sgRNA1_End) %>%
    dplyr::mutate(sgRNA2_BP = sgRNA2_Start + (sgRNA2_End - sgRNA2_Start)/2, .after = sgRNA2_End)
  
  return(combined)
  
}

#' Create pseudo-single logFCs from dual KO logFCs in CRISPRcleanR^2 analysis
#'
#' This function creates pseudo-single logFCs from dual KO logFCs. Considering a sgRNA in position i (1 or 2), 
#' it assign the pseudo single value as the average logFC across all the dual guide pairs having that sgRNA in position i. 
#' The resulting data frame contains details about the sgRNAs, their genomic locations, 
#' and log fold changes from both the pseudo single and the single KO, if the sgRNA sequence matches (if not replaced by NAs)
#'
#' @param dual_FC Data frame with dual KO logFCs (already averaged per replications, logFCs stored in "avgFC" column)
#' @param single_FC Data frame with single KO logFCs.
#' @param match_dual_single_seq Data frame containing matched sequences between the dual and single libraries (output of ccr2.matchDualandSingleSeq).
#' @param guide_id Identifier for the sgRNA position to be considered (1 or 2).
#'
#' @return A data frame containing details about the sgRNAs, their genomic locations, 
#' and log fold changes from both the pseudo single and the single KO, if the sgRNA sequence matches (if not replaced by NAs).
#'
#' @examples
#' @export
ccr2.createPseudoSingle <- function(
  dual_FC, 
  single_FC, 
  match_dual_single_seq, 
  guide_id
) { 
  
  # filter library seq based on single_FC match
  #match_dual_single_seq <- match_dual_single_seq %>% 
  #  dplyr::filter(ID_single %in% rownames(single_FC))
  
  # filter dual per guides with a single counterpart
  guide <- dual_FC 
  #%>% 
   # dplyr::filter((sgRNA1_WGE_ID %in% match_dual_single_seq$ID &   
  #                   sgRNA2_WGE_ID %in% match_dual_single_seq$ID) |  
  #                  ((grepl("NONTARGET", Gene1) | grepl("CTRL", Gene1)) & sgRNA2_WGE_ID %in% match_dual_single_seq$ID) | 
  #                  ((grepl("NONTARGET", Gene2) | grepl("CTRL", Gene2)) & sgRNA1_WGE_ID %in% match_dual_single_seq$ID))
  
  var <- sprintf("sgRNA%i", guide_id)
  other_var <- setdiff(c("sgRNA1", "sgRNA2"), sprintf("sgRNA%i", guide_id))
  
  guide <- guide %>% 
    dplyr::filter(!is.na(!!sym(sprintf("%s_BP", var)))) %>% 
    dplyr::group_by(!!sym(sprintf("%s_WGE_ID", var))) %>% 
    dplyr::summarise(sgRNA_Library = unique(!!sym(sprintf("%s_Library", var))), 
                     info_subtype = paste0(unique(info_subtype), collapse = ","), 
                     genes = unique(!!sym(sprintf("Gene%i", guide_id))), 
                     avgFC = mean(avgFC), 
                     n = dplyr::n(), 
                     matched_sgRNA_ID = paste0(sort(!!sym(sprintf("%s_WGE_ID", other_var))), collapse = ","),
                     #lib = paste0(unique(lib), collapse = ","),
                     CHR =  unique(!!sym(sprintf("%s_Chr", var))), 
                     startp = unique(!!sym(sprintf("%s_Start", var))), 
                     endp = unique(!!sym(sprintf("%s_End", var))),
                     BP = unique(!!sym(sprintf("%s_BP", var))), 
                     SEQ = unique(!!sym(sprintf("%s_WGE_Sequence", var)))) %>%
    dplyr::arrange(CHR, BP) %>%
    dplyr::ungroup() %>%
    dplyr::rename(sgRNA_ID = !!sym(sprintf("%s_WGE_ID", var)))
  
  # add matching single ID
  guide <- dplyr::left_join(
    guide, 
    match_dual_single_seq %>% 
      dplyr::select(ID, ID_single, SEQ_single) %>%
      dplyr::rename(sgRNA_ID = ID), by = "sgRNA_ID"
  )
  
  # add logFC from single
  guide <- dplyr::left_join(
    guide, 
    single_FC %>% 
      tibble::rownames_to_column(var = "ID_single") %>%
      dplyr::select(ID_single, avgFC, correction, correctedFC) %>%
      dplyr::rename(avgFC_single = avgFC, 
                    correction_single = correction, 
                    correctedFC_single = correctedFC), 
    by = "ID_single")
  
  return(guide)
  
}

#' Create pseudo-single logFCs from dual KO logFCs in CRISPRcleanR^2 analysis for both positions
#'
#' This function creates pseudo-single logFCs from dual KO logFCs, both for position 1 and position 2 (see ccr2.createPseudoSingle)
#' 
#' @param dual_FC Data frame with dual KO logFCs (already averaged per replications, logFCs stored in "avgFC" column).
#' @param single_FC Data frame with single KO logFCs.
#' @param match_dual_single_seq Data frame containing matched sequences between the dual and single libraries (output of ccr2.matchDualandSingleSeq).
#' @param saveToFig Logical indicating whether to save the plot to a file.
#' @param display Logical indicating whether to display the plot.
#' @param saveFormat The format in which to save the plot (e.g., "png", "pdf"). Default is NULL.
#' @param EXPname Name for the experiment. Default is an empty string.
#' @param outdir Directory where the plot file will be saved. Default is "./".
#'
#' @return A list containing data frames for pseudo-single guides in position 1 (`sgRNA1`) and position 2 (`sgRNA2`).
#'
#' @examples
#' @export
ccr2.createPseudoSingle_combine <- function(
  dual_FC, 
  single_FC, 
  match_dual_single_seq, 
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = NULL, 
  EXPname = "", 
  outdir = "./"
) { 
  
  guide1 <- ccr2.createPseudoSingle(
    dual_FC = dual_FC,
    single_FC = single_FC,
    match_dual_single_seq = match_dual_single_seq,
    guide_id = 1
  )
  
  guide2 <- ccr2.createPseudoSingle( 
    dual_FC = dual_FC,
    single_FC = single_FC,
    match_dual_single_seq = match_dual_single_seq,
    guide_id = 2
  )
  
  if (saveToFig) {
    display <- TRUE
    file_name <- sprintf("%s%s_PseudoGuideCommon.%s", outdir, EXPname, saveFormat)
  }
  
  
  if (display) {
    
    common_guides <- dplyr::inner_join(guide1, guide2, by = "sgRNA_ID")
    
    if (nrow(common_guides) > 0) {
      
      matched_ID_1 <- lapply(common_guides$matched_sgRNA_ID.x, 
                             function(x) str_split(x, pattern = ",")[[1]])
      matched_ID_2 <- lapply(common_guides$matched_sgRNA_ID.y, 
                             function(x) str_split(x, pattern = ",")[[1]])
      common_guides$jaccard_sim <- mapply(
        function(x,y) length(intersect(x,y))/length(union(x,y)), 
        x = matched_ID_1, y = matched_ID_2, SIMPLIFY = TRUE
        )
      breaks_size <- round(seq(min(common_guides$jaccard_sim), 
                                max(common_guides$jaccard_sim), length.out = 6), 
                            digits = 2)
      
      #breaks_size1 <- round(seq(min(common_guides$n.x), 
      #                          max(common_guides$n.x), length.out = 6))
      #breaks_size2 <- round(seq(min(common_guides$n.y), 
      #                          max(common_guides$n.y), length.out = 6))
      text_corr <- data.frame(
        label = sprintf("Pear. corr. = %.3f", 
                        cor(common_guides$avgFC.x, common_guides$avgFC.y)), 
        xpos = -Inf, 
        ypos = Inf)
      
      pl <- ggplot(common_guides, aes(x = avgFC.x, y = avgFC.y, 
                                      color = jaccard_sim)) + 
        geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
        geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
        geom_point(alpha = 0.5, size = 2) +
        theme_bw() + 
        theme(legend.position = "right", plot.title = element_text(hjust = 0.5)) + 
        scale_colour_viridis_c(breaks = breaks_size) +
        # scale_size_continuous(breaks = breaks_size1) + 
        geom_text(data = text_corr, 
                  aes(label = label, x = xpos, y = ypos), size = 5, 
                  hjust = -0.1, vjust = 1.1, inherit.aes = F) +
        xlab("Psuedo single guide position 1") + 
        ylab("Pseudo single guide position 2") + 
        labs(color = "Jaccard Sim.\n(matched guides)") + 
        ggtitle(sprintf("Guides in common (N=%i)", nrow(common_guides)))
      print(pl)
    }
    
    if (saveToFig) {
      ggsave(filename = file_name, plot = pl, width = 6, height = 5)
    }
  }
  return(list(sgRNA1 = guide1, sgRNA2 = guide2))
}


#' Model comparison between pseudo-single and single guides logFCs for a specific position
#'
#' This function performs a linear regression to model the relationship between 
#' average log fold changes of pseudo-single guides and their corresponding single guides 
#' for a specific guide position. 
#' The regression also include the number of guide pairs used to compute the pseudo-single logFCs ("n"). 
#' In addition, it removed outliers based on extreme cooks distance from the regression (> 99% quantile) 
#' and it recomputes the linear regression after those outliers removal.
#'
#' @param pseudo_single_FC Data frame containing logFCs for pseudo-single guides and logFCs for single KO in matching guides (output of ccr2.createPseudoSingle_combine).
#' @param guide_id Identifier for the sgRNA position to be considered (1 or 2).
#' @param display A logical indicating whether to display plots. Default is TRUE.
#' @param saveToFig A logical indicating whether to save plots to files. Default is FALSE.
#' @param saveFormat The format in which to save the plots (e.g., "pdf", "png"). Default is "pdf".
#' @param outdir Directory where the plot files will be saved. Default is "./".
#' @param EXPname Name for the experiment. Default is an empty string.
#'
#' @return A list containing the model fit, summary statistics of the model, and guides removed due to as outliers from Cook's distance
#'
#' @examples
#' @export
ccr2.modelSingleVSPseudoSingle <- function(
  pseudo_single_FC, 
  guide_id, 
  display = TRUE, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = ""
) { 
  
  # remove NA in single screens
  matched_df <- pseudo_single_FC[!is.na(pseudo_single_FC$ID_single), ]
  # matched_df$n_bin <- cut(
  #  matched_df$n, 
  #  breaks = c(0, 1, seq(10, round(max(matched_df$n), -1) + 10, 10))
  #  )
  # matched_df$n <- factor(matched_df$n,levels = sort(unique(matched_df$n)))
  fmla <- "avgFC ~ 0 + avgFC_single + n"
  fit_model <- glm(formula = as.formula(fmla), 
                   data = matched_df,
                   family = gaussian(link = "identity"))
  # remove outliers
  cooksD <- cooks.distance(fit_model) 
  cooksD <- cooksD[!is.na(cooksD)]
  # remove extreme cases (quatile > 0.99)
  thr_influential <- quantile(cooksD, prob = 0.99, na.rm = T)
  influential <- cooksD[cooksD > thr_influential]
  # influential <- cooksD[(cooksD > (1.5 * mean(cooksD, na.rm = TRUE)))]

  if (length(influential) > 0) {
    # repeat:
    matched_df_filt <- matched_df[!rownames(matched_df) %in% names(influential), ] %>%
      as.data.frame()
    rownames(matched_df_filt) <- paste0(matched_df_filt$genes, "_", 
                                        matched_df_filt$sgRNA_ID)
    fit_model <- glm(formula = as.formula(fmla), 
                     data = matched_df_filt,
                     family = gaussian(link = "identity"))  
    rm_guides_cooks <- matched_df[names(influential), ]
    
  }else{
    matched_df_filt <- matched_df
    rm_guides_cooks <- NA
  }
  
  df_plot <- data.frame(id = rownames(fit_model$model),
                        avgFC = fit_model$model$avgFC,
                        avgFC_single = fit_model$model$avgFC_single,
                        fitted_avgFC = fit_model$fitted.values,
                        residuals_avgFC = residuals(fit_model),
                        n_guides_in_mean = as.numeric(as.character(fit_model$model$n)))
                        # n_guides_in_mean = as.numeric(fit_model$weights))
                        # lib = dual_collapsed_filt$lib)
  
  if (saveToFig) {
    display <- TRUE
    file_name_comp <- sprintf("%s%s_PseudoGuide%i_vs_single.%s", 
                              outdir, EXPname, guide_id, saveFormat)
    file_name_comp_pred <- sprintf("%s%s_PseudoGuide%i_vs_FittedPseudoGuide%i.%s", 
                              outdir, EXPname, guide_id, guide_id, saveFormat)
    file_name_resid <- sprintf("%s%s_PseudoGuide%i_vs_single_resid.%s", 
                               outdir, EXPname, guide_id, saveFormat)
    file_name_qqplot <- sprintf("%s%s_PseudoGuide%i_vs_single_qqplot.%s", 
                                outdir, EXPname, guide_id, saveFormat)
  }
  
  if (display) {
    
    n_min <- min(df_plot$n_guides_in_mean)
    n_max <- max(df_plot$n_guides_in_mean)
    breaks_size <- round(seq(n_min, 
                             n_max, length.out = 6))
    text_corr <- data.frame(
      label = sprintf("Pear. corr. = %.3f", cor(df_plot$avgFC, df_plot$fitted_avgFC)), 
      xpos = -Inf, 
      ypos = Inf)
    
    pl_comp_pred <- ggplot(df_plot, aes(x = avgFC, y =  fitted_avgFC)) + 
      geom_point(aes(y = fitted_avgFC, 
                     # size = n_guides_in_mean, 
                     color = n_guides_in_mean), 
                 alpha = 0.5) +
      # geom_line(aes(y =  fitted_avgFC), colour = "red", linewidth = 1) + 
      geom_abline(intercept = 0, slope = 1, 
                  colour = "black", size = 1, linetype = "dashed") + 
      #geom_smooth(formula = y ~ 0 + x, method = "loess", 
      #            span = 0.3, 
      #            method.args = list(degree = 2)) + 
      theme_bw() + 
      theme(axis.title = element_text(size = 12), 
            axis.text = element_text(size = 11), 
            legend.position = "bottom", 
            plot.title = element_text(hjust = 0.5)) + 
      scale_colour_viridis_c(breaks = breaks_size) +
      geom_text(data = text_corr, 
                aes(label = label, x = xpos, y = ypos), size = 5, 
                hjust = -0.1, vjust = 1.1, inherit.aes = F) +
      xlab("Pseudo Single avgFC") + 
      ylab("Fitted Pseudo Single avgFC") + 
      ggtitle(sprintf("Guide position %i", guide_id))
    print(pl_comp_pred)
    
    text_corr <- data.frame(
      label = sprintf("Pear. corr. = %.3f", cor(df_plot$avgFC, df_plot$avgFC_single)), 
      xpos = -Inf, 
      ypos = Inf)
    
    pl_comp <- ggplot(df_plot, aes(x = avgFC_single, y =  avgFC)) + 
      geom_line(aes(y =  fitted_avgFC), colour = "red", linewidth = 1) + 
      geom_point(aes(y =  avgFC, 
                     # size = n_guides_in_mean, 
                     color = n_guides_in_mean), 
                     alpha = 0.5) +
      geom_abline(intercept = 0, slope = 1, 
                  colour = "black", size = 1, linetype = "dashed") + 
      #geom_smooth(formula = y ~ 0 + x, method = "loess", 
      #            span = 0.3, 
      #            method.args = list(degree = 2)) + 
      theme_bw() + 
      theme(axis.title = element_text(size = 12), 
            axis.text = element_text(size = 11), 
            legend.position = "bottom", 
            plot.title = element_text(hjust = 0.5)) + 
      scale_colour_viridis_c(breaks = breaks_size) +
      geom_text(data = text_corr, 
                aes(label = label, x = xpos, y = ypos), size = 5, 
                hjust = -0.1, vjust = 1.1, inherit.aes = F) +
      xlab("Single avgFC") + 
      ylab("Pseudo Single avgFC") + 
      ggtitle(sprintf("Guide position %i", guide_id))
    print(pl_comp)
    
    # residuals VS fitted 
    df_plot$outliers <- ""
    df_plot$outliers[abs(df_plot$residuals_avgFC) > 2] <- df_plot$id[abs(df_plot$residuals_avgFC) > 2]
    pl_pred_vs_res <- ggplot(df_plot, aes(x = fitted_avgFC, 
                                          y = residuals_avgFC)) + 
      geom_text_repel(color = "black", aes(label = outliers), 
                      min.segment.length = 0, size = 3) + 
      geom_hline(yintercept = 0, colour = "black", size = 1, linetype = "dashed") + 
      geom_smooth(formula = y ~ x, method = "loess") +
      geom_point(aes(y = residuals_avgFC, 
                     #size = n_guides_in_mean, 
                     color = n_guides_in_mean), 
                 alpha = 0.5) +
      theme_bw() + 
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5)) + 
      scale_colour_viridis_c(breaks = breaks_size) +
      #scale_size_continuous(breaks = breaks_size) + 
      ylab("Residuals") + 
      xlab("Fitted avgFC") + 
      ggtitle(sprintf("Guide position %i", guide_id))
    print(pl_pred_vs_res)
    
    # qqplot 
    pl_qq <- ggplot(df_plot, aes(sample = residuals_avgFC)) + 
      geom_hline(yintercept = 0, colour = "red", size = 0.5, linetype = "dashed") + 
      geom_vline(xintercept = 0, colour = "red", size = 0.5, linetype = "dashed") + 
      stat_qq(alpha = 0.5, size = 2) + 
      stat_qq_line(size = 1) +
      theme_bw() + 
      theme(legend.position = "right", 
            plot.title = element_text(hjust = 0.5)) + 
      ggtitle(sprintf("Guide position %i", guide_id)) + 
      ylab("Observed") + 
      xlab("Theoretical") 
    
    x.pnts <- ggplot_build(pl_qq)$data[[3]]$x
    y.pnts <- ggplot_build(pl_qq)$data[[3]]$y
    
    pl_qq <- pl_qq +  
      geom_text_repel(color = "black", 
                      label = df_plot$outliers[order(df_plot$residuals_avgFC)], 
                      x = x.pnts, 
                      y = y.pnts, 
                      min.segment.length = 0, size = 3, 
                      max.overlaps = Inf)
    print(pl_qq)  
  }
  
  if (saveToFig) {
    ggsave(filename = file_name_comp, plot = pl_comp, width = 4.5, height = 5.5)
    ggsave(filename = file_name_comp_pred, plot = pl_comp_pred, width = 4.5, height = 5.5)
    ggsave(filename = file_name_resid, plot = pl_pred_vs_res, width = 6, height = 5)
    ggsave(filename = file_name_qqplot, plot = pl_qq, width = 5, height = 5)
  }
  
  return(list(model = fit_model,
              model_summary = summary(fit_model), 
              rm_guides_cooks = rm_guides_cooks))
}

#' Inject data from single KO space to pseudo-single space (for a specific guide position) at the genome-wide level
#'
#' This function predicts pseudo-single logFCs values at the genome-wide level 
#' based on the model fitted to the relationship between matching pseudo-single and single guides. 
#' The predicted values are filling the genome-wide space based on single KO guides (data injection). 
#' The actual pseudo-single guides instead of the predicted ones are used where available.
#'
#' @param model_single_to_pseudo Fitted model object representing the relationship between single and pseudo-single logFCs (output of ccr2.modelSingleVSPseudoSingle).
#' @param pseudo_single_FC Data frame containing logFCs for pseudo-single guides and logFCs for single KO in matching guides (output of ccr2.createPseudoSingle_combine).
#' @param single_FC Data frame containing logFCs for single KO.
#' @param guide_id Identifier for the guide position (1 or 2).
#' @param correctGW A character specifying the type of pseudo-single logFCs mean centering after data-injection ("CHR": per chromosome, "GW": genome-wide, or NULL for no centering). Default is NULL. #TODO: USED FOR TESTING, REMOVE!
#'
#' @return A data frame containing information for pseudo-single guides with injected values.
#'
#' @examples
#' @export
ccr2.injectData <- function(
  model_single_to_pseudo,
  pseudo_single_FC,
  single_FC, 
  guide_id, 
  correctGW = NULL
) { 
  
  # set n equal to the value for library combinations, rationale: class with the majority of common pairs!
  #id_libcomb <- grepl("LibraryCombinations", pseudo_single_FC$info_subtype)
  # n_fixed <- names(which.max(table(pseudo_single_FC$n[id_libcomb])))
  n_fixed <- 1
  print(paste0("Fixed n_guides_in_mean for library combinations: ", n_fixed))
  
  ## predict on single guides (make a new function) ##
  single_FC <- single_FC %>% 
    dplyr::rename(avgFC_single = avgFC, correction_single = correction) %>%
    dplyr::mutate(n = n_fixed)
  
  pseudo_single_GW <- predict(model_single_to_pseudo, newdata = single_FC)
  pseudo_single_GW <- single_FC[, 1:7] %>% 
    dplyr::mutate(sgRNA_ID_single = rownames(single_FC), 
                  avgFC = unname(pseudo_single_GW), 
                  guideID = guide_id, 
                  sgRNA_ID_dual = NA) %>%
    dplyr::rename(avgFC_original = avgFC_single)
  
  # replace fitted with original values when available
  matched_df <- pseudo_single_FC[!is.na(pseudo_single_FC$ID_single), ]
  pseudo_single_GW$avgFC[match(matched_df$ID_single, 
                                pseudo_single_GW$sgRNA_ID_single)] <- matched_df$avgFC
  pseudo_single_GW$sgRNA_ID_dual[match(matched_df$ID_single, 
                               pseudo_single_GW$sgRNA_ID_single)] <- matched_df$sgRNA_ID
  
  # add values not matched by position
  not_matched_to_add <- pseudo_single_FC %>%
    dplyr::filter(is.na(ID_single)) %>%
    dplyr::mutate(avgFC_original = NA, 
                  SEQ = NA, 
                  sgRNA_ID_single = NA, 
                  guideID = guide_id) %>%
    dplyr::rename(sgRNA_ID_dual = sgRNA_ID)
  
  col_names <- colnames(pseudo_single_GW)
  pseudo_single_GW <- dplyr::bind_rows(
    pseudo_single_GW, 
    not_matched_to_add[, col_names]) %>% 
    dplyr::arrange(CHR, BP)  
  
  if (!is.null(correctGW)) {
    
    if (correctGW == "CHR") {
      pseudo_single_GW <- pseudo_single_GW %>% 
        dplyr::group_by(CHR) %>%
        dplyr::rename(avgFC_uncorr = avgFC) %>%
        dplyr::mutate(avgFC = avgFC_uncorr - mean(avgFC_uncorr)) %>%
        dplyr::ungroup()
    } else {
      if (correctGW == "GW") {
        pseudo_single_GW <- pseudo_single_GW %>% 
          dplyr::rename(avgFC_uncorr = avgFC) %>%
          dplyr::mutate(avgFC = avgFC_uncorr - mean(avgFC_uncorr)) %>%
          dplyr::ungroup()
      } else {
        print("correctGW MUST be either CHR (per chromosome) or GW (genome-wide) or NULL (no correction)")   
      }
    }
  } else {
    pseudo_single_GW <- pseudo_single_GW %>% 
      dplyr::mutate(avgFC_uncorr = NA)
  }
  
  return(pseudo_single_GW)
  
}

#' Post-processing of CRISPRcleanR corrected psuedo-single logFCs 
#'
#' This function filters and post-process psuedo-single logFCs after genome-wide 
#' CRISPRcleanR correction (output of ccr.GWclean). It filters for those guides originally available in the pseudo-single logFCs (i.e. present in dual KO)
#' It provides information about the pseudo-single correction values and gene segments in a specific guide position (1 or 2).
#'
#' @param dataInjection_correctedFCs Data frame containing corrected injected pseudo-single logFCs via CRISPRcleanR (output of ccr.GWclean).
#' @param dataInjection_segments Data frame containing information about identified gene segments from injected pseudo-single logFCs via CRISPRcleanR (output of ccr.GWclean).
#' @param pseudo_single Data frame containing information for pseudo-single logFCs (before CRISPRcleanR correction, output of ccr2.createPseudoSingle_combine)
#' @param guide_id Identifier for the guide position (1 or 2).
#' @param saveToFig Logical, indicating whether to save the generated plots to files. Default is FALSE.
#' @param display Logical, indicating whether to display the plots. Default is TRUE.
#' @param saveFormat File format for saving the plots (e.g., ".pdf", ".png"). Default is ".pdf".
#' @param outdir Directory to save the plots. Default is "./".
#' @param EXPname Experiment name to include in the plot filenames.
#'
#' @return A list containing data frames with information about corrected pseudo-single logFCs 
#' and gene segments identified by CRISPRcleanR, for those guides originally available in the dual KO screen.
#'
#' @examples
#' @export
ccr2.filterGWclean <- function(
  dataInjection_correctedFCs, 
  dataInjection_segments,
  pseudo_single, 
  guide_id, 
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = ".pdf",
  outdir = "./", 
  EXPname = ""
) {
  
  
  pseudo_single_correctedFCs <- dataInjection_correctedFCs %>% 
    dplyr::mutate(correction = correctedFC - avgFC) %>% 
    dplyr::select(sgRNA_ID_dual,
                  guideIdx, guideID, 
                  correction, correctedFC) %>%
    dplyr::rename(sgRNA_ID = sgRNA_ID_dual)
  
  complete_out <- dplyr::left_join(
    pseudo_single, 
    pseudo_single_correctedFCs, 
    by = "sgRNA_ID") %>%
    dplyr::mutate(correction_scaled = n*correction)
  
  # get segments
  tmp <- dataInjection_correctedFCs %>% 
    dplyr::mutate(GW_guideIdx = 1:length(guideIdx)) %>%
    dplyr::filter(!is.na(sgRNA_ID_dual)) %>%
    dplyr::mutate(correction = correctedFC - avgFC)
  
  gw_pseudo_single_segments <- dataInjection_segments %>%
    dplyr::mutate(
      seg_id = sprintf("chr%i_%i_%i", CHR, startp, endp), 
      guide_start = as.numeric(str_trim(str_split_i(string = guideIdx, pattern = "[,]", i = 1))), 
      guide_end = as.numeric(str_trim(str_split_i(string = guideIdx, pattern = "[,]", i = 2)))
    )
  
  tmp$segment_id <- sapply(tmp$GW_guideIdx,
                           function(x) 
                             gw_pseudo_single_segments$seg_id[x <= gw_pseudo_single_segments$guide_end & 
                                                                x >= gw_pseudo_single_segments$guide_start])
  
  pseudo_single_segments <- tmp %>% 
    dplyr::group_by(segment_id) %>%
    dplyr::summarise(
      CHR = unique(CHR), 
      startp = min(startp), 
      endp = max(endp),
      correction = mean(correction), # small fluctuation, unique does not work
      n_guides = length(sgRNA_ID_dual), 
      all_guides = paste0(sgRNA_ID_dual, collapse = ","), 
      n_genes = length(unique(genes)),
      all_genes = paste0(unique(genes), collapse = ",")) %>%
    arrange(CHR, startp) %>%
    ungroup()
  
  
  if (saveToFig) {
    display <- TRUE
    file1_name <- sprintf("%s%s_correction_singleVSpseudosingle_guide%i.%s", 
                         outdir, EXPname, guide_id, saveFormat)
    file2_name <- sprintf("%s%s_hist_genes_per_segment_guide%i.%s", 
                          outdir, EXPname, guide_id, saveFormat)
  }
  
  if (display) {
    
    pl1 <- ggplot(complete_out, 
                 aes(x = correction_single, y = correction)) + 
      geom_point(alpha = 0.7, size = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
      theme_bw() + 
      theme(legend.position = "right") + 
      ggtitle(sprintf("Guide position %i", guide_id)) +
      xlab("Single correction") + 
      ylab("Pseudo single correction")
    print(pl1)
    
    pl2 <- ggplot(pseudo_single_segments, 
                 aes(x = n_genes)) + 
      geom_bar(width = 0.5) +
      theme_bw() + 
      theme(legend.position = "right") + 
      ggtitle(sprintf("Guide position %i", guide_id)) +
      xlab("N. of genes in a segment") + 
      ylab("N. of segments")
    print(pl2)
    
  }
  
  if (saveToFig) {
    ggsave(filename = file1_name, plot = pl1, width = 4.5, height = 4.5)
    ggsave(filename = file2_name, plot = pl2, width = 3, height = 4.5)
  }
  
  return(list(pseudo_single_FC = complete_out, 
              pseudo_single_seg = pseudo_single_segments))
  
}

#' Approximate the solutions of linear system to retrieve dual KO CRISPRcleanR^2 correction
#'
#' This function solves a linear system (approximated solution, underdetermined system plus box constraints) to identify dual KO CRISPRcleanR^2 correction, 
#' taking into account pseudo-single correction and how pseudo-single guides are built from dual KO guides.
#' It provides information about the dual KO corrected logFCs as well as the correction values for pseudo-single logFCs.
#' The function also plots the actual pseudo-single correction versus the approximated solution, original VS corrected dual logFCs and psuedo-single correction (position1 and position2) compared to retrieved dual correction.
#'
#' @param dual_FC Data frame with dual KO logFCs (already averaged per replications, logFCs stored in "avgFC" column).
#' @param single_FC_corrected Data frame containing CRISPRcleanR corrected logFCs for single KO.
#' @param match_dual_single_seq Data frame containing matched sequences between the dual and single libraries (output of ccr2.matchDualandSingleSeq).
#' @param pseudo_single_p1_correctedFCs Data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides in position 1 (output of ccr2.filterGWclean).
#' @param pseudo_single_p2_correctedFCs Data frame containing cCRISPRcleanR corrected logFCs for pseudo-single guides in position 2 (output of ccr2.filterGWclean).
#' @param split_graph Logical, indicating whether to split the graph representing the coefficient matrix of the system into largest connected component plus everything else. It depends on the position 1 VS position 2 pairs configurations. When possible, split reduces the computational complexity. Default is TRUE.
#' @param display Logical, indicating whether to display the plots. Default is TRUE.
#' @param saveToFig Logical, indicating whether to save the generated plots to files. Default is FALSE.
#' @param saveFormat File format for saving the plots (e.g., ".pdf", ".png"). Default is ".pdf".
#' @param outdir Directory to save the plots. Default is "./".
#' @param EXPname Experiment name to include in the plot filenames.
#'
#' @return A list containing 
#' -  max_correction: maximum dual correction allowed based on single KO correction (used as upper bound for the solution).
#' -  min_correction: minimum dual correction allowed based on single KO correction (used as upper bound for the solution).
#' -  dual_FC: Data frame of dual KO including CRISPRcleanR^2 correction.
#' -  pseudo_single: Data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides (both positions), it also includes the fitted pseudo-single correction obtained as system solution.
#' -  matrix_interactions: matrix with logical entries, rows = sgRNAs in position 1 and cols =  sgRNAs in position 2. Each entry indicates if the pair sgRNA pos1 ~ sgRNA pos 2 exists in the given dual KO dataset.
#' -  sys_solution: system solution by the chosen solver (OSQP solver from CVXR, output of \code{\link{get_pairwise_correction}})
#' 
#' @examples
#' @seealso 
#' \code{\link{get_pairwise_correction}}
#' 
#' @export
ccr2.solveLinearSys <- function(
  dual_FC, 
  single_FC_corrected, 
  match_dual_single_seq, 
  pseudo_single_p1_correctedFCs, 
  pseudo_single_p2_correctedFCs, 
  split_graph = TRUE,
  display = TRUE, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = ""
) {
  
  # get maximum of correction from single screens
  unique_ID <- unique(c(dual_FC$sgRNA1_WGE_ID, 
                        dual_FC$sgRNA2_WGE_ID))
  id_keep <-  match_dual_single_seq$ID_single[match_dual_single_seq$ID %in% unique_ID]
  id_keep <- id_keep[!is.na(id_keep)]
  single_FC_corrected_filt <- single_FC_corrected[unique(id_keep), ]
  single_FC_corrected_filt <- single_FC_corrected_filt[!is.na(single_FC_corrected_filt$CHR), ]
  max_corr <- max(max(single_FC_corrected_filt$correction)*2, 1) # linear relationship
  min_corr <- min(min(single_FC_corrected_filt$correction), 0) # single boundary
  
  # 1. create matrix of interactions 
  dual_FC <- dual_FC %>% 
    dplyr::mutate(sgRNA_ID_pair = paste0(sgRNA1_WGE_ID, "~", sgRNA2_WGE_ID)) 

  combinations <- t(sapply(pseudo_single_p1_correctedFCs$sgRNA_ID, 
                           function(x) paste0(x, "~", pseudo_single_p2_correctedFCs$sgRNA_ID)))
  matrix_interactions <- apply(combinations, 2, 
                               function(x) as.numeric(x %in% dual_FC$sgRNA_ID_pair))
  rownames(matrix_interactions) <- pseudo_single_p1_correctedFCs$sgRNA_ID
  colnames(matrix_interactions) <- pseudo_single_p2_correctedFCs$sgRNA_ID
  
  # can be split in multiple componets?
  combinations_edges <- as.vector(combinations) 
  combinations_edges <- combinations_edges[combinations_edges %in% dual_FC$sgRNA_ID_pair]
  combinations_edges <- do.call(rbind, str_split(combinations_edges, pattern = "~"))
  
  g <- igraph::graph_from_edgelist(el = combinations_edges, 
                                   directed = FALSE)
  groups_g <-  igraph::clusters(g)
  max_groups <- which.max(groups_g$csize)
  
  sub_nodes_small <- names(which(groups_g$membership != max_groups))
  sub_nodes_large <- names(which(groups_g$membership == max_groups))
  
  if (split_graph & length(sub_nodes_small) > 0) {
    
    print("Combinations split, largest connected components solved separately")
    
    # solve large system:
    print(paste("Largest connected component with sgRNA n =", length(sub_nodes_large)))
    id_large_row <- rownames(matrix_interactions) %in% sub_nodes_large
    id_large_col <- colnames(matrix_interactions) %in% sub_nodes_large
    
    res_large <- get_pairwise_correction(
      matrix_interactions = matrix_interactions[id_large_row, id_large_col], 
      combinations = combinations[id_large_row, id_large_col], 
      dual_FC = dual_FC, 
      pseudo_single_p1 = pseudo_single_p1_correctedFCs[id_large_row, ], 
      pseudo_single_p2 = pseudo_single_p2_correctedFCs[id_large_col, ], 
      max_corr = max_corr, 
      min_corr = min_corr
    )
    
    # solve small system:
    print(paste("Remaining sgRNA n =", length(sub_nodes_small)))
    id_small_row <- rownames(matrix_interactions) %in% sub_nodes_small
    id_small_col <- colnames(matrix_interactions) %in% sub_nodes_small
    res_small <- get_pairwise_correction(
        matrix_interactions = matrix_interactions[id_small_row, id_small_col], 
        combinations = combinations[id_small_row, id_small_col], 
        dual_FC = dual_FC, 
        pseudo_single_p1 = pseudo_single_p1_correctedFCs[id_small_row, ], 
        pseudo_single_p2 = pseudo_single_p2_correctedFCs[id_small_col, ], 
        max_corr = max_corr, 
        min_corr = min_corr
    )
    
    res_all <- list(
      sys_solution = list(large = res_large$sys_solution, 
                          small = res_small$sys_solution),
      df_corr = rbind(res_large$df_corr, res_small$df_corr), 
      pseudo_single_correction = rbind(res_large$pseudo_single_correction, 
                                       res_small$pseudo_single_correction)
    )
    
  }else{
    
    print("Linear system solved without split of combinations")
    
    res_all <- get_pairwise_correction(
      matrix_interactions = matrix_interactions, 
      combinations = combinations, 
      dual_FC = dual_FC, 
      pseudo_single_p1 = pseudo_single_p1_correctedFCs, 
      pseudo_single_p2 = pseudo_single_p2_correctedFCs, 
      max_corr = max_corr, 
      min_corr = min_corr
    )
    
  }
  
  sys_solution <- res_all$sys_solution
  df_corr <- res_all$df_corr
  pseudo_single_correction <- res_all$pseudo_single_correction
  dual_FC_correctedFC <- dplyr::left_join(dual_FC, df_corr, by = "sgRNA_ID_pair") %>%
    dplyr::mutate(correctedFC = avgFC + correction)
  
  if (saveToFig) {
    display <- TRUE
    file_name_sys <- sprintf("%s%s_solutionSystem.%s", outdir, EXPname, saveFormat)
    file_name_out <- sprintf("%s%s_FC_vs_correctedFC.%s", outdir, EXPname, saveFormat)
    file_name_corr <- sprintf("%s%s_corrections_pseudosingle_vs_pair.%s", outdir, EXPname, saveFormat)
  }
  
  if (display) {
    
    pl_sys <- ggplot(pseudo_single_correction, aes(x = pseudo_single, 
                                                   y = pseudo_single_fitted, 
                                                   color = guide_pos, 
                                                   label = outliers, 
                                                   size = n)) + 
      geom_point(alpha = 0.5) +
      geom_text_repel(size = 3, min.segment.length = 0, color = "black") + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
      theme_bw() + 
      theme(legend.position = "right") + 
      xlab("Pseudo single correction") + 
      ylab("Fitted correction") +
      scale_size(breaks = c(1, 3, seq(5, max(pseudo_single_correction$n), length.out = 4)),
        range = c(1, 6))
    print(pl_sys)
    
    pl_out <- ggplot(dual_FC_correctedFC, aes(x = avgFC, 
                                              y = correctedFC)) + 
      geom_point(alpha = 0.5, size = 1) +
      geom_abline(slope = 1, intercept = 0, 
                  linetype = "dashed", color = "red") + 
      theme_bw() + 
      theme(legend.position = "right") + 
      xlab("Uncorrected logFC") + ylab("Corrected logFC")
    print(pl_out)
    
    ### plot pseudo_single and pairwise correction ###
    # dual_FC_corr_plot <- dual_FC_correctedFC[!is.na(dual_FC_correctedFC$correction), ]
    dual_FC_corr_plot <- dual_FC_correctedFC
    dual_FC_corr_plot$correction_1 <- pseudo_single_correction$pseudo_single[match(dual_FC_corr_plot$sgRNA1_WGE_ID, pseudo_single_correction$sgRNA_ID)]
    dual_FC_corr_plot$correction_2 <- pseudo_single_correction$pseudo_single[match(dual_FC_corr_plot$sgRNA2_WGE_ID, pseudo_single_correction$sgRNA_ID)]
    
    pl_corr <- ggplot(dual_FC_corr_plot, 
                      mapping = aes(x = correction_1, 
                                    y = correction_2, 
                                    color = correction)) + 
      geom_point() + 
      scale_color_viridis_c() + 
      theme_bw() + 
      xlab("Correction pseudo single pos 1") + 
      ylab("Correction pseudo single pos 2")
    print(pl_corr)
    
  }
  
  if (saveToFig) {
    ggsave(filename = file_name_sys, plot = pl_sys, width = 4.7, height = 4.5)
    ggsave(filename = file_name_out, plot = pl_out, width = 4, height = 4)
    ggsave(filename = file_name_corr, plot = pl_corr, width = 5, height = 4)
  }
  
  return(list(
    max_correction = max_corr, 
    min_correction = min_corr,
    dual_FC = dual_FC_correctedFC, 
    pseudo_single = pseudo_single_correction, 
    matrix_interactions = matrix_interactions, 
    sys_solution = sys_solution))
}

#' Plot logFC distributions for before and after CRISPRcleanR^2 correction, separating per dual guides classes
#' 
#' This function generates boxplots and density plots for the distribution of logFC values (before and after CRISPRcleanR^2 correction), 
#' comparing dual guides classes.
#'
#' @param dual_FC_correctedFC Data frame containing information about corrected logFC values for dual KO.
#' @param saveToFig Logical, indicating whether to save the generated plot to a file. Default is FALSE.
#' @param saveFormat File format for saving the plot (e.g., ".pdf", ".png"). Default is ".pdf".
#' @param outdir Directory to save the plot. Default is "./".
#' @param EXPname Experiment name to include in the plot filename.
#'
#' @return NULL
#'
#' @examples
#' @export
ccr2.plotClasses <- function(
  dual_FC_correctedFC, 
  saveToFig = F, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname =""
) {
  
  df_plot <- dual_FC_correctedFC %>% dplyr::filter(!is.na(correction))
  df_plot <- data.frame(avgFC = c(df_plot$avgFC, df_plot$correctedFC), 
                        type = c(rep("pre-CRISPRCleanR^2", nrow(df_plot)), 
                                 rep("post-CRISPRCleanR^2", nrow(df_plot))),
                        info_subtype = rep(df_plot$info_subtype, 2), 
                        info = rep(df_plot$info, 2))
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  
  pl1 <- ggplot(df_plot, mapping = aes(x = info, y = avgFC, fill = type)) + 
    geom_boxplot(alpha = 0.7) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
    ylab("logFC") + 
    theme_bw() + 
    theme(axis.text = element_text(size = 12), 
          axis.title.y = element_blank(), 
          legend.title = element_blank(), 
          legend.position = "bottom") +
    scale_fill_brewer(palette = "Paired") +
    coord_flip() 
  
  pl2 <- ggplot(df_plot, mapping = aes(x = avgFC, fill = type)) + 
    geom_density(alpha = 0.7) +
    geom_vline(xintercept = 0, color = "red",
               linetype = "dashed", size = 0.8) +
    ylab("density") + 
    theme_classic() + 
    theme(axis.title.x = element_blank(), 
          legend.title = element_blank(), 
          legend.position = "none") +
    scale_fill_brewer(palette = "Paired")
  pl <- ggpubr::ggarrange(plotlist = list(pl2, pl1), ncol = 1, 
                          heights = c(0.2, 1), align = "v")
  
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_PreANDPostFC_vs_Classes.%s", 
                         outdir, EXPname, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 6, height = 6)
  }
  
}

#' Add Copy Number Alteration (CNA) information to dual KO data.
#'
#' This function takes a data frame containing dual KO information and a data frame with copy number alteration (CNA) information. 
#' It adds CNA information to the dual KO data based on genomic positions. 
#' It gives info on CN for gene in position 1 and position 2, 
#' as well as the sum, maximum and product of the CN of genes in position 1 and 2.
#'
#' @param CNA A data frame containing copy number alteration information. It must include the following columns:
#'  - CHROM: chromosome id in "chri" format   
#'  - start: start of the CN region, genomic location in base pair
#'  - end: end of the CN region, genomic location in base pair
#'  - C: copy number (absolute) value 
#' @param dual_FC A data frame containing dual KO data.
#' 
#' @return A data frame similar to dual_FC with added CNA information. The dual guides with no CNA info in at least one position are removed.
#' 
#' @examples
#' @export
ccr2.add_CNA <- function(
  CNA, 
  dual_FC
) { 
  
  # remove pairs with NA chr (nontarget)
  dual_FC <- dual_FC %>%
    dplyr::filter(!is.na(sgRNA1_Chr) & !is.na(sgRNA2_Chr))
  
  # use GenomicRanges package to match
  CNA <- CNA %>%
    dplyr::mutate(dplyr::across("CHROM", .fns = \(x) 
                                str_replace(x, pattern = "chr", replacement = "")))
  
  CNA_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    CNA[, c("CHROM", "start", "end")],
    seqnames.field = "CHROM",
    start.field = "start",
    end.field = "end")

  pos1_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    dual_FC[, c("sgRNA1_Chr", "sgRNA1_Start", "sgRNA1_End")],
    seqnames.field = "sgRNA1_Chr",
    start.field = "sgRNA1_Start",
    end.field = "sgRNA1_End")
  
  pos2_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    dual_FC[, c("sgRNA2_Chr", "sgRNA2_Start", "sgRNA2_End")],
    seqnames.field = "sgRNA2_Chr",
    start.field = "sgRNA2_Start",
    end.field = "sgRNA2_End")
  
  pos1_overlap <- GenomicRanges::findOverlaps(
    pos1_genom_range, 
    CNA_genom_range)
  
  pos2_overlap <- GenomicRanges::findOverlaps(
    pos2_genom_range, 
    CNA_genom_range)
  
  dual_FC_withCNA <- dual_FC %>%
    dplyr::mutate(Gene1_CN = NA, Gene2_CN = NA) 
  dual_FC_withCNA$Gene1_CN[queryHits(pos1_overlap)] <- CNA$C[subjectHits(pos1_overlap)]
  dual_FC_withCNA$Gene2_CN[queryHits(pos2_overlap)] <- CNA$C[subjectHits(pos2_overlap)]
  
  dual_FC_withCNA <- dual_FC_withCNA %>% 
    dplyr::mutate(
      Gene1_CN = round(Gene1_CN),
      Gene2_CN = round(Gene2_CN),
      Prod_CN = Gene1_CN * Gene2_CN, 
      Sum_CN = Gene1_CN + Gene2_CN) %>%
    dplyr::filter(!is.na(Prod_CN))
  dual_FC_withCNA$Max_CN <- mapply(function(x, y) max(x,y), 
                                   x = dual_FC_withCNA$Gene1_CN, 
                                   y = dual_FC_withCNA$Gene2_CN, SIMPLIFY = T)
  
  dual_FC_withCNA$Gene1_CN <- factor(dual_FC_withCNA$Gene1_CN, 
                                     levels = sort(unique(dual_FC_withCNA$Gene1_CN)))
  dual_FC_withCNA$Gene2_CN <- factor(dual_FC_withCNA$Gene2_CN, 
                                     levels = sort(unique(dual_FC_withCNA$Gene2_CN)))
  dual_FC_withCNA$Sum_CN <- factor(dual_FC_withCNA$Sum_CN, 
                                   levels = sort(unique(dual_FC_withCNA$Sum_CN)))
  dual_FC_withCNA$Prod_CN <- factor(dual_FC_withCNA$Prod_CN, 
                                    levels = sort(unique(dual_FC_withCNA$Prod_CN)))
  dual_FC_withCNA$Max_CN  <- factor(dual_FC_withCNA$Max_CN, 
                                    levels = sort(unique(dual_FC_withCNA$Max_CN)))  
  
  return(dual_FC_withCNA)
  
}


#' Plot Copy Number Alteration (CNA) information for dual KO data before and after CRISPRcleanr^2 correction.
#'
#' This function creates boxplots to visualize the relationship between CNA and dual KO logFCs. 
#' It outputs two types of boxplots: the first plots the Maximum CNA between position 1 and position 2 versus logFCs
#' It provides insights into how CNAs configurations may affect the logFC of dual guides.
#'
#' @param dual_FC_correctedFC A data frame containing CRISPRcleanR^2 corrected logFCs for dual KO.
#' @param CNA A data frame containing copy number alteration information.
#' @param saveToFig Logical, indicating whether to save the plots to files (default is FALSE).
#' @param saveFormat Character, the format for saving the plots (default is "pdf").
#' @param outdir Character, the directory to save the plots (default is "./").
#' @param EXPname Character, a prefix to add to the saved plot filenames (default is "").
#' @param excludeGene Character vector, genes to be excluded from dual_FC_correctedFC (default is NULL).
#' @param var_to_plot Character, the variable to plot. If "observed" (default), plots logFCs. Otherwise a user-defined variable, such as bliss z-scores.
#'
#' @return NULL
#'
#' @examples
#' @seealso
#' \code{\link{ccr2.add_CNA}}
#'
#' @export
ccr2.plotCNA <- function(
  dual_FC_correctedFC, 
  CNA, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = "", 
  excludeGene = NULL, 
  var_to_plot = "observed"
) {
  
  dual_FC_withCNA <- ccr2.add_CNA(CNA = CNA, dual_FC = dual_FC_correctedFC) %>%
    dplyr::filter(!is.na(correction))
 
  if (var_to_plot == "observed") {
    var_name <- "avgFC"
    var_name_corrected <- "correctedFC"
    ylab_name <- "logFC"
  }else{
    var_name <- var_to_plot
    var_name_corrected <- paste(var_to_plot, "corrected", sep = "_")
    ylab_name <- var_to_plot
  }
  
  df_plot <- data.frame(info_subtype = rep(dual_FC_withCNA$info_subtype, 2),
                        logFC = c(dual_FC_withCNA[, var_name, drop = TRUE], 
                                  dual_FC_withCNA[, var_name_corrected, drop = TRUE]), 
                        Sum_CN = rep(dual_FC_withCNA$Sum_CN, 2), 
                        Gene1_CN = rep(dual_FC_withCNA$Gene1_CN, 2), 
                        Gene2_CN = rep(dual_FC_withCNA$Gene2_CN, 2), 
                        Prod_CN = rep(dual_FC_withCNA$Prod_CN, 2), 
                        Max_CN = rep(dual_FC_withCNA$Max_CN, 2), 
                        type = c(rep("pre-CRISPRCleanR^2", nrow(dual_FC_withCNA)), 
                                 rep("post-CRISPRCleanR^2", nrow(dual_FC_withCNA))), 
                        Gene1 = rep(dual_FC_withCNA$Gene1, 2), 
                        Gene2 = rep(dual_FC_withCNA$Gene2, 2))
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  df_plot$comb_CN <- sprintf("%s ~ %s", df_plot$Gene1_CN, df_plot$Gene2_CN)
  
  if (!is.null(excludeGene)) {
    df_plot <-  df_plot %>% 
      dplyr::filter(!Gene1 %in% excludeGene, !Gene2 %in% excludeGene)
    title_plot <- paste("Exclude", paste(excludeGene, collapse = ","))
    add_to_name <- paste0("_rm", paste(excludeGene, collapse = "_"))
  } else {
    title_plot = NULL
    add_to_name = "" 
  }
  
  pl_CN <- ggplot(df_plot, aes(x = Max_CN, y = logFC, 
                               fill = type)) + 
    geom_boxplot(outlier.size = 0.5) + 
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
    theme_bw() + 
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom", 
          axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("Max CN guide1 & guide2") + 
    ylab(ylab_name) + 
    ggtitle(title_plot)
  
  df_plot$Gene2_CN <- sprintf("CN guide2: %s", df_plot$Gene2_CN)
  df_plot$Gene2_CN <- factor(
    df_plot$Gene2_CN, 
    levels = sprintf("CN guide2: %s", sort(as.numeric(levels(dual_FC_withCNA$Gene2_CN)))))
  
  pl_CN_comb <- ggplot(df_plot, aes(x = Gene1_CN, y = logFC, fill = type)) + 
    geom_boxplot(outlier.size = 0.5) + 
    theme_bw() + 
    facet_wrap(.~Gene2_CN) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", linewidth = 0.5) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", 
          legend.title = element_blank(), 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("CN guide1") + ylab(ylab_name) + ggtitle(title_plot)
  
  print(pl_CN)
  print(pl_CN_comb)
  
  if (saveToFig) {
    file_name_maxCN <- sprintf("%s%s_MaxCN_vs_%s%s.%s", 
                               outdir, EXPname, ylab_name, add_to_name, saveFormat)
    file_name_combCN <- sprintf("%s%s_CombCN_vs_%s%s.%s", 
                                outdir, EXPname, ylab_name, add_to_name, saveFormat)
    ggsave(filename = file_name_maxCN, plot = pl_CN, width = 5, height = 4.5)
    ggsave(filename = file_name_combCN, plot = pl_CN_comb, width = 9, height = 7)
  }
  
}

#' Plot density of dual KO logFCs, separating for extreme Copy Number Alteration (CNA) events.
#'
#' This function creates density plots to visualize the distribution of dual KO logFCs 
#' separating for extreme CNA (higher than a certain threshold). The density is shown before and after CRISRPcleanR^2 correction.
#' It provides insights into how CNAs may influence the logFC distribution.
#'
#' @param dual_FC_correctedFC A data frame containing CRISPRcleanR^2 corrected logFCs for dual KO.
#' @param CNA A data frame containing copy number alteration information.
#' @param saveToFig Logical, indicating whether to save the plots to files (default is FALSE).
#' @param saveFormat Character, the format for saving the plots (default is "pdf").
#' @param outdir Character, the directory to save the plots (default is "./").
#' @param EXPname Character, a prefix to add to the saved plot filenames (default is "").
#' @param excludeGene Character vector, genes to be excluded from dual_FC_correctedFC (default is NULL).
#' @param CN_thr Numeric (absolute), a threshold for defining extreme CNA.
#' @param var_to_plot Character, the variable to plot. If "observed" (default), plots logFCs. Otherwise a user-defined variable, such as bliss z-scores.
#'
#' @return NULL
#' 
#' @examples
#' @seealso
#' \code{\link{ccr2.add_CNA}}
#'
#' @export
ccr2.plotCNAdensity <- function(
  dual_FC_correctedFC, 
  CNA, 
  saveToFig = F, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = "", 
  excludeGene = NULL, 
  CN_thr, 
  var_to_plot = "observed") {
  
  if (var_to_plot == "observed") {
    var_name <- "avgFC"
    var_name_corrected <- "correctedFC"
    xlab_name <- "logFC"
  }else{
    var_name <- var_to_plot
    var_name_corrected <- paste(var_to_plot, "corrected", sep = "_")
    xlab_name <- var_to_plot
  }
  
  dual_FC_withCNA <- ccr2.add_CNA(CNA = CNA, dual_FC = dual_FC_correctedFC) %>%
    dplyr::filter(!is.na(correction)) %>%
    dplyr::mutate(CN_class = dplyr::case_when(
      as.numeric(as.character(Gene1_CN)) >= !!(CN_thr) | 
      as.numeric(as.character(Gene2_CN)) >= !!(CN_thr) ~ 
        sprintf("CN guide1 or guide2 >= %i", CN_thr), 
      TRUE ~ "Others"))

  df_plot <- data.frame(
    logFC = c(dual_FC_withCNA[, var_name, drop = T], 
              dual_FC_withCNA[, var_name_corrected, drop = T]),
    CN_class = rep(dual_FC_withCNA$CN_class, 2), 
    type = c(rep("pre-CRISPRCleanR^2", nrow(dual_FC_withCNA)), 
             rep("post-CRISPRCleanR^2", nrow(dual_FC_withCNA))), 
    Gene1 = rep(dual_FC_withCNA$Gene1, 2), 
    Gene2 = rep(dual_FC_withCNA$Gene2, 2))
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  
  if (!is.null(excludeGene)) {
    df_plot <-  df_plot %>% 
      dplyr::filter(!Gene1 %in% excludeGene, !Gene2 %in% excludeGene)
    title_plot <- paste("Exclude", paste(excludeGene, collapse = ","))
    add_to_name <- paste0("_rm", paste(excludeGene, collapse = "_"))
  } else {
    title_plot = NULL
    add_to_name = ""
  }
  
  pl <- ggplot(df_plot, aes(x = logFC, fill = CN_class)) + 
    geom_density(alpha = 0.7) + 
    theme_bw() + 
    facet_wrap(.~type, nrow = 2) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", size = 0.8) +
    theme(axis.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          axis.title = element_text(size = 13),
          legend.position = "bottom", legend.title = element_blank()) +
    scale_fill_manual(values = c("red", "grey70")) +
    xlab(xlab_name) + 
    ggtitle(title_plot)
  
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_density%s_highCN%s.%s", 
                         outdir, EXPname, xlab_name, add_to_name, saveFormat)
    ggplot2::ggsave(filename =  file_name, plot = pl, width = 5, height = 5)
  }
  
}

#' Plot Copy Number Alteration (CNA) information for single KO data considering only guides matching with dual KO data, before and after CRISPRcleanR correction.
#'
#' This function creates two plots:
#' 1. A boxplot showing the distribution of logFC for single KO guides in common with dual screen, stratified by copy number (CN).
#' 2. A scatter plot comparing uncorrected logFC with corrected logFC for single KO guides in common with dual screen.
#'
#' @param dual_FC_correctedFC A data frame containing CRISPRcleanR^2 corrected logFCs for dual KO.
#' @param match_dual_single_seq Data frame containing matched sequences between the dual and single libraries (output of ccr2.matchDualandSingleSeq).
#' @param single_correctedFCs Data frame containing single screen data with CRISPRcleanR corrected logFC.
#' @param CNA A data frame containing copy number alteration information.
#' @param saveToFig Logical, whether to save the plots to files (default is FALSE).
#' @param saveFormat Character, the format for saving the plots (default is "pdf").
#' @param outdir Character, the directory path for saving the plots (default is "./").
#' @param EXPname Character, an identifier for the experiment (default is "").
#'
#' @return A data frame containing single screen data for matching guides in dual screen data with annotated copy number information.
#' 
#' @example
#' @seealso
#' \code{\link{ccr2.add_CNA}}
#'  
#' @export
ccr2.plotMatchingSingle <- function(
  dual_FC_correctedFC,  
  match_dual_single_seq, 
  single_correctedFCs, 
  CNA, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = ""
) {
  
  id_corr <- !is.na(dual_FC_correctedFC$correction)
  unique_ID <- unique(c(dual_FC_correctedFC$sgRNA1_WGE_ID[id_corr], 
                        dual_FC_correctedFC$sgRNA2_WGE_ID[id_corr]))
  
  id_keep <-  match_dual_single_seq$ID_single[match_dual_single_seq$ID %in% unique_ID]
  id_keep <- id_keep[!is.na(id_keep)]
  single_correctedFCs_filt <- single_correctedFCs[unique(id_keep), ]
  single_correctedFCs_filt <- single_correctedFCs_filt[!is.na(single_correctedFCs_filt$CHR), ]
  
  # match with CNA info
  CNA <- CNA %>%
    dplyr::mutate(dplyr::across("CHROM", .fns = \(x) 
                                str_replace(x, pattern = "chr", replacement = "")))
  
  CNA_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    CNA[, c("CHROM", "start", "end")],
    seqnames.field = "CHROM",
    start.field = "start",
    end.field = "end")
  
  single_genom_range <- GenomicRanges::makeGRangesFromDataFrame(
    single_correctedFCs_filt[, c("CHR", "startp", "endp")],
    seqnames.field = "CHR",
    start.field = "startp",
    end.field = "endp")
  
  overlap <- GenomicRanges::findOverlaps(
    single_genom_range, 
    CNA_genom_range)
  
  FC_withCNA <- single_correctedFCs_filt %>%
    dplyr::mutate(Gene_CN = NA) 
  FC_withCNA$Gene_CN[queryHits(overlap)] <- CNA$C[subjectHits(overlap)]
  FC_withCNA <- FC_withCNA %>%
    dplyr::mutate(Gene_CN = round(Gene_CN))
  
  df_plot <- data.frame(
    logFC = c(FC_withCNA$avgFC, FC_withCNA$correctedFC), 
    CN = rep(FC_withCNA$Gene_CN, 2), 
    type = c(rep("pre-CRISPRCleanR^2", nrow(FC_withCNA)), 
             rep("post-CRISPRCleanR^2", nrow(FC_withCNA))), 
    Gene = rep(FC_withCNA$genes, 2))
  df_plot$type <- factor(df_plot$type, 
                         levels = c("pre-CRISPRCleanR^2", "post-CRISPRCleanR^2"))
  
  pl_CN <- ggplot(df_plot, aes(x = as.factor(CN), y = logFC, fill = type)) + 
    geom_boxplot() +
    theme_bw() + 
    # geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
    theme(axis.text = element_text(size = 12),
          legend.position = "bottom", 
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("CN guide") + 
    ylab("logFC") +
    ggtitle("Single screen", 
            subtitle = "guides in common with dual screen") 
  print(pl_CN)
  
  pl_out <- ggplot(FC_withCNA , aes(x = avgFC, y = correctedFC)) + 
    geom_point(alpha = 0.7, size = 1) +
    geom_abline(slope = 1, intercept = 0, 
                linetype = "dashed", color = "red") + 
    theme_bw() + 
    theme(legend.position = "right") + 
    xlab("Uncorrected logFC") + 
    ylab("Corrected logFC") + 
    ggtitle("Single screen", 
            subtitle = "guides in common with dual screen") 
  print(pl_out)
  
  if (saveToFig) {
    file_name_out <- sprintf("%s%s_singlescreen_FC_vs_correctedFC.%s", 
                             outdir, EXPname, saveFormat)
    ggsave(filename = file_name_out, plot = pl_out, width = 4, height = 4)
    
    file_name_CN <- sprintf("%s%s_singlescreen_CN_vs_FC.%s", 
                            outdir, EXPname, saveFormat)
    ggsave(filename =  file_name_CN, plot = pl_CN, width = 4, height = 4)
  }
  
  return(FC_withCNA)
}

#' Plot Bliss synergy z-scores before and after CRISPRcleanR^2 correction.
#'
#' This function generates a scatter plot comparing Bliss z-scores before and
#' after correction CRISPRcleanR^2. The points are colored according dual guide class (e.g. "PositiveControls").
#'
#' @param df A data frame containing Bliss z-scores before and after correction (output of ccr2.compute_bliss).
#'        It should have columns 'bliss_zscore', 'bliss_zscore_corrected', and 'info'.
#' @param saveToFig Logical, indicating whether to save the plot to a file (default is FALSE)
#' @param saveFormat Character, the format in which to save the plot (default is "pdf").
#' @param outdir Character, the directory where the plot file will be saved (default is "./").
#' @param EXPname Character, an additional name to include in the plot file name (default is "").
#'
#' @return NULL
#'
#' @examples
#' @seealso 
#' \code{\link{ccr2.compute_bliss}}
#' 
#' @export
ccr2.plot_correction <- function(df, 
                                 saveToFig = FALSE, 
                                 saveFormat = "pdf",
                                 outdir = "./", 
                                 EXPname = "") {
  
  # plot synergy before and after correction:
  pl <- ggplot(df, aes(x = bliss_zscore, 
                       y = bliss_zscore_corrected, 
                       color = info)) + 
    geom_abline(slope = 1, intercept = 0, 
                linetype = "dashed", color = "black") + 
    geom_point(alpha = 0.5, size = 1) +
    scale_color_brewer(palette = "Set3" ) +
    theme_bw() + 
    theme(legend.position = "right") + 
    xlab("Uncorrected Bliss z-score") + 
    ylab("Corrected Bliss z-score")
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_bliss_original_vs_corrected.%s", 
                             outdir, EXPname, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 6, height = 4)
  }
  
}

#' Plot Bliss fit for expected vs observed dual KO logFC.
#'
#' This function generates a scatter plot comparing expected vs observed logFC in dual KO.
#' Two plots are shown, one before and one after CRISPRcleanR^2 correction. 
#'
#' @param dual_FC A data frame containing information about dual KO logFCs, with observed and expected logFCs in dual guide pairs based on bliss model (output of ccr2.compute_bliss).
#' @param saveToFig Logical, indicating whether to save the plot to a file .
#' @param saveFormat Character, the format in which to save the plot (default is "pdf", NOTE: the generated output will be heavy in size).
#' @param outdir Character, the directory where the plot file will be saved (default is "./").
#' @param EXPname Character, an additional name to include in the plot file name (default is "").
#'
#' @return NULL
#'
#' @examples
#' @seealso 
#' \code{\link{ccr2.compute_bliss}}
#' 
#' @export
ccr2.plot_bliss_fit <- function(dual_FC, 
                                saveToFig = FALSE, 
                                saveFormat = "pdf",
                                outdir = "./", 
                                EXPname = "") {
  
  df_plot_1 <- dual_FC %>%
    dplyr::filter(!is.na(correction)) %>%
    dplyr::select(sgRNA_ID_pair, info, info_subtype, avgFC, expected_pairFC) %>%
    dplyr::mutate(type = "Uncorrected") %>%
    dplyr::rename(observed = avgFC, expected = expected_pairFC)
  
  df_plot_2 <- dual_FC %>%
    dplyr::filter(!is.na(correction)) %>%
    dplyr::select(sgRNA_ID_pair, info, info_subtype, correctedFC, expected_pairFC_corrected) %>%
    dplyr::mutate(type = "Corrected") %>%
    dplyr::rename(observed = correctedFC, expected = expected_pairFC_corrected)
  
  df_plot <- dplyr::bind_rows(df_plot_1, df_plot_2)
  
  pl <- ggplot(df_plot , aes(x = expected, y = observed, 
                             color = info)) + 
    geom_point(alpha = 0.5, size = 0.5) +
    facet_wrap(.~type) + 
    geom_smooth(method = "loess", 
                inherit.aes = FALSE,
                se = FALSE, 
                aes(x = expected, y = observed)) + 
    theme_bw() + 
    geom_abline(intercept = 0, slope = 1, 
                linetype = "dashed", color = "red") +
    theme(legend.position = "bottom", 
          plot.title = element_text(hjust = 0.5)) + 
    xlab("Expected logFC") + ylab("Observed logFC")
  print(pl)
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_blissfit_expected_vs_observed.%s", 
                         outdir, EXPname, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 9, height = 5)
  }
  
}

#' Plot Bliss Z-score vs logFC to identify synergistic and lethal pairs.
#'
#' This function generates a scatter plot comparing Bliss Z-score with logFC in dual KO screen. 
#' Here, scaled logFCs are used, centering positive controls to -1 and negative controls to 0 (see ccr2.scale_pos_neg).
#' It can be used to plot results both before and after CRISPRcleanR^2 correction. 
#' Dual guides with bliss z-score and logFCs passing pre-specified thresholds are stored.
#'
#' @param dual_FC A data frame containing information about dual KO logFCs.
#'        It should have columns 'avgFC_scaled', 'correctedFC_scaled' (see ccr2.scale_pos_neg), and 'bliss_zscore', 
#'        'bliss_zscore_corrected' (see ccr2.compute_bliss), and 'info'. 
#' @param corrected Logical, indicating whether to use CRISPRcleanR^2 corrected values for logFC and Bliss Z-score (default is FALSE).
#' @param THR_FC Numeric, the threshold for logFC below which points will be considered lethal and stored (default is -1).
#' @param THR_BLISS Numeric, the threshold for Bliss Z-score below which points will be considered synergistic and stored (default is -1).
#' @param saveToFig Logical, indicating whether to save the plot to a file (default is FALSE).
#' @param saveFormat Character, the format in which to save the plot (default is "pdf").
#' @param outdir Character, the directory where the plot file will be saved (default is "./").
#' @param EXPname Character, an additional name to include in the plot file name (default is "").
#'
#' @return A data frame containing dual guide pairs classified as lethal and synergistic based on the specified thresholds. Same structure as input dual_FC data frame.
#'
#' @examples
#' @seealso 
#' \code{\link{ccr2.compute_bliss}},
#' \code{\link{ccr2.scale_pos_neg}}
#' 
#' @export
ccr2.plot_bliss_vs_FC <- function(dual_FC, 
                                  corrected = FALSE, 
                                  THR_FC = -1,
                                  THR_BLISS = -1, 
                                  saveToFig = FALSE, 
                                  saveFormat = "pdf",
                                  outdir = "./", 
                                  EXPname = ""){
  
  dual_FC <- dual_FC[!is.na(dual_FC$correction),]
  
  name_var_y <- "avgFC_scaled" 
  # name_var_y <- "avgFC" # use the actual values, not scaled
  name_var_x <- "bliss_zscore"
  title_pl <- "Uncorrected"
  if (corrected) {
    name_var_y <- "correctedFC_scaled" 
    # name_var_y <- "correctedFC" # use the actual values, not scaled
    name_var_x <- "bliss_zscore_corrected"
    title_pl <- "Corrected"
  }
  
  pl <- ggplot(dual_FC, 
         aes(x = get(name_var_x) , 
             y = get(name_var_y), color = info)) + 
    geom_point(alpha = 0.5) + 
    theme_bw() + 
    geom_vline(xintercept = THR_BLISS, color = "black", linetype = "dashed") + 
    geom_hline(yintercept = THR_FC, color = "black", linetype = "dashed") + 
    theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5)) + 
    xlab("Bliss Z-score") + ylab("logFC") + 
    ggtitle(title_pl)
  print(pl)
  
  top_res <- dual_FC %>%
    dplyr::filter(get(name_var_x) < THR_BLISS & get(name_var_y) < THR_FC 
                  & !is.na(get(name_var_x))) %>%
    dplyr::arrange(get(name_var_x))
  
  if (saveToFig) {
    file_name <- sprintf("%s%s_blisszscore_vs_FC_%s.%s", 
                         outdir, EXPname, title_pl, saveFormat)
    ggsave(filename = file_name, plot = pl, width = 5, height = 6)
  }
  
  return(top_res)
  
}


#' Run CRISPRcleanR^2 analysis for dual KO screen ("double-cut" pairs).
#'
#' This function performs a comprehensive CRISPRcleanR^2 analysis in dual KO screen for "double-cut" pairs,
#' meaning both guides in position1 and position 2 target an existing genomic region. The following steps are performed:
#' - creation of pseudo single scores from dual KO logFCs, 
#' - modeling of single KO logFCs vs. pseudo single logFCs for matching guides, 
#' - injection of data into pseudo single logFCs space from genome-wide single KO screen,
#' - application of CRISPRcleanR to pseudo-single logFCs (genome-wide from original plus injected data), 
#' - approximate linear systems with box constraints to obtain pairwise correction from pseudo-single correction.
#' It returns information such as the maximum and minimum corrections allowed, 
#' model performance for single VS pseudo-single, 
#' identified pseudo-single segments and correction, and
#' final dual logFC data frame with CRISPRcleanR^2 correction.
#' 
#' 
#' @param dual_FC The data frame containing dual KO logFCs (output of ccr2.logFCs2chromPos).
#' @param single_correctedFCs The data frame containing corrected single KO logFCs already corrected via CRISPRcleanR (output of ccr.run_complete).
#' @param libraryAnnotation_dual Data frame containing information from the dual library (one row per guide pair).
#' @param match_dual_single_seq A data frame with one row per sgRNA in dual library (union of position 1 and position 2), including matched single library IDs and sequences (output of ccr2.matchDualandSingleSeq).
#' @param EXPname A character string specifying the experiment name (default is "").
#' @param saveToFig Logical, whether to save figures (default is FALSE).
#' @param display Logical, whether to display figures (default is TRUE).
#' @param saveFormat The format in which to save figures (default is NULL).
#' @param outdir The directory to save figures (default is "./").
#' @param correctGW A character specifying the type of pseudo-single logFCs mean centering after data-injection ("CHR": per chromosome, "GW": genome-wide, or NULL for no centering). See ccr2.injectData.
#' @param ... Additional parameters to be passed.
#'
#' @return A list containing maximum and minimum corrections allowed, model performance and estimated 
#'         for single VS pseudo-single, identified pseudo-single segments and correction, and final dual logFC data frame with CRISPRcleanR^2 correction.
#'
#' @examples
#' @seealso
#' \code{\link{ccr2.createPseudoSingle_combine}},
#' \code{\link{ccr2.modelSingleVSPseudoSingle}},
#' \code{\link{ccr2.injectData}},
#' \code{\link{ccr.GWclean}},
#' \code{\link{ccr2.filterGWclean}},
#' \code{\link{ccr2.solveLinearSys}},
#' \code{\link{ccr2.logFCs2chromPos}},
#' \code{\link{ccr.run_complete}},
#' \code{\link{ccr2.matchDualandSingleSeq}}
#'
#' @export
ccr2.run <- function(
  dual_FC,
  single_correctedFCs, 
  libraryAnnotation_dual, 
  match_dual_single_seq, 
  EXPname = "",
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = NULL, 
  outdir = "./", 
  correctGW, 
  ...
) {
  
  # exclude non-target pairs
  # id_nontarget <- grepl("NONTARGET", dual_FC$info_subtype)
  id_nontarget <- (is.na(dual_FC$sgRNA1_Chr) | is.na(dual_FC$sgRNA2_Chr))
  dual_FC <- dual_FC[!id_nontarget,]
  
  ### create pseudo single scores ###
  # collapse to single info, add single logFC results baesd on sequence matching
  dual_pseudo_single_FC <- ccr2.createPseudoSingle_combine(
    dual_FC = dual_FC, 
    single_FC = single_correctedFCs, 
    match_dual_single_seq = match_dual_single_seq, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
  print("Computed pseudo single")
  
  ##########################################################
  ### create model to convert single to multiple screens ###
  model_guide1 <- ccr2.modelSingleVSPseudoSingle(
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    display = display,
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
  model_guide2 <- ccr2.modelSingleVSPseudoSingle(
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    display = display, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
  models_performance <- data.frame(
    position = c("guide1", "guide2"), 
    n = c(nrow(model_guide1$model$model), nrow(model_guide2$model$model)), 
    cor = c(cor(model_guide1$model$fitted.values, model_guide1$model$model$avgFC), 
            cor(model_guide2$model$fitted.values, model_guide2$model$model$avgFC)), 
    R2 = c(with(model_guide1$model_summary, 1 - deviance/null.deviance),
           with(model_guide2$model_summary, 1 - deviance/null.deviance)))
  
  models_estimates <- rbind(as.data.frame(coef(summary(model_guide1$model))) %>% 
                        dplyr::mutate(position = "guide1") %>%
                        dplyr::mutate(feature = rownames(coef(summary(model_guide1$model))), 
                                      .before = "Estimate"), 
                      as.data.frame(coef(summary(model_guide2$model))) %>% 
                        dplyr::mutate(position = "guide2") %>%
                        dplyr::mutate(feature = rownames(coef(summary(model_guide2$model))), 
                                      .before = "Estimate"))
  
  ### inject data ###
  dataInjection_guide1 <- ccr2.injectData(
    model_single_to_pseudo = model_guide1$model, 
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA1, 
    single_FC = single_correctedFCs, 
    guide_id = 1, 
    correctGW = correctGW
  )
  
  dataInjection_guide2 <- ccr2.injectData(
    model_single_to_pseudo = model_guide1$model, 
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA2, 
    single_FC = single_correctedFCs, 
    guide_id = 2, 
    correctGW = correctGW
  )
  
  print("Injected data on pseudo single position 1 and 2 spaces")
  
  ##########################################
  ### apply CRISPRCleanR to single lines ###
  dataInjection_guide1_correctedFCs <- ccr.GWclean(
    dataInjection_guide1,
    display = display, 
    label = sprintf("guide1_%s", EXPname), 
    saveTO = outdir)
  
  dataInjection_guide2_correctedFCs <- ccr.GWclean(
    dataInjection_guide2,
    display = display, 
    label = sprintf("guide2_%s", EXPname), 
    saveTO = outdir)
  
  # consider only guides in dual screen
  pseudo_single_p1_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide1_correctedFCs$corrected_logFCs,
    dataInjection_segments = dataInjection_guide1_correctedFCs$segments,
    pseudo_single = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display)
  
  pseudo_single_p2_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide2_correctedFCs$corrected_logFCs,
    dataInjection_segments = dataInjection_guide2_correctedFCs$segments,
    pseudo_single = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display)
  
  # save identified segments
  tmp1 <- pseudo_single_p1_correctedFCs$pseudo_single_seg %>%
    dplyr::mutate(guide_pos = "1")
  tmp2 <- pseudo_single_p2_correctedFCs$pseudo_single_seg %>%
    dplyr::mutate(guide_pos = "2")
  pseudo_single_segments <- rbind(tmp1, tmp2)
  
  print("CRISPRcleanR applied to GW pseudo singles (position 1 and 2)")
  
  ##################################################
  ### solve linear system to get pair correction ###
  sys_solution <- ccr2.solveLinearSys(
    single_FC_corrected = single_correctedFCs,
    match_dual_single_seq = match_dual_single_seq,
    dual_FC = dual_FC, 
    pseudo_single_p1_correctedFCs = pseudo_single_p1_correctedFCs$pseudo_single_FC, 
    pseudo_single_p2_correctedFCs = pseudo_single_p2_correctedFCs$pseudo_single_FC, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display, 
    ...)
  print("System solved, collpased correction converted to pairwise correction")
  
  pseudo_single_correction <- sys_solution$pseudo_single
  dual_FC_correctedFC <- sys_solution$dual_FC
  max_correction <- sys_solution$max_correction
  min_correction <- sys_solution$min_correction
  
  return(
    list(max_correction = max_correction, 
         min_correction = min_correction, 
         model_perf = models_performance,
         model_est = models_estimates,
         pseudo_single_segments = pseudo_single_segments,
         pseudo_single = pseudo_single_correction, 
         dual_FC = dual_FC_correctedFC)
    )
}



#' Run CRISPRcleanR^2 analysis for dual KO screen ("non-target" pairs).
#'
#' This function performs a comprehensive CRISPRcleanR^2 analysis in dual KO screen for "non-target" pairs,
#' meaning one between the guides in position1 and position 2 is targeting a NON existing genomic-region. 
#' It consider only guide pairs having a non-target guide and performs the following steps:
#' - creation of pseudo single scores from dual KO logFCs (only matched with non-targets), 
#' - modeling of single KO logFCs vs. pseudo single logFCs for matching guides, 
#' - injection of data into pseudo single logFCs (only matched with non-targets) space from genome-wide single KO screen,
#' - application of CRISPRcleanR to pseudo-single logFCs (genome-wide from original plus injected data), 
#' - assignment of dual KO correction based on the correction in the new space.
#' Differently from \code{\link{ccr2.run}}, there is no need to approximate the linear system to retrieve the correction for dual KO.
#' Since each guide pair is actually targeting only one genomic-region, we can assign the correction obtained from CRISPRcleanR directly.
#' It returns information such as the model performance for single VS pseudo-single, identified pseudo-single segments, and
#' final dual logFC data frame with CRISPRcleanR^2 correction (only for non-target pairs).
#'
#' @param dual_FC The data frame containing dual KO logFCs (output of ccr2.logFCs2chromPos).
#' @param single_correctedFCs The data frame containing corrected single KO logFCs already corrected via CRISPRcleanR (output of ccr.run_complete).
#' @param libraryAnnotation_dual Data frame containing information from the dual library (one row per guide pair).
#' @param match_dual_single_seq TA data frame with one row per sgRNA in dual library (union of position 1 and position 2), including matched single library IDs and sequences (output of ccr2.matchDualandSingleSeq).
#' @param EXPname A character string specifying the experiment name (default is "").
#' @param saveToFig Logical, whether to save figures (default is FALSE).
#' @param display Logical, whether to display figures (default is TRUE).
#' @param saveFormat The format in which to save figures (default is NULL).
#' @param outdir The directory to save figures (default is "./").
#' @param correctGW A character specifying the type of pseudo-single logFCs mean centering after data-injection ("CHR": per chromosome, "GW": genome-wide, or NULL for no centering). See ccr2.injectData.
#' 
#' @return A list containing model performance and estimated 
#'         for single VS pseudo-single, identified pseudo-single segments, 
#'         and final dual logFC data frame with CRISPRcleanR^2 correction.
#'
#' @examples
#' @seealso
#' \code{\link{ccr2.createPseudoSingle_combine}},
#' \code{\link{ccr2.modelSingleVSPseudoSingle}},
#' \code{\link{ccr2.injectData}},
#' \code{\link{ccr.GWclean}},
#' \code{\link{ccr2.filterGWclean}},
#' \code{\link{ccr2.logFCs2chromPos}},
#' \code{\link{ccr.run_complete}},
#' \code{\link{ccr2.matchDualandSingleSeq}}
#'
#' @export
ccr2.run_nontarget <- function(
  dual_FC,
  single_correctedFCs, 
  libraryAnnotation_dual, 
  match_dual_single_seq, 
  EXPname = "",
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = NULL, 
  outdir = "./", 
  correctGW
) {
  
  # process only gene_nontarget or nontarget_gene
  id_nontarget <- (is.na(dual_FC$sgRNA1_Chr) | is.na(dual_FC$sgRNA2_Chr)) & 
    !(is.na(dual_FC$sgRNA1_Chr) & is.na(dual_FC$sgRNA2_Chr))
    
  dual_FC_nt <- dual_FC[id_nontarget,]
  
  dual_pseudo_single_FC <- ccr2.createPseudoSingle_combine(
    dual_FC = dual_FC_nt, 
    single_FC = single_correctedFCs, 
    match_dual_single_seq = match_dual_single_seq,
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  model_guide1 <- ccr2.modelSingleVSPseudoSingle(
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  model_guide2 <- ccr2.modelSingleVSPseudoSingle(
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  models_performance <- data.frame(
    position = c("guide1", "guide2"), 
    n = c(nrow(model_guide1$model$model), nrow(model_guide2$model$model)), 
    cor = c(cor(model_guide1$model$fitted.value, model_guide1$model$model$avgFC), 
            cor(model_guide2$model$fitted.value, model_guide2$model$model$avgFC)), 
    R2 = c(with(model_guide1$model_summary, 1 - deviance/null.deviance),
           with(model_guide2$model_summary, 1 - deviance/null.deviance)))
  
  models_estimates <- rbind(as.data.frame(coef(summary(model_guide1$model))) %>% 
                        dplyr::mutate(position = "guide1") %>%
                        dplyr::mutate(feature = rownames(coef(summary(model_guide1$model))), 
                                      .before = "Estimate"), 
                      as.data.frame(coef(summary(model_guide2$model))) %>% 
                        dplyr::mutate(position = "guide2") %>%
                        dplyr::mutate(feature = rownames(coef(summary(model_guide2$model))), 
                                      .before = "Estimate"))
  
  ### inject data ###
  dataInjection_guide1 <- ccr2.injectData(
    model_single_to_pseudo = model_guide1$model, 
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA1, 
    single_FC = single_correctedFCs, 
    guide_id = 1, 
    correctGW = correctGW
  )
  
  dataInjection_guide2 <- ccr2.injectData(
    model_single_to_pseudo = model_guide1$model, 
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA2, 
    single_FC = single_correctedFCs, 
    guide_id = 2, 
    correctGW = correctGW
  )
  
  ### apply CRISPRCleanR to single lines ###
  dataInjection_guide1_correctedFCs <- ccr.GWclean(
    dataInjection_guide1,
    display = display, 
    label = sprintf("guide1_%s", EXPname), 
    saveTO = sprintf("%sNONTARGET_PAIR_", outdir))
  
  dataInjection_guide2_correctedFCs <- ccr.GWclean(
    dataInjection_guide2,
    display = display, 
    label = sprintf("guide2_%s", EXPname), 
    saveTO = sprintf("%sNONTARGET_PAIR_", outdir))
  
  # consider only guides in dual screen
  pseudo_single_p1_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide1_correctedFCs$corrected_logFCs,
    dataInjection_segments = dataInjection_guide1_correctedFCs$segments,
    pseudo_single = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  pseudo_single_p2_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide2_correctedFCs$corrected_logFCs,
    dataInjection_segments = dataInjection_guide2_correctedFCs$segments,
    pseudo_single = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  # assign results (logFC)
  dual_FC_nt <- dual_FC_nt %>% 
    dplyr::mutate(sgRNA_ID_pair = paste0(sgRNA1_WGE_ID, "~", sgRNA2_WGE_ID))
  
  tmp1 <- pseudo_single_p1_correctedFCs$pseudo_single_FC %>% 
    dplyr::rename(sgRNA1_WGE_ID = sgRNA_ID) %>%
    dplyr::select(sgRNA1_WGE_ID, correction) 
  tmp1 <- dplyr::inner_join(dual_FC_nt, tmp1, by = "sgRNA1_WGE_ID")
  
  tmp2 <- pseudo_single_p2_correctedFCs$pseudo_single_FC %>% 
    dplyr::rename(sgRNA2_WGE_ID = sgRNA_ID) %>%
    dplyr::select(sgRNA2_WGE_ID, correction) 
  tmp2 <- dplyr::inner_join(dual_FC_nt, tmp2, by = "sgRNA2_WGE_ID")
  
  dual_FC_correctedFC <- rbind(tmp1, tmp2) %>%
    dplyr::mutate(correctedFC = avgFC + correction)
  
  # assign results (segments)
  tmp1 <- pseudo_single_p1_correctedFCs$pseudo_single_seg %>%
    dplyr::mutate(guide_pos = "1")
  
  tmp2 <- pseudo_single_p2_correctedFCs$pseudo_single_seg %>%
    dplyr::mutate(guide_pos = "2")
  
  pseudo_single_segments <- rbind(tmp1, tmp2)
  
  if (!all(dual_FC_nt$ID %in% dual_FC_correctedFC$ID)) {
    stop("In NONTARGET correction some guides are not corrected")
  }
 
  return(list(model_perf = models_performance,
              model_est = models_estimates,
              pseudo_single_segments = pseudo_single_segments,
              dual_FC = dual_FC_correctedFC))
  
}


#' Run complete CRISPRcleanR^2 analysis for dual screens.
#'
#' This function is the complete wrap-up for the entire pipeline, from raw count of both single and dual KO screens to CRISPRcleanR^2 correction, synergy estimation and results visualization.
#' In particular, this function
#' - runs CRISPRcleanR pipeline on single KO raw counts (genome-wide) (see \code{\link{ccr.run_complete}}),
#' - matches single KO and dual KO libraries based on sequence (see \code{\link{ccr2.matchDualandSingleSeq}}),
#' - pre-processes dual KO screen, converting raw counts into logFCs (if necessary, see \code{\link{ccr2.NormfoldChanges}}), averages logFCs across replicates and annotate sgRNA pairs with chromosomal positions (see \code{\link{ccr2.logFCs2chromPos}}),
#' - computes singletons (guide pairs matched with non-essential gene or intergenic region) summary statistics to be used in bliss model (see \code{\link{ccr2.get_summary_singletons}}),
#' - computes dual KO logFCs correction for "non-target" pairs (see \code{\link{ccr2.run_nontarget}}),
#' - computes dual KO logFCs correction for "double-cuts" dual pairs, targeting two existing genomic regions (see \code{\link{ccr2.run}}),
#' - combines the results from both non-target and double-cuts pairs,
#' - computes synergy based on the bliss model for both uncorrected and corrected dual KO logFCs (see \code{\link{ccr2.compute_bliss}}),
#' - visualize results such as relationship with CNA status (see \code{\link{ccr2.plot_correction}}, \code{\link{ccr2.plotClasses}}, \code{\link{ccr2.plotCNA}}, \code{\link{ccr2.plotCNAdensity}}),
#' - visualize results for single KO for guides available also in with dual KO library (see \code{\link{ccr2.plotMatchingSingle}}),
#' - centers logFCs of positive controls to -1 and negative controls to 0 for both uncorrected and corrected logFCs (see \code{\link{ccr2.scale_pos_neg}}),
#' - selects synergistic and lethal guide pairs based on logFCs threshold -1 and bliss z-score threshold 0.5 (see \code{\link{ccr2.plot_bliss_vs_FC}}), 
#' - compute logFCs at the gene pair level (median) and repeats the last 2 steps to obtain synergistic and lethal gene pairs.
#'
#' @param filename_single A string specifying the path of a tsv file containing the raw sgRNA counts. This must be a tab delimited file with one row per sgRNA and the following columns/headers: 
#'  - sgRNA: containing alphanumerical identifiers of the sgRNA under consideration;
#'  - gene: containing HGNC symbols of the genes targeted by the sgRNA under consideration;
#' followed by the columns containing the sgRNAs' counts for the controls and columns for library trasfected samples.
#' @param min_reads_single This parameter defines a filter threshold value for sgRNAs, based on their average counts in the control sample. 
#'  Specifically, it indicates the minimal number of counts that each individual sgRNA needs to have in the controls (on average) in order to be included in the output.
#' @param libraryAnnotation_single A data frame containing the sgRNA annotations, with a named row for each sgRNA, and columns for targeted genes, genomic coordinates and possibly other information. 
#' @param min_reads Minimum number of reads required in the control sample (e.g. plasmid) for a pair to be retained (default is 30).
#' @param libraryAnnotation_dual Data frame containing information from the dual library (one row per guide pair).
#' @param dual_count Data frame containing the raw counts from dual KO screen, raw counts start from the 9th column. The 9th column include the control raw counts (default is NULL).
#' @param dual_logFC Data frame containing dual KO log fold changes (default is NULL).
#' @param EXPname Name of the experiment.
#' @param display Logical, whether to display plots (default is FALSE).
#' @param outdir Directory to save the output files (default is "./").
#' @param saveToFig Logical, whether to save figures.
#' @param saveFormat File format for saving figures (default is "pdf").
#' @param correctGW A character specifying the type of pseudo-single logFCs mean centering after data-injection ("CHR": per chromosome, "GW": genome-wide, or NULL for no centering). Default is NULL. #TODO: USED FOR TESTING, REMOVE!
#' @param excludeGene Character vector, genes to be excluded from dual_FC_correctedFC (default is NULL).
#' @param CNA A data frame containing copy number alteration information. It must include the following columns:
#'  - CHROM: chromosome id in "chri" format   
#'  - start: start of the CN region, genomic location in base pair
#'  - end: end of the CN region, genomic location in base pair
#'  - C: copy number (absolute) value 
#' @param CN_thr Numeric (absolute), a threshold for defining extreme CNA.
#' @param ... Additional parameters to be passed to other functions.
#'
#' @return A list containing:
#'  - dual: data frame of dual KO logFCs including CRISPRcleanR^2 corrected values, bliss z-score synergy, scaled logFCs per positive and negative controls (one row per guide pair),
#'  - dual_gene: data frame same as dual but at the gene pair level (scores collapsed using median),
#'  - single_gw: data frame of single KO logFCs including CRISPRcleanR corrected values at the genome-wide level,
#'  - single: data frame of single KO logFCs including CRISPRcleanR corrected values only for those sgRNAs also present in dual KO library,
#'  - pseudo_single_segments: identified genomic regions (segments) by CRISPRcleanR applied to pseudo single logFCs,
#'  - system_solition: data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides (both positions), it also includes the fitted pseudo-single correction obtained as system solution,
#'  - model_perf:  model performance for single VS pseudo-single, including both double-cut pairs and non-target pairs,
#'  - model_est: model estimates (coefficients) for single VS pseudo-single, including both double-cut pairs and non-target pairs,
#'  - top_corrected: using CRISPRcleanR^2 corrected results, dual guide pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
#'  - top_gene_corrected: using CRISPRcleanR^2 corrected results, dual gene pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
#'  - top_uncorrected: using uncorrected (original) results, dual guide pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
#'  - top_gene_uncorrected: using uncorrected (original) results, dual gene pairs classified as lethal and synergistic based on scaled logFCS < -1 and bliss z-score < -0.5,
#'  - max_correction: maximum dual correction allowed based on single KO correction (used as upper bound for the optimization problem solution).,
#'  - min_correction: minimum dual correction allowed based on single KO correction (used as upper bound for the optimization problem solution).
#'
#' @examples
#' @seealso
#' \code{\link{ccr.run_complete}},
#' \code{\link{ccr2.matchDualandSingleSeq}},
#' \code{\link{ccr2.NormfoldChanges}},
#' \code{\link{ccr2.logFCs2chromPos}},
#' \code{\link{ccr2.get_summary_singletons}},
#' \code{\link{ccr2.match_singletons_singleFC}},
#' \code{\link{ccr2.run_nontarget}},
#' \code{\link{ccr2.run}},
#' \code{\link{ccr2.compute_bliss}},
#' \code{\link{ccr2.plot_correction}},
#' \code{\link{ccr2.plotClasses}},
#' \code{\link{ccr2.plotCNA}},
#' \code{\link{ccr2.plotCNAdensity}},
#' \code{\link{ccr2.plotMatchingSingle}},
#' \code{\link{ccr2.scale_pos_neg}},
#' \code{\link{ccr2.plot_bliss_vs_FC}}
#' 
#' @export
ccr2.run_complete <- function(
  filename_single, 
  min_reads_single = 30,
  libraryAnnotation_single,  
  min_reads = 30, 
  libraryAnnotation_dual, 
  dual_count = NULL,
  dual_logFC = NULL,
  EXPname,
  display = FALSE, 
  outdir = "./", 
  saveToFig, 
  saveFormat = "pdf", 
  correctGW, 
  excludeGene_plot = NULL, 
  CNA, 
  CN_thr = 8,
  ...
) {
  if ((is.null(dual_count) & is.null(dual_logFC)) | 
      (!is.null(dual_count) & !is.null(dual_logFC))) {
    stop("ONLY ONE between dual_count and dual_logFC MUST be not NULL")
  }
  
  #################################
  ### get single screens output ###
  
  single_screen <- ccr.run_complete(
    filename_single = filename_single, 
    min_reads_single = min_reads_single,
    EXPname = EXPname, 
    libraryAnnotation_single = libraryAnnotation_single, 
    display = display, 
    outdir = outdir)
  
  single_correctedFCs <- single_screen$FC
  libraryAnnotation_single <- single_screen$library
  
  print("Processed single screen")
  
  #########################################################
  ### match single and dual libraries based on sequence ###
  library_matched_seq <- ccr2.matchDualandSingleSeq(
    dual_library = libraryAnnotation_dual, 
    single_library = libraryAnnotation_single)
  
  print("Matched single and dual libraries")
  
  ###############################      
  ### get dual screen output ###
  if (!is.null(dual_count)) {
    print("Convert count to logFCs")
    dual_FC <- ccr2.NormfoldChanges(
      Dframe = dual_count, 
      min_reads = min_reads)  
    dual_logFC <- dual_FC$logFCs
  }
  
  dual_FC <- ccr2.logFCs2chromPos(
    dual_FC = dual_logFC, 
    dual_library = libraryAnnotation_dual)
  
  # get singletons averages:
  singletons_summary <- ccr2.get_summary_singletons(
    dual_FC = dual_FC)
  
  # plot comparison with single screen
  merged_singletons_single <- ccr2.match_singletons_singleFC(
    match_dual_single_seq = library_matched_seq, 
    single_FC = single_correctedFCs, 
    singletons_summary_FC = singletons_summary,
    saveToFig = saveToFig, 
    display = display, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
  #############################################################################
  ### get correction for pairs with NONTARGET, use pseudo single correction ###
  tmp <- ccr2.run_nontarget(
    dual_FC = dual_FC, 
    match_dual_single_seq = library_matched_seq, 
    single_correctedFCs = single_correctedFCs, 
    libraryAnnotation_dual = libraryAnnotation_dual, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname, 
    correctGW = correctGW
  )
  dual_FC_nt_correctedFC <- tmp$dual_FC
  model_perf_nt <-  tmp$model_perf %>% 
    mutate(type = "NONTARGET_PAIR")
  model_est_nt <-  tmp$model_est %>% 
    mutate(type = "NONTARGET_PAIR")
  
  pseudo_single_segments_nt <- tmp$pseudo_single_segments %>% 
    mutate(type = "NONTARGET_PAIR")
  
  ##########################################
  ### get correction for pairs targeting ###
  tmp <- ccr2.run(
    dual_FC = dual_FC, 
    match_dual_single_seq = library_matched_seq, 
    single_correctedFCs = single_correctedFCs, 
    libraryAnnotation_dual = libraryAnnotation_dual, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname, 
    correctGW = correctGW,
    ...
  )
  
  dual_FC_correctedFC <- tmp$dual_FC
  pseudo_single_correction <- tmp$pseudo_single
  max_correction <- tmp$max_correction
  min_correction <- tmp$min_correction
  
  model_perf <-  tmp$model_perf %>% 
    mutate(type = "DOUBLE_CUT_PAIR") %>% 
    bind_rows(model_perf_nt)
  model_est <-  tmp$model_est %>% 
    mutate(type = "DOUBLE_CUT_PAIR") %>% 
    bind_rows(model_est_nt)
  
  pseudo_single_segments <- tmp$pseudo_single_segments %>% 
    mutate(type = "DOUBLE_CUT_PAIR") %>% 
    bind_rows(pseudo_single_segments_nt)
  
  # combine all res 
  dual_FC_correctedFC <- rbind(dual_FC_correctedFC, dual_FC_nt_correctedFC)
  # add non-target / non-target
  dual_FC_missing <- dual_FC %>% 
    dplyr::filter(!ID %in% dual_FC_correctedFC$ID) %>%
    dplyr::mutate(sgRNA_ID_pair = paste0(sgRNA1_WGE_ID, "~", sgRNA2_WGE_ID),
                  correction = NA, correctedFC = NA)
  dual_FC_correctedFC <- rbind(dual_FC_correctedFC, dual_FC_missing)
  
  # add synergy from bliss model
  dual_FC_correctedFC <- ccr2.compute_bliss(
    dual_FC = dual_FC_correctedFC)
  
  dual_FC_correctedFC <- ccr2.compute_bliss(
    dual_FC = dual_FC_correctedFC,
    corrected = TRUE)
  print("Computed synergy via bliss model")
  
  ########################
  ### visualize output ###
  # plot bliss before and after
  ccr2.plot_correction(df = dual_FC_correctedFC, 
                       EXPname = EXPname, 
                       saveToFig = saveToFig, 
                       saveFormat = saveFormat,
                       outdir = outdir)
  # plot bliss fit
  #ccr2.plot_bliss_fit(dual_FC = dual_FC_correctedFC, 
  #                    EXPname = EXPname, 
  #                    saveToFig = saveToFig, 
  #                    saveFormat = saveFormat,
  #                    outdir = outdir)
  
  # divide by class
  ccr2.plotClasses(
    dual_FC_correctedFC = dual_FC_correctedFC,
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir)
  
  # distribution wrt CN
  ccr2.plotCNA(
    dual_FC_correctedFC = dual_FC_correctedFC, 
    CNA = CNA, 
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir)
  
  ccr2.plotCNA(
    dual_FC_correctedFC = dual_FC_correctedFC, 
    CNA = CNA, 
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    var_to_plot = "bliss_zscore")
  
  # density distribution for high CN for any guide1 or guide2
  ccr2.plotCNAdensity(
    dual_FC_correctedFC = dual_FC_correctedFC, 
    CNA = CNA, 
    CN_thr = CN_thr,  
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir)
  
  ccr2.plotCNAdensity(
    dual_FC_correctedFC = dual_FC_correctedFC, 
    CNA = CNA, 
    CN_thr = CN_thr,  
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    var_to_plot = "bliss_zscore")
  
  
  if (!is.null(excludeGene_plot)) {
    
    # exclude pairs including a gene
    ccr2.plotCNA(dual_FC_correctedFC = dual_FC_correctedFC, 
                 CNA = CNA, 
                 excludeGene = excludeGene_plot,  
                 EXPname = EXPname, 
                 saveToFig = saveToFig,
                 saveFormat = saveFormat,
                 outdir = outdir)
    
    ccr2.plotCNA(dual_FC_correctedFC = dual_FC_correctedFC, 
                 CNA = CNA, 
                 excludeGene = excludeGene_plot,  
                 EXPname = EXPname, 
                 saveToFig = saveToFig, 
                 saveFormat = saveFormat,
                 outdir = outdir, 
                 var_to_plot = "bliss_zscore")
    
    ccr2.plotCNAdensity(dual_FC_correctedFC = dual_FC_correctedFC, 
                        CNA = CNA, 
                        excludeGene = excludeGene_plot, 
                        CN_thr = CN_thr, 
                        EXPname = EXPname, 
                        saveToFig = saveToFig, 
                        saveFormat = saveFormat,
                        outdir = outdir)
    
    ccr2.plotCNAdensity(dual_FC_correctedFC = dual_FC_correctedFC, 
                        CNA = CNA, 
                        excludeGene = excludeGene_plot, 
                        CN_thr = CN_thr, 
                        EXPname = EXPname, 
                        saveToFig = saveToFig, 
                        saveFormat = saveFormat,
                        outdir = outdir, 
                        var_to_plot = "bliss_zscore")
  }
  
  
  # compare with correction on single genes
  single_correctedFCs_filt <- ccr2.plotMatchingSingle(
    dual_FC_correctedFC = dual_FC_correctedFC, 
    match_dual_single_seq = library_matched_seq, 
    single_correctedFCs = single_correctedFCs, 
    CNA = CNA,  
    EXPname = EXPname, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir)
  
  # center on positive and negative classes
  dual_FC_correctedFC <- ccr2.scale_pos_neg(dual_FC = dual_FC_correctedFC)
  dual_FC_correctedFC <- ccr2.scale_pos_neg(dual_FC = dual_FC_correctedFC, 
                                            corrected = TRUE)
  
  # plot bliss vs avgFC
  top_corrected <- ccr2.plot_bliss_vs_FC(dual_FC = dual_FC_correctedFC, 
                                         corrected = TRUE, 
                                         THR_FC = -1, 
                                         THR_BLISS = -0.5,
                                         EXPname = EXPname, 
                                         saveToFig = saveToFig, 
                                         saveFormat = saveFormat,
                                         outdir = outdir)
  
  top_uncorrected <- ccr2.plot_bliss_vs_FC(dual_FC = dual_FC_correctedFC, 
                                           corrected = FALSE, 
                                           THR_FC = -1, 
                                           THR_BLISS = -0.5,
                                           EXPname = EXPname, 
                                           saveToFig = saveToFig, 
                                           saveFormat = saveFormat,
                                           outdir = outdir)
  
  
  # get gene level
  dual_FC_gene_correctedFC <- dual_FC_correctedFC %>%
    dplyr::group_by(Gene_Pair) %>%
    dplyr::summarise(
      n_guides = dplyr::n(),
      ID = paste0(ID, collapse = ","), 
      ID_lib = paste0(ID_lib, collapse = ","), 
      lib = paste0(lib, collapse = ","), 
      info_subtype = paste0(unique(info_subtype), collapse = ","), 
      info = paste0(unique(info), collapse = ","), 
      Gene1 = unique(Gene1), 
      Gene1_Chr = paste0(unique(sgRNA1_Chr), collapse = ","), 
      Gene1_BP = mean(sgRNA1_BP), 
      Gene2 = unique(Gene2), 
      Gene2_BP = mean(sgRNA2_BP), 
      Gene2_Chr = paste0(unique(sgRNA2_Chr), collapse = ","), 
      correction = median(correction),
      avgFC = median(avgFC),
      bliss_zscore = median(bliss_zscore), 
      correctedFC = median(correctedFC),
      bliss_zscore_corrected = median(bliss_zscore_corrected)) %>%
    dplyr::ungroup()
  
  # center on positive and negative classes
  dual_FC_gene_correctedFC <- ccr2.scale_pos_neg(dual_FC = dual_FC_gene_correctedFC)
  dual_FC_gene_correctedFC <- ccr2.scale_pos_neg(dual_FC = dual_FC_gene_correctedFC, 
                                                 corrected = TRUE)
  # plot bliss vs avgFC
  top_gene_corrected <- ccr2.plot_bliss_vs_FC(dual_FC = dual_FC_gene_correctedFC,
                                             corrected = TRUE, 
                                             THR_FC = -1, 
                                             THR_BLISS = -0.5,
                                             EXPname = EXPname, 
                                             saveToFig = saveToFig, 
                                             saveFormat = saveFormat,
                                             outdir = sprintf("%sGENELEVEL_", outdir))
  
  top_gene_uncorrected <- ccr2.plot_bliss_vs_FC(dual_FC = dual_FC_gene_correctedFC, 
                                              THR_FC = -1, 
                                              THR_BLISS = -0.5,
                                              corrected = FALSE, 
                                              EXPname = EXPname, 
                                              saveToFig = saveToFig, 
                                              saveFormat = saveFormat,
                                              outdir = sprintf("%sGENELEVEL_", outdir))
  
  return(list(dual = dual_FC_correctedFC, 
              dual_gene = dual_FC_gene_correctedFC,
              single_gw = single_correctedFCs,  
              single = single_correctedFCs_filt,
              pseudo_single_segments = pseudo_single_segments,
              system_solition = pseudo_single_correction,
              model_perf = model_perf,
              model_est = model_est,
              top_corrected = top_corrected, 
              top_gene_corrected = top_gene_corrected, 
              top_uncorrected = top_uncorrected, 
              top_gene_uncorrected = top_gene_uncorrected, 
              max_correction = max_correction, 
              min_correction = min_correction))
  
}

#' Get Pairwise Correction approximating linear system based on pseudo-single construction.
#'
#' This function calculates pairwise correction using a system solver.
#' The correction is found as the optimal solution of a convex problem (approximation of the linear system) with constraints on the solution (maximum correction and minimum correction allowed).
#' The CVXR package is used with default solver OSQP.
#' The function creates a sparse matrix representing the system of equations,
#' solves the problem approximating linear system solution with constraints, and returns the pairwise correction. 
#' NOTE: the box constraints prevent the solution from "exploding" introducing highly negative correction and hence false positives (avoid induced highly lethal pairs). 
#'
#' @param matrix_interactions matrix with logical entries, rows = sgRNAs in position 1 and cols =  sgRNAs in position 2. Each entry indicates if the pair sgRNA pos1 ~ sgRNA pos 2 exists in the given dual KO dataset.
#' @param combinations matrix with same size as matrix_interactions. Each entry is a character formed by sgRNA pos 1 ~ sgRNA pos 2 expliciting the guide pair combination (it could not exist in dual KO screen).
#' @param dual_FC Data frame with dual KO logFCs (already averaged per replications, logFCs stored in "avgFC" column).
#' @param pseudo_single_p1 Data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides in position 1 (output of \code{\link{ccr2.filterGWclean}}).
#' @param pseudo_single_p2 Data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides in position 2 (output of \code{\link{ccr2.filterGWclean}}).
#' @param max_corr maximum dual correction allowed based on single KO correction (used as upper bound for the solution).
#' @param min_corr minimum dual correction allowed based on single KO correction (used as upper bound for the solution).
#'
#' @return A list containing 
#'  - sys_solution: system solution by the chosen solver (OSQP solver from CVXR),
#'  - df_corr: data frame of sgRNA ID pair and pairwise correction, 
#'  - pseudo_single_correction: Data frame containing CRISPRcleanR corrected logFCs for pseudo-single guides (both positions), it also includes the fitted pseudo-single correction obtained as system solution.
#' 
#' @examples
#' @export
#' @seealso 
#' \code{\link{ccr2.solveLinearSys}}
get_pairwise_correction <- function(matrix_interactions, 
                                    combinations, 
                                    dual_FC, 
                                    pseudo_single_p1, 
                                    pseudo_single_p2, 
                                    max_corr, 
                                    min_corr){
  
  # 2. create matrix of system
  # n.eq=nrow(int)+ncol(int) times n.unknowns=nrow(int)*ncol(int)
  sparse_mat_j <- c()
  sparse_mat_i <- c()
  
  for (idx_row in 1:nrow(matrix_interactions)) {
    tmp <- matrix_interactions
    tmp[-idx_row, ] <- 0
    tmp_vec <- as.vector(t(tmp))
    tmp_j <- which(tmp_vec != 0)
    tmp_i <- rep(idx_row, length(tmp_j))
    sparse_mat_j <- c(sparse_mat_j, tmp_j)
    sparse_mat_i <- c(sparse_mat_i, tmp_i)
  }
  
  for (idx_col in 1:ncol(matrix_interactions)) {
    tmp <- matrix_interactions
    tmp[, -idx_col] <- 0
    tmp_vec <- as.vector(t(tmp))
    tmp_j <- which(tmp_vec != 0)
    tmp_i <- rep(idx_col + nrow(matrix_interactions), length(tmp_j))
    sparse_mat_j <- c(sparse_mat_j, tmp_j)
    sparse_mat_i <- c(sparse_mat_i, tmp_i)
  }
  
  nrow_sys <- as.integer(nrow(matrix_interactions) + ncol(matrix_interactions))
  ncol_sys <- as.integer(nrow(matrix_interactions) * ncol(matrix_interactions))
  matrix_sys <- sparseMatrix(
    i = sparse_mat_i, 
    j = sparse_mat_j,
    x = 1, 
    dims = c(nrow_sys, ncol_sys), 
    dimnames = list(c(rownames(matrix_interactions), colnames(matrix_interactions)), 
                    as.vector(t(combinations)))) 
  
  # matrix_sys include all possible combinations, 
  # remove those that are not present in dual_FC
  matrix_sys <- matrix_sys[, colnames(matrix_sys) %in% dual_FC$sgRNA_ID_pair]
  print(paste("dimension system matrix:", dim(matrix_sys)[1], dim(matrix_sys)[2]))
  print(paste("n. entries not zero:", length(matrix_sys@i)))
  
  # 3. create coefficient vector
  correction_vect <- c(pseudo_single_p1$correction_scaled, 
                       pseudo_single_p2$correction_scaled)
  
  
  # 4. solve system, (constrained)
  x <- CVXR::Variable(ncol(matrix_sys))
  objective <- CVXR::Minimize(sum((matrix_sys %*% x - correction_vect)^2))
  problem <- CVXR::Problem(objective, constraints = list(x >= min_corr, x <= max_corr))
  sys_solution <- CVXR::solve(problem)
  correction_pair <- sys_solution$getValue(x)[,1]
  print(paste("Solution status is", sys_solution$status))
  print(paste("The solver used is", sys_solution$solver))
  print(paste("Problem solved in", sys_solution$solve_time, "sec"))
  
  # TODO: removed tested methods
  # # try with glmnet
  # time_solve <- system.time(sys_solution <- glmnet(
  #   x = matrix_sys, 
  #   y = correction_vect, 
  #   # lower.limits = 0, 
  #   # upper.limits = 2,
  #   alpha = 1, lambda = 0, # no penalization
  #   intercept = FALSE, 
  #   standardize = FALSE)
  # )
  # print(sprintf("system solved after: %.2fs", time_solve[3]))
  # correction_pair <- as.vector(sys_solution$beta[,1])
  # # try vith bvls
  # time_solve <- system.time(sys_solution <- bvls::bvls(
  #            as.matrix(matrix_sys), 
  #            correction_vect, 
  #            bl = rep(0, ncol(matrix_sys)), # how to set up bounds?
  #            bu = rep(1, ncol(matrix_sys)))
  # )
  #print(sprintf("system solved after: %.2fs", time_solve[3]))
  #correction_pair <- sys_solution$x
  # # try with nnls
  # a <- Matrix::crossprod(matrix_sys)
  # b <- Matrix::crossprod(matrix_sys, correction_vect)
  # print("matrix computed")
  # time_solve <- system.time(sys_solution <- RcppML::nnls(
  #   fast_nnls = TRUE,
  #   a = as.matrix(a), 
  #   b = as.matrix(b))
  # )
  # print(sprintf("system solved after: %.2fs", time_solve[3]))
  # correction_pair <- sys_solution
  # try with lsei
  # time_solve <- system.time(sys_solution <- limSolve::lsei(
  #   A = matrix_sys, 
  #   B = correction_vect, 
  #   G = diag(x = 1, nrow = ncol(matrix_sys), ncol = ncol(matrix_sys)), 
  #   H = rep(0, ncol(matrix_sys)))
  # )
  # print(sprintf("system solved after: %.2fs", time_solve[3]))
  # correction_pair <- sys_solution$X
  # # try with nnls
  # time_solve <- system.time(sys_solution <- nnls::nnls(
  #   matrix_sys, correction_vect
  #   )
  # )
  # print(sprintf("system solved after: %.2fs", time_solve[3]))
  # correction_pair <- sys_solution$x
  # use moore-penrose approximation
  # time_solve <- system.time(correction_pair <- solve_system_fun( # function removed, if interested check v.0.5.1)
  #   matrix_sys = matrix_sys,
  #   b_coeff = correction_vect)
  # )
  # sys_solution <- NULL
  
  # save output
  df_corr <- data.frame(sgRNA_ID_pair = colnames(matrix_sys), 
                        correction = correction_pair) 
  
  pseudo_single_correction <- data.frame(
    sgRNA_ID = c(pseudo_single_p1$sgRNA_ID, 
                 pseudo_single_p2$sgRNA_ID), 
    gene = c(pseudo_single_p1$genes, 
             pseudo_single_p2$genes),
    pseudo_single_fitted = as.vector(matrix_sys %*% correction_pair), 
    pseudo_single = correction_vect, 
    guide_pos = c(rep(1, nrow(pseudo_single_p1)), 
                  rep(2, nrow(pseudo_single_p2))), 
    n = c(pseudo_single_p1$n, pseudo_single_p2$n)
  )
  
  # for plot
  pseudo_single_correction$pseudo_single_fitted <- pseudo_single_correction$pseudo_single_fitted/c(pseudo_single_p1$n,  
                                                                                                   pseudo_single_p2$n)
  pseudo_single_correction$pseudo_single <- pseudo_single_correction$pseudo_single/c(pseudo_single_p1$n,  
                                                                                     pseudo_single_p2$n)
  pseudo_single_correction$guide_pos <- factor(pseudo_single_correction$guide_pos)
  pseudo_single_correction$ID <- paste0(pseudo_single_correction$gene, "_", 
                                        pseudo_single_correction$sgRNA_ID)
  pseudo_single_correction$outliers <- ""
  id_out <- abs(pseudo_single_correction$pseudo_single - pseudo_single_correction$pseudo_single_fitted) > 0.5
  pseudo_single_correction$outliers[id_out] <- pseudo_single_correction$ID[id_out]
  
  return(
    list(
      sys_solution = sys_solution,
      df_corr = df_corr, 
      pseudo_single_correction = pseudo_single_correction
    )
  )
  
}

