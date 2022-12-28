#' @import ggplot2
NULL

#' Load matched files for dual KO
#'
#' get_input_data() load files as indicated in param_file. Must include count,
#' library and copy_number fields.
#'
#' @param param_file .tsv file, the first column indicates the type of data, the
#'   second the file location. Can include multiple files that contain the name
#'   count or library BUT the lenght must coincide.
#'
#' @return a list * CNA: copy number, one row per gene * count:  dual KO count,
#'   one row per guides combination * library: metadata for dual KO library *
#'   CL_name: cell line name * out_fold: location to store results
#' @export
#' 
get_input_data <- function(param_file) {
  
  input_info <- suppressWarnings(readr::read_table(
    param_file, 
    col_names = FALSE, 
    show_col_types = FALSE))
  
  for (i in 1:nrow(input_info)) {
    assign(input_info$X1[i], input_info$X2[i])    
  }
  
  count_files <- input_info$X1[grepl("count", input_info$X1)]
  library_files <- input_info$X1[grepl("library", input_info$X1)]
  if (length(count_files) != length(library_files)) {
    stop("library and count files MUST match")
  }
  
  n_files <- length(count_files)
  dual_count_list <- lapply(count_files, function(x)
    suppressWarnings(readr::read_table(sprintf('%s%s', input_fold, get(x)), 
                                show_col_types = FALSE)))
  
  dual_library_list <- lapply(library_files , function(x)
    suppressWarnings(readxl::read_xlsx(sprintf('%s%s', input_fold, get(x)))))
  
  for (idx_file in seq_len(n_files)) {
    
    dual_count_list[[idx_file]] <- dual_count_list[[idx_file]] %>%
      dplyr::mutate(ID = paste0(ID, "_lib", idx_file))
    
    dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>%
      dplyr::filter(!is.na(ID)) %>%
      dplyr::mutate(ID = paste0(ID, "_lib", idx_file))
    
    if (!identical(dual_count_list[[idx_file]]$ID, dual_library_list[[idx_file]]$ID)) {
      stop("count and library rows must match")
    }
    
    id_change <- grep("CTRL", dual_library_list[[idx_file]]$sgRNA1_WGE_ID)
    if (length(id_change) > 0) {
      dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_change] <- str_replace(
        string =  dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_change], 
        pattern = "CTRL", 
        replacement = "NONTARGET")
    }
    
    id_change <- grep("CTRL", dual_library_list[[idx_file]]$sgRNA2_WGE_ID)
    if (length(id_change) > 0) {
      dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_change] <- str_replace(
        string =  dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_change], 
        pattern = "CTRL", 
        replacement = "NONTARGET")
    }
    
    id_nontarget <- grep("NONTARGET", dual_count_list[[idx_file]]$Gene1)
    if (length(id_nontarget) > 0) {
      dual_count_list[[idx_file]]$sgRNA1_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene1[id_nontarget]
      dual_library_list[[idx_file]]$sgRNA1_WGE_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene1[id_nontarget]
    }
    
    id_nontarget <- grep("NONTARGET", dual_count_list[[idx_file]]$Gene2)
    if (length(id_nontarget) > 0) {
      dual_count_list[[idx_file]]$sgRNA2_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene2[id_nontarget]
      dual_library_list[[idx_file]]$sgRNA2_WGE_ID[id_nontarget] <- dual_count_list[[idx_file]]$Gene2[id_nontarget]
    }
    
    dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>% 
      dplyr::mutate(COMB_ID =  paste(sgRNA1_WGE_Sequence, sgRNA2_WGE_Sequence, sep = "_"))
    
    if ("MyNotes" %in% colnames(dual_library_list[[idx_file]])) {
      dual_library_list[[idx_file]] <- dual_library_list[[idx_file]] %>% 
        dplyr::select(-MyNotes)
    }
    
  }
  
  dual_CNA <- readr::read_table(sprintf('%s%s', input_fold,  copy_number_file), 
                         show_col_types = FALSE)
  
  return(list(CNA = dual_CNA, 
              count = dual_count_list, 
              library =  dual_library_list, 
              CL_name = CL_name, 
              out_fold = out_fold))
  
}

#' Title
#'
#' @param filename_single 
#' @param min_reads_single 
#' @param EXPname 
#' @param libraryAnnotation_single 
#' @param display 
#' @param outdir 
#'
#' @return
#' @export
#'
#' @examples
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
                                                           rownames(libraryAnnotation_single))])
  # get CRISPRCleanR single correction
  single_correctedFCs <- ccr.GWclean(gwSortedFCs,
                                     saveTO = outdir, 
                                     display = display,
                                     label = sprintf("single_%s", EXPname)) 
  
  single_correctedFCs <- single_correctedFCs$corrected_logFCs %>% 
    dplyr::mutate(correction = correctedFC - avgFC)
  
  libraryAnnotation_single <- libraryAnnotation_single %>% 
    dplyr::filter(CODE %in% rownames(single_correctedFCs)) 
  
  return(list(
    FC = single_correctedFCs, 
    library = libraryAnnotation_single
    ))
  
}

#' Title
#'
#' @param dual_library 
#' @param single_library 
#'
#' @return
#' @export
#'
#' @examples
ccr2.matchDualandSingleSeq <- function(dual_library, single_library) {
  
  dual_library_seq <- data.frame(
    GENES = c(dual_library$sgRNA1_Approved_Symbol, dual_library$sgRNA2_Approved_Symbol), 
    ID = c(dual_library$sgRNA1_WGE_ID, dual_library$sgRNA2_WGE_ID), 
    LIBRARY = c(dual_library$sgRNA1_Library, dual_library$sgRNA2_Library), 
    SEQ = c(dual_library$sgRNA1_WGE_Sequence, dual_library$sgRNA2_WGE_Sequence)
    )
  
  dual_library_seq <- dual_library_seq %>% 
    dplyr::filter(!duplicated(SEQ)) %>%
    dplyr::mutate(ID_single = NA, SEQ_single = NA)
  # filter(!duplicated(ID), LIBRARY == single_library_name) # from input consider all the libraries
  
  single_library_seq <- data.frame(
    GENES = single_library$GENES, 
    ID = rownames(single_library),
    SEQ = single_library$seq
    )
  
  # filter for genes and non-targeting
  single_library_seq <- single_library_seq %>%
    dplyr::filter(!duplicated(SEQ), 
                  GENES %in% dual_library_seq$GENES | GENES == "NON-TARGETING")
  
  #dual_library_seq <- dual_library_seq %>%
  #    filter(GENES %in% single_library_seq$GENES | is.na(GENES))
  
  # NOTE: RACK1, MARCHF5 and INTS6L removed because they have another name in single, how to solve?
  
  single_nchar <- unique(nchar(single_library_seq$SEQ))
  dual_nchar <- unique(nchar(dual_library_seq$SEQ))
  
  if (length(single_nchar) > 1) {
    stop("single screen library MUST have the same length for each guide")
  }
  
  if (all(single_nchar <= dual_nchar)) {
    
    match_id_dual  <-  lapply(single_library_seq$SEQ, function(x) 
      grep(pattern = x, x = dual_library_seq$SEQ))
    len_match <- sapply(match_id_dual, length)
    to_keep <- which(len_match == 1) # keep onle guides with a unique association
    match_id_dual <- unlist(match_id_dual[to_keep])
    
    dual_library_seq$ID_single[match_id_dual] <- single_library_seq$ID[to_keep]
    dual_library_seq$SEQ_single[match_id_dual] <- single_library_seq$SEQ[to_keep]
    
  } else {
    
    if (all(dual_nchar <= single_nchar)) {
      
      match_id_single  <-  lapply(dual_library_seq$SEQ, function(x) 
        grep(pattern = x, x = single_library_seq$SEQ))
      len_match <- sapply(match_id_single, length)
      to_keep <- which(len_match == 1) # keep only guides with a unique association
      match_id_single <- unlist(match_id_single[to_keep])
      
      dual_library_seq$ID_single[to_keep] <- single_library_seq$ID[match_id_single]
      dual_library_seq$SEQ_single[to_keep] <- single_library_seq$SEQ[match_id_single]
      
    } else {
      
      stop("dual screen guides have variable length, > and < than single screen")
      
    }
    
  }
  
  return(dual_library_seq)
}


#' Title
#'
#' @param Dframe 
#' @param min_reads 
#'
#' @return
#' @export
#'
#' @examples
ccr2.NormfoldChanges <- function(
  Dframe, 
  min_reads = 30
) {
  
  # remove column already normalized and logFC
  counts <- Dframe %>% 
    dplyr::select(-dplyr::ends_with("_norm"), -dplyr::ends_with("_logFC"))
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
  
  norm_fact <- t(matrix(rep(colSums(counts), nrow(counts)), 
                        ncol(counts), nrow(counts)))
  norm_counts <- counts/norm_fact*10000000
  
  norm_controls <- norm_counts[, 1]
  norm_counts <- norm_counts[, -1]
  
  log_fc <- apply(norm_counts, 2, function(x) log2((x + 0.5)/(norm_controls + 0.5)))
  colnames(log_fc) <- paste(colnames(log_fc), "logFC", sep = "_")
  log_fc <- cbind(info, log_fc)
  
  return(log_fc)
}

#' Title
#'
#' @param dual_FC 
#' @param dual_library 
#'
#' @return
#' @export
#'
#' @examples
ccr2.logFCs2chromPos <- function(
  dual_FC, 
  dual_library
) {
  
  dual_library <-  dual_library %>% 
    dplyr::select(ID, sgRNA1_Library, 
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
    dplyr::select(ID, MyNote, Note, Gene_Pair, Gene1, Gene2, 
                  # dplyr::ends_with("_Scaled_logFC")) %>%
                  dplyr::ends_with("_logFC")) %>%
    dplyr::rename(info = Note, info_subtype = MyNote) %>% 
    dplyr::mutate(avgFC = rowMeans(dplyr::select(.,ends_with("_logFC"))))
  #dplyr::mutate(avgFC = rowMeans(dplyr::select(., ends_with("_Scaled_logFC"))))
  avg_FC <-  avg_FC %>% dplyr::select(-dplyr::ends_with("_logFC"))
  
  # merge with position info
  combined <- dplyr::left_join(avg_FC, dual_library, by = "ID") %>%
    dplyr::mutate(sgRNA1_BP = sgRNA1_Start + (sgRNA1_End - sgRNA1_Start)/2, .after = sgRNA1_End) %>%
    dplyr::mutate(sgRNA2_BP = sgRNA2_Start + (sgRNA2_End - sgRNA2_Start)/2, .after = sgRNA2_End)
  
  return(combined)
  
}
