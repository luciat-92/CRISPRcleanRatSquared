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
  dual_library_list <- list()
  dual_count_list <- lapply(count_files, function(x)
    suppressWarnings(readr::read_table(sprintf('%s%s', input_fold, get(x)), 
                                show_col_types = FALSE)))
  
  for (idx_file in seq_len(n_files)) {
    
    # load library and specify col types:
    names_col_lib <- names(readxl::read_xlsx(
      path = sprintf('%s%s', input_fold, get(library_files[idx_file])), 
      n_max = 0))
    col_types_lib <- ifelse(grepl("Chr", names_col_lib) | grepl("WGE_ID", names_col_lib), 
                            "text", "guess")
    dual_library_list[[idx_file]] <- readxl::read_xlsx(
      path = sprintf('%s%s', input_fold, get(library_files[idx_file])),
      col_types = col_types_lib)
    
    print(sprintf("############ load file n. %i ############", idx_file))
    
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
          .default = as.character(sgRNA1_Chr)), 
        sgRNA2_Chr = dplyr::case_when(
          sgRNA2_Chr == "X" ~ "23",
          sgRNA2_Chr == "Y" ~ "24",
          .default = as.character(sgRNA2_Chr))
        ) %>%
      dplyr::mutate(
        sgRNA1_Chr = as.numeric(sgRNA1_Chr), 
        sgRNA2_Chr = as.numeric(sgRNA2_Chr)
        )
  }
  
  CNA <- readr::read_table(sprintf('%s%s', input_fold,  copy_number_file), 
                         show_col_types = FALSE) %>%
    dplyr::mutate(
      CHROM = dplyr::case_when(
      CHROM == "chrX" ~ "chr23",
      CHROM == "chrY" ~ "chr24",
      .default = as.character(CHROM))
    )
  
  if ("Sampleid" %in% colnames(CNA)) {
    CNA <- CNA %>% 
      dplyr::filter(Sampleid %in% CL_name)
  }
    
  return(list(CNA = CNA, 
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
        .default = as.character(CHRM))
      ) %>%
    dplyr::mutate(CHRM = as.numeric(CHRM))
  
  return(list(
    FC = single_correctedFCs, 
    segment = single_ccr_segments, 
    library = libraryAnnotation_single
    ))
  
}


#' Title
#'
#' @param dual_FC 
#' @param corrected 
#'
#' @return
#' @export
#'
#' @examples
ccr2.get_summary_singletons <- function(dual_FC, 
                                       # singletons_n_guides, 
                                       corrected = F){
  
  name_var <- "avgFC"
  if (corrected) {
    name_var <- "correctedFC"
  }
  
  singletons_pos1 <- dual_FC %>%
    dplyr::filter(grepl("Singletons", info_subtype)) %>%
    dplyr::group_by(sgRNA1_WGE_ID) %>%
    dplyr::summarise(n_ID = dplyr::n()) %>%
    dplyr::arrange(desc(n_ID)) %>%
    dplyr::filter(n_ID > 30) %>% # how to set this??
    dplyr::pull(sgRNA1_WGE_ID) %>%
    sort()
  
  singletons_pos2 <- dual_FC %>%
    dplyr::filter(grepl("Singletons", info_subtype)) %>%
    dplyr::group_by(sgRNA2_WGE_ID) %>%
    dplyr::summarise(n_ID = dplyr::n()) %>%
    dplyr::arrange(desc(n_ID)) %>%
    dplyr::filter(n_ID > 30) %>% # how to set this??
    #dplyr::slice_max(n_ID) %>% ## only working if the number is always the same! adjust
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

#' Title
#'
#' @param dual_FC 
#' @param corrected 
#'
#' @return
#' @export
#'
#' @examples
ccr2.scale_pos_neg <- function(dual_FC, corrected = FALSE){
  
  dual_FC_original <- dual_FC
  dual_FC <- dual_FC[!is.na(dual_FC$correction),]
  
  name_var <- "avgFC"
  if (corrected) {
    name_var <- "correctedFC"
  }
  
  median_neg <-  median(dual_FC[dual_FC$info == "NegativeControls", name_var], na.rm = TRUE)
  median_pos <-  median(dual_FC[dual_FC$info == "PositiveControls", name_var], na.rm = TRUE)
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


#' Title
#'
#' @param match_dual_single_seq 
#' @param single_FC 
#' @param singletons_summary_FC 
#' @param saveToFig 
#' @param display 
#' @param saveFormat 
#' @param EXPname 
#' @param outdir 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param dual_FC 
#' @param corrected 
#'
#' @return
#' @export
#'
#' @examples
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


#' Title
#'
#' @param dual_seq 
#' @param single_library_seq 
#'
#' @return
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
  
  # NOTE: RACK1, MARCHF5 and INTS6L removed because they have another name in single,
  # how to solve? (partially solved with hg38 matching)
  
  # Create matrix to match single x dual SEQ
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

#' Title
#'
#' psuedo single avgFC via mean for each gene in guide1 or guide2 across all the other genes
#' store also number of genes from which mean is computed
#'
#' @param dual_FC 
#' @param guide_id 
#' @param single_FC 
#' @param match_dual_single_seq 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param dual_FC 
#' @param saveToFig 
#' @param display 
#' @param saveFormat 
#' @param EXPname 
#' @param outdir 
#' @param single_FC 
#' @param match_dual_single_seq 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' Model to convert combined dual to single screens
#' that depends on the specific gene selection
#'
#' @param guide_id 
#' @param correctGW 
#' @param display 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#' @param pseudo_single_FC 
#'
#' @return
#' @export
#'
#' @examples
ccr2.modelSingleVSPseudoSingle <- function(
  pseudo_single_FC, 
  guide_id, 
  correctGW = NULL, 
  display=TRUE, 
  saveToFig=FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = ""
) { 
  
  # remove NA in single screens
  matched_df <- pseudo_single_FC[!is.na(pseudo_single_FC$ID_single), ]
  # fmla <- as.formula("avgFC ~ avgFC_single*correction_single")
  fmla <- "avgFC ~ 0 + avgFC_single"

  fit_model <- glm(formula = as.formula(fmla), 
                   data = matched_df,
                   weights = matched_df$n, 
                   family = gaussian(link = "identity"))
  # remove outliers
  cooksD <- cooks.distance(fit_model)
  # remove extreme cases (quatile > 0.99)
  thr_influential <- quantile(cooksD, prob = 0.99)
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
                     weights = matched_df_filt$n, 
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
                        n_guides_in_mean = as.numeric(fit_model$weights))
                        # lib = dual_collapsed_filt$lib)
  
  if (saveToFig) {
    display <- TRUE
    file_name_comp <- sprintf("%s%s_PseudoGuide%i_vs_single.%s", 
                              outdir, EXPname, guide_id, saveFormat)
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
      label = sprintf("Pear. corr. = %.3f", cor(df_plot$avgFC_single, df_plot$avgFC)), 
      xpos = -Inf, 
      ypos = Inf)
    
    pl_comp <- ggplot(df_plot, aes(x = avgFC_single, y = avgFC)) + 
      geom_line(aes(y =  fitted_avgFC), colour = "red", linewidth = 1) + 
      geom_abline(intercept = 0, slope = 1, 
                  colour = "black", size = 1, linetype = "dashed") + 
      geom_point(aes(y = avgFC, 
                     # size = n_guides_in_mean, 
                     color = n_guides_in_mean), 
                 alpha = 0.5) +
      #geom_smooth(formula = y ~ 0 + x, method = "loess", 
      #            span = 0.3, 
      #            method.args = list(degree = 2)) + 
      theme_bw() + 
      theme(axis.title = element_text(size = 12), 
            axis.text = element_text(size = 11), 
            legend.position = "bottom", 
            plot.title = element_text(hjust = 0.5)) + 
      scale_colour_viridis_c(breaks = breaks_size) +
      # scale_size_continuous(breaks = breaks_size) + 
      geom_text(data = text_corr, 
                aes(label = label, x = xpos, y = ypos), size = 5, 
                hjust = -0.1, vjust = 1.1, inherit.aes = F) +
      ylab("Pseudo Single avgFC") + 
      xlab("Single avgFC") + 
      ggtitle(sprintf("Guide position %i", guide_id))
    print(pl_comp)
    
    # residuals VS fitted 
    df_plot$outliers <- ""
    df_plot$outliers[abs(df_plot$residuals_avgFC) > 10] <- df_plot$id[abs(df_plot$residuals_avgFC) > 10]
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
    ggsave(filename = file_name_resid, plot = pl_pred_vs_res, width = 6, height = 5)
    ggsave(filename = file_name_qqplot, plot = pl_qq, width = 5, height = 5)
  }
  
  return(list(model = fit_model,
              model_summary = summary(fit_model), 
              rm_guides_cooks = rm_guides_cooks))
}

#' Title
#'
#' @param model_single_to_pseudo 
#' @param single_FC 
#' @param guide_id 
#' @param correctGW 
#' @param pseudo_single_FC 
#'
#' @return
#' @export
#'
#' @examples
ccr2.injectData <- function(
  model_single_to_pseudo,
  pseudo_single_FC,
  single_FC, 
  guide_id, 
  correctGW = NULL
) { 
  
  ## predict on single guides (make a new function) ##
  single_FC <- single_FC %>% 
    dplyr::rename(avgFC_single = avgFC, correction_single = correction)
  
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

#' Title
#'
#' adjust ccr output applied to pseudo single and filter per considered guide in dual
#'
#' @param dataInjection_correctedFCs 
#' @param guide_id 
#' @param saveToFig 
#' @param display 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#' @param pseudo_single 
#'
#' @return
#' @export
#'
#' @examples
ccr2.filterGWclean <- function(
  dataInjection_correctedFCs, 
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
  
  if (saveToFig) {
    display <- TRUE
    file_name <- sprintf("%s%s_correction_singleVSpseudosingle_guide%i.%s", 
                         outdir, EXPname, guide_id, saveFormat)
  }
  
  if (display) {
    
    pl <- ggplot(complete_out, 
                 aes(x = correction_single, y = correction)) + 
      geom_point(alpha = 0.7, size = 1.5) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "red") + 
      geom_hline(yintercept = 0, linetype = "dashed", color = "red") + 
      theme_bw() + 
      theme(legend.position = "right") + 
      ggtitle(sprintf("Guide position %i", guide_id)) +
      xlab("Single correction") + 
      ylab("Pseudo single correction")
    print(pl)
    
  }
  
  if (saveToFig) {
    ggsave(filename = file_name, plot = pl, width = 4.5, height = 4.5)
  }
  
  return(complete_out)
  
}

#' Title
#'
#' @param dual_FC 
#' @param display 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#' @param pseudo_single_p1_correctedFCs 
#' @param pseudo_single_p2_correctedFCs 
#'
#' @return
#' @export
#'
#' @examples
ccr2.solveLinearSys <- function(
  dual_FC, 
  pseudo_single_p1_correctedFCs, 
  pseudo_single_p2_correctedFCs, 
  display = TRUE, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname = ""
) {
  
  ### 1. create matrix of interactions ###
  dual_FC <- dual_FC %>% 
    dplyr::mutate(sgRNA_ID_pair = paste0(sgRNA1_WGE_ID, "~", sgRNA2_WGE_ID)) 

  combinations <- t(sapply(pseudo_single_p1_correctedFCs$sgRNA_ID, 
                           function(x) paste0(x, "~", pseudo_single_p2_correctedFCs$sgRNA_ID)))
  matrix_interactions <- apply(combinations, 2, 
                               function(x) as.numeric(x %in% dual_FC$sgRNA_ID_pair))
  rownames(matrix_interactions) <- pseudo_single_p1_correctedFCs$sgRNA_ID
  colnames(matrix_interactions) <- pseudo_single_p2_correctedFCs$sgRNA_ID
  # remove rows and columns with no index, correspond to pairs including "NONTARGET"
  id_rm_row <- rowSums(matrix_interactions) == 0
  id_rm_col <- colSums(matrix_interactions) == 0
  matrix_interactions <- matrix_interactions[!id_rm_row, ]
  matrix_interactions <- matrix_interactions[,!id_rm_col]
  combinations <- combinations[!id_rm_row, ]
  combinations <- combinations[, !id_rm_col]
  
  ### 2. create matrix of system: ### 
  # n.eq=nrow(int)+ncol(int) times n.unknowns=nrow(int)*ncol(int)
  tmp_G1_vect <- list()
  for (idx_row in 1:nrow(matrix_interactions)) {
    tmp <- matrix_interactions
    tmp[-idx_row, ] <- 0
    tmp_G1_vect[[idx_row]] <- as.vector(t(tmp))
  }
  
  tmp_G2_vect <- list()
  for (idx_col in 1:ncol(matrix_interactions)) {
    tmp <- matrix_interactions
    tmp[, -idx_col] <- 0
    tmp_G2_vect[[idx_col]] <- as.vector(t(tmp))
  }
  
  matrix_sys <- rbind(do.call(rbind, tmp_G1_vect), do.call(rbind, tmp_G2_vect))
  rownames(matrix_sys) <- c(rownames(matrix_interactions), colnames(matrix_interactions))
  colnames(matrix_sys) <- as.vector(t(combinations))
  print(dim(matrix_sys))
  
  ### 3. create coefficient matrix ###
  correction_vect <- c(pseudo_single_p1_correctedFCs$correction_scaled[!id_rm_row], 
                       pseudo_single_p2_correctedFCs$correction_scaled[!id_rm_col])
  
  # the matrix is very sparse,
  # time_solve <- system.time(correction_pair <- MASS::ginv(matrix_sys) %*% correction_vect)
  sparse_matrix_sys <- Matrix(matrix_sys, sparse = TRUE)
  print(length(sparse_matrix_sys@i))
  time_solve <- system.time(sparse_inv <- spginv(sparse_matrix_sys))
  correction_pair <- sparse_inv %*% correction_vect
  print(sprintf("system solved after: %.2fs", time_solve[3]))
  
  # save output
  df_corr <- data.frame(sgRNA_ID_pair = colnames(matrix_sys), 
                        correction = correction_pair) %>%
    dplyr::filter(sgRNA_ID_pair %in% dual_FC$sgRNA_ID_pair)
  # add back pairs including non target sgrna
  tmp1 <-  pseudo_single_p1_correctedFCs %>% 
    dplyr::filter(id_rm_row) %>% 
    dplyr::select(sgRNA_ID, correction) %>%
    dplyr::mutate(sgRNA_ID_pair = dual_FC$sgRNA_ID_pair[match(sgRNA_ID, dual_FC$sgRNA1_WGE_ID)]) %>%
    dplyr::select(-sgRNA_ID)
  
  tmp2 <-  pseudo_single_p2_correctedFCs %>% 
    dplyr::filter(id_rm_col) %>% 
    dplyr::select(sgRNA_ID, correction) %>%
    dplyr::mutate(sgRNA_ID_pair = dual_FC$sgRNA_ID_pair[match(sgRNA_ID, dual_FC$sgRNA2_WGE_ID)]) %>%
    dplyr::select(-sgRNA_ID)
  
  df_corr <- rbind(df_corr, tmp1, tmp2)
  
  dual_FC_correctedFC <- dplyr::left_join(dual_FC, df_corr) %>%
    dplyr::mutate(correctedFC = avgFC + correction)
  
  pseudo_single_correction <- data.frame(
    sgRNA_ID = c(pseudo_single_p1_correctedFCs$sgRNA_ID[!id_rm_row], 
                 pseudo_single_p2_correctedFCs$sgRNA_ID[!id_rm_col]), 
    gene = c(pseudo_single_p1_correctedFCs$genes[!id_rm_row], 
             pseudo_single_p2_correctedFCs$genes[!id_rm_col]),
    pseudo_single_fitted = matrix_sys %*% correction_pair, 
    pseudo_single = correction_vect, 
    guide_pos = c(rep(1, nrow(pseudo_single_p1_correctedFCs[!id_rm_row, ])), 
                  rep(2, nrow(pseudo_single_p2_correctedFCs[!id_rm_col, ]))))
  
  # plot
  pseudo_single_correction$pseudo_single_fitted <- pseudo_single_correction$pseudo_single_fitted/c(pseudo_single_p1_correctedFCs$n[!id_rm_row],  
                                                                                   pseudo_single_p2_correctedFCs$n[!id_rm_col])
  pseudo_single_correction$pseudo_single <- pseudo_single_correction$pseudo_single/c(pseudo_single_p1_correctedFCs$n[!id_rm_row],  
                                                                     pseudo_single_p2_correctedFCs$n[!id_rm_col])
  pseudo_single_correction$guide_pos <- factor(pseudo_single_correction$guide_pos)
  pseudo_single_correction$ID <- paste0(pseudo_single_correction$gene, "_", 
                                        pseudo_single_correction$sgRNA_ID)
  pseudo_single_correction$outliers <- ""
  id_out <- abs(pseudo_single_correction$pseudo_single - pseudo_single_correction$pseudo_single_fitted) > 0.5
  pseudo_single_correction$outliers[id_out] <- pseudo_single_correction$ID[id_out]
  
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
                                  label = outliers)) + 
      geom_point(alpha = 0.7, size = 2) +
      geom_text_repel(size = 3, min.segment.length = 0, color = "black") + 
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") + 
      theme_bw() + 
      theme(legend.position = "right") + 
      xlab("Pseudo single correction") + ylab("Fitted correction")
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
    dual_FC_corr_plot <- dual_FC_correctedFC[!is.na(dual_FC_correctedFC$correction), ]
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
    matrix_system = matrix_sys, 
    dual_FC = dual_FC_correctedFC, 
    pseudo_single = pseudo_single_correction))
}

#' Title
#'
#' plot: distribution before and after based on classes
#'
#' @param dual_FC_correctedFC 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' combine copy number info with dual FC data.frame
#'
#' @param CNA 
#' @param dual_FC 
#'
#' @return
#'
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

  pos1_genom_range <- makeGRangesFromDataFrame(
    dual_FC[, c("sgRNA1_Chr", "sgRNA1_Start", "sgRNA1_End")],
    seqnames.field = "sgRNA1_Chr",
    start.field = "sgRNA1_Start",
    end.field = "sgRNA1_End")
  
  pos2_genom_range <- makeGRangesFromDataFrame(
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
      Gene1_CN = dplyr::case_when(Gene1_CN > 9 ~ 8, TRUE ~ round(Gene1_CN)), 
      Gene2_CN = dplyr::case_when(Gene2_CN > 9 ~ 8, TRUE ~ round(Gene2_CN)), 
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


#' Title
#'
#' plot: distribution with respect to CN
#'
#' @param dual_FC_correctedFC 
#' @param CNA 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#' @param excludeGene 
#' @param var_to_plot 
#'
#' @return
#' @export
#'
#' @examples
ccr2.plotCNA <- function(
  dual_FC_correctedFC, 
  CNA, 
  saveToFig = FALSE, 
  saveFormat = "pdf",
  outdir = "./", 
  EXPname ="", 
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
                        logFC = c(dual_FC_withCNA[, var_name], dual_FC_withCNA[, var_name_corrected]), 
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
    geom_boxplot(outlier.size = 1) + 
    # geom_jitter(position = position_jitter(width = 0.05), size = 0.5) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
    theme_bw() + 
    #stat_summary(fun = mean, geom = "line", aes(group = type), color = "red")  + 
    #stat_summary(fun = mean, geom = "point", aes(group = type), color = "red") +
    # geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
    theme(axis.text = element_text(size = 12),
          legend.title = element_blank(),
          legend.position = "bottom", 
          axis.text.x = element_text(angle = 0, hjust = 1)) +
    scale_fill_brewer(palette = "Paired") +
    xlab("Max CN guide1 & guide2") + 
    ylab(ylab_name) + 
    ggtitle(title_plot)
  
  df_plot$Gene2_CN <- sprintf("CN guide2: %s", df_plot$Gene2_CN)
  pl_CN_comb <- ggplot(df_plot, aes(x = Gene1_CN, y = logFC, fill = type)) + 
    geom_boxplot() + 
    #geom_jitter(position = position_jitter(width = 0.05), size = 0.5)+
    theme_bw() + 
    facet_wrap(.~Gene2_CN) +
    geom_hline(yintercept = 0, color = "red", linetype = "dashed", size = 0.8) +
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

#' Title
#'
#' plot: density divided per max CNA threshold
#'
#' @param dual_FC_correctedFC 
#' @param CNA 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#' @param excludeGene 
#' @param CN_thr 
#' @param var_to_plot 
#'
#' @return
#' @export
#'
#' @examples
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
    logFC = c(dual_FC_withCNA[, var_name], dual_FC_withCNA[, var_name_corrected]),
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

#' Title
#'
#' plot: before and after correction of matching single dual gene
#'
#' @param dual_FC_correctedFC 
#' @param match_dual_single_seq 
#' @param single_correctedFCs 
#' @param CNA 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#'
#' @return
#' @export
#'
#' @examples
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
  
  single_genom_range <- makeGRangesFromDataFrame(
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
    dplyr::mutate(Gene_CN = dplyr::case_when(
      Gene_CN > 9 ~ 8, TRUE ~ round(Gene_CN)
      ))
  
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

#' Title
#'
#' @param df 
#' @param outdir 
#' @param EXPname 
#' @param saveFormat 
#' @param saveToFig 
#'
#' @return
#' @export
#'
#' @examples
ccr2.plot_correction <- function(df, 
                                 saveToFig = FALSE, 
                                 saveFormat = "pdf",
                                 outdir = "./", 
                                 EXPname = "") {
  
  # plot synergy before and after correction:
  pl <- ggplot(df, aes(x = bliss_zscore, 
                       y = bliss_zscore_corrected, 
                       color = info_subtype)) + 
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


#' Title
#'
#' @param dual_FC 
#' @param corrected 
#' @param THR_FC 
#' @param THR_BLISS 
#' @param saveToFig 
#' @param saveFormat 
#' @param outdir 
#' @param EXPname 
#'
#' @return
#' @export
#'
#' @examples
ccr2.plot_bliss_vs_FC <- function(dual_FC, corrected = FALSE, 
                                  THR_FC = -1,
                                  THR_BLISS = -1, 
                                  saveToFig = FALSE, 
                                  saveFormat = "pdf",
                                  outdir = "./", 
                                  EXPname = ""){
  
  dual_FC <- dual_FC[!is.na(dual_FC$correction),]
  
  name_var_y <- "avgFC_scaled"
  name_var_x <- "bliss_zscore"
  title_pl <- "Uncorrected"
  if (corrected) {
    name_var_y <- "correctedFC_scaled"
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


#' Title
#'
#' @param EXPname 
#' @param display 
#' @param outdir 
#' @param saveToFig 
#' @param saveFormat 
#' @param libraryAnnotation_dual 
#' @param dual_FC
#' @param single_correctedFCs 
#' @param match_dual_single_seq 
#'
#' @return
#' @export
#'
#' @examples
ccr2.run <- function(
  single_correctedFCs, 
  libraryAnnotation_dual, 
  dual_FC,
  match_dual_single_seq, 
  EXPname = "",
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = NULL, 
  outdir = "./", 
  correctGW
) {
  
  # exclude non-target pairs
  id_nontarget <- grepl("NONTARGET", dual_FC$info_subtype)
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
    correctGW = correctGW, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
  model_guide2 <- ccr2.modelSingleVSPseudoSingle(
    pseudo_single_FC = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    display = display, 
    correctGW = correctGW, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname)
  
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
    pseudo_single = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display)
  
  pseudo_single_p2_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide2_correctedFCs$corrected_logFCs,
    pseudo_single = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display)
  
  print("CRISPRcleanR applied to GW pseudo singles (position 1 and 2)")
  
  ##################################################
  ### solve linear system to get pair correction ###
  sys_solution <- ccr2.solveLinearSys(
    dual_FC = dual_FC, 
    pseudo_single_p1_correctedFCs = pseudo_single_p1_correctedFCs, 
    pseudo_single_p2_correctedFCs = pseudo_single_p2_correctedFCs, 
    saveToFig = saveToFig, 
    saveFormat = saveFormat,
    outdir = outdir, 
    EXPname = EXPname,
    display = display)
  print("System solved, collpased correction converted to pairwise correction")
  
  pseudo_single_correction <- sys_solution$pseudo_single
  dual_FC_correctedFC <- sys_solution$dual_FC
  
  return(
    list(pseudo_single = pseudo_single_correction, 
        dual_FC = dual_FC_correctedFC)
    )
}


#' Title
#'
#' @param EXPname 
#' @param display 
#' @param outdir 
#' @param saveToFig 
#' @param saveFormat 
#' @param libraryAnnotation_dual 
#' @param dual_FC
#' @param single_correctedFCs 
#' @param match_dual_single_seq 
#'
#' @return
#' @export
#'
#' @examples
ccr2.run_nontarget <- function(
  single_correctedFCs, 
  libraryAnnotation_dual, 
  dual_FC,
  match_dual_single_seq, 
  EXPname = "",
  saveToFig = FALSE, 
  display = TRUE, 
  saveFormat = NULL, 
  outdir = "./", 
  correctGW
) {
  
  # process only gene_nontarget or nontarget_gene
  id_nontarget <- grepl("NONTARGET", dual_FC$info_subtype) & 
    dual_FC$info_subtype != "NONTARGET-NONTARGET"
  
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
    pseudo_single = dual_pseudo_single_FC$sgRNA1,
    guide_id = 1, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  pseudo_single_p2_correctedFCs <- ccr2.filterGWclean(
    dataInjection_correctedFCs = dataInjection_guide2_correctedFCs$corrected_logFCs,
    pseudo_single = dual_pseudo_single_FC$sgRNA2,
    guide_id = 2, 
    saveToFig = saveToFig, 
    display = display,  
    saveFormat = saveFormat,
    outdir = sprintf("%sNONTARGET_PAIR_", outdir), 
    EXPname = EXPname)
  
  # assign results
  dual_FC_nt <- dual_FC_nt %>% 
    dplyr::mutate(sgRNA_ID_pair = paste0(sgRNA1_WGE_ID, "~", sgRNA2_WGE_ID))
  
  tmp1 <- pseudo_single_p1_correctedFCs %>% 
    dplyr::select(sgRNA_ID, correction) %>%
    dplyr::mutate(sgRNA_ID_pair = dual_FC_nt$sgRNA_ID_pair[match(sgRNA_ID, dual_FC_nt$sgRNA1_WGE_ID)]) %>%
    dplyr::select(-sgRNA_ID)
  
  tmp2 <- pseudo_single_p2_correctedFCs %>% 
    dplyr::select(sgRNA_ID, correction) %>%
    dplyr::mutate(sgRNA_ID_pair = dual_FC_nt$sgRNA_ID_pair[match(sgRNA_ID, dual_FC_nt$sgRNA2_WGE_ID)]) %>%
    dplyr::select(-sgRNA_ID)
  
  df_corr <- rbind(tmp1, tmp2)
  dual_FC_correctedFC <- dplyr::left_join(dual_FC_nt, df_corr) %>%
    dplyr::mutate(correctedFC = avgFC + correction)
  
  return(dual_FC_correctedFC)
  
}


#' Title
#'
#' @param filename_single 
#' @param min_reads_single 
#' @param EXPname 
#' @param libraryAnnotation_single 
#' @param display 
#' @param outdir 
#' @param min_reads 
#' @param saveToFig 
#' @param saveFormat 
#' @param libraryAnnotation_dual 
#' @param dual_count 
#' @param correctGW 
#' @param excludeGene_plot 
#' @param CNA 
#' @param CN_thr 
#'
#' @return
#' @export
#'
#' @examples
ccr2.run_complete <- function(
  filename_single, 
  min_reads_single = 30,
  libraryAnnotation_single,  
  min_reads = 30, 
  libraryAnnotation_dual, 
  dual_count, 
  EXPname,
  display = FALSE, 
  outdir = "./", 
  saveToFig, 
  saveFormat = "pdf", 
  correctGW, 
  excludeGene_plot = NULL, 
  CNA, 
  CN_thr = 8
) {
  
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
  dual_FC <- ccr2.NormfoldChanges(
    Dframe = dual_count, 
    min_reads = min_reads)
  
  dual_FC <- ccr2.logFCs2chromPos(
    dual_FC = dual_FC$logFCs, 
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
  dual_FC_nt_correctedFC <- ccr2.run_nontarget(
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
    correctGW = correctGW
  )
  dual_FC_correctedFC <- tmp$dual_FC
  pseudo_single_correction <- tmp$pseudo_single
  
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
  ccr2.plot_bliss_fit(dual_FC = dual_FC_correctedFC, 
                      EXPname = EXPname, 
                      saveToFig = saveToFig, 
                      saveFormat = saveFormat,
                      outdir = outdir)
  
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
  
  # density distribution for CN > 8 for any guide1 or guide2
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
                                         EXPname = EXPname, 
                                         saveToFig = saveToFig, 
                                         saveFormat = saveFormat,
                                         outdir = outdir)
  
  top_uncorrected <- ccr2.plot_bliss_vs_FC(dual_FC = dual_FC_correctedFC, 
                                           corrected = FALSE, 
                                           EXPname = EXPname, 
                                           saveToFig = saveToFig, 
                                           saveFormat = saveFormat,
                                           outdir = outdir)
  
  return(list(dual = dual_FC_correctedFC, 
              single = single_correctedFCs_filt,
              top_corrected = top_corrected, 
              top_uncorrected = top_uncorrected, 
              system_solition = pseudo_single_correction))
  
}

#' Title
#'
#' @param x 
#'
#' @return
#'
spginv <- function(x) {
  Xsvd <- sparsesvd::sparsesvd(x)
  Positive <- Xsvd$d > max(sqrt(.Machine$double.eps) * Xsvd$d[1L], 0)
  if (all(Positive))
    Xsvd$v %*% (1/Xsvd$d * t(Xsvd$u))
  else if (!any(Positive))
    array(0, dim(x)[2L:1L])
  else Xsvd$v[, Positive, drop = FALSE] %*% ((1/Xsvd$d[Positive]) * t(Xsvd$u[, Positive, drop = FALSE]))
}


