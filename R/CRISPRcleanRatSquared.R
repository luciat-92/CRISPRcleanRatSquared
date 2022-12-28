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