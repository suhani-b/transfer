#' @description
#'  `test_interaction.poisson.single` fits a gamma-poisson glm for each gene with respect to an interaction between two spectra factors
#'
#'  @details Fits poisson glm for all genes in anndata expressed in some number of cells for a given cell type. Regression is tailored to understand main effects of two cell type specific factors and their interaction.
#'  Additional factors and covariates can be controlled as well as gene detection rate.
#' @param adata anndata object; same number of obs as cell scores
#' @param cell_scores data.frame; spectra cell scores matrix
#' @param factor_pairs 2xn matrix;  each column is a pair of factors to compute interaction
#' @param cell_type character string; cell type label to filter for in adata$obs$obs_key
#' @param obs_key character string; key for cell type label column in adata.obs
#' @param layer_key character string; key for layer of adata object to use for expression values
#' @param cell_thresh integer; minimum number of cells a transcript has to be detected to be included in the test
#' @param use_gdr boolean;  whether or not to use the gene detection rate (number of genes per cell) as a covariate
#' @param additional_covariates character vector; additional covariates to control for in regression model in adata.obs
#' @param control_other_factors character vector; all factors to control for in regression
#'
#' @return list of glm fit and DE results for each term in the model

test_interaction.poisson.single <- function(dat, obs,
                                            cell_scores,
                                            factor_pairs,
                                            cell_type,
                                            obs_key = "subtype_broad",
                                            layer_key = "counts",
                                            cell_thresh = 50,
                                            use_gdr = FALSE,
                                            # rm_gene_regex = "^AB_",
                                            additional_covariates = NULL,
                                            control_factors = NULL,
                                            seed_use = 1234) {
    set.seed(seed_use)

    # set up expression object
    # filter genes before converting to dense to preserve memory
    dat_b <- dat > 0
    gdr <- Matrix::rowSums(dat_b)
    gene_count <- Matrix::colSums(dat_b)
    gene_keep <- gene_count >= cell_thresh
    dat_f <- dat[, gene_keep]

    # free up memory
    rm(dat_b, dat)
    gc()

    # convert to dense matrix
    dat_f <- Matrix::as.matrix(dat_f)
    # scale expression values
    ## remain unscaled for now
    # dat_f = dat_f %>% base::scale()

    # if some cells are removed from cell scores that are present in anndata then I need to make sure those are accounted for in the expression object
    dat_f <- dat_f[rownames(dat_f) %in% rownames(cell_scores), ]

    # combine with cell scores of factors and gene detection rate
    # dat_f = cbind(dat_f, cell_scores[rownames(dat_f),], list(GDR = gdr), adata$obs[rownames(dat_f),])
    # create metadata df
    meta <- cbind(cell_scores[rownames(dat_f), ], list(GDR = gdr), obs[rownames(dat_f), ])

    # loop through each column of factor matrix and compute interaction score
    res_i <- list()

    # get test genes which are all that are included in the expression data
    # gene_test = adata$var_names[gene_keep]
    ## for troubleshoot
    # gene_test = adata$var_names[1:10]
    # gene_test = gene_test[1:10]

    int_term <- paste(apply(factor_pairs, 2, function(x) {
        paste(paste0("`", x, "`"), collapse = "*")
    }), collapse = " + ")
    int_factors <- unique(factor_pairs)
    # make model formula with all terms
    # the response terms cannot be all specified at once with cbind in a glm so there needs to be a for loop for each DV
    if (use_gdr == TRUE) {
        model_formula <- paste0(" ~ 1 + ", int_term, " + GDR")
    }
    if (use_gdr == FALSE) {
        model_formula <- paste0(" ~ 1 + ", int_term)
    }
    if (!(is.null(control_factors))) {
        add_factors <- paste(lapply(control_factors[!(control_factors %in% int_factors)], function(x) {
            paste0("`", x, "`")
        }) %>% unlist(), collapse = " + ")

        model_formula <- paste(model_formula, add_factors, sep = " + ")
    }
    if (!(is.null(additional_covariates))) {
        model_formula <- paste(model_formula, paste(additional_covariates, collapse = " + "), sep = " + ")
    }

    print(model_formula)
    print("fitting glm")

    fit.gp <- glm_gp(dat_f,
        design = as.formula(model_formula), col_data = meta,
        # size_factors = 'deconvolution',
        size_factors = "normed_sum",
        overdispersion = T,
        overdispersion_shrinkage = T,
        on_disk = F,
        verbose = T
    )

    # modify names in fit.gp to allow for de testing
    colnames(fit.gp$model_matrix) <- make.names(colnames(fit.gp$model_matrix))

    # test DE for each term that I want to in the model
    res_i <- list()
    res_i[["model"]] <- fit.gp
    res_i[["de"]] <- list()

    # get all cytokine factors and their interaction terms for testing
    int_terms_test <- gsub("[`\\+]", "", int_term) %>%
        gsub("\\*", ":", .) %>%
        strsplit(., " ") %>%
        unlist()
    int_terms_test <- int_terms_test[int_terms_test != ""]
    indiv_terms_test <- gsub("[`\\+]", "", int_term) %>%
        gsub("\\*", ":", .) %>%
        strsplit(., " ") %>%
        unlist() %>%
        strsplit(., ":") %>%
        unlist() %>%
        unique()
    test_terms <- c(test_terms, indiv_terms_test)

    print("testing DE")
    for (i in test_terms) {
        print(i)
        res_i[["de"]][[i]] <- test_de(fit.gp, i)
    }

    return(res_i)
}
