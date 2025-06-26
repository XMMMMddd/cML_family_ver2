library(tibble)
generate_mr_trio_data <- function(
    # 基本参数
    n = 1000, seed = NULL, maf = 0.3,
    # --- 暴露(Exposure)的效应系数 ---
    beta_os_to_oe_exp = c(0.2, 0.2), # 子女自身SNP
    beta_fs_to_oe_exp = c(0.05, 0.05), # 父亲SNP (遗传叠加)
    beta_ms_to_oe_exp = c(0.05, 0.05), # 母亲SNP (遗传叠加)
    sigma_confounding_exp = 0.1, # 混杂U的方差
    sigma_correlation_exp = 0.1, # 共享环境C的方差
    # --- 结局(Outcome)的效应系数 ---
    beta_os_to_oe_out = c(0.1, 0.1), # 水平多效性
    beta_fs_to_oe_out = c(0.02, 0.1),
    beta_ms_to_oe_out = c(0.02, 0.1),
    sigma_confounding_out = 0.1, # 混杂U对结局的方差
    sigma_correlation_out = 0.1, # 共享环境C对结局的方差
    # --- 直接设定的参数 ---
    beta_exp_to_out = 0.4, # 因果效应
    sigma_lon_sq = 1, # 暴露的独立残差方差 (基准)
    sigma_psi_sq = 1, # 结局的独立残差方差 (基准)
    # --- 模拟流程控制参数 ---
    assortative_mating_strength = 1000) {
    # ===================================================================
    # 1. 初始化和内置辅助函数
    # ===================================================================
    set.seed(seed)
    n_snps <- length(beta_os_to_oe_exp)
    beta_zero_vector <- rep(0, n_snps) # 用于无父代效应的场景

    # --- 内置函数 1: 生成SNP矩阵 (修正了之前缺失的问题) ---
    generate_snp_hwe_matrix <- function(n_samples, maf, n_snps) {
        matrix(rbinom(n_samples * n_snps, size = 2, prob = maf),
            nrow = n_samples, ncol = n_snps
        )
    }

    # --- 内置函数 2: 模拟减数分裂传递等位基因 ---
    get_transmitted_allele <- function(parent_snps) {
        transmitted_alleles <- parent_snps / 2
        het_positions <- which(parent_snps == 1)
        if (length(het_positions) > 0) {
            transmitted_alleles[het_positions] <- rbinom(length(het_positions), 1, 0.5)
        }
        return(transmitted_alleles)
    }

    # --- 内置函数 3: 生成表型 (暴露/结局) ---
    generate_phenotype <- function(snps, snps_f, snps_m,
                                   beta_o, beta_f, beta_m,
                                   confounder, residual_var,
                                   causal_term = 0) {
        genetic_effect <- snps %*% beta_o + snps_f %*% beta_f + snps_m %*% beta_m
        residual_effect <- rnorm(nrow(snps), mean = 0, sd = sqrt(residual_var))

        phenotype <- causal_term + genetic_effect + confounder + residual_effect
        return(as.vector(phenotype))
    }

    # --- 内置函数 4: 选型婚配 (修正了对混杂U的处理) ---
    assortative_mating_function <- function(d1, d2, am_strength) {
        n <- nrow(d1$snps)
        # 为排序添加随机噪声
        d1$sort_pheno <- d1$expose + rnorm(n, 0, sd = sqrt(am_strength))
        d2$sort_pheno <- d2$expose + rnorm(n, 0, sd = sqrt(am_strength))

        # 根据加噪的表型排序
        order1 <- order(d1$sort_pheno)
        order2 <- order(d2$sort_pheno)

        # 重新排列d1和d2的所有数据，确保所有个体信息(包括混杂U)被一同洗牌
        d1_sorted <- lapply(d1, function(x) if (is.matrix(x)) x[order1, , drop = FALSE] else x[order1])
        d2_sorted <- lapply(d2, function(x) if (is.matrix(x)) x[order2, , drop = FALSE] else x[order2])

        return(list(d1_sorted, d2_sorted))
    }

    # ===================================================================
    # 2. 生成初始人群和混杂因素
    # ===================================================================
    # 生成n个个体的稳定混杂因素U (会被代际传递和选型婚配影响)
    confounding_exp_pop <- rnorm(n, 0, sd = sqrt(sigma_confounding_exp))
    confounding_out_pop <- rnorm(n, 0, sd = sqrt(sigma_confounding_out))

    # --- 祖父母代 ---
    gpf_data <- list(
        snps = generate_snp_hwe_matrix(n, maf, n_snps),
        confound_exp = confounding_exp_pop, confound_out = confounding_out_pop
    )
    gmf_data <- list(
        snps = generate_snp_hwe_matrix(n, maf, n_snps),
        confound_exp = confounding_exp_pop, confound_out = confounding_out_pop
    )
    gpm_data <- list(
        snps = generate_snp_hwe_matrix(n, maf, n_snps),
        confound_exp = confounding_exp_pop, confound_out = confounding_out_pop
    )
    gmm_data <- list(
        snps = generate_snp_hwe_matrix(n, maf, n_snps),
        confound_exp = confounding_exp_pop, confound_out = confounding_out_pop
    )

    # 为祖父母生成表型 (他们没有父代，所以父代效应为0)
    gpf_data$expose <- generate_phenotype(
        gpf_data$snps,
        gpf_data$snps, gpf_data$snps, beta_os_to_oe_exp,
        beta_zero_vector, beta_zero_vector, gpf_data$confound_exp, sigma_lon_sq
    )
    gmf_data$expose <- generate_phenotype(
        gmf_data$snps,
        gmf_data$snps, gmf_data$snps, beta_os_to_oe_exp,
        beta_zero_vector, beta_zero_vector, gmf_data$confound_exp, sigma_lon_sq
    )
    gpm_data$expose <- generate_phenotype(
        gpm_data$snps,
        gpm_data$snps, gpm_data$snps, beta_os_to_oe_exp,
        beta_zero_vector, beta_zero_vector, gpm_data$confound_exp, sigma_lon_sq
    )
    gmm_data$expose <- generate_phenotype(gmm_data$snps, gmm_data$snps, gmm_data$snps, beta_os_to_oe_exp, beta_zero_vector, beta_zero_vector, gmm_data$confound_exp, sigma_lon_sq)

    # ===================================================================
    # 3. 第一轮选型婚配 (祖父母 -> 父母)
    # ===================================================================
    gpf_gmf_mated <- assortative_mating_function(gpf_data, gmf_data, assortative_mating_strength)
    gpm_gmm_mated <- assortative_mating_function(gpm_data, gmm_data, assortative_mating_strength)

    # --- 父母代 ---
    # 父母继承了祖父母的基因和混杂U (经过了洗牌)
    f_data <- list(
        snps = get_transmitted_allele(gpf_gmf_mated[[1]]$snps) + get_transmitted_allele(gpf_gmf_mated[[2]]$snps),
        confound_exp = gpf_gmf_mated[[1]]$confound_exp, # 继承父系混杂
        confound_out = gpf_gmf_mated[[1]]$confound_out
    )
    m_data <- list(
        snps = get_transmitted_allele(gpm_gmm_mated[[1]]$snps) + get_transmitted_allele(gpm_gmm_mated[[2]]$snps),
        confound_exp = gpm_gmm_mated[[1]]$confound_exp, # 继承母系混杂
        confound_out = gpm_gmm_mated[[1]]$confound_out
    )

    # 为父母生成暴露 (用于他们之间的婚配排序)
    f_data$expose <- generate_phenotype(f_data$snps, gpf_gmf_mated[[1]]$snps, gpf_gmf_mated[[2]]$snps, beta_os_to_oe_exp, beta_fs_to_oe_exp, beta_ms_to_oe_exp, f_data$confound_exp, sigma_lon_sq)
    m_data$expose <- generate_phenotype(m_data$snps, gpm_gmm_mated[[1]]$snps, gpm_gmm_mated[[2]]$snps, beta_os_to_oe_exp, beta_fs_to_oe_exp, beta_ms_to_oe_exp, m_data$confound_exp, sigma_lon_sq)

    # ===================================================================
    # 4. 第二轮选型婚配 (父母 -> 子女)
    # ===================================================================
    f_m_mated <- assortative_mating_function(f_data, m_data, assortative_mating_strength)
    father_final <- f_m_mated[[1]]
    mother_final <- f_m_mated[[2]]

    # --- 子女代 ---
    offspring_snps <- get_transmitted_allele(father_final$snps) + get_transmitted_allele(mother_final$snps)

    # !!! 逻辑修正: 在父母配对完成后，为新家庭生成共享环境C !!!
    correlation_exp_final <- rnorm(n, 0, sd = sqrt(sigma_correlation_exp))
    correlation_out_final <- rnorm(n, 0, sd = sqrt(sigma_correlation_out))

    # ===================================================================
    # 5. 生成最终的暴露和结局数据
    # ===================================================================
    # --- 暴露 ---
    offspring_expose_no_corr <- generate_phenotype(offspring_snps, father_final$snps, mother_final$snps, beta_os_to_oe_exp, beta_fs_to_oe_exp, beta_ms_to_oe_exp, father_final$confound_exp, sigma_lon_sq)
    offspring_expose <- offspring_expose_no_corr + correlation_exp_final

    father_expose <- father_final$expose + correlation_exp_final
    mother_expose <- mother_final$expose + correlation_exp_final

    # --- 结局 ---
    offspring_outcome_no_corr <- generate_phenotype(offspring_snps, father_final$snps, mother_final$snps, beta_os_to_oe_out, beta_fs_to_oe_out, beta_ms_to_oe_out, father_final$confound_out, sigma_psi_sq, causal_term = beta_exp_to_out * offspring_expose)
    offspring_outcome <- offspring_outcome_no_corr + correlation_out_final

    father_outcome <- generate_phenotype(father_final$snps, gpf_gmf_mated[[1]]$snps, gpf_gmf_mated[[2]]$snps, beta_os_to_oe_out, beta_fs_to_oe_out, beta_ms_to_oe_out, father_final$confound_out, sigma_psi_sq, causal_term = beta_exp_to_out * father_expose) + correlation_out_final
    mother_outcome <- generate_phenotype(mother_final$snps, gpm_gmm_mated[[1]]$snps, gpm_gmm_mated[[2]]$snps, beta_os_to_oe_out, beta_fs_to_oe_out, beta_ms_to_oe_out, mother_final$confound_out, sigma_psi_sq, causal_term = beta_exp_to_out * mother_expose) + correlation_out_final

    # ===================================================================
    # 6. 计算各组分的方差解释比例
    # ===================================================================
    # --- 暴露X的方差分解 ---
    var_x_total <- var(offspring_expose)
    var_x_components <- tibble(
        Variable = "Exposure",
        Component = c(
            "Offspring SNPs", "Father SNPs (Nurture)", "Mother SNPs (Nurture)",
            "Confounding (U)", "Shared Environment (C)", "Residual (epsilon)"
        ),
        Variance = c(
            var(offspring_snps %*% beta_os_to_oe_exp),
            var(father_final$snps %*% beta_fs_to_oe_exp),
            var(mother_final$snps %*% beta_ms_to_oe_exp),
            var(father_final$confound_exp),
            var(correlation_exp_final),
            sigma_lon_sq # 使用理论值，因经验值不易直接获得
        )
    )
    var_x_components$Proportion <- var_x_components$Variance / var_x_total

    # --- 结局Y的方差分解 ---
    var_y_total <- var(offspring_outcome)
    var_y_components <- tibble(
        Variable = "Outcome",
        Component = c(
            "Causal Effect of X", "Offspring SNPs (Pleiotropy)", "Father SNPs (Pleiotropy)", "Mother SNPs (Pleiotropy)",
            "Confounding (U)", "Shared Environment (C)", "Residual (psi)"
        ),
        Variance = c(
            var(beta_exp_to_out * offspring_expose),
            var(offspring_snps %*% beta_os_to_oe_out),
            var(father_final$snps %*% beta_fs_to_oe_out),
            var(mother_final$snps %*% beta_ms_to_oe_out),
            var(father_final$confound_out),
            var(correlation_out_final),
            sigma_psi_sq # 使用理论值
        )
    )
    var_y_components$Proportion <- var_y_components$Variance / var_y_total

    variance_explained <- rbind(var_x_components, var_y_components)

    # ===================================================================
    # 7. 最终返回结果
    # ===================================================================
    # --- 创建最终的数据集 (修正了之前data.frame的问题) ---
    final_data <- tibble(
        father_snps = I(father_final$snps),
        mother_snps = I(mother_final$snps),
        offspring_snps = I(offspring_snps),
        father_expose = father_expose,
        mother_expose = mother_expose,
        offspring_expose = offspring_expose,
        father_outcome = father_outcome,
        mother_outcome = mother_outcome,
        offspring_outcome = offspring_outcome
    )

    return(list(
        data = final_data,
        variance_explained = variance_explained
    ))
}

# %% 系数计算函数

generate_mr_trio_data_matrix_ultra <- function(
    # 基本参数
    n = 1000, seed = NULL, maf = 0.3,
    # --- 暴露(Exposure)的效应系数 ---
    beta_os_to_oe_exp = c(0.2, 0.2), # 子女自身SNP
    beta_fs_to_oe_exp = c(0.05, 0.05), # 父亲SNP (遗传叠加)
    beta_ms_to_oe_exp = c(0.05, 0.05), # 母亲SNP (遗传叠加)
    sigma_confounding_exp = 0.1, # 混杂U的方差
    sigma_correlation_exp = 0.1, # 共享环境C的方差
    # --- 暴露(Exposure)的效应系数（异质性） ---
    h_beta_os_to_oe_exp = c(0.2, 0.2), # 子女自身SNP
    h_beta_fs_to_oe_exp = c(0.05, 0.05), # 父亲SNP (遗传叠加)
    h_beta_ms_to_oe_exp = c(0.05, 0.05), # 母亲SNP (遗传叠加)
    # --- 结局(Outcome)的效应系数 ---
    beta_os_to_oe_out = c(0.1, 0.1), # 水平多效性
    beta_fs_to_oe_out = c(0.02, 0.1),
    beta_ms_to_oe_out = c(0.02, 0.1),
    sigma_confounding_out = 0.1, # 混杂U对结局的方差
    sigma_correlation_out = 0.1, # 共享环境C对结局的方差
    # --- 直接设定的参数 ---
    beta_exp_to_out = 0.4, # 因果效应
    sigma_lon_sq = 1, # 暴露的独立残差方差 (基准)
    sigma_psi_sq = 1, # 结局的独立残差方差 (基准)
    # --- 模拟流程控制参数 ---
    assortative_mating_strength = 1000) {
    # --- 1. 生成两个数据 ---

    data_exp <- generate_mr_trio_data(
        # 基本参数
        n = n, seed = seed, maf = maf,
        # --- 暴露(Exposure)的效应系数 ---
        beta_os_to_oe_exp = beta_os_to_oe_exp, # 子女自身SNP
        beta_fs_to_oe_exp = beta_fs_to_oe_exp, # 父亲SNP (遗传叠加)
        beta_ms_to_oe_exp = beta_ms_to_oe_exp, # 母亲SNP (遗传叠加)
        sigma_confounding_exp = sigma_confounding_exp, # 混杂U的方差
        sigma_correlation_exp = sigma_correlation_exp, # 共享环境C的方差
        # --- 结局(Outcome)的效应系数 ---
        beta_os_to_oe_out = beta_os_to_oe_out, # 水平多效性
        beta_fs_to_oe_out = beta_fs_to_oe_out,
        beta_ms_to_oe_out = beta_ms_to_oe_out,
        sigma_confounding_out = sigma_confounding_out, # 混杂U对结局的方差
        sigma_correlation_out = sigma_correlation_out, # 共享环境C对结局的方差
        # --- 直接设定的参数 ---
        beta_exp_to_out = beta_exp_to_out, # 因果效应
        sigma_lon_sq = sigma_lon_sq, # 暴露的独立残差方差 (基准)
        sigma_psi_sq = sigma_psi_sq, # 结局的独立残差方差 (基准)
        # --- 模拟流程控制参数 ---
        assortative_mating_strength = assortative_mating_strength
    )

    data_out <- generate_mr_trio_data(
        # 基本参数
        n = n, seed = seed, maf = maf,
        # --- 暴露(Exposure)的效应系数 ---
        beta_os_to_oe_exp = h_beta_os_to_oe_exp, # 子女自身SNP
        beta_fs_to_oe_exp = h_beta_fs_to_oe_exp, # 父亲SNP (遗传叠加)
        beta_ms_to_oe_exp = h_beta_ms_to_oe_exp, # 母亲SNP (遗传叠加)
        sigma_confounding_exp = sigma_confounding_exp, # 混杂U的方差
        sigma_correlation_exp = sigma_correlation_exp, # 共享环境C的方差
        # --- 结局(Outcome)的效应系数 ---
        beta_os_to_oe_out = beta_os_to_oe_out, # 水平多效性
        beta_fs_to_oe_out = beta_fs_to_oe_out,
        beta_ms_to_oe_out = beta_ms_to_oe_out,
        sigma_confounding_out = sigma_confounding_out, # 混杂U对结局的方差
        sigma_correlation_out = sigma_correlation_out, # 共享环境C对结局的方差
        # --- 直接设定的参数 ---
        beta_exp_to_out = beta_exp_to_out, # 因果效应
        sigma_lon_sq = sigma_lon_sq, # 暴露的独立残差方差 (基准)
        sigma_psi_sq = sigma_psi_sq, # 结局的独立残差方差 (基准)
        # --- 模拟流程控制参数 ---
        assortative_mating_strength = assortative_mating_strength
    )
    results <- list(data_exp = data_exp, data_out = data_out)
    return(results)
}
if (FALSE) {
    test_1 <- generate_mr_trio_data_matrix_ultra()
    test_1
}
