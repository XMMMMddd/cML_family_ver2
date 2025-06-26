# %% 加载必要的函数
library(dplyr)
library(MASS)
library(nlme)
library(parallel)
library(pbapply)
# install.packages("devtools") # 安装其他代码的包
# library(devtools)
# devtools::install_github("xue-hr/MRcML")
library(MRcML) # cML的实现
# devtools::install_github("rondolab/MR-PRESSO") # MR-PRSSO实现
library(MRPRESSO) # MR_PRESSO 实现
# install.packages("MendelianRandomization")
library(MendelianRandomization)
# devtools::install_github("gqi/MRMix")
library(MRMix)
# devtools::install_github("noahlorinczcomi/MRBEE")
library(MRBEE) # 暂时没用到
# %% 加载必要的包
source("三联体数据生成函数/三联体数据生成函数.r") # 加载数据生成函数 library(dplyr) library(MASS) 里面加载了这些包
source("FGWAS处理函数/FGWAS.r") # 加载 FGWAS 函数1
# source("cml家庭函数/cml家庭函数二代/cml_family_version2.r")
source("cML-family/cml_family_2_b.r")
# %% 辅助函数

cor2cov_big_matrix <- function(cor_matrix, param_se) {


}
## 计算se和估计值的函数
calculate_mr_stats <- function(estimate, se) {
    # 计算 Z 分数
    z_score <- estimate / se
    # 计算 P 值 (双侧检验)
    p_val <- 2 * pnorm(abs(z_score), lower.tail = FALSE)
    # 返回包含 Z 分数和 P 值的列表
    return(list(z = z_score, p_value = p_val))
}

triplet_family_simulation_once_robust <- function(
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
    # ======================== 核心修改点：1. 定义辅助函数 ========================
    # 此函数确保任何赋值给 results 列表的值，如果是 NULL 或长度为0，都会被转换为 NA
    sanitize_value <- function(val) {
        if (is.null(val) || length(val) == 0) {
            return(NA_real_) # 使用 NA_real_ 保证数值类型
        }
        return(val)
    }
    # =========================================================================

    # 首先声明结果
    results <- list(
        a_theta_point = NA, a_theta_se = NA, a_z = NA,
        a_p_value = NA, a_duration = NA,
        b_theta_point = NA, b_theta_se = NA, b_z = NA,
        b_p_value = NA, b_duration = NA,
        c_theta_point = NA, c_theta_se = NA, c_z = NA,
        c_p_value = NA, c_duration = NA,
        d_theta_point = NA, d_theta_se = NA, d_z = NA,
        d_p_value = NA, d_duration = NA,
        e_theta_point = NA, e_theta_se = NA, e_z = NA,
        e_p_value = NA, e_duration = NA,
        f_theta_point = NA, f_theta_se = NA, f_z = NA,
        f_p_value = NA, f_duration = NA,
        g_theta_point = NA, g_theta_se = NA, g_z = NA,
        g_p_value = NA, g_duration = NA,
        h_theta_point = NA, h_theta_se = NA, h_z = NA,
        h_p_value = NA, h_duration = NA,
        i_theta_point = NA, i_theta_se = NA, i_z = NA,
        i_p_value = NA, i_duration = NA
    ) # 初始化结果列表

    # 生成数据 (这部分代码不变)
    phase_two_data_full <- generate_mr_trio_data_matrix_ultra(
        n = n, seed = seed, maf = maf,
        beta_os_to_oe_exp = beta_os_to_oe_exp,
        beta_fs_to_oe_exp = beta_fs_to_oe_exp,
        beta_ms_to_oe_exp = beta_ms_to_oe_exp,
        sigma_confounding_exp = sigma_confounding_exp,
        sigma_correlation_exp = sigma_correlation_exp,
        h_beta_os_to_oe_exp = h_beta_os_to_oe_exp,
        h_beta_fs_to_oe_exp = h_beta_fs_to_oe_exp,
        h_beta_ms_to_oe_exp = h_beta_ms_to_oe_exp,
        beta_os_to_oe_out = beta_os_to_oe_out,
        beta_fs_to_oe_out = beta_fs_to_oe_out,
        beta_ms_to_oe_out = beta_ms_to_oe_out,
        sigma_confounding_out = sigma_confounding_out,
        sigma_correlation_out = sigma_correlation_out,
        beta_exp_to_out = beta_exp_to_out,
        sigma_lon_sq = sigma_lon_sq,
        sigma_psi_sq = sigma_psi_sq,
        assortative_mating_strength = assortative_mating_strength
    )
    data_exp <- phase_two_data_full$data_exp
    data_out <- phase_two_data_full$data_out
    phase_two_data_analysis_exp <- fgwas_for_data_matrix(data_exp$data, processing_func = FMR_trio_IFGLS, predicted_outcome = "expose")
    phase_two_data_analysis_out <- fgwas_for_data_matrix(data_out$data, processing_func = FMR_trio_IFGLS, predicted_outcome = "outcome")
    phase_two_data_analysis_exp_gls <- fgwas_for_data_matrix_single(data_exp$data, predicted_outcome = "expose")
    phase_two_data_analysis_out_gls <- fgwas_for_data_matrix_single(data_out$data, predicted_outcome = "outcome")

    # -------- FGWAS为基础 --------

    # ===================== 核心修改点：2. 包裹每一个赋值操作 =====================

    # a: cml_overlap
    tryCatch(
        {
            res_a_time_1 <- Sys.time()
            res_a <- mr_cML_Overlap(phase_two_data_analysis_exp$beta_hat[, 1], phase_two_data_analysis_out$beta_hat[, 1], phase_two_data_analysis_exp$beta_hat_se, phase_two_data_analysis_out$beta_hat_se, n = n, rho = 0)
            res_a_time_2 <- Sys.time()
            res_a_estimtor <- calculate_mr_stats(res_a$MA_BIC_theta, res_a$MA_BIC_se)
            results$a_theta_point <- sanitize_value(res_a$MA_BIC_theta)
            results$a_theta_se <- sanitize_value(res_a$MA_BIC_se)
            results$a_z <- sanitize_value(res_a_estimtor$z)
            results$a_p_value <- sanitize_value(res_a$MA_BIC_p)
            results$a_duration <- sanitize_value(as.numeric(res_a_time_2 - res_a_time_1))
        },
        error = function(e) {
            results$a_duration <<- NA
        }
    )

    # b: cml_family_ver2 方法
    tryCatch(
        {
            res_b_time_1 <- Sys.time()
            res_b <- cml_family_ultral_b(beta_hat_exp = phase_two_data_analysis_exp$beta_hat, beta_hat_out = phase_two_data_analysis_out$beta_hat, beta_sigma_exp = phase_two_data_analysis_exp$Sigma_inv, beta_sigma_out = phase_two_data_analysis_out$Sigma_inv, n = n, cov_index = "normal", alpha_range = c(-0.000001, 0.000001))
            res_b_time_2 <- Sys.time()
            res_b_estimtor <- calculate_mr_stats(res_b[[2]]$weighted_alpha, res_b[[2]]$weighted_se)
            results$b_theta_point <- sanitize_value(res_b[[2]]$weighted_alpha)
            results$b_theta_se <- sanitize_value(res_b[[2]]$weighted_se)
            results$b_z <- sanitize_value(res_b_estimtor$z)
            results$b_p_value <- sanitize_value(res_b_estimtor$p_value)
            results$b_duration <- sanitize_value(as.numeric(res_b_time_2 - res_b_time_1))
        },
        error = function(e) {
            results$b_duration <<- NA
        }
    )

    # c: IVW
    tryCatch(
        {
            mr_input_obj <- MendelianRandomization::mr_input(bx = as.vector(phase_two_data_analysis_exp$beta_hat[, 1]), bxse = as.vector(phase_two_data_analysis_exp$beta_hat_se), by = as.vector(phase_two_data_analysis_out$beta_hat[, 1]), byse = as.vector(phase_two_data_analysis_out$beta_hat_se))
            res_c_time_1 <- Sys.time()
            res_c <- MendelianRandomization::mr_ivw(mr_input_obj)
            res_c_time_2 <- Sys.time()
            res_c_estimtor <- calculate_mr_stats(res_c$Estimate, res_c$StdError)
            results$c_theta_point <- sanitize_value(res_c$Estimate)
            results$c_theta_se <- sanitize_value(res_c$StdError)
            results$c_z <- sanitize_value(res_c_estimtor$z)
            results$c_p_value <- sanitize_value(res_c$Pvalue)
            results$c_duration <- sanitize_value(as.numeric(res_c_time_2 - res_c_time_1))
        },
        error = function(e) {
            results$c_duration <<- NA
        }
    )

    # d: MR-Egger
    tryCatch(
        {
            res_d_time_1 <- Sys.time()
            res_d <- MendelianRandomization::mr_egger(mr_input_obj)
            res_d_time_2 <- Sys.time()
            res_d_estimtor <- calculate_mr_stats(res_d$Estimate, res_d$StdError.Est)
            results$d_theta_point <- sanitize_value(res_d$Estimate)
            results$d_theta_se <- sanitize_value(res_d$StdError.Est)
            results$d_z <- sanitize_value(res_d_estimtor$z)
            results$d_p_value <- sanitize_value(res_d$Causal.pval)
            results$d_duration <- sanitize_value(as.numeric(res_d_time_2 - res_d_time_1))
        },
        error = function(e) {
            results$d_duration <<- NA
        }
    )

    # e: MR-Lasso
    tryCatch(
        {
            res_e_time_1 <- Sys.time()
            res_e <- MendelianRandomization::mr_lasso(mr_input_obj)
            res_e_time_2 <- Sys.time()
            res_e_estimtor <- calculate_mr_stats(res_e$Estimate, res_e$StdError)
            results$e_theta_point <- sanitize_value(res_e$Estimate)
            results$e_theta_se <- sanitize_value(res_e$StdError)
            results$e_z <- sanitize_value(res_e_estimtor$z)
            results$e_p_value <- sanitize_value(res_e$Pvalue)
            results$e_duration <- sanitize_value(as.numeric(res_e_time_2 - res_e_time_1))
        },
        error = function(e) {
            results$e_duration <<- NA
        }
    )

    # -------- 基础LMM为基础 --------

    # f: cml_overlap
    tryCatch(
        {
            res_f_time_1 <- Sys.time()
            res_f <- mr_cML_Overlap(phase_two_data_analysis_exp_gls$beta_hat, phase_two_data_analysis_out_gls$beta_hat, phase_two_data_analysis_exp_gls$beta_hat_se, phase_two_data_analysis_out_gls$beta_hat_se, n = n, rho = 0)
            res_f_time_2 <- Sys.time()
            res_f_estimtor <- calculate_mr_stats(res_f$MA_BIC_theta, res_f$MA_BIC_se)
            results$f_theta_point <- sanitize_value(res_f$MA_BIC_theta)
            results$f_theta_se <- sanitize_value(res_f$MA_BIC_se)
            results$f_z <- sanitize_value(res_f_estimtor$z)
            results$f_p_value <- sanitize_value(res_f$MA_BIC_p)
            results$f_duration <- sanitize_value(as.numeric(res_f_time_2 - res_f_time_1))
        },
        error = function(e) {
            results$f_duration <<- NA
        }
    )

    # g: IVW
    tryCatch(
        {
            mr_input_obj_gls <- MendelianRandomization::mr_input(bx = as.vector(phase_two_data_analysis_exp_gls$beta_hat), bxse = as.vector(phase_two_data_analysis_exp_gls$beta_hat_se), by = as.vector(phase_two_data_analysis_out_gls$beta_hat), byse = as.vector(phase_two_data_analysis_out_gls$beta_hat_se))
            res_g_time_1 <- Sys.time()
            res_g <- MendelianRandomization::mr_ivw(mr_input_obj_gls)
            res_g_time_2 <- Sys.time()
            res_g_estimtor <- calculate_mr_stats(res_g$Estimate, res_g$StdError)
            results$g_theta_point <- sanitize_value(res_g$Estimate)
            results$g_theta_se <- sanitize_value(res_g$StdError)
            results$g_z <- sanitize_value(res_g_estimtor$z)
            results$g_p_value <- sanitize_value(res_g$Pvalue)
            results$g_duration <- sanitize_value(as.numeric(res_g_time_2 - res_g_time_1))
        },
        error = function(e) {
            results$g_duration <<- NA
        }
    )

    # h: MR-Egger
    tryCatch(
        {
            res_h_time_1 <- Sys.time()
            res_h <- MendelianRandomization::mr_egger(mr_input_obj_gls)
            res_h_time_2 <- Sys.time()
            res_h_estimtor <- calculate_mr_stats(res_h$Estimate, res_h$StdError.Est)
            results$h_theta_point <- sanitize_value(res_h$Estimate)
            results$h_theta_se <- sanitize_value(res_h$StdError.Est)
            results$h_z <- sanitize_value(res_h_estimtor$z)
            results$h_p_value <- sanitize_value(res_h$Causal.pval)
            results$h_duration <- sanitize_value(as.numeric(res_h_time_2 - res_h_time_1))
        },
        error = function(e) {
            results$h_duration <<- NA
        }
    )

    # i: MR-Lasso
    tryCatch(
        {
            res_i_time_1 <- Sys.time()
            res_i <- MendelianRandomization::mr_lasso(mr_input_obj_gls)
            res_i_time_2 <- Sys.time()
            res_i_estimtor <- calculate_mr_stats(res_i$Estimate, res_i$StdError)
            results$i_theta_point <- sanitize_value(res_i$Estimate)
            results$i_theta_se <- sanitize_value(res_i$StdError)
            results$i_z <- sanitize_value(res_i_estimtor$z)
            results$i_p_value <- sanitize_value(res_i$Pvalue)
            results$i_duration <- sanitize_value(as.numeric(res_i_time_2 - res_i_time_1))
        },
        error = function(e) {
            results$i_duration <<- NA
        }
    )

    return(results)
}

# %% 模拟
if (FALSE) {
    source("统计模拟/MR 模拟（单次）.R")
    a <- triplet_family_simulation_once_robust(
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
        assortative_mating_strength = 1000
    )
    a
    phase_two_data_full <- generate_mr_trio_data_matrix_ultra(
        n = 1000, n_snps = 10, n_pleiotropic = 0, n_expose_heterogeneity = 0,
        p_f = 0.3, p_m = 0.3, beta_fs_to_oe_exp = 0.1, beta_ms_to_oe_exp = 0.1,
        beta_os_to_oe_exp = 0.3,
        h_beta_fs_to_oe_exp = 0.2,
        h_beta_ms_to_oe_exp = 0.1,
        h_beta_os_to_oe_exp = 0.3,
        mean_beta_fs_to_oe_out = 0.1,
        sd_beta_fs_to_oe_out = 0.05,
        mean_beta_ms_to_oe_out = 0.1,
        sd_beta_ms_to_oe_out = 0.05,
        mean_beta_os_to_oe_out = 0.1,
        sd_beta_os_to_oe_out = 0.05,
        p_neg_pleiotropy = 0.5,
        assortative_mating_strength = 0,
        crowd_differences = 0,
        beta_exp_to_out = 0, confounding_exp = 0.2,
        confounding_out = 0.2,
        correlation = 0.2, seed = NULL
    )
    # 常规纳入分析的数据（FGWAS）
    data_exp <- phase_two_data_full$data_exp
    data_out <- phase_two_data_full$data_out

    phase_two_data_analysis_exp <- fgwas_for_data_matrix(data_exp,
        predicted_outcome = "expose"
    )
    phase_two_data_analysis_out <- fgwas_for_data_matrix(data_out,
        predicted_outcome = "outcome"
    )
    # 基于普通的lmm的数据
    phase_two_data_analysis_exp_gls <- fgwas_for_data_matrix_single(
        data_exp,
        predicted_outcome = "expose"
    )
    phase_two_data_analysis_out_gls <- fgwas_for_data_matrix_single(
        data_out,
        predicted_outcome = "outcome"
    )

    res_a <- mr_cML_Overlap(
        phase_two_data_analysis_exp$beta_hat[, 1],
        phase_two_data_analysis_out$beta_hat[, 1],
        phase_two_data_analysis_exp$beta_hat_se,
        phase_two_data_analysis_out$beta_hat_se,
        n = n, # 确认样本量参数
        rho = 0
    )

    res_b <- cml_family_ultral_b(
        beta_hat_exp = phase_two_data_analysis_exp$beta_hat,
        beta_hat_out = phase_two_data_analysis_out$beta_hat,
        beta_sigma_exp = phase_two_data_analysis_exp$Sigma_inv,
        beta_sigma_out = phase_two_data_analysis_out$Sigma_inv,
        n = n
    )
    res_b[[2]]$weighted_se

    arrange(res_b[[1]], bic)
    res_c <- cml_family_ultral(
        phase_two_data_analysis_exp$beta_hat,
        phase_two_data_analysis_out$beta_hat,
        phase_two_data_analysis_exp$Sigma_inv,
        phase_two_data_analysis_out$Sigma_inv,
        n = n
    )
    res_c[[2]]$weighted_se
}

if (FALSE) {
    a <- triplet_family_simulation_once_robust(
        n = 10, # 要生成的总数据集数量
        num_pleiotropic = 0,
        n_expose_heterogeneity = 0,
        N_exp = 1000, N_out = 1000,
        p_f = 0.3, p_m = 0.3,
        # --- 暴露效应 (当非零时的大小) ---
        beta_FStoOE_exp = 0, beta_MStoOE_exp = 0,
        beta_OStoOE_exp = 0.3,
        # --- 暴露异质性 ---
        h_beta_FStoOE_exp = 0.5, h_beta_MStoOE_exp = 0.5,
        h_beta_OStoOE_exp = 0.3,
        # --- 结局效应 / 水平多效性 (当非零时的分布参数) ---
        # 子代基因型 -> 结局 (遗传叠加效应或基因多效性)
        mean_beta_FStoOE_out = 0.1, sd_beta_FStoOE_out = 0.05,
        mean_beta_MStoOE_out = 0.1, sd_beta_MStoOE_out = 0.05,
        # 父母基因型 -> 结局 (应为子代自身SNP对结局的效应，即基因多效性; 参数名可能需对应调整)
        mean_beta_OStoOE_out = 0, sd_beta_OStoOE_out = 0,
        prop_negative_pleiotropy = 0, # 在指定为多效性的SNP中，其效应为负的比例 (0 到 1)
        # 选型婚配
        assortative_mating_prob = 0,
        assortative_mating_strength = 1000, # 选型婚配对结局的影响因子
        # 人群分层
        ## 定义人群分层的差异(次等位基因频率差异)
        crowd_stratification_differences = 0, # 用于模拟两个具有不同等位基因频率的亚群
        # --- 其他效应 ---
        beta_exp_to_out = 0, # 暴露对结局的真实因果效应
        beta_confounding_exp = 0.2, # 影响暴露的混杂因素的方差 (效应大小为1)
        beta_confounding_out = 0.2, # 影响结局的混杂因素的方差 (效应大小为1)
        correlation = 0.2, # 共享环境因素的方差
        seed = NULL
    )

    a$b_p_value
    a
}
