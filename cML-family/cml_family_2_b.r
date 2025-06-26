# %% 载入要用到的包
library(microbenchmark)
# %% 定义函数
# beta_sigma_exp 是指逆矩阵
run_iterative_mle_b <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0) {
    num_snps <- nrow(beta_hat_exp)

    # --- 1. 数据预处理 ---
    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)

    for (m in 1:num_snps) {
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]
        Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
        Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
    }

    # --- 2. 初始化参数 ---
    alpha_current <- alpha_init
    gamma_current <- gamma_hat

    converged <- FALSE
    final_sets <- list()

    # ✨ 新增：初始化向量来追踪 alpha 和 gamma 的历史值
    alpha_history <- rep(NA, max_iter)
    gamma_history <- array(NA, dim = c(2, num_snps, max_iter))

    # --- 3. 迭代主循环 ---
    for (iter in 1:max_iter) {
        alpha_prev <- alpha_current
        gamma_prev <- gamma_current # ✨ 新增：保存前一次的gamma值

        # ==========================================================
        # 步骤 1: 选择 SNPs 集合 A 与 B
        # ==========================================================
        var_beta_o <- sapply(Sigma_beta_list, function(mat) mat[1, 1])
        var_beta_p <- sapply(Sigma_beta_list, function(mat) mat[2, 2])

        t_m <- (beta_hat[1, ] - alpha_current * gamma_current[1, ])^2 / var_beta_o
        f_m <- (beta_hat[2, ] - alpha_current * gamma_current[2, ])^2 / var_beta_p

        a <- min(a, num_snps)
        b <- min(b, num_snps)

        set_A_indices <- order(t_m)[1:a]
        set_B_indices <- order(f_m)[1:b]

        A_intersect_B <- intersect(set_A_indices, set_B_indices)
        A_diff_B <- setdiff(set_A_indices, set_B_indices)
        B_diff_A <- setdiff(set_B_indices, set_A_indices)
        others <- c(setdiff(
            1:num_snps,
            union(set_A_indices, set_B_indices)
        ), B_diff_A)

        # ==========================================================
        # 步骤 2: 更新 gamma
        # ==========================================================
        gamma_new <- matrix(0, nrow = 2, ncol = num_snps)

        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                term1_inv <- Sigma_gamma_inv_list[[m]]
                term2_inv <- alpha_current^2 * Sigma_beta_inv_list[[m]]
                combined_inv <- term1_inv + term2_inv
                term1_val <- term1_inv %*% gamma_hat[, m]
                term2_val <- alpha_current * Sigma_beta_inv_list[[m]] %*%
                    beta_hat[, m]
                gamma_new[, m] <- solve(combined_inv) %*% (term1_val + term2_val)
            }
        }

        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_go <- gamma_hat[1, m] * var_beta_o[m] - alpha_current * Sigma_gamma_list[[m]][1, 1] * beta_hat[1, m]
                den_go <- var_beta_o[m] - alpha_current^2 * Sigma_gamma_list[[m]][1, 1]
                gamma_o_new <- num_go / den_go
                term_gp <- Sigma_gamma_list[[m]][1, 2] * (-alpha_current / var_beta_o[m]) * (beta_hat[1, m] - alpha_current * gamma_o_new)
                gamma_p_new <- gamma_hat[2, m] + term_gp
                gamma_new[, m] <- c(gamma_o_new, gamma_p_new)
            }
        }

        if (length(others) > 0) {
            gamma_new[, others] <- gamma_hat[, others]
        }

        gamma_current <- gamma_new

        # ==========================================================
        # 步骤 3: 更新 alpha
        # ==========================================================
        num_alpha <- 0
        den_alpha <- 0

        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                num_alpha <- num_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% beta_hat[, m]
                den_alpha <- den_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% gamma_current[, m]
            }
        }
        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_alpha <- num_alpha + (gamma_current[1, m] * beta_hat[1, m] / var_beta_o[m])
                den_alpha <- den_alpha + (gamma_current[1, m]^2 / var_beta_o[m])
            }
        }

        if (abs(den_alpha) > 1e-9) {
            alpha_current <- as.numeric(num_alpha / den_alpha)
        } else {
            warning(paste("Iteration", iter, ": Denominator for alpha is near zero. Halting."))
            break
        }

        # ✨ 新增：在每次迭代后记录 alpha 和 gamma 的值
        alpha_history[iter] <- alpha_current
        gamma_history[, , iter] <- gamma_current

        # --- 4. 检查收敛（同时检验alpha和gamma） ---
        # ✨ 修改：计算alpha和gamma的变化量
        alpha_change <- abs(alpha_current - alpha_prev)
        gamma_change <- max(abs(gamma_current - gamma_prev)) # 使用最大绝对变化
        # 也可以使用Frobenius范数：gamma_change <- norm(gamma_current - gamma_prev, "F")

        # if (alpha_change < tol && gamma_change < tol) {
        #     converged <- TRUE
        #     final_sets <- list(A = set_A_indices, B = set_B_indices)
        #     message(paste("✅ Algorithm converged after", iter, "iterations."))
        #     message(paste("   Final alpha change:", round(alpha_change, 8)))
        #     message(paste("   Final gamma change:", round(gamma_change, 8)))
        #     break
        # }
    }

    # if (!converged) {
    #     warning(paste("⚠️ Algorithm did not converge after", max_iter, "iterations."))
    #     warning(paste("   Final alpha change:", round(alpha_change, 8)))
    #     warning(paste("   Final gamma change:", round(gamma_change, 8)))
    #     final_sets <- list(A = set_A_indices, B = set_B_indices)
    # }

    # ✨ 新增：清理历史记录，移除未使用的 NA 值
    alpha_history <- alpha_history[!is.na(alpha_history)]
    gamma_history <- gamma_history[, , !is.na(alpha_history), drop = FALSE]

    return(list(
        alpha_final = as.numeric(alpha_current),
        gamma_final = t(gamma_current),
        iterations = iter,
        converged = converged,
        sets = final_sets,
        alpha_history = alpha_history,
        gamma_history = gamma_history, # ✨ 新增：返回gamma的历史记录
        # ✨ 将所有预处理的列表一并返回
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        )
    ))
}
# b可以为零完全退化
run_iterative_mle_b_zero <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0) {
    num_snps <- nrow(beta_hat_exp)

    # --- 1. 数据预处理 ---
    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)

    for (m in 1:num_snps) {
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]
        Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
        Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
    }

    # --- 2. 初始化参数 ---
    alpha_current <- alpha_init
    gamma_current <- gamma_hat

    converged <- FALSE
    final_sets <- list()

    # ✨ 新增：初始化向量来追踪 alpha 和 gamma 的历史值
    alpha_history <- rep(NA, max_iter)
    gamma_history <- array(NA, dim = c(2, num_snps, max_iter))

    # --- 3. 迭代主循环 ---
    for (iter in 1:max_iter) {
        alpha_prev <- alpha_current
        gamma_prev <- gamma_current # ✨ 新增：保存前一次的gamma值

        # ==========================================================
        # 步骤 1: 选择 SNPs 集合 A 与 B
        # ==========================================================
        var_beta_o <- sapply(Sigma_beta_list, function(mat) mat[1, 1])
        var_beta_p <- sapply(Sigma_beta_list, function(mat) mat[2, 2])

        t_m <- (beta_hat[1, ] - alpha_current * gamma_current[1, ])^2 / var_beta_o

        # 🔧 修改：调整集合大小上限
        a <- min(a, num_snps)
        b <- min(b, num_snps) # b 可以为 0

        # 选择集合 A
        set_A_indices <- if (a > 0) order(t_m)[1:a] else integer(0)

        # 🔧 修改：处理 b = 0 的情况
        if (b > 0) {
            f_m <- (beta_hat[2, ] - alpha_current * gamma_current[2, ])^2 / var_beta_p
            set_B_indices <- order(f_m)[1:b]
        } else {
            set_B_indices <- integer(0) # 空集合
        }

        # 🔧 修改：重新计算集合关系，考虑 B 可能为空的情况
        if (length(set_B_indices) > 0) {
            A_intersect_B <- intersect(set_A_indices, set_B_indices)
            A_diff_B <- setdiff(set_A_indices, set_B_indices)
            B_diff_A <- setdiff(set_B_indices, set_A_indices)
            others <- c(setdiff(
                1:num_snps,
                union(set_A_indices, set_B_indices)
            ), B_diff_A)
        } else {
            # 当 B 为空时
            A_intersect_B <- integer(0)
            A_diff_B <- set_A_indices # A 的所有元素都在 A \ B 中
            B_diff_A <- integer(0)
            others <- setdiff(1:num_snps, set_A_indices)
        }

        # ==========================================================
        # 步骤 2: 更新 gamma
        # ==========================================================
        gamma_new <- matrix(0, nrow = 2, ncol = num_snps)

        # 🔧 修改：A ∩ B 的处理（当 b=0 时，此集合为空）
        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                term1_inv <- Sigma_gamma_inv_list[[m]]
                term2_inv <- alpha_current^2 * Sigma_beta_inv_list[[m]]
                combined_inv <- term1_inv + term2_inv
                term1_val <- term1_inv %*% gamma_hat[, m]
                term2_val <- alpha_current * Sigma_beta_inv_list[[m]] %*%
                    beta_hat[, m]
                gamma_new[, m] <- solve(combined_inv) %*% (term1_val + term2_val)
            }
        }

        # 🔧 修改：A \ B 的处理（当 b=0 时，这就是整个集合 A）
        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_go <- gamma_hat[1, m] * var_beta_o[m] - alpha_current * Sigma_gamma_list[[m]][1, 1] * beta_hat[1, m]
                den_go <- var_beta_o[m] - alpha_current^2 * Sigma_gamma_list[[m]][1, 1]

                # 🔧 新增：检查分母是否接近零
                if (abs(den_go) > 1e-9) {
                    gamma_o_new <- num_go / den_go
                } else {
                    warning(paste("Iteration", iter, "SNP", m, ": Denominator for gamma_o is near zero. Using previous value."))
                    gamma_o_new <- gamma_current[1, m]
                }

                term_gp <- Sigma_gamma_list[[m]][1, 2] * (-alpha_current / var_beta_o[m]) * (beta_hat[1, m] - alpha_current * gamma_o_new)
                gamma_p_new <- gamma_hat[2, m] + term_gp
                gamma_new[, m] <- c(gamma_o_new, gamma_p_new)
            }
        }

        # 其他 SNPs 保持原值
        if (length(others) > 0) {
            gamma_new[, others] <- gamma_hat[, others]
        }

        gamma_current <- gamma_new

        # ==========================================================
        # 步骤 3: 更新 alpha
        # ==========================================================
        num_alpha <- 0
        den_alpha <- 0

        # 🔧 修改：A ∩ B 的贡献（当 b=0 时为空，不会执行）
        if (length(A_intersect_B) > 0) {
            for (m in A_intersect_B) {
                num_alpha <- num_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% beta_hat[, m]
                den_alpha <- den_alpha + t(gamma_current[, m]) %*% Sigma_beta_inv_list[[m]] %*% gamma_current[, m]
            }
        }

        # 🔧 修改：A \ B 的贡献（当 b=0 时，这就是整个集合 A 的贡献）
        if (length(A_diff_B) > 0) {
            for (m in A_diff_B) {
                num_alpha <- num_alpha + (gamma_current[1, m] * beta_hat[1, m] / var_beta_o[m])
                den_alpha <- den_alpha + (gamma_current[1, m]^2 / var_beta_o[m])
            }
        }

        # 🔧 新增：更robust的alpha更新
        if (abs(den_alpha) > 1e-9) {
            alpha_current <- as.numeric(num_alpha / den_alpha)
        } else {
            warning(paste("Iteration", iter, ": Denominator for alpha is near zero. Using previous value."))
            alpha_current <- alpha_prev
        }

        # ✨ 新增：在每次迭代后记录 alpha 和 gamma 的值
        alpha_history[iter] <- alpha_current
        gamma_history[, , iter] <- gamma_current

        # --- 4. 检查收敛（同时检验alpha和gamma） ---
        # ✨ 修改：计算alpha和gamma的变化量
        alpha_change <- abs(alpha_current - alpha_prev)
        gamma_change <- max(abs(gamma_current - gamma_prev)) # 使用最大绝对变化

        # 🔧 新增：可选的收敛检查（您可以取消注释来启用）
        if (alpha_change < tol && gamma_change < tol) {
            converged <- TRUE
            final_sets <- list(A = set_A_indices, B = set_B_indices)
            message(paste("✅ Algorithm converged after", iter, "iterations."))
            message(paste("   Final alpha change:", round(alpha_change, 8)))
            message(paste("   Final gamma change:", round(gamma_change, 8)))
            message(paste("   |A| =", length(set_A_indices), ", |B| =", length(set_B_indices)))
            break
        }
    }

    if (!converged) {
        warning(paste("⚠️ Algorithm did not converge after", max_iter, "iterations."))
        warning(paste("   Final alpha change:", round(alpha_change, 8)))
        warning(paste("   Final gamma change:", round(gamma_change, 8)))
        final_sets <- list(A = set_A_indices, B = set_B_indices)
    }

    # ✨ 新增：清理历史记录，移除未使用的 NA 值
    alpha_history <- alpha_history[!is.na(alpha_history)]
    gamma_history <- gamma_history[, , !is.na(alpha_history), drop = FALSE]

    # 🔧 新增：输出集合大小信息
    message(paste("📊 Final sets: |A| =", length(set_A_indices), ", |B| =", length(set_B_indices)))
    message(paste("📊 A ∩ B:", length(A_intersect_B), ", A \\ B:", length(A_diff_B), ", B \\ A:", length(B_diff_A)))

    return(list(
        alpha_final = as.numeric(alpha_current),
        gamma_final = t(gamma_current),
        iterations = iter,
        converged = converged,
        sets = final_sets,
        alpha_history = alpha_history,
        gamma_history = gamma_history, # ✨ 新增：返回gamma的历史记录
        # ✨ 将所有预处理的列表一并返回
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        )
    ))
}
# 计算方差
calculate_alpha_variance_b <- function(estimation_result) {
    # --- 1. 从结果对象中优雅地提取所有所需变量 ---
    alpha <- estimation_result$alpha_final
    gamma <- t(estimation_result$gamma_final)
    sets <- estimation_result$sets

    # 从嵌套的 data 列表中提取
    beta_hat <- t(estimation_result$data$beta_hat_out)
    Sigma_gamma_list <- estimation_result$data$Sigma_gamma_list
    Sigma_beta_list <- estimation_result$data$Sigma_beta_list
    Sigma_gamma_inv_list <- estimation_result$data$Sigma_gamma_inv_list
    Sigma_beta_inv_list <- estimation_result$data$Sigma_beta_inv_list

    all_selected_snps <- sort(sets$A)
    K <- length(all_selected_snps)

    if (K == 0) {
        warning("没有被选中的SNPs，无法计算方差。")
        return(list(variance = NA, std_error = NA, hessian_eigenvalues = NA))
    }

    total_dim <- 1 + 2 * K
    H <- matrix(0, nrow = total_dim, ncol = total_dim)

    # --- 2. 填充矩阵 H ---
    H_alpha_alpha_total <- 0

    for (i in 1:K) {
        m <- all_selected_snps[i]
        row_start <- 2 + 2 * (i - 1)
        solve(Sigma_beta_list[[m]])
        # H_αα 和 H_γm,α 的计算不变
        if (m %in% intersect(sets$A, sets$B)) {
            H_alpha_alpha_total <- H_alpha_alpha_total +
                t(gamma[, m]) %*% solve(Sigma_beta_list[[m]]) %*% gamma[, m]
            H_gmgm <- alpha^2 * Sigma_beta_inv_list[[m]] +
                Sigma_gamma_inv_list[[m]]
            H_gma <- -Sigma_beta_inv_list[[m]] %*%
                (beta_hat[, m] - 2 * alpha * gamma[, m])
        } else if (m %in% setdiff(sets$A, sets$B)) {
            var_beta_o_m <- Sigma_beta_list[[m]][1, 1]
            var_gamma_o_m <- Sigma_gamma_list[[m]][1, 1] # 获取gamma_o的方差

            H_alpha_alpha_total <- H_alpha_alpha_total +
                gamma[1, m]^2 / var_beta_o_m
            # ✨ 核心修正：根据最终似然函数更新 H_gmgm
            H_gmgm <- matrix(c(alpha^2 / var_beta_o_m + 1 / var_gamma_o_m, 0, 0, 0), nrow = 2)
            H_gma <- c((-beta_hat[1, m] + 2 * alpha * gamma[1, m]) / var_beta_o_m, 0)
        }

        # 填充矩阵
        H[
            row_start:(row_start + 1),
            row_start:(row_start + 1)
        ] <- H_gmgm
        H[row_start:(row_start + 1), 1] <- H_gma
        H[1, row_start:(row_start + 1)] <- t(H_gma)
    }

    H[1, 1] <- H_alpha_alpha_total
    H <- H[rowSums(H != 0) > 0, colSums(H != 0) > 0]
    # --- 3. 求逆并提取结果（带正则化） ---
    eigenvals <- eigen(H, only.values = TRUE)$values
    min_eigenval <- min(eigenvals)
    regularization_applied <- FALSE

    # 检查特征值并应用正则化
    if (any(eigenvals <= 1e-8)) { # 使用更严格的阈值
        warning("Hessian 矩阵非正定或接近奇异。正在应用正则化。")

        # 计算需要添加的正则化项
        reg_value <- max(1e-6, abs(min_eigenval) + 1e-6)

        # 应用正则化：添加到对角线
        H_regularized <- H + reg_value * diag(nrow(H))

        # 验证正则化后的矩阵
        eigenvals_reg <- eigen(H_regularized, only.values = TRUE)$values

        if (any(eigenvals_reg <= 0)) {
            # 如果还是有问题，使用更强的正则化
            reg_value <- max(1e-4, abs(min(eigenvals_reg)) + 1e-4)
            H_regularized <- H + reg_value * diag(nrow(H))
            warning(paste("应用强正则化，正则化值:", reg_value))
        }

        H <- H_regularized
    }
    H_inv <- tryCatch(solve(H), error = function(e) NULL)

    if (is.null(H_inv)) {
        return(list(variance = NA, std_error = NA, hessian_eigenvalues = eigenvals))
    }

    var_alpha <- H_inv[1, 1]

    if (var_alpha < 0) {
        warning("计算出的 Alpha 方差为负数。")
        return(list(variance = var_alpha, std_error = NA, hessian_eigenvalues = eigenvals))
    }

    se_alpha <- sqrt(var_alpha)
    return(list(variance = var_alpha, std_error = se_alpha, hessian = H))
}

# 首先编译Rcpp函数
Rcpp::sourceCpp("cML-family/cml_family_2_b_rcpp.cpp")

# rcpp优化版本
run_iterative_mle_b_optimized <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b,
    max_iter = 100, tol = 1e-6, alpha_init = 0, Sigma_beta_list,
    Sigma_gamma_list,
    Sigma_beta_inv_list,
    Sigma_gamma_inv_list) {
    # 检查是否已编译Rcpp函数
    if (!exists("compute_statistics_cpp")) {
        stop("请先编译Rcpp函数: Rcpp::sourceCpp('mle_optimization.cpp')")
    }

    num_snps <- nrow(beta_hat_exp)

    # --- 1. 数据预处理 ---
    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)


    # --- 2. 初始化参数 ---
    alpha_current <- alpha_init
    gamma_current <- gamma_hat

    converged <- FALSE
    final_sets <- list()

    # 初始化历史记录
    alpha_history <- rep(NA, max_iter)
    gamma_history <- array(NA, dim = c(2, num_snps, max_iter))

    # --- 3. 迭代主循环 ---
    for (iter in 1:max_iter) {
        alpha_prev <- alpha_current
        gamma_prev <- gamma_current

        # ==========================================================
        # 步骤 1: 使用Rcpp计算统计量和选择集合
        # ==========================================================
        stats_result <- compute_statistics_cpp(
            beta_hat, gamma_current,
            alpha_current, Sigma_beta_list
        )
        t_m <- stats_result$t_m
        f_m <- stats_result$f_m
        var_beta_o <- stats_result$var_beta_o
        var_beta_p <- stats_result$var_beta_p

        a <- min(a, num_snps)
        b <- min(b, num_snps)

        set_A_indices <- order(t_m)[1:a]
        set_B_indices <- order(f_m)[1:b]

        A_intersect_B <- intersect(set_A_indices, set_B_indices)
        A_diff_B <- setdiff(set_A_indices, set_B_indices)
        B_diff_A <- setdiff(set_B_indices, set_A_indices)
        others <- c(setdiff(1:num_snps, union(set_A_indices, set_B_indices)), B_diff_A)

        # ==========================================================
        # 步骤 2: 使用Rcpp更新gamma
        # ==========================================================
        gamma_new <- matrix(0, nrow = 2, ncol = num_snps)

        # 更新A∩B集合
        if (length(A_intersect_B) > 0) {
            gamma_new <- update_gamma_intersect_cpp(
                A_intersect_B, gamma_hat, beta_hat, alpha_current,
                Sigma_gamma_inv_list, Sigma_beta_inv_list, gamma_new
            )
        }

        # 更新A-B集合
        if (length(A_diff_B) > 0) {
            gamma_new <- update_gamma_diff_cpp(
                A_diff_B, gamma_hat, beta_hat, alpha_current,
                Sigma_gamma_list, var_beta_o, gamma_new
            )
        }

        # 更新其他集合
        if (length(others) > 0) {
            gamma_new[, others] <- gamma_hat[, others]
        }

        gamma_current <- gamma_new

        # ==========================================================
        # 步骤 3: 使用Rcpp更新alpha
        # ==========================================================
        alpha_terms <- update_alpha_terms_cpp(
            A_intersect_B, A_diff_B, gamma_current, beta_hat,
            Sigma_beta_inv_list, var_beta_o
        )

        if (abs(alpha_terms$den_alpha) > 1e-9) {
            alpha_current <- alpha_terms$num_alpha / alpha_terms$den_alpha
        } else {
            warning(paste("Iteration", iter, ": Denominator for alpha is near zero. Halting."))
            break
        }

        # 记录历史值
        alpha_history[iter] <- alpha_current
        gamma_history[, , iter] <- gamma_current

        # --- 4. 检查收敛（使用Rcpp计算变化量） ---
        alpha_change <- abs(alpha_current - alpha_prev)
        gamma_change <- compute_max_change_cpp(gamma_current, gamma_prev)

        if (alpha_change < tol && gamma_change < tol) {
            converged <- TRUE
            final_sets <- list(A = set_A_indices, B = set_B_indices)
            message(paste("✅ Algorithm converged after", iter, "iterations."))
            message(paste("   Final alpha change:", round(alpha_change, 8)))
            message(paste("   Final gamma change:", round(gamma_change, 8)))
            break
        }
    }

    if (!converged) {
        warning(paste("⚠️ Algorithm did not converge after", max_iter, "iterations."))
        if (exists("alpha_change") && exists("gamma_change")) {
            warning(paste("   Final alpha change:", round(alpha_change, 8)))
            warning(paste("   Final gamma change:", round(gamma_change, 8)))
        }
        final_sets <- list(A = set_A_indices, B = set_B_indices)
    }

    # 清理历史记录
    alpha_history <- alpha_history[!is.na(alpha_history)]
    gamma_history <- gamma_history[, , !is.na(alpha_history), drop = FALSE]


    return(list(
        alpha_final = as.numeric(alpha_current),
        gamma_final = t(gamma_current),
        iterations = iter,
        converged = converged,
        sets = final_sets,
        alpha_history = alpha_history,
        gamma_history = gamma_history,
        data = list(
            beta_hat_exp = beta_hat_exp,
            beta_hat_out = beta_hat_out,
            Sigma_gamma_list = Sigma_gamma_list,
            Sigma_beta_list = Sigma_beta_list,
            Sigma_gamma_inv_list = Sigma_gamma_inv_list,
            Sigma_beta_inv_list = Sigma_beta_inv_list
        )
    ))
}
run_multi_start_mle_b_optimized <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, a, b, n,
    n_starts = 10,
    alpha_range = c(-1, 1),
    max_iter = 100,
    tol = 1e-6,
    seed = NULL,
    verbose = FALSE) {
    # 设置随机种子
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # 检查是否已编译Rcpp函数
    if (!exists("compute_statistics_cpp")) {
        stop("请先编译Rcpp函数: Rcpp::sourceCpp('mle_optimization.cpp')")
    }

    # 获取SNP数量和索引
    num_snps <- nrow(beta_hat_exp)
    all_snps <- 1:num_snps

    num_snps <- nrow(beta_hat_exp)

    # --- 1. 数据预处理 ---
    gamma_hat <- t(beta_hat_exp)
    beta_hat <- t(beta_hat_out)

    Sigma_gamma_inv_list <- vector("list", num_snps)
    Sigma_beta_inv_list <- vector("list", num_snps)
    Sigma_gamma_list <- vector("list", num_snps)
    Sigma_beta_list <- vector("list", num_snps)

    for (m in 1:num_snps) {
        Sigma_gamma_inv_list[[m]] <- beta_sigma_exp[(2 * m - 1):(2 * m), ]
        Sigma_beta_inv_list[[m]] <- beta_sigma_out[(2 * m - 1):(2 * m), ]
        Sigma_gamma_list[[m]] <- solve(Sigma_gamma_inv_list[[m]])
        Sigma_beta_list[[m]] <- solve(Sigma_beta_inv_list[[m]])
    }

    # 生成多个初始值
    alpha_inits <- if (n_starts == 1) {
        0 # 单一起点使用0
    } else {
        c(0, runif(n_starts - 1, min = alpha_range[1], max = alpha_range[2]))
    }

    # 存储所有结果
    all_results <- vector("list", n_starts)
    all_bic_values <- numeric(n_starts)
    convergence_status <- logical(n_starts)

    if (verbose) {
        cat("开始多起点优化...\n")
        cat(sprintf("起点数量: %d\n", n_starts))
        cat(sprintf("Alpha初始值范围: [%.2f, %.2f]\n", alpha_range[1], alpha_range[2]))
        cat(paste(rep("=", 50), collapse = ""), "\n")
    }

    # 对每个起点进行优化
    for (start_idx in 1:n_starts) {
        if (verbose) {
            cat(sprintf(
                "起点 %d/%d: alpha_init = %.4f\n",
                start_idx, n_starts, alpha_inits[start_idx]
            ))
        }



        tryCatch(
            {
                # 运行单次优化
                result <- run_iterative_mle_b_optimized(
                    beta_hat_exp = beta_hat_exp,
                    beta_hat_out = beta_hat_out,
                    beta_sigma_exp = beta_sigma_exp,
                    beta_sigma_out = beta_sigma_out,
                    a = a,
                    b = b,
                    max_iter = max_iter,
                    tol = tol,
                    alpha_init = alpha_inits[start_idx],
                    Sigma_beta_list,
                    Sigma_gamma_list,
                    Sigma_beta_inv_list,
                    Sigma_gamma_inv_list
                )
                result$data <- list(
                    beta_hat_exp = beta_hat_exp,
                    beta_hat_out = beta_hat_out,
                    Sigma_gamma_list = Sigma_gamma_list,
                    Sigma_beta_list = Sigma_beta_list,
                    Sigma_gamma_inv_list = Sigma_gamma_inv_list,
                    Sigma_beta_inv_list = Sigma_beta_inv_list
                )
                # 计算BIC值
                bic_value <- calculate_bic_cml_b(
                    input_model = result,
                    all_snps = all_snps,
                    n = n
                )

                # 存储结果
                all_results[[start_idx]] <- result
                all_bic_values[start_idx] <- bic_value
                convergence_status[start_idx] <- result$converged

                if (verbose) {
                    cat(sprintf(
                        "  -> 收敛: %s, 迭代次数: %d, BIC: %.4f, 最终alpha: %.6f\n",
                        ifelse(result$converged, "是", "否"),
                        result$iterations,
                        bic_value,
                        result$alpha_final
                    ))
                }
            },
            error = function(e) {
                if (verbose) {
                    cat(sprintf("  -> 错误: %s\n", e$message))
                }
                all_results[[start_idx]] <- NULL
                all_bic_values[start_idx] <- Inf
                convergence_status[start_idx] <- FALSE
            }
        )
    }

    # 找到有效结果
    valid_indices <- which(!is.infinite(all_bic_values) & !sapply(all_results, is.null))

    if (length(valid_indices) == 0) {
        stop("所有起点都失败了，请检查输入数据或调整参数")
    }

    # 选择BIC最小的结果
    best_idx <- valid_indices[which.min(all_bic_values[valid_indices])]
    best_result <- all_results[[best_idx]]
    best_bic <- all_bic_values[best_idx]

    if (verbose) {
        cat("\n", paste(rep("=", 50), collapse = ""), "\n")
        cat("优化完成!\n")
        cat(sprintf("有效结果数量: %d/%d\n", length(valid_indices), n_starts))
        cat(sprintf("收敛结果数量: %d/%d\n", sum(convergence_status[valid_indices]), length(valid_indices)))
        cat(sprintf("最佳起点: %d (alpha_init = %.4f)\n", best_idx, alpha_inits[best_idx]))
        cat(sprintf("最佳BIC: %.4f\n", best_bic))
        cat(sprintf("最佳alpha: %.6f\n", best_result$alpha_final))

        # 显示BIC排名前3的结果
        if (length(valid_indices) > 1) {
            sorted_indices <- valid_indices[order(all_bic_values[valid_indices])]
            cat("\nBIC排名前3的结果:\n")
            for (i in 1:min(3, length(sorted_indices))) {
                idx <- sorted_indices[i]
                cat(sprintf(
                    "  %d. 起点%d: BIC=%.4f, alpha=%.6f, 收敛=%s\n",
                    i, idx, all_bic_values[idx],
                    all_results[[idx]]$alpha_final,
                    ifelse(convergence_status[idx], "是", "否")
                ))
            }
        }
    }

    # 创建增强的返回结果
    enhanced_result <- best_result
    enhanced_result$multi_start_info <- list(
        n_starts = n_starts,
        best_start_index = best_idx,
        best_alpha_init = alpha_inits[best_idx],
        best_bic = best_bic,
        all_bic_values = all_bic_values,
        all_alpha_finals = sapply(all_results, function(x) if (is.null(x)) NA else x$alpha_final),
        convergence_status = convergence_status,
        valid_starts = length(valid_indices),
        alpha_inits = alpha_inits
    )
    enhanced_result$data <- list(
        beta_hat_exp = beta_hat_exp,
        beta_hat_out = beta_hat_out,
        Sigma_gamma_list = Sigma_gamma_list,
        Sigma_beta_list = Sigma_beta_list,
        Sigma_gamma_inv_list = Sigma_gamma_inv_list,
        Sigma_beta_inv_list = Sigma_beta_inv_list
    )
    return(enhanced_result)
}
calculate_bic_cml_b <- function(input_model, all_snps, n, penalty_factor = 1 / 4) {
    # 获取 SNP 的总数 m
    m <- length(all_snps)
    A_hat <- input_model$sets$A
    B_hat <- input_model$sets$B
    alpha <- input_model$alpha_final
    gamma_model <- input_model$gamma_final
    beta_exp_sigma_inv <- input_model$data$Sigma_gamma_inv_list
    beta_out_sigma_inv <- input_model$data$Sigma_beta_inv_list
    beta_exp_hat <- input_model$data$beta_hat_exp
    beta_out_hat <- input_model$data$beta_hat_out
    # 1. 确定各个 SNP 集合的索引
    # -------------------------------------
    snps_intersect <- intersect(A_hat, B_hat)
    snps_A_only <- setdiff(A_hat, B_hat)
    snps_other <- setdiff(all_snps, A_hat)

    # 将 SNP ID 转换为在矩阵/列表中的行索引
    idx_intersect <- match(snps_intersect, all_snps)
    idx_A_only <- match(snps_A_only, all_snps)
    idx_other <- match(snps_other, all_snps)

    # 初始化各部分计算结果
    term1 <- 0
    term2 <- 0
    term3 <- 0


    # 2. 计算第一项 (m ∈ Â ∩ B̂)
    # -------------------------------------
    if (length(idx_intersect) > 0) {
        term1_values <- sapply(idx_intersect, function(i) {
            # 提取对应 SNP 的数据
            beta_hat_i <- beta_out_hat[i, ]
            gamma_hat_i <- beta_exp_hat[i, ]
            gamma_mod_i <- gamma_model[i, ]
            Sigma_beta_inv_i <- beta_out_sigma_inv[[i]]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            # 计算差值向量
            diff_beta <- beta_hat_i - alpha * gamma_mod_i
            diff_gamma <- gamma_hat_i - gamma_mod_i

            # 计算二次型 (t(x) %*% M %*% x)
            val_beta <- t(diff_beta) %*% Sigma_beta_inv_i %*% diff_beta
            val_gamma <- t(diff_gamma) %*% Sigma_gamma_inv_i %*% diff_gamma

            return(val_beta + val_gamma)
        })
        term1 <- sum(term1_values)
    }

    # 3. 计算第二项 (m ∈ Â \ B̂)
    # -------------------------------------
    if (length(idx_A_only) > 0) {
        term2_values <- sapply(idx_A_only, function(i) {
            # 提取 'o' 分量 (第一列)
            beta_o_hat <- beta_out_hat[i, 1]
            gamma_o_hat <- beta_exp_hat[i, 1]
            gamma_o_mod <- gamma_model[i, 1]

            # 从逆协方差矩阵计算方差 σ²
            # σ²_o = (Σ⁻¹)_pp / det(Σ⁻¹)
            Sigma_beta_inv_i <- beta_out_sigma_inv[[i]]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            sigma2_beta_o <- Sigma_beta_inv_i[2, 2] / det(Sigma_beta_inv_i)
            sigma2_gamma_o <- Sigma_gamma_inv_i[2, 2] / det(Sigma_gamma_inv_i)

            # 计算该项的值
            val_beta <- (beta_o_hat - alpha * gamma_o_mod)^2 / sigma2_beta_o
            val_gamma <- (gamma_o_hat - gamma_o_mod)^2 / sigma2_gamma_o

            return(val_beta + val_gamma)
        })
        term2 <- sum(term2_values)
    }
    # 5. 计算第三项 (其他 m)
    # -------------------------------------
    if (length(idx_other) > 0) {
        term4_values <- sapply(idx_other, function(i) {
            gamma_hat_i <- beta_exp_hat[i, ]
            gamma_mod_i <- gamma_model[i, ]
            Sigma_gamma_inv_i <- beta_exp_sigma_inv[[i]]

            diff_gamma <- gamma_hat_i - gamma_mod_i
            val_gamma <- t(diff_gamma) %*% Sigma_gamma_inv_i %*% diff_gamma

            return(val_gamma)
        })
        term3 <- sum(term4_values)
    }

    # 6. 计算第四项 (惩罚项)
    # -------------------------------------

    double_penalty <- setdiff(all_snps, intersect(A_hat, B_hat))
    single_penalty <- setdiff(all_snps, A_hat)
    single_penalty <- setdiff(single_penalty, double_penalty)
    penalty <- penalty_factor * log(n) *
        (2 * length(double_penalty) + length(single_penalty))

    # 7. 加总所有项得到最终 BIC 值
    # -------------------------------------
    bic_total <- term1 + term2 + term3 + penalty

    return(bic_total)
}

cml_family_ultral_b <- function(
    beta_hat_exp, beta_hat_out, beta_sigma_exp,
    beta_sigma_out, n = 1000, cov_index = "normal",
    n_starts = 5, verbose = FALSE,
    alpha_range = c(-0.001, 0.001),
    max_iter = 100, tol = 1e-6, penalty_factor = 1 / 4) {
    all_snps <- 1:nrow(beta_hat_exp)
    data_table_full <- expand.grid(a = all_snps, b = all_snps)

    # 只保留 b <= a 的组合
    data_table <- data_table_full %>% filter(b <= a)

    # 初始化结果列
    data_table$alpha <- NA
    data_table$alpha_se <- NA
    data_table$bic <- NA
    data_table$converged <- NA
    data_table$iterations <- NA
    data_table$valid_starts <- NA
    data_table$error_message <- NA

    # 统计变量
    total_combinations <- nrow(data_table)
    successful_combinations <- 0
    failed_combinations <- 0
    failed_indices <- c()

    if (verbose) {
        cat("开始CML Family分析...\n")
        cat(sprintf("总组合数: %d\n", total_combinations))
        cat(sprintf("每个组合的起点数: %d\n", n_starts))
        cat(paste(rep("=", 60), collapse = ""), "\n")
    }

    # 存储最后一个成功的结果用于返回
    last_successful_result <- NULL

    for (i in 1:nrow(data_table)) {
        a_val <- data_table[i, 1]
        b_val <- data_table[i, 2]

        if (verbose) {
            cat(sprintf("组合 %d/%d: a=%d, b=%d", i, total_combinations, a_val, b_val))
        }

        # 尝试运行优化
        tryCatch(
            {
                # 运行多起点优化
                data_any <- run_multi_start_mle_b_optimized(
                    beta_hat_exp = beta_hat_exp,
                    beta_hat_out = beta_hat_out,
                    beta_sigma_exp = beta_sigma_exp,
                    beta_sigma_out = beta_sigma_out,
                    a = a_val,
                    b = b_val,
                    n = n,
                    n_starts = n_starts,
                    alpha_range = alpha_range,
                    max_iter = max_iter,
                    tol = tol,
                    verbose = FALSE # 关闭详细输出避免过多信息
                )

                # 检查是否有有效的结果
                if (is.null(data_any) || is.null(data_any$alpha_final)) {
                    stop("优化返回了NULL结果")
                }

                # 检查多起点信息
                if (is.null(data_any$multi_start_info) ||
                    data_any$multi_start_info$valid_starts == 0) {
                    stop("所有起点都失败了")
                }

                # 提取基本结果
                data_table$alpha[i] <- data_any$alpha_final
                data_table$converged[i] <- data_any$converged
                data_table$iterations[i] <- data_any$iterations
                data_table$valid_starts[i] <- data_any$multi_start_info$valid_starts

                # 计算标准误
                tryCatch(
                    {
                        if (cov_index == "normal") {
                            alpha_se <- calculate_alpha_variance_b(data_any)
                        } else {
                            alpha_se <- calculate_alpha_variance_b(data_any)
                        }
                        data_table$alpha_se[i] <- alpha_se$std_error
                    },
                    error = function(e) {
                        if (verbose) {
                            cat(" [SE计算失败]")
                        }
                        data_table$alpha_se[i] <<- NA
                        data_table$error_message[i] <<- paste("SE计算失败:", e$message)
                    }
                )

                # 计算BIC
                tryCatch(
                    {
                        bic_value <- calculate_bic_cml_b(
                            input_model = data_any,
                            all_snps = all_snps,
                            n = n, penalty_factor = penalty_factor
                        )
                        data_table$bic[i] <- bic_value
                    },
                    error = function(e) {
                        if (verbose) {
                            cat(" [BIC计算失败]")
                        }
                        data_table$bic[i] <<- NA
                        if (is.na(data_table$error_message[i])) {
                            data_table$error_message[i] <<- paste("BIC计算失败:", e$message)
                        } else {
                            data_table$error_message[i] <<- paste(
                                data_table$error_message[i],
                                "; BIC计算失败:", e$message
                            )
                        }
                    }
                )

                # 标记为成功
                successful_combinations <- successful_combinations + 1
                last_successful_result <- data_any

                if (verbose) {
                    cat(sprintf(
                        " ✓ [alpha=%.6f, BIC=%.2f, 有效起点=%d/%d]\n",
                        data_any$alpha_final,
                        ifelse(is.na(data_table$bic[i]), -999, data_table$bic[i]),
                        data_any$multi_start_info$valid_starts,
                        n_starts
                    ))
                }
            },
            error = function(e) {
                # 处理完全失败的情况
                failed_combinations <- failed_combinations + 1
                failed_indices <- c(failed_indices, i)
                data_table$error_message[i] <<- e$message

                if (verbose) {
                    cat(sprintf(" ✗ [失败: %s]\n", e$message))
                }
            }
        )
    }

    if (verbose) {
        cat(paste(rep("=", 60), collapse = ""), "\n")
        cat("CML Family分析完成!\n")
        cat(sprintf(
            "成功组合: %d/%d (%.1f%%)\n",
            successful_combinations, total_combinations,
            100 * successful_combinations / total_combinations
        ))
        cat(sprintf(
            "失败组合: %d/%d (%.1f%%)\n",
            failed_combinations, total_combinations,
            100 * failed_combinations / total_combinations
        ))

        if (failed_combinations > 0) {
            cat("失败的组合 (a,b):\n")
            for (idx in failed_indices) {
                cat(sprintf(
                    "  (%d,%d): %s\n",
                    data_table$a[idx], data_table$b[idx],
                    data_table$error_message[idx]
                ))
            }
        }
    }

    # 检查是否有足够的成功结果进行加权估计
    valid_results <- !is.na(data_table$alpha) & !is.na(data_table$alpha_se) & !is.na(data_table$bic)

    if (sum(valid_results) == 0) {
        warning("没有任何组合成功完成，无法进行加权估计")
        result <- NULL
    } else {
        if (verbose && sum(valid_results) < successful_combinations) {
            cat(sprintf(
                "注意: 只有 %d/%d 个成功组合有完整的结果用于加权估计\n",
                sum(valid_results), successful_combinations
            ))
        }

        # 只对有效结果进行加权估计
        tryCatch(
            {
                result <- weighted_estimation_robust(data_table[valid_results, ])
            },
            error = function(e) {
                warning(paste("加权估计失败:", e$message))
                result <<- NULL
            }
        )
    }

    # 添加汇总信息
    summary_info <- list(
        total_combinations = total_combinations,
        successful_combinations = successful_combinations,
        failed_combinations = failed_combinations,
        valid_for_weighting = sum(valid_results),
        failed_indices = failed_indices,
        success_rate = successful_combinations / total_combinations
    )

    return(list(
        data_table = data_table,
        weighted_result = result,
        last_successful_model = last_successful_result,
        summary = summary_info
    ))
}
weighted_estimation_robust <- function(data_table) {
    # 调试信息

    # 确保数据是data.frame格式
    if (!is.data.frame(data_table)) {
        cat("转换为data.frame...\n")
        data_table <- as.data.frame(data_table)
    }

    # 检查必要的列是否存在
    required_cols <- c("alpha", "alpha_se", "bic")
    missing_cols <- required_cols[!required_cols %in% names(data_table)]
    if (length(missing_cols) > 0) {
        stop("缺少必要的列: ", paste(missing_cols, collapse = ", "))
    }



    # 手动检查NA值
    alpha_na <- is.na(data_table$alpha)
    alpha_se_na <- is.na(data_table$alpha_se)
    bic_na <- is.na(data_table$bic)

    # 找出完整的行
    valid_rows <- !alpha_na & !alpha_se_na & !bic_na



    if (sum(valid_rows) == 0) {
        warning("没有完整的数据行用于计算")
        return(list(
            weighted_alpha = NA,
            weighted_variance = NA,
            weighted_se = NA,
            n_valid = 0
        ))
    }

    # 提取有效数据
    alpha <- data_table$alpha[valid_rows]
    alpha_se <- data_table$alpha_se[valid_rows]
    bic <- data_table$bic[valid_rows]

    # 确保数据是数值型
    alpha <- as.numeric(alpha)
    alpha_se <- as.numeric(alpha_se)
    bic <- as.numeric(bic)

    # 检查是否还有NA值
    if (any(is.na(alpha)) || any(is.na(alpha_se)) || any(is.na(bic))) {
        stop("数据转换后仍有NA值")
    }


    # 计算基于BIC的权重
    min_bic <- min(bic)
    weights <- exp(-0.5 * (bic - min_bic))
    weights <- weights / sum(weights)



    # 计算加权点估计
    weighted_alpha <- sum(weights * alpha)

    # 计算加权标准误
    se_components <- weights * sqrt(alpha_se^2 + (alpha - weighted_alpha)^2)
    weighted_se <- sum(se_components)

    # 加权方差
    weighted_variance <- weighted_se^2

    # 返回结果
    return(list(
        weighted_alpha = weighted_alpha,
        weighted_variance = weighted_variance,
        weighted_se = weighted_se,
        weights = weights,
        n_valid = sum(valid_rows),
        valid_indices = which(valid_rows)
    ))
}
# %% 实验
if (FALSE) {
    source("样本生成函数/三联体家庭结构生成函数.R")
    source("FGWAS 优化版本/FGWAS 函数.R")
    n_sample <- 2000
    test <- generate_mr_trio_data_matrix_ultra(
        n = n_sample, n_snps = 10, n_pleiotropic = 0,
        n_expose_heterogeneity = 0,
        p_f = 0.3, p_m = 0.3,
        beta_fs_to_oe_exp = 0.1,
        beta_ms_to_oe_exp = 0.1,
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
        crowd_differences = 0, beta_exp_to_out = 0,
        confounding_exp = 0, confounding_out = 0,
        correlation = 0.7, seed = NULL
    )
    data_df_exp <- fgwas_for_data_matrix(test$data_exp,
        processing_func = FMR_trio_IFGLS,
        predicted_outcome = "expose"
    )
    data_df_out <- fgwas_for_data_matrix(test$data_out,
        processing_func = FMR_trio_IFGLS,
        predicted_outcome = "outcome"
    )

    test_3 <- cml_family_ultral_b(
        beta_hat_exp = data_df_exp$beta_hat, beta_hat_out = data_df_out$beta_hat,
        beta_sigma_exp = data_df_exp$Sigma_inv,
        beta_sigma_out = data_df_out$Sigma_inv, n = n_sample, cov_index = "normal"
    )
    head(arrange(test_3[[1]], bic))




    alpha_se_2 <- calculate_alpha_variance_b(test_3)$std_error

    alpha <- test_3$alpha_final
    alpha_se_1 <- 0
    for (i in 1:nrow(data_df_exp$beta_hat)) {
        matrix_index <- (2 * i - 1):(2 * i)
        current_exp_se <- solve(data_df_exp$Sigma_inv[matrix_index, ])[1, 1]
        current_out_se <- solve(data_df_out$Sigma_inv[matrix_index, ])[1, 1]
        current_beta_exp <- data_df_exp$beta_hat[i, 1]
        current_beta_out <- data_df_out$beta_hat[i, 1]
        fenzi <- current_beta_exp^2 * current_out_se +
            current_beta_out^2 * current_exp_se
        fenmu <- (alpha^2 * current_exp_se + current_out_se)^2

        alpha_se_1 <- fenzi / fenmu + alpha_se_1
        print(fenzi)
        print(fenmu)
        print(alpha_se_1)
    }
    sqrt(1 / alpha_se_1)
    alpha_se_2$std_error
    test_4$MA_BIC_se
    mr_input_obj <- MendelianRandomization::mr_input(
        bx = as.vector(phase_two_data_analysis_exp$beta_hat[, 1]),
        bxse = as.vector(phase_two_data_analysis_exp$beta_hat_se),
        by = as.vector(phase_two_data_analysis_out$beta_hat[, 1]),
        byse = as.vector(phase_two_data_analysis_out$beta_hat_se)
    )
    res_c <- MendelianRandomization::mr_ivw(mr_input_obj)
    res_c$StdError

    p_value <- 2 * pnorm(abs(test_3$alpha_final / alpha_se$std_error), lower.tail = FALSE)
    p_value
}

# %% 分析方差(cML_family_se)
if (FALSE) {
    n_simulations <- 1000 # 模拟次数，可以根据需要调整
    results <- data.frame(
        alpha = numeric(n_simulations),
        alpha_se = numeric(n_simulations),
        z_stat = numeric(n_simulations),
        p_value = numeric(n_simulations)
    )

    # 设置进度条
    cat("开始模拟分析...\n")
    pb <- txtProgressBar(min = 0, max = n_simulations, style = 3)

    # 进行模拟
    for (i in 1:n_simulations) {
        # 生成数据
        test <- generate_mr_trio_data_matrix_ultra(
            n = 1000, n_snps = 10, n_pleiotropic = 0, n_expose_heterogeneity = 0,
            p_f = 0.3, p_m = 0.3,
            beta_fs_to_oe_exp = 0.3,
            beta_ms_to_oe_exp = 0.3,
            beta_os_to_oe_exp = 0.3,
            h_beta_fs_to_oe_exp = 0.2, h_beta_ms_to_oe_exp = 0.1,
            h_beta_os_to_oe_exp = 0.3,
            mean_beta_fs_to_oe_out = 0.1, sd_beta_fs_to_oe_out = 0.05,
            mean_beta_ms_to_oe_out = 0.1, sd_beta_ms_to_oe_out = 0.05,
            mean_beta_os_to_oe_out = 0.1, sd_beta_os_to_oe_out = 0.05,
            p_neg_pleiotropy = 0.5, assortative_mating_strength = 0,
            crowd_differences = 0, beta_exp_to_out = 0,
            confounding_exp = 0.2, confounding_out = 0.2, correlation = 0.2,
            seed = NULL
        )

        # 处理暴露数据
        data_df_exp <- fgwas_for_data_matrix(test$data_exp,
            processing_func = FMR_trio_ols,
            predicted_outcome = "expose"
        )

        # 处理结局数据
        data_df_out <- fgwas_for_data_matrix(test$data_out,
            processing_func = FMR_trio_ols,
            predicted_outcome = "outcome"
        )

        # 运行MLE优化
        test_result <- run_iterative_mle_b_zero(
            beta_hat_exp = data_df_exp$beta_hat, beta_hat_out = data_df_out$beta_hat,
            beta_sigma_exp = data_df_exp$Sigma_inv,
            beta_sigma_out = data_df_out$Sigma_inv,
            a = 10, b = 10
        )

        # 计算标准误
        alpha_se <- calculate_alpha_variance_b(test_result)$std_error


        # 计算z统计量和p值
        z_stat <- test_result$alpha_final / alpha_se
        p_value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

        # 存储结果
        results[i, ] <- c(test_result$alpha_final, alpha_se, z_stat, p_value)

        # 更新进度条
        setTxtProgressBar(pb, i)
    }

    close(pb)
    cat("\n模拟完成！\n")

    # 计算描述性统计
    cat("=== 描述性统计 ===\n")
    cat("Alpha估计值:\n")
    print(summary(results$alpha))
    cat("\nZ统计量:\n")
    print(summary(results$z_stat))
    cat("\nP值:\n")
    print(summary(results$p_value))

    # 计算显著性比例
    sig_005 <- mean(results$p_value < 0.05)
    sig_001 <- mean(results$p_value < 0.01)
    cat(sprintf("\nP < 0.05的比例: %.3f\n", sig_005))
    cat(sprintf("P < 0.01的比例: %.3f\n", sig_001))

    # 可视化结果
    library(ggplot2)
    library(gridExtra)

    # 1. Alpha估计值的分布
    p1 <- ggplot(results, aes(x = alpha)) +
        geom_histogram(bins = 30, fill = "skyblue", color = "black", alpha = 0.7) +
        geom_vline(xintercept = mean(results$alpha), color = "red", linetype = "dashed", size = 1) +
        labs(title = "Alpha估计值分布", x = "Alpha", y = "频数") +
        theme_minimal()

    # 2. Z统计量的分布
    p2 <- ggplot(results, aes(x = z_stat)) +
        geom_histogram(bins = 30, fill = "lightgreen", color = "black", alpha = 0.7) +
        geom_vline(xintercept = 0, color = "red", linetype = "dashed", size = 1) +
        geom_vline(xintercept = c(-1.96, 1.96), color = "orange", linetype = "dashed", size = 0.8) +
        labs(title = "Z统计量分布", x = "Z统计量", y = "频数") +
        theme_minimal()

    # 3. P值的分布
    p3 <- ggplot(results, aes(x = p_value)) +
        geom_histogram(bins = 30, fill = "lightcoral", color = "black", alpha = 0.7) +
        geom_vline(xintercept = 0.05, color = "red", linetype = "dashed", size = 1) +
        geom_vline(xintercept = 0.01, color = "orange", linetype = "dashed", size = 1) +
        labs(title = "P值分布", x = "P值", y = "频数") +
        theme_minimal()

    # 4. Q-Q图检验Z统计量是否符合标准正态分布
    p4 <- ggplot(results, aes(sample = z_stat)) +
        stat_qq() +
        stat_qq_line(color = "red") +
        labs(title = "Z统计量Q-Q图", x = "理论分位数", y = "样本分位数") +
        theme_minimal()

    # 组合图形
    grid.arrange(p1, p2, p3, p4, ncol = 2)

    # 额外分析：检验统计性质
    cat("\n=== 统计性质检验 ===\n")
    cat(sprintf("Z统计量均值: %.4f (理论值: 0)\n", mean(results$z_stat)))
    cat(sprintf("Z统计量标准差: %.4f (理论值: 1)\n", sd(results$z_stat)))

    # Shapiro-Wilk正态性检验（如果样本量不太大）
    if (n_simulations <= 5000) {
        shapiro_test <- shapiro.test(results$z_stat)
        cat(sprintf("Z统计量Shapiro-Wilk正态性检验 p值: %.4f\n", shapiro_test$p.value))
    }

    # Kolmogorov-Smirnov检验
    ks_test <- ks.test(results$z_stat, "pnorm", 0, 1)
    cat(sprintf("Z统计量K-S正态性检验 p值: %.4f\n", ks_test$p.value))
}

# %% 似然函数可视化
# 假设您已经运行了您的算法
# estimation_result <- run_iterative_mle_b(...)
if (FALSE) {
    estimation_result <- test_3
    # --- 从结果中提取对象 ---
    alpha_hat <- estimation_result$alpha_final
    gamma_hat_final <- t(estimation_result$gamma_final) # 转置回 2xN 格式
    all_data <- estimation_result$data

    # --- 挑选一个你想可视化的、高度相关的SNP的索引 ---
    # 根据您的截图，我们选第一个SNP，因为您已经验证了它的高相关性
    m <- 1

    # --- 提取这个特定SNP m的相关数据 ---
    beta_hat_m <- all_data$beta_hat_out[m, ]
    gamma_hat_m <- all_data$beta_hat_exp[m, ]

    Sigma_beta_m <- all_data$Sigma_beta_list[[m]]
    Sigma_gamma_m <- all_data$Sigma_gamma_list[[m]]

    Sigma_beta_inv_m <- solve(Sigma_beta_m)
    Sigma_gamma_inv_m <- solve(Sigma_gamma_m)

    # 找到这个SNP最终的gamma估计值，这将是我们的绘图中心
    gamma_m_final <- gamma_hat_final[, m]

    # 定义一个函数，用于计算在给定gamma_m值下的（负）对数似然 l
    # 注意：我们只计算与gamma_m相关的部分，因为其他部分在绘图时是常数
    calculate_likelihood_for_snp_m <- function(gamma_vec,
                                               alpha_val,
                                               beta_hat_vec,
                                               gamma_hat_vec,
                                               Sigma_beta_inv,
                                               Sigma_gamma_inv) {
        term1 <- t(gamma_hat_vec - gamma_vec) %*% Sigma_gamma_inv %*% (gamma_hat_vec - gamma_vec)
        term2 <- t(beta_hat_vec - alpha_val * gamma_vec) %*% Sigma_beta_inv %*% (beta_hat_vec - alpha_val * gamma_vec)

        # 乘以0.5是标准形式，但对于可视化形状而言不是必须的
        return(0.5 * (term1 + term2))
    }

    # 以最终估计值为中心，创建一个gamma_o和gamma_p的值的序列
    # grid_range可以调整，它决定了我们观察的“视野”大小
    grid_range <- 0.1
    grid_points <- 100 # 网格密度

    go_seq <- seq(gamma_m_final[1] - grid_range, gamma_m_final[1] + grid_range, length.out = grid_points)
    gp_seq <- seq(gamma_m_final[2] - grid_range, gamma_m_final[2] + grid_range, length.out = grid_points)

    # 初始化一个矩阵来存储每个网格点的似然值
    likelihood_matrix <- matrix(NA, nrow = grid_points, ncol = grid_points)

    # 遍历网格，填充似然矩阵
    for (i in 1:grid_points) {
        for (j in 1:grid_points) {
            current_gamma <- c(go_seq[i], gp_seq[j])
            likelihood_matrix[i, j] <- calculate_likelihood_for_snp_m(
                gamma_vec = current_gamma,
                alpha_val = alpha_hat,
                beta_hat_vec = beta_hat_m,
                gamma_hat_vec = gamma_hat_m,
                Sigma_beta_inv = Sigma_beta_inv_m,
                Sigma_gamma_inv = Sigma_gamma_inv_m
            )
        }
    }

    # --- 绘制等高线图 ---
    par(pty = "s") # 设置绘图区域为正方形，以正确显示形状

    contour(
        x = go_seq, y = gp_seq, z = likelihood_matrix,
        levels = quantile(likelihood_matrix, probs = seq(0, 1, 0.05)), # 选择一些等高线水平
        xlab = expression(gamma[o]), # 使用表达式显示希腊字母
        ylab = expression(gamma[p]),
        main = paste("Likelihood Contour for SNP", m, "(Corr =", round(cov2cor(Sigma_beta_m)[1, 2], 2), ")"),
        drawlabels = FALSE, # 不在图上显示等高线数值
        lwd = 1.5
    )

    # 在图上标记出最终的估计点（“谷底”）
    points(gamma_m_final[1], gamma_m_final[2], pch = 4, col = "red", cex = 2, lwd = 3)

    legend("topright", legend = "MLE Estimate", pch = 4, col = "red", pt.cex = 2, pt.lwd = 2, bty = "n")
}
