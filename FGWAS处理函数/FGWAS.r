# %% 需要用到的包
library(lme4)
library(dplyr)
library(MASS)
library(nlme) # <-- 已添加
library(parallel)
library(pbapply)
# %% FGWAS 多个数据

FMR_trio_optimized_intercept <- function(data_test, y_index = "expose") {
    n <- nrow(data_test) # 获取样本量

    # --- 1. 计算协方差矩阵 Omega (此部分无变化) ---
    if (y_index == "expose") {
        exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
    } else {
        exposures <- data_test[, c("Father_outcome", "Mother_outcome", "Offspring_outcome")]
    }
    father_exp <- exposures[[1]]
    mother_exp <- exposures[[2]]
    offspring_exp <- exposures[[3]]

    sa2 <- (cov(father_exp, offspring_exp) + cov(mother_exp, offspring_exp)) / 2
    ss2 <- cov(father_exp, mother_exp)
    se2_F <- var(father_exp)
    se2_M <- var(mother_exp)
    se2_O <- var(offspring_exp)

    Sigma <- matrix(c(
        se2_F, ss2,   sa2,
        ss2,   se2_M, sa2,
        sa2,   sa2,   se2_O
    ), nrow = 3, byrow = TRUE)

    Omega_inv <- tryCatch(
        {
            chol2inv(chol(Sigma))
        },
        error = function(e) {
            warning("协方差矩阵 Sigma 求逆失败 (可能奇异): ", e$message)
            matrix(NaN, 3, 3)
        }
    )
    if (any(is.nan(Omega_inv))) {
        ## <<< CHANGED: 返回值的维度增加
        return(list(beta = matrix(NaN, 1, 3), informatrix = matrix(NaN, 3, 3), beta_se = rep(NaN, 3)))
    }


    # --- 2. 准备基因型矩阵 X, X_par, 以及截距项矩阵 X_int ---
    snps <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps[[1]]
    mother_snps <- snps[[2]]
    offspring_snps <- snps[[3]]

    p <- mean(c(mother_snps, offspring_snps)) / 2

    G1 <- ifelse(father_snps == 0, 2 * p, ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p)))
    G2 <- ifelse(mother_snps == 0, 2 * p, ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p)))
    G3 <- father_snps + mother_snps

    # 构建基因型矩阵
    X <- as.matrix(snps)
    X_par <- cbind(G1, G2, G3)
    Y <- as.matrix(exposures)

    ## <<< NEW: 为截距项创建设计矩阵
    # 每一行都是 [1, 1, 1]，代表为父、母、子都设置一个共同的截距
    X_int <- matrix(1, nrow = n, ncol = 3)


    # --- 3. 向量化计算 fenzi 和 informatrix (包含截距项) ---
    # 预计算常用乘积
    Y_OmInv <- Y %*% Omega_inv
    X_OmInv <- X %*% Omega_inv
    Xp_OmInv <- X_par %*% Omega_inv
    ## <<< NEW: 计算截距项的常用乘积
    X_int_OmInv <- X_int %*% Omega_inv

    # 计算 fenzi (现在是 1x3 矩阵) 的三个元素
    ## <<< NEW: 增加截距项的 fenzi
    fenzi_0 <- sum(Y_OmInv * X_int) # 截距项
    fenzi_1 <- sum(Y_OmInv * X) # 原 beta_G 项
    fenzi_2 <- sum(Y_OmInv * X_par) # 原 beta_G_par 项

    ## <<< CHANGED: fenzi 现在是 1x3
    fenzi <- matrix(c(fenzi_0, fenzi_1, fenzi_2), nrow = 1)

    # 计算 informatrix (现在是 3x3 矩阵) 的元素
    ## <<< NEW: 增加与截距项相关的元素
    info_00 <- sum(X_int_OmInv * X_int)
    info_01 <- sum(X_int_OmInv * X)
    info_02 <- sum(X_int_OmInv * X_par)

    info_11 <- sum(X_OmInv * X)
    info_12 <- sum(X_OmInv * X_par)
    info_22 <- sum(Xp_OmInv * X_par)

    ## <<< CHANGED: informatrix 现在是 3x3 的对称矩阵
    informatrix <- matrix(c(
        info_00, info_01, info_02,
        info_01, info_11, info_12,
        info_02, info_12, info_22
    ), nrow = 3, byrow = TRUE)


    # --- 4. 计算 Beta 和标准误 ---
    informatrix_inv <- tryCatch(
        {
            chol2inv(chol(informatrix))
        },
        error = function(e) {
            warning("信息矩阵 informatrix 求逆失败 (可能奇异): ", e$message)
            matrix(NaN, 3, 3)
        }
    )

    ## <<< CHANGED: 调整 beta 和 se 的维度
    if (any(is.nan(informatrix_inv))) {
        beta <- matrix(NaN, 1, 3)
        beta_se <- rep(NaN, 3)
    } else {
        # 计算 beta (现在是 1x3 矩阵)
        beta <- fenzi %*% informatrix_inv
        # 为beta的三个系数添加名称，方便解读
        colnames(beta) <- c("Intercept", "beta_G", "beta_G_par")

        # 计算所有 beta 的标准误
        # 标准误是方差的平方根，方差是协方差矩阵（即 informatrix_inv）的对角元素
        beta_se <- sqrt(pmax(diag(informatrix_inv), 0))
    }

    # --- 5. 返回结果 ---
    # 返回与原函数结构相似但维度已更新的列表
    A <- list(
        beta = beta[, 2:3], # 1x3 矩阵，包含截距, beta_G, beta_G_par
        informatrix = informatrix[2:3, 2:3], # 3x3 信息矩阵
        beta_se = beta_se[2] # 包含3个标准误的向量
    )

    return(A)
}
# lmm两个变量(FGWAS)
FMR_trio_optimized <- function(data_test, y_index = "expose") {
    # --- 1. 计算协方差矩阵 Omega ---
    # 使用 colMeans/cov/var 提高效率和简洁性

    if (y_index == "expose") {
        exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
        father_exp <- exposures[[1]]
        mother_exp <- exposures[[2]]
        offspring_exp <- exposures[[3]]
    } else {
        exposures <- data_test[, c(
            "Father_outcome",
            "Mother_outcome", "Offspring_outcome"
        )]
        father_exp <- exposures[[1]]
        mother_exp <- exposures[[2]]
        offspring_exp <- exposures[[3]]
    }


    # 计算协方差/方差 （注意：cov计算的是样本协方差，除以 n-1）
    # 如果原始代码的意图是总体协方差（除以 n），需要调整
    sa2 <- (cov(father_exp, offspring_exp) + cov(mother_exp, offspring_exp)) / 2
    ss2 <- cov(father_exp, mother_exp)
    se2_F <- var(father_exp)
    se2_M <- var(mother_exp)
    se2_O <- var(offspring_exp)

    Sigma <- matrix(c(
        se2_F, ss2,   sa2,
        ss2,   se2_M, sa2,
        sa2,   sa2,   se2_O
    ), nrow = 3, byrow = TRUE)

    # 计算 Omega 的逆
    # 使用 tryCatch 处理可能的奇异矩阵
    Omega_inv <- tryCatch(
        {
            chol2inv(chol(Sigma))
        },
        error = function(e) {
            warning("协方差矩阵 Sigma 求逆失败 (可能奇异): ", e$message)
            # 返回一个 NaN 矩阵或根据需要停止
            matrix(NaN, 3, 3)
        }
    )
    if (any(is.nan(Omega_inv))) {
        return(list(beta = matrix(NaN, 1, 2), informatrix = matrix(NaN, 2, 2), beta_se = NaN))
    }


    # --- 2. 准备基因型矩阵 X 和 X_par ---
    # 提取 SNP 数据
    snps <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps[[1]]
    mother_snps <- snps[[2]]
    offspring_snps <- snps[[3]]

    # 计算等位基因频率 p
    # 使用 mean 直接计算所有 SNP 的均值，然后除以 2
    ## 父母只用一个防止人群分层
    p <- mean(c(mother_snps, offspring_snps)) / 2

    # 计算 G1, G2, G3 (向量化)
    # 注意: 原始 G1/G2 公式的含义可能需要确认，这里按原样实现
    G1 <- ifelse(father_snps == 0, 2 * p,
        ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p))
    )
    G2 <- ifelse(mother_snps == 0, 2 * p,
        ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p))
    )
    G3 <- father_snps + mother_snps

    # 构建基因型矩阵
    X <- as.matrix(snps) # n x 3
    X_par <- cbind(G1, G2, G3) # n x 3
    Y <- as.matrix(exposures) # n x 3


    # --- 3. 向量化计算 fenzi 和 informatrix (原 fenmu) ---
    # 预计算常用乘积
    Y_OmInv <- Y %*% Omega_inv # n x 3
    X_OmInv <- X %*% Omega_inv # n x 3
    Xp_OmInv <- X_par %*% Omega_inv # n x 3

    # 计算 fenzi (1x2 矩阵) 的两个元素
    # sum(Y_OmInv * X) 等价于 sum(diag(t(Y) %*% Omega_inv %*% X))
    fenzi_1 <- sum(Y_OmInv * X)
    fenzi_2 <- sum(Y_OmInv * X_par)
    fenzi <- matrix(c(fenzi_1, fenzi_2), nrow = 1)

    # 计算 informatrix (2x2 矩阵) 的元素 (原 fenmu)
    # sum(X_OmInv * X) 等价于 sum(diag(t(X) %*% Omega_inv %*% X))
    info_11 <- sum(X_OmInv * X) # 原 beta_se 累加值
    info_12 <- sum(X_OmInv * X_par) # 原 beta_cov 累加值 (info_21 = info_12)
    info_22 <- sum(Xp_OmInv * X_par) # 原 beta_par_se 累加值

    informatrix <- matrix(c(
        info_11, info_12,
        info_12, info_22
    ), nrow = 2, byrow = TRUE)


    # --- 4. 计算 Beta 和标准误 ---
    # 计算 informatrix 的逆
    informatrix_inv <- tryCatch(
        {
            chol2inv(chol(informatrix))
        },
        error = function(e) {
            warning("信息矩阵 informatrix 求逆失败 (可能奇异): ", e$message)
            matrix(NaN, 2, 2)
        }
    )

    if (any(is.nan(informatrix_inv))) {
        beta <- matrix(NaN, 1, 2)
        beta_se1 <- NaN
    } else {
        # 计算 beta (1x2 矩阵)
        beta <- fenzi %*% informatrix_inv

        # 计算第一个 beta 的标准误
        # 标准误是方差的平方根，方差是协方差矩阵（即 informatrix_inv）的对角元素
        beta_se1 <- sqrt(pmax(informatrix_inv[1, 1], 0)) # 使用 pmax 避免对负数开方
    }

    # --- 5. 返回结果 ---
    # 返回与原函数相同的结构
    A <- list(
        beta = beta, # 1x2 矩阵
        informatrix = informatrix, # 2x2 信息矩阵
        beta_se = beta_se1 # 第一个 beta 的标准误 (标量)
        # 如果需要第二个beta的SE: sqrt(pmax(informatrix_inv[2, 2], 0))
        # 如果需要完整的beta协方差矩阵: informatrix_inv
    )

    return(A)
}
FMR_trio_FGLS <- function(data_test, y_index = "expose") {
    n <- nrow(data_test)

    # --- 准备 Y 和 X 矩阵 (这部分提前，为计算残差做准备) ---
    if (y_index == "expose") {
        exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
    } else {
        exposures <- data_test[, c("Father_outcome", "Mother_outcome", "Offspring_outcome")]
    }

    snps_data <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps_data[[1]]
    mother_snps <- snps_data[[2]]
    offspring_snps <- snps_data[[3]]

    p <- mean(c(mother_snps, offspring_snps)) / 2

    G1 <- ifelse(father_snps == 0, 2 * p, ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p)))
    G2 <- ifelse(mother_snps == 0, 2 * p, ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p)))
    G3 <- father_snps + mother_snps

    Y <- as.matrix(exposures)
    X_G <- as.matrix(snps_data)
    X_par <- cbind(G1, G2, G3)
    X_int <- matrix(1, nrow = n, ncol = 3)


    # --- 1. <<< NEW: 基于初步OLS模型的残差来估算协方差矩阵 Omega ---

    # 为了得到初步的残差，我们可以对每个家庭成员（父、母、子）分别运行OLS
    # 这是估算残差协方差结构的一种直接且有效的方法

    # 将宽数据转换为长数据以方便运行OLS
    # 当然也可以用矩阵运算，但这样更直观
    data_long <- data.frame(
        y = c(Y[, 1], Y[, 2], Y[, 3]),
        x_g = c(X_G[, 1], X_G[, 2], X_G[, 3]),
        x_par = c(X_par[, 1], X_par[, 2], X_par[, 3])
    )

    # 运行初步的OLS模型
    prelim_ols <- lm(y ~ x_g + x_par, data = data_long)

    # 获取残差并将其重塑回 n x 3 的矩阵形式
    residuals_matrix <- matrix(residuals(prelim_ols), nrow = n, ncol = 3)

    # 现在，使用残差来计算协方差矩阵 Sigma
    resid_F <- residuals_matrix[, 1]
    resid_M <- residuals_matrix[, 2]
    resid_O <- residuals_matrix[, 3]

    # 计算方法与之前相同，但输入的是残差
    sa2 <- (cov(resid_F, resid_O) + cov(resid_M, resid_O)) / 2
    ss2 <- cov(resid_F, resid_M)
    se2_F <- var(resid_F)
    se2_M <- var(resid_M)
    se2_O <- var(resid_O)

    Sigma <- matrix(c(
        se2_F, ss2,   sa2,
        ss2,   se2_M, sa2,
        sa2,   sa2,   se2_O
    ), nrow = 3, byrow = TRUE)

    # --- 后续步骤与 FMR_trio_optimized_intercept 完全相同 ---

    Omega_inv <- tryCatch(chol2inv(chol(Sigma)), error = function(e) matrix(NaN, 3, 3))
    if (any(is.nan(Omega_inv))) {
        return(list(beta = matrix(NaN, 1, 3), informatrix = matrix(NaN, 3, 3), beta_se = rep(NaN, 3)))
    }

    # --- 向量化计算 fenzi 和 informatrix ---
    Y_OmInv <- Y %*% Omega_inv
    X_int_OmInv <- X_int %*% Omega_inv
    X_G_OmInv <- X_G %*% Omega_inv
    Xp_OmInv <- X_par %*% Omega_inv

    fenzi_0 <- sum(Y_OmInv * X_int)
    fenzi_1 <- sum(Y_OmInv * X_G)
    fenzi_2 <- sum(Y_OmInv * X_par)
    fenzi <- matrix(c(fenzi_0, fenzi_1, fenzi_2), nrow = 1)

    info_00 <- sum(X_int_OmInv * X_int)
    info_01 <- sum(X_int_OmInv * X_G)
    info_02 <- sum(X_int_OmInv * X_par)
    info_11 <- sum(X_G_OmInv * X_G)
    info_12 <- sum(X_G_OmInv * X_par)
    info_22 <- sum(Xp_OmInv * X_par)

    informatrix <- matrix(c(
        info_00, info_01, info_02,
        info_01, info_11, info_12,
        info_02, info_12, info_22
    ), nrow = 3, byrow = TRUE)

    # --- 计算 Beta 和标准误 ---
    informatrix_inv <- tryCatch(chol2inv(chol(informatrix)), error = function(e) matrix(NaN, 3, 3))

    if (any(is.nan(informatrix_inv))) {
        beta <- matrix(NaN, 1, 3)
        beta_se <- rep(NaN, 3)
    } else {
        beta <- fenzi %*% informatrix_inv
        colnames(beta) <- c("Intercept", "beta_G", "beta_G_par")
        beta_se <- sqrt(pmax(diag(informatrix_inv), 0))
    }

    # --- 返回结果 ---
    A <- list(
        beta = beta[, 2:3], # 1x2 矩阵，包含 beta_G 和 beta_G_par
        informatrix = informatrix[2:3, 2:3],
        beta_se = beta_se[2]
    )

    return(A)
}
FMR_trio_IFGLS <- function(data_test, y_index = "expose",
                           max_iter = 20, tolerance = 1e-7) {
    n <- nrow(data_test)
    cat("Starting IFGLS procedure...\n")

    # --- 1. 准备 Y 和 X 矩阵 (仅执行一次) ---
    if (y_index == "expose") {
        exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
    } else {
        exposures <- data_test[, c("Father_outcome", "Mother_outcome", "Offspring_outcome")]
    }

    snps_data <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps_data[[1]]
    mother_snps <- snps_data[[2]]
    offspring_snps <- snps_data[[3]]

    p <- mean(c(mother_snps, offspring_snps)) / 2

    G1 <- ifelse(father_snps == 0, 2 * p,
        ifelse(father_snps == 1, 1 + 2 * p, 2 * (1 + p))
    )
    G2 <- ifelse(mother_snps == 0, 2 * p,
        ifelse(mother_snps == 1, 1 + 2 * p, 2 * (1 + p))
    )
    G3 <- father_snps + mother_snps

    Y <- as.matrix(exposures)
    X_G <- as.matrix(snps_data)
    X_par <- cbind(G1, G2, G3)
    X_int <- matrix(1, nrow = n, ncol = 3)


    # --- 2. 初始化：使用OLS获得beta的初始估计值 ---
    data_long <- data.frame(
        y = c(Y[, 1], Y[, 2], Y[, 3]),
        x_g = c(X_G[, 1], X_G[, 2], X_G[, 3]),
        x_par = c(X_par[, 1], X_par[, 2], X_par[, 3])
    )
    prelim_ols <- lm(y ~ x_g + x_par, data = data_long)
    beta_current <- matrix(coef(prelim_ols), nrow = 1)
    colnames(beta_current) <- c("Intercept", "beta_G", "beta_G_par")

    cat("Initial OLS estimates (beta):", round(beta_current, 5), "\n")

    # --- 3. 迭代循环 ---
    converged <- FALSE
    for (i in 1:max_iter) {
        # 存储上一次的beta以供比较
        beta_old <- beta_current

        # --- 步骤 A: 基于当前的beta估计值计算残差 ---
        predicted_Y <- X_int * beta_current[1] + X_G * beta_current[2] + X_par * beta_current[3]
        residuals_matrix <- Y - predicted_Y

        # --- 步骤 B: 使用残差更新协方差矩阵 Sigma ---
        resid_F <- residuals_matrix[, 1]
        resid_M <- residuals_matrix[, 2]
        resid_O <- residuals_matrix[, 3]

        sa2 <- (cov(resid_F, resid_O) + cov(resid_M, resid_O)) / 2
        ss2 <- cov(resid_F, resid_M)
        se2_F <- var(resid_F)
        se2_M <- var(resid_M)
        se2_O <- var(resid_O)

        Sigma <- matrix(c(se2_F, ss2, sa2, ss2, se2_M, sa2, sa2, sa2, se2_O), nrow = 3, byrow = TRUE)

        Omega_inv <- tryCatch(chol2inv(chol(Sigma)), error = function(e) matrix(NaN, 3, 3))
        if (any(is.nan(Omega_inv))) {
            warning("Sigma matrix became singular during iteration. Stopping.")
            break
        }

        # --- 步骤 C: 使用更新后的Sigma(Omega)执行一步GLS，得到新的beta ---
        Y_OmInv <- Y %*% Omega_inv
        X_int_OmInv <- X_int %*% Omega_inv
        X_G_OmInv <- X_G %*% Omega_inv
        Xp_OmInv <- X_par %*% Omega_inv

        fenzi <- matrix(c(sum(Y_OmInv * X_int), sum(Y_OmInv * X_G), sum(Y_OmInv * X_par)), nrow = 1)

        info_00 <- sum(X_int_OmInv * X_int)
        info_01 <- sum(X_int_OmInv * X_G)
        info_02 <- sum(X_int_OmInv * X_par)
        info_11 <- sum(X_G_OmInv * X_G)
        info_12 <- sum(X_G_OmInv * X_par)
        info_22 <- sum(Xp_OmInv * X_par)

        informatrix <- matrix(c(info_00, info_01, info_02, info_01, info_11, info_12, info_02, info_12, info_22), nrow = 3, byrow = TRUE)

        informatrix_inv <- tryCatch(chol2inv(chol(informatrix)), error = function(e) matrix(NaN, 3, 3))

        if (any(is.nan(informatrix_inv))) {
            warning("Information matrix became singular during iteration. Stopping.")
            beta_current <- matrix(NaN, 1, 3)
            break
        }

        beta_current <- fenzi %*% informatrix_inv
        colnames(beta_current) <- c("Intercept", "beta_G", "beta_G_par")

        # --- 步骤 D: 检查收敛 ---
        # 计算新旧beta估计值之间的变化（平方差之和）
        change <- sum((beta_current - beta_old)^2)
        cat("Iteration", i, "| Change in beta:", sprintf("%.4g", change), "\n")

        if (change < tolerance) {
            cat("Convergence reached after", i, "iterations.\n")
            converged <- TRUE
            break
        }
    } # 迭代循环结束

    if (!converged) {
        warning("IFGLS did not converge after ", max_iter, " iterations.")
    }

    # --- 4. 计算最终的标准误并返回结果 ---
    informatrix_inv_small <- informatrix_inv[2:3, 2:3]
    informatrix_small <- solve(informatrix_inv_small)
    beta_se <- sqrt(informatrix_inv[2, 2])
    A <- list(
        beta = beta_current[, 2:3],
        informatrix = informatrix_small,
        beta_se = beta_se
    )

    return(A)
}

FMR_trio_ols <- function(data_test, y_index = "expose") {
    n <- nrow(data_test)
    # --- 1. 准备 Y 和 X 矩阵 (仅执行一次) ---
    if (y_index == "expose") {
        exposures <- data_test[, c("Father_expose", "Mother_expose", "Offspring_expose")]
    } else {
        exposures <- data_test[, c("Father_outcome", "Mother_outcome", "Offspring_outcome")]
    }

    snps_data <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps_data[[1]]
    mother_snps <- snps_data[[2]]
    offspring_snps <- snps_data[[3]]

    offspring_expose <- exposures[, 3]

    G3 <- father_snps + mother_snps
    lm_model <- summary(lm(offspring_expose$Offspring_expose ~ offspring_snps +
        father_snps + mother_snps))
    beta_current <- coef(lm_model)
    informatrix <- solve(vcov(lm_model)[2:3, 2:3])
    beta_se <- lm_model$coefficients["offspring_snps", "Std. Error"]
    A <- list(
        beta = beta_current[2:3],
        informatrix = informatrix,
        beta_se = beta_se
    )

    return(A)
}

FMR_single_ols <- function(data_test, y_index = "expose") {
    n <- nrow(data_test)

    # --- 1. 准备 Y 和 X 矩阵 (仅执行一次) ---
    if (y_index == "expose") {
        offspring_expose <- data_test$Offspring_expose
    } else {
        offspring_expose <- data_test$Offspring_outcome
    }

    snps_data <- data_test[, c("Father_SNPs", "Mother_SNPs", "Offspring_SNPs")]
    father_snps <- snps_data[[1]]
    mother_snps <- snps_data[[2]]
    offspring_snps <- snps_data[[3]]


    G3 <- father_snps + mother_snps

    lm_model <- summary(lm(offspring_expose ~ offspring_snps))
    beta_current <- coef(lm_model)
    beta_se <- lm_model$coefficients["offspring_snps", "Std. Error"]
    A <- list(
        beta = beta_current[2],
        beta_se = beta_se
    )

    return(A)
}

extract_snp <- function(data_df, snp_number) {
    # 创建新的数据框
    Father_SNPs <- data_df$father_snps[, snp_number]

    result_df <- tibble(
        Father_SNPs = data_df$father_snps[, snp_number],
        Father_expose = data_df$father_expose,
        Father_outcome = data_df$father_outcome,
        Mother_SNPs = data_df$mother_snps[, snp_number],
        Mother_expose = data_df$mother_expose,
        Mother_outcome = data_df$mother_outcome,
        Offspring_SNPs = data_df$offspring_snps[, snp_number],
        Offspring_expose = data_df$offspring_expose,
        Offspring_outcome = data_df$offspring_outcome
    )

    return(result_df)
}
fgwas_for_data_matrix <- function(data_df,
                                  predicted_outcome = "expose",
                                  processing_func = FMR_trio_IFGLS) { # 假设 FMR_trio 存在
    n_snps <- ncol(data_df$father_snps)




    beta_hat_matrix <- matrix(0, nrow = n_snps, ncol = 2)
    Sigma_inv_matrix <- matrix(0, nrow = 2 * n_snps, ncol = 2)
    beta_hat_se <- matrix(0, nrow = n_snps, ncol = 1)
    # --- 2， 遍历数据集并应用处理函数 ---
    for (i in 1:n_snps) {
        subset_data <- extract_snp(data_df, i)

        processed_result <- processing_func(
            subset_data,
            predicted_outcome
        )
        # 存储 beta 估计值
        beta_hat_matrix[i, ] <- processed_result[[1]]
        # 存储 2x2 的 Sigma_inv 块
        row_indices <- (2 * i - 1):(2 * i)
        Sigma_inv_matrix[row_indices, ] <- processed_result[[2]]
        beta_hat_se[i] <- processed_result[[3]]
        # beta_hat_matrix[i,] 和 Sigma_inv_matrix[row_indices,] 保持 NA
    } # 结束 for 循环

    cat("数据集处理完成。\n")

    # --- 返回结果 ---
    # 检查是否有任何数据集处理成功

    final_result <- list(
        beta_hat = beta_hat_matrix, Sigma_inv = Sigma_inv_matrix,
        beta_hat_se = beta_hat_se
    )
    return(final_result)
}

fgwas_for_data_matrix_single <- function(data_df,
                                         predicted_outcome = "expose",
                                         processing_func = FMR_single_ols) {
    n_snps <- ncol(data_df$father_snps)
    # --- 1. 构建提取子数据集的函数 ---


    beta_hat <- matrix(0, nrow = n_snps, ncol = 1)
    beta_hat_se <- matrix(0, nrow = n_snps, ncol = 1)
    # --- 2， 遍历数据集并应用处理函数 ---
    for (i in 1:n_snps) {
        subset_data <- extract_snp(data_df, i)

        processed_result <- processing_func(
            subset_data,
            predicted_outcome
        )
        # 存储 beta 估计值
        beta_hat[i] <- processed_result[[1]]
        beta_hat_se[i] <- processed_result[[2]]
        # beta_hat_matrix[i,] 和 Sigma_inv_matrix[row_indices,] 保持 NA
    } # 结束 for 循环

    cat("数据集处理完成。\n")

    # --- 返回结果 ---
    # 检查是否有任何数据集处理成功

    final_result <- list(
        beta_hat = beta_hat,
        beta_hat_se = beta_hat_se
    )
    return(final_result)
}
# %% 测试


if (FALSE) {
    test <- generate_mr_trio_data_matrix_ultra()
    data_exp <- test$data_exp$data
    test_1 <- fgwas_for_data_matrix_single(data_exp)
    test_2 <- fgwas_for_data_matrix(data_exp)
    test_2
    test_1
}
