# %%
# 多核模拟函数
library(parallel)
library(doParallel)
library(foreach)
# 载入必要的R包
# 请确保已经安装了这些包: install.packages(c("parallel", "dplyr"))
library(parallel)
library(dplyr)
library(future)
# install.packages("furrr", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
library(furrr)
library(dplyr)
library(progressr)

# %% 函数构建

# 多核MR模拟实现
# 适用于Windows平台

library(parallel)
library(dplyr)
library(tibble)


run_parallel_simulation <- function(n_simulations,
                                    script_path,
                                    n_cores = NULL,
                                    ...) {
    # --- 1. 检查和加载依赖包 ---
    stopifnot(
        "请先安装 'parallel' 包 (通常是R自带的)" = requireNamespace("parallel", quietly = TRUE),
        "请先安装 'dplyr' 包: install.packages('dplyr')" = requireNamespace("dplyr", quietly = TRUE),
        "请先安装 'tibble' 包: install.packages('tibble')" = requireNamespace("tibble", quietly = TRUE),
        "请先安装 'data.table' 包: install.packages('data.table')" = requireNamespace("data.table", quietly = TRUE)
    )

    # --- 2. 并行环境设置 ---
    if (!file.exists(script_path)) stop(paste("错误: 依赖脚本未找到，请检查路径:", script_path))
    if (is.null(n_cores)) {
        n_cores <- parallel::detectCores() - 1
        if (n_cores < 1) n_cores <- 1
    }
    message(paste("检测到", parallel::detectCores(), "个核心，将使用", n_cores, "个核心进行并行计算..."))
    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    on.exit(message("并行集群已自动关闭。"), add = TRUE)

    # --- 3. 导出依赖 (核心修正点) ---
    # 将所有 ... 参数打包成一个列表
    dot_args <- list(...)

    # 导出依赖脚本和打包好的参数列表
    tryCatch(
        {
            # 导出脚本路径和参数列表这两个对象
            parallel::clusterExport(cl, c("script_path", "dot_args"), envir = environment())
            # 在每个核心上加载脚本
            parallel::clusterEvalQ(cl, {
                source(script_path)
                library(MendelianRandomization)
            })
        },
        error = function(e) stop(paste("在核心上加载依赖脚本时出错: ", e$message))
    )

    # --- 4. 运行并行模拟 (核心修正点) ---
    message(paste("开始运行", n_simulations, "次模拟..."))
    start_time <- Sys.time()

    results_list <- parallel::parLapply(cl, 1:n_simulations, function(i) {
        # 在每个核心上，合并参数。用户传入的 dot_args 优先级更高。
        # 如果用户传入了 seed，则使用用户的 seed；否则，使用循环的 i 作为 seed。
        final_args <- dot_args
        if (!"seed" %in% names(final_args)) {
            final_args$seed <- i
        }

        tryCatch(
            {
                # 关键修复：使用 do.call 而不是 !!! 操作符
                res <- do.call(triplet_family_simulation_once_robust, final_args)
                list(success = TRUE, result = res)
            },
            error = function(e) {
                list(success = FALSE, error = paste0("模拟 #", i, " 失败: ", e$message))
            }
        )
    })

    end_time <- Sys.time()
    message(paste("并行模拟完成！总耗时:", round(end_time - start_time, 2), "秒"))

    # --- 5. 数据整合 (使用最稳健的版本) ---
    message("正在汇总和检查结果...")
    success_flags <- sapply(results_list, `[[`, "success")
    successful_results <- results_list[success_flags]
    failed_results <- results_list[!success_flags]

    if (length(failed_results) > 0) {
        warning(paste(length(failed_results), "次模拟中有", length(failed_results), "次失败了。"))
        # 打印前5个错误
        lapply(failed_results[1:min(5, length(failed_results))], function(x) message(x$error))
    }

    if (length(successful_results) == 0) {
        message("所有模拟都失败了，无法生成结果 tibble。")
        return(list(results = tibble::tibble(), errors = failed_results))
    }

    # 使用之前版本中被证明最稳健的清洗和合并逻辑
    successful_data_list <- lapply(successful_results, `[[`, "result")
    sanitized_list <- lapply(successful_data_list, function(single_run_list) {
        lapply(single_run_list, function(element) {
            if (is.null(element) || length(element) == 0) {
                return(NA_real_)
            }
            return(element)
        })
    })

    final_dt <- data.table::rbindlist(sanitized_list, use.names = TRUE, fill = TRUE)
    final_tibble <- tibble::as_tibble(final_dt) %>%
        tibble::add_column(simulation_id = which(success_flags), .before = 1)

    message("所有任务完成！")
    return(list(results = final_tibble, errors = failed_results))
}
# 将仿真结果转换为画图函数所需的汇总统计量格式
calculate_summary_stats <- function(simulation_results,
                                    true_value = 0,
                                    alpha_levels = c(0.05, 0.01),
                                    confidence_levels = c(0.95, 0.99)) {
    # 检查输入
    if (is.null(simulation_results) || nrow(simulation_results) == 0) {
        stop("simulation_results 为空或无数据")
    }

    # 识别方法列（以_theta_point结尾的列）
    theta_cols <- grep("_theta_point$", names(simulation_results), value = TRUE)
    if (length(theta_cols) == 0) {
        stop("在 simulation_results 中未找到 _theta_point 列")
    }

    # 提取方法标识符（去掉_theta_point后缀）
    methods <- gsub("_theta_point$", "", theta_cols)

    # 初始化结果列表
    summary_stats <- list()

    # 为每个方法计算汇总统计量
    for (method in methods) {
        # 提取该方法的相关列
        theta_col <- paste0(method, "_theta_point")
        se_col <- paste0(method, "_theta_se")
        z_col <- paste0(method, "_z")
        p_col <- paste0(method, "_p_value")
        duration_col <- paste0(method, "_duration")

        # 检查必需的列是否存在
        required_cols <- c(theta_col, se_col, p_col)
        missing_cols <- setdiff(required_cols, names(simulation_results))
        if (length(missing_cols) > 0) {
            warning(paste("方法", method, "缺少必需的列:", paste(missing_cols, collapse = ", ")))
            next
        }

        # 提取数据（移除NA值）
        theta_vals <- simulation_results[[theta_col]]
        se_vals <- simulation_results[[se_col]]
        p_vals <- simulation_results[[p_col]]

        # 创建完整数据的逻辑向量
        complete_cases <- !is.na(theta_vals) & !is.na(se_vals) & !is.na(p_vals)

        if (sum(complete_cases) == 0) {
            warning(paste("方法", method, "没有完整的数据"))
            next
        }

        # 过滤到完整数据
        theta_complete <- theta_vals[complete_cases]
        se_complete <- se_vals[complete_cases]
        p_complete <- p_vals[complete_cases]

        # 初始化该方法的统计量
        method_stats <- list()

        # 1. 基本统计量
        method_stats$n_simulations <- sum(complete_cases)
        method_stats$mean_estimate <- mean(theta_complete)
        method_stats$median_estimate <- median(theta_complete)
        method_stats$sd_estimate <- sd(theta_complete)
        method_stats$mean_se <- mean(se_complete)

        # 2. Bias
        method_stats$bias <- mean(theta_complete) - true_value

        # 3. Power（对于不同的显著性水平）
        for (alpha in alpha_levels) {
            power_name <- paste0("power_", gsub("\\.", "", sprintf("%.3f", alpha)))
            method_stats[[power_name]] <- mean(p_complete < alpha, na.rm = TRUE)
        }

        # 4. Coverage（对于不同的置信水平）
        for (conf_level in confidence_levels) {
            # 计算置信区间
            alpha_ci <- 1 - conf_level
            z_critical <- qnorm(1 - alpha_ci / 2)

            # 置信区间的下界和上界
            ci_lower <- theta_complete - z_critical * se_complete
            ci_upper <- theta_complete + z_critical * se_complete

            # 计算覆盖率（真实值是否在置信区间内）
            coverage <- mean(ci_lower <= true_value & true_value <= ci_upper, na.rm = TRUE)

            coverage_name <- paste0("coverage_", gsub("\\.", "", sprintf("%.2f", conf_level * 100)))
            method_stats[[coverage_name]] <- coverage
        }

        # 5. 运行时间统计（如果有的话）
        if (duration_col %in% names(simulation_results)) {
            duration_vals <- simulation_results[[duration_col]][complete_cases]
            method_stats$mean_duration <- mean(duration_vals, na.rm = TRUE)
            method_stats$median_duration <- median(duration_vals, na.rm = TRUE)
        }

        # 6. 其他有用的统计量
        method_stats$min_estimate <- min(theta_complete)
        method_stats$max_estimate <- max(theta_complete)
        method_stats$q25_estimate <- quantile(theta_complete, 0.25)
        method_stats$q75_estimate <- quantile(theta_complete, 0.75)

        # 7. 相对偏差和相对标准误
        if (true_value != 0) {
            method_stats$relative_bias <- method_stats$bias / abs(true_value)
        }
        method_stats$monte_carlo_se <- method_stats$sd_estimate / sqrt(method_stats$n_simulations)

        # 将该方法的统计量添加到结果中
        summary_stats[[method]] <- method_stats
    }

    # 检查是否有任何成功处理的方法
    if (length(summary_stats) == 0) {
        stop("没有成功处理任何方法的数据")
    }

    # 添加元信息
    attr(summary_stats, "true_value") <- true_value
    attr(summary_stats, "alpha_levels") <- alpha_levels
    attr(summary_stats, "confidence_levels") <- confidence_levels
    attr(summary_stats, "total_simulations") <- nrow(simulation_results)

    return(summary_stats)
}

# 辅助函数：打印汇总统计量的概览
print_summary_stats <- function(summary_stats, digits = 4) {
    cat("仿真结果汇总统计量\n")
    cat("===================\n\n")

    # 获取元信息
    true_value <- attr(summary_stats, "true_value")
    total_sims <- attr(summary_stats, "total_simulations")

    cat("真实值:", true_value, "\n")
    cat("总仿真次数:", total_sims, "\n")
    cat("方法数量:", length(summary_stats), "\n\n")

    # 为每个方法打印关键统计量
    for (method in names(summary_stats)) {
        stats <- summary_stats[[method]]
        cat("方法", method, ":\n")
        cat("  有效仿真次数:", stats$n_simulations, "\n")
        cat("  平均估计值:", round(stats$mean_estimate, digits), "\n")
        cat("  偏差:", round(stats$bias, digits), "\n")

        # Power统计
        power_stats <- names(stats)[grepl("^power_", names(stats))]
        for (p_stat in power_stats) {
            alpha_val <- gsub("power_", "", p_stat)
            alpha_val <- as.numeric(paste0("0.", alpha_val))
            cat("  ", p_stat, " (α=", alpha_val, "):", round(stats[[p_stat]], digits), "\n")
        }

        # Coverage统计
        coverage_stats <- names(stats)[grepl("^coverage_", names(stats))]
        for (c_stat in coverage_stats) {
            conf_val <- gsub("coverage_", "", c_stat)
            conf_val <- as.numeric(conf_val) / 100
            cat("  ", c_stat, " (", conf_val * 100, "% CI):", round(stats[[c_stat]], digits), "\n")
        }

        cat("\n")
    }
}

# 使用示例
if (FALSE) {
    # 运行100次模拟，使用4个核心

    test_1 <- generate_mr_trio_data_matrix_ultra(
        n = 1000, seed = NULL, maf = 0.3,
        # --- 暴露(Exposure)的效应系数 ---
        beta_os_to_oe_exp = rep(0.1, 10), # 子女自身SNP
        beta_fs_to_oe_exp = rep(0.03, 10), # 父亲SNP (遗传叠加)
        beta_ms_to_oe_exp = rep(0.03, 10), # 母亲SNP (遗传叠加)
        sigma_confounding_exp = 0.1, # 混杂U的方差
        sigma_correlation_exp = 0.1, # 共享环境C的方差
        # --- 暴露(Exposure)的效应系数（异质性） ---
        h_beta_os_to_oe_exp = rep(0.1, 10), # 子女自身SNP
        h_beta_fs_to_oe_exp = rep(0.03, 10), # 父亲SNP (遗传叠加)
        h_beta_ms_to_oe_exp = rep(0.03, 10), # 母亲SNP (遗传叠加)
        # --- 结局(Outcome)的效应系数 ---
        beta_os_to_oe_out = c(rep(0, 3), rep(0, 7)), # 水平多效性
        beta_fs_to_oe_out = c(rep(0, 3), rep(0, 7)),
        beta_ms_to_oe_out = c(rep(0, 3), rep(0, 7)),
        sigma_confounding_out = 0.1, # 混杂U对结局的方差
        sigma_correlation_out = 0.1, # 共享环境C对结局的方差
        # --- 直接设定的参数 ---
        beta_exp_to_out = 0, # 因果效应
        sigma_lon_sq = 1, # 暴露的独立残差方差 (基准)
        sigma_psi_sq = 1, # 结局的独立残差方差 (基准)
        # --- 模拟流程控制参数 ---
        assortative_mating_strength = 1000
    )
    test_1$data_exp$variance_explained

    test_1 <- run_parallel_simulation(
        n_simulations = 1000,
        script_path = "统计模拟/MR 模拟（单次）.R",
        n_cores = NULL,
        n = 1000, seed = NULL, maf = 0.3,
        # --- 暴露(Exposure)的效应系数 ---
        beta_os_to_oe_exp = rep(0.1, 10), # 子女自身SNP
        beta_fs_to_oe_exp = rep(0.03, 10), # 父亲SNP (遗传叠加)
        beta_ms_to_oe_exp = rep(0.03, 10), # 母亲SNP (遗传叠加)
        sigma_confounding_exp = 0.1, # 混杂U的方差
        sigma_correlation_exp = 0.1, # 共享环境C的方差
        # --- 暴露(Exposure)的效应系数（异质性） ---
        h_beta_os_to_oe_exp = rep(0.1, 10), # 子女自身SNP
        h_beta_fs_to_oe_exp = rep(0.03, 10), # 父亲SNP (遗传叠加)
        h_beta_ms_to_oe_exp = rep(0.03, 10), # 母亲SNP (遗传叠加)
        # --- 结局(Outcome)的效应系数 ---
        beta_os_to_oe_out = c(rep(0, 3), rep(0, 7)), # 水平多效性
        beta_fs_to_oe_out = c(rep(0, 3), rep(0, 7)),
        beta_ms_to_oe_out = c(rep(0, 3), rep(0, 7)),
        sigma_confounding_out = 0.1, # 混杂U对结局的方差
        sigma_correlation_out = 0.1, # 共享环境C对结局的方差
        # --- 直接设定的参数 ---
        beta_exp_to_out = 0.1, # 因果效应
        sigma_lon_sq = 1, # 暴露的独立残差方差 (基准)
        sigma_psi_sq = 1, # 结局的独立残差方差 (基准)
        # --- 模拟流程控制参数 ---
        assortative_mating_strength = 1000
    )
    test_1_1 <- calculate_summary_stats(test_1$results, true_value = 0.1)
    print_summary_stats(test_1_1)
}

# %% batch 运行总控函数
batch_simulation_controller <- function(
    param_combinations_list,
    script_path,
    output_dir = "统计模拟/统计模拟结果",
    save_individual = TRUE,
    save_combined = TRUE,
    progress_interval = 1,
    timestamp = TRUE,
    verbose = TRUE) {
    # --- 1. 输入验证和初始化 ---
    stopifnot(
        "param_combinations_list 必须是一个列表" = is.list(param_combinations_list),
        "param_combinations_list 不能为空" = length(param_combinations_list) > 0,
        "script_path 必须存在" = file.exists(script_path)
    )

    # 确保输出目录存在
    if (!dir.exists(output_dir)) {
        dir.create(output_dir, recursive = TRUE)
        if (verbose) message(paste("创建输出目录:", output_dir))
    }

    # 生成时间戳
    time_stamp <- if (timestamp) format(Sys.time(), "%Y%m%d_%H%M%S") else ""

    # 初始化结果存储
    n_combinations <- length(param_combinations_list)
    all_results <- vector("list", n_combinations)
    all_metadata <- vector("list", n_combinations)

    if (verbose) {
        message("=== 批量模拟开始 ===")
        message(paste("总共", n_combinations, "个参数组合"))
        message(paste("输出目录:", output_dir))
        message(paste("时间戳:", time_stamp))
    }

    # --- 2. 主循环：逐个运行参数组合 ---
    overall_start_time <- Sys.time()

    for (i in seq_along(param_combinations_list)) {
        if (verbose && (i %% progress_interval == 0 || i == 1)) {
            message(paste(">>> 开始第", i, "/", n_combinations, "个参数组合"))
        }

        combination_start_time <- Sys.time()
        current_params <- param_combinations_list[[i]]

        # 确保每个参数组合包含script_path
        if (!"script_path" %in% names(current_params)) {
            current_params$script_path <- script_path
        }

        tryCatch(
            {
                # 运行单个参数组合的并行模拟
                result <- do.call(run_parallel_simulation, current_params)

                combination_end_time <- Sys.time()
                combination_duration <- as.numeric(combination_end_time - combination_start_time)

                # 存储结果
                all_results[[i]] <- result

                # 存储元数据
                all_metadata[[i]] <- list(
                    combination_id = i,
                    parameters = current_params,
                    start_time = combination_start_time,
                    end_time = combination_end_time,
                    duration_seconds = combination_duration,
                    success_rate = nrow(result$results) / (nrow(result$results) + length(result$errors)),
                    total_simulations = nrow(result$results) + length(result$errors),
                    successful_simulations = nrow(result$results),
                    failed_simulations = length(result$errors)
                )

                # 保存单个结果（如果需要）
                if (save_individual) {
                    individual_filename <- if (timestamp) {
                        paste0("simulation_", time_stamp, "_combination_", sprintf("%03d", i), ".rds")
                    } else {
                        paste0("simulation_combination_", sprintf("%03d", i), ".rds")
                    }
                    individual_filepath <- file.path(output_dir, individual_filename)

                    saveRDS(list(
                        results = result,
                        metadata = all_metadata[[i]]
                    ), individual_filepath)
                }

                if (verbose) {
                    message(paste(
                        "    第", i, "个组合完成 - 耗时:",
                        round(combination_duration, 2), "秒",
                        "- 成功率:", round(all_metadata[[i]]$success_rate * 100, 1), "%"
                    ))
                }
            },
            error = function(e) {
                combination_end_time <- Sys.time()
                combination_duration <- as.numeric(combination_end_time - combination_start_time)

                # 记录失败的组合
                all_results[[i]] <- list(
                    results = tibble::tibble(),
                    errors = list(paste("整个参数组合失败:", e$message))
                )

                all_metadata[[i]] <- list(
                    combination_id = i,
                    parameters = current_params,
                    start_time = combination_start_time,
                    end_time = combination_end_time,
                    duration_seconds = combination_duration,
                    success_rate = 0,
                    total_simulations = 0,
                    successful_simulations = 0,
                    failed_simulations = 1,
                    error_message = e$message
                )

                if (verbose) {
                    message(paste("    第", i, "个组合失败:", e$message))
                }
            }
        )
    }

    overall_end_time <- Sys.time()
    overall_duration <- as.numeric(overall_end_time - overall_start_time)

    # --- 3. 结果汇总和保存 ---
    if (verbose) message("正在汇总所有结果...")

    # 合并所有成功的结果
    all_successful_results <- list()
    for (i in seq_along(all_results)) {
        if (nrow(all_results[[i]]$results) > 0) {
            temp_result <- all_results[[i]]$results %>%
                tibble::add_column(combination_id = i, .before = 1)
            all_successful_results[[i]] <- temp_result
        }
    }

    if (length(all_successful_results) > 0) {
        combined_results <- dplyr::bind_rows(all_successful_results)
    } else {
        combined_results <- tibble::tibble()
    }

    # 创建元数据汇总
    metadata_summary <- dplyr::bind_rows(lapply(all_metadata, function(x) {
        tibble::tibble(
            combination_id = x$combination_id,
            duration_seconds = x$duration_seconds,
            success_rate = x$success_rate,
            total_simulations = x$total_simulations,
            successful_simulations = x$successful_simulations,
            failed_simulations = x$failed_simulations,
            has_error = "error_message" %in% names(x)
        )
    }))

    # 创建最终汇总对象
    final_summary <- list(
        combined_results = combined_results,
        metadata_summary = metadata_summary,
        individual_results = all_results,
        individual_metadata = all_metadata,
        overall_summary = list(
            total_combinations = n_combinations,
            successful_combinations = sum(metadata_summary$successful_simulations > 0),
            total_runtime_seconds = overall_duration,
            total_successful_simulations = sum(metadata_summary$successful_simulations),
            total_failed_simulations = sum(metadata_summary$failed_simulations),
            overall_success_rate = sum(metadata_summary$successful_simulations) /
                (sum(metadata_summary$successful_simulations) + sum(metadata_summary$failed_simulations)),
            start_time = overall_start_time,
            end_time = overall_end_time,
            timestamp = time_stamp
        )
    )

    # 保存合并结果（如果需要）
    if (save_combined) {
        combined_filename <- if (timestamp) {
            paste0("batch_simulation_", time_stamp, "_complete.rds")
        } else {
            "batch_simulation_complete.rds"
        }
        combined_filepath <- file.path(output_dir, combined_filename)
        saveRDS(final_summary, combined_filepath)

        if (verbose) message(paste("合并结果已保存到:", combined_filepath))
    }

    # --- 4. 最终报告 ---
    if (verbose) {
        message("=== 批量模拟完成 ===")
        message(paste("总运行时间:", round(overall_duration / 60, 2), "分钟"))
        message(paste("成功的参数组合:", final_summary$overall_summary$successful_combinations, "/", n_combinations))
        message(paste("总成功模拟次数:", final_summary$overall_summary$total_successful_simulations))
        message(paste("总失败模拟次数:", final_summary$overall_summary$total_failed_simulations))
        message(paste("总体成功率:", round(final_summary$overall_summary$overall_success_rate * 100, 2), "%"))
    }

    return(final_summary)
}

# 辅助函数：创建参数组合列表的便捷函数
create_param_combinations <- function(base_params, varying_params) {
    # base_params: 所有组合共享的基础参数列表
    # varying_params: 命名列表，每个元素是一个向量，表示该参数的不同取值

    # 创建所有参数组合的网格
    param_grid <- expand.grid(varying_params, stringsAsFactors = FALSE)

    # 为每个组合创建完整的参数列表
    combinations <- list()
    for (i in 1:nrow(param_grid)) {
        combination <- base_params
        for (param_name in names(varying_params)) {
            combination[[param_name]] <- param_grid[i, param_name]
        }
        combinations[[i]] <- combination
    }

    return(combinations)
}

if (FALSE) {
    # === 批量模拟使用示例 ===

    # 方法1: 手动创建参数组合列表
    # 这种方法适合完全自定义的参数组合

    param_combinations_manual <- list(
        # 组合1: 低因果效应
        list(
            n_simulations = 1000,
            script_path = "统计模拟/MR 模拟（单次）.R",
            n_cores = NULL,
            n = 1000, maf = 0.3,
            beta_os_to_oe_exp = rep(0.1, 10),
            beta_fs_to_oe_exp = rep(0.03, 10),
            beta_ms_to_oe_exp = rep(0.03, 10),
            sigma_confounding_exp = 0.1,
            sigma_correlation_exp = 0.1,
            h_beta_os_to_oe_exp = rep(0.1, 10),
            h_beta_fs_to_oe_exp = rep(0.03, 10),
            h_beta_ms_to_oe_exp = rep(0.03, 10),
            beta_os_to_oe_out = rep(0, 10),
            beta_fs_to_oe_out = rep(0, 10),
            beta_ms_to_oe_out = rep(0, 10),
            sigma_confounding_out = 0.1,
            sigma_correlation_out = 0.1,
            beta_exp_to_out = 0.1, # 低因果效应
            sigma_lon_sq = 1,
            sigma_psi_sq = 1,
            assortative_mating_strength = 1000
        ),

        # 组合2: 中等因果效应
        list(
            n_simulations = 1000,
            script_path = "统计模拟/MR 模拟（单次）.R",
            n_cores = NULL,
            n = 1000, maf = 0.3,
            beta_os_to_oe_exp = rep(0.1, 10),
            beta_fs_to_oe_exp = rep(0.03, 10),
            beta_ms_to_oe_exp = rep(0.03, 10),
            sigma_confounding_exp = 0.1,
            sigma_correlation_exp = 0.1,
            h_beta_os_to_oe_exp = rep(0.1, 10),
            h_beta_fs_to_oe_exp = rep(0.03, 10),
            h_beta_ms_to_oe_exp = rep(0.03, 10),
            beta_os_to_oe_out = rep(0, 10),
            beta_fs_to_oe_out = rep(0, 10),
            beta_ms_to_oe_out = rep(0, 10),
            sigma_confounding_out = 0.1,
            sigma_correlation_out = 0.1,
            beta_exp_to_out = 0.3, # 中等因果效应
            sigma_lon_sq = 1,
            sigma_psi_sq = 1,
            assortative_mating_strength = 1000
        ),

        # 组合3: 高因果效应
        list(
            n_simulations = 1000,
            script_path = "统计模拟/MR 模拟（单次）.R",
            n_cores = NULL,
            n = 1000, maf = 0.3,
            beta_os_to_oe_exp = rep(0.1, 10),
            beta_fs_to_oe_exp = rep(0.03, 10),
            beta_ms_to_oe_exp = rep(0.03, 10),
            sigma_confounding_exp = 0.1,
            sigma_correlation_exp = 0.1,
            h_beta_os_to_oe_exp = rep(0.1, 10),
            h_beta_fs_to_oe_exp = rep(0.03, 10),
            h_beta_ms_to_oe_exp = rep(0.03, 10),
            beta_os_to_oe_out = rep(0, 10),
            beta_fs_to_oe_out = rep(0, 10),
            beta_ms_to_oe_out = rep(0, 10),
            sigma_confounding_out = 0.1,
            sigma_correlation_out = 0.1,
            beta_exp_to_out = 0.5, # 高因果效应
            sigma_lon_sq = 1,
            sigma_psi_sq = 1,
            assortative_mating_strength = 1000
        )
    )

    # 方法2: 使用辅助函数自动生成参数组合
    # 这种方法适合系统性的参数扫描

    # 定义基础参数（所有组合共享）
    base_params <- list(
        n_simulations = 1000,
        script_path = "统计模拟/MR 模拟（单次）.R",
        n = 1000,
        maf = 0.3,
        beta_os_to_oe_exp = rep(0.1, 10),
        beta_fs_to_oe_exp = rep(0.03, 10),
        beta_ms_to_oe_exp = rep(0.03, 10),
        sigma_confounding_exp = 0.1,
        sigma_correlation_exp = 0.1,
        h_beta_os_to_oe_exp = rep(0.1, 10),
        h_beta_fs_to_oe_exp = rep(0.03, 10),
        h_beta_ms_to_oe_exp = rep(0.03, 10),
        beta_os_to_oe_out = rep(0, 10),
        beta_fs_to_oe_out = rep(0, 10),
        beta_ms_to_oe_out = rep(0, 10),
        sigma_confounding_out = 0.1,
        sigma_correlation_out = 0.1,
        sigma_lon_sq = 1,
        sigma_psi_sq = 1,
        assortative_mating_strength = 1000
    )

    # 定义变化的参数
    varying_params <- list(
        beta_exp_to_out = c(0.1, 0.2, 0.3, 0.4, 0.5), # 因果效应
        n = c(500, 1000, 2000), # 样本量
        maf = c(0.2, 0.3, 0.4) # 次等位基因频率
    )

    # 自动生成所有参数组合
    param_combinations_auto <- create_param_combinations(base_params, varying_params)

    # 这将创建 5 × 3 × 3 = 45 个参数组合

    # === 运行批量模拟 ===

    # 选择要使用的参数组合（这里用手动创建的3个组合作为示例）
    results <- batch_simulation_controller(
        param_combinations_list = param_combinations_manual,
        script_path = "统计模拟/MR 模拟（单次）.R",
        save_individual = TRUE, # 保存每个组合的单独结果
        save_combined = TRUE, # 保存合并的结果
        progress_interval = 1, # 每个组合都显示进度
        timestamp = TRUE, # 添加时间戳
        verbose = TRUE # 显示详细信息
    )

    # === 结果分析示例 ===

    # 1. 查看整体摘要
    print(results$overall_summary)

    # 2. 查看每个组合的元数据
    print(results$metadata_summary)

    # 3. 分析合并后的结果
    if (nrow(results$combined_results) > 0) {
        library(dplyr)

        # 按组合计算各方法的平均估计值
        method_performance <- results$combined_results %>%
            group_by(combination_id) %>%
            summarise(
                mean_a_theta = mean(a_theta_point, na.rm = TRUE),
                mean_b_theta = mean(b_theta_point, na.rm = TRUE),
                mean_c_theta = mean(c_theta_point, na.rm = TRUE),
                .groups = "drop"
            )

        print("各组合的方法性能:")
        print(method_performance)

        # 计算估计偏差（假设真实值已知）
        true_effects <- c(0.1, 0.3, 0.5) # 对应三个组合的真实因果效应

        bias_analysis <- results$combined_results %>%
            group_by(combination_id) %>%
            summarise(
                bias_a = mean(a_theta_point - true_effects[combination_id], na.rm = TRUE),
                bias_b = mean(b_theta_point - true_effects[combination_id], na.rm = TRUE),
                bias_c = mean(c_theta_point - true_effects[combination_id], na.rm = TRUE),
                .groups = "drop"
            )

        print("各方法的偏差分析:")
        print(bias_analysis)
    }

    # === 加载已保存的结果 ===

    # 如果需要重新加载结果
    # loaded_results <- readRDS("batch_simulation_results/batch_simulation_YYYYMMDD_HHMMSS_complete.rds")

    # 加载单个组合的结果
    # individual_result <- readRDS("batch_simulation_results/simulation_YYYYMMDD_HHMMSS_combination_001.rds")

    # === 高级使用：大规模参数扫描 ===

    # 如果要运行大规模参数扫描（比如上面的45个组合），可以这样做：
    # large_scale_results <- batch_simulation_controller(
    #     param_combinations_list = param_combinations_auto[1:10],  # 先运行前10个作为测试
    #     script_path = "统计模拟/MR 模拟（单次）.R",
    #     output_dir = "large_scale_simulation",
    #     save_individual = TRUE,
    #     save_combined = TRUE,
    #     progress_interval = 2,  # 每2个组合显示一次进度
    #     verbose = TRUE
    # )
}

# %% 运行
