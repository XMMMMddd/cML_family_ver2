library(data.table)
library(ggplot2)
library(ggdist)

# %% 可视化函数

# MR方法性能指标可视化函数
plot_mr_performance <- function(summary_stats,
                                method_names = NULL,
                                colors = NULL,
                                save_plot = FALSE,
                                output_dir = "plots",
                                width = 12,
                                height = 8) {
    # 加载必需的包
    required_packages <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)

    # 获取方法列表
    methods <- names(summary_stats)
    n_methods <- length(methods)

    # 设置方法名称
    if (is.null(method_names)) {
        method_labels <- setNames(methods, methods) # 直接使用a, b, c...
    } else {
        # 检查命名向量
        if (is.null(names(method_names))) {
            stop("method_names 必须是一个命名向量，名称对应方法标识符（a, b, c, ...）")
        }
        method_labels <- method_names
        # 确保所有方法都有标签
        missing_methods <- setdiff(methods, names(method_labels))
        if (length(missing_methods) > 0) {
            warning("以下方法缺少标签，将使用默认标签: ", paste(missing_methods, collapse = ", "))
            for (m in missing_methods) {
                method_labels[m] <- m # 直接使用字母标识符
            }
        }
    }

    # 设置颜色
    if (is.null(colors)) {
        if (n_methods <= 9) {
            colors <- RColorBrewer::brewer.pal(max(3, n_methods), "Set1")[1:n_methods]
        } else {
            colors <- rainbow(n_methods)
        }
    } else {
        if (length(colors) < n_methods) {
            warning("颜色数量不足，将重复使用颜色")
            colors <- rep(colors, length.out = n_methods)
        }
    }

    # 提取数据并转换为长格式
    extract_metric <- function(metric_name) {
        values <- sapply(methods, function(m) {
            val <- summary_stats[[m]][[metric_name]]
            if (is.null(val) || is.na(val)) {
                return(NA)
            }
            return(val)
        })

        data.frame(
            Method = methods,
            Method_Label = method_labels[methods],
            Metric = metric_name,
            Value = values,
            stringsAsFactors = FALSE
        )
    }

    # 提取所有指标
    metrics_to_plot <- c("power_005", "power_001", "coverage_95", "coverage_99", "bias", "rmse")
    metric_labels <- c(
        "power_005" = "Power (α=0.05)",
        "power_001" = "Power (α=0.01)",
        "coverage_95" = "Coverage (95% CI)",
        "coverage_99" = "Coverage (99% CI)",
        "bias" = "Bias",
        "rmse" = "RMSE"
    )

    plot_data <- do.call(rbind, lapply(metrics_to_plot, extract_metric))
    plot_data$Metric_Label <- metric_labels[plot_data$Metric]

    # 移除NA值
    plot_data <- plot_data[!is.na(plot_data$Value), ]

    if (nrow(plot_data) == 0) {
        stop("没有有效的数据可以绘图")
    }

    # 创建绘图函数
    create_plot <- function(data, title_suffix = "") {
        # 基本的条形图 - 调整柱子宽度和透明度
        p <- ggplot(data, aes(x = Method_Label, y = Value, fill = Method)) +
            geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) + # 宽度调回0.9
            scale_fill_manual(values = setNames(colors, methods))

        # 添加方法名称标签在柱子上方（水平显示，调小字体）
        p <- p + geom_text(
            aes(label = Method_Label),
            vjust = -0.2,
            size = 2.8, # 从3.2调小到2.8
            fontface = "bold",
            color = "black"
        )

        p <- p + facet_wrap(~Metric_Label, scales = "free_y", ncol = 3) +
            labs(
                title = paste("MR Methods Performance Comparison", title_suffix),
                x = "Method",
                y = "Value",
                fill = "Method"
            ) +
            theme_minimal() +
            theme(
                axis.text.x = element_blank(), # 隐藏x轴标签，因为方法名在柱子上方
                axis.ticks.x = element_blank(), # 隐藏x轴刻度
                axis.text.y = element_text(size = 9), # 从10调小到9
                axis.title = element_text(size = 11, face = "bold"), # 从12调小到11
                plot.title = element_text(size = 13, face = "bold", hjust = 0.5), # 从14调小到13
                strip.text = element_text(size = 10, face = "bold"), # 从11调小到10
                legend.position = "bottom",
                legend.title = element_text(size = 10, face = "bold"), # 从11调小到10
                legend.text = element_text(size = 9), # 从10调小到9
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = "white", color = "grey90"),
                strip.background = element_rect(fill = "grey95", color = "grey90")
            ) +
            guides(fill = guide_legend(nrow = 2))

        # 调整y轴范围以容纳方法名称标签
        p <- p + expand_limits(y = max(data$Value, na.rm = TRUE) * 1.2)

        return(p)
    }

    # 保存结果的列表
    plot_list <- list()

    # 为每个指标单独创建图
    for (metric in metrics_to_plot) {
        metric_data <- plot_data[plot_data$Metric == metric, ]

        if (nrow(metric_data) > 0) {
            metric_label <- metric_labels[metric]

            p_single <- ggplot(metric_data, aes(x = Method_Label, y = Value, fill = Method)) +
                geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
                scale_fill_manual(values = setNames(colors, methods)) +
                geom_text(
                    aes(label = Method_Label),
                    vjust = -0.2,
                    size = 2.8,
                    fontface = "bold",
                    color = "black"
                ) +
                labs(
                    title = metric_label,
                    x = "Method",
                    y = "Value",
                    fill = "Method"
                ) +
                theme_minimal() +
                theme(
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.text.y = element_text(size = 9),
                    axis.title = element_text(size = 11, face = "bold"),
                    plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
                    legend.position = "bottom",
                    legend.title = element_text(size = 10, face = "bold"),
                    legend.text = element_text(size = 9),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "white", color = "grey90")
                ) +
                guides(fill = guide_legend(nrow = 2)) +
                expand_limits(y = max(metric_data$Value, na.rm = TRUE) * 1.2)

            print(p_single)
            plot_list[[metric]] <- p_single
        }
    }

    # 创建包含所有指标的综合图
    p_all <- create_plot(plot_data, "- All Metrics")
    print(p_all)
    plot_list[["all_metrics"]] <- p_all

    # 保存图片
    if (save_plot) {
        if (!dir.exists(output_dir)) {
            dir.create(output_dir, recursive = TRUE)
        }

        # 保存每个指标的单独图片
        for (metric in names(plot_list)) {
            if (metric != "all_metrics") {
                metric_label_clean <- gsub("[^A-Za-z0-9_]", "_", metric_labels[metric])
                filename <- paste0("mr_performance_", metric, ".png")
                ggsave(file.path(output_dir, filename),
                    plot_list[[metric]],
                    width = width, height = height, dpi = 300
                )
            }
        }

        # 保存综合图
        if ("all_metrics" %in% names(plot_list)) {
            ggsave(file.path(output_dir, "mr_performance_all_metrics.png"),
                plot_list[["all_metrics"]],
                width = width, height = height, dpi = 300
            )
        }

        cat("图片已保存到:", output_dir, "\n")
        cat("单独指标图片:", paste0("mr_performance_", metrics_to_plot, ".png"), "\n")
        cat("综合图片: mr_performance_all_metrics.png\n")
    }

    # 返回绘图数据和图片对象，以便进一步分析
    return(list(
        plot_data = plot_data,
        method_labels = method_labels,
        colors_used = setNames(colors, methods),
        plots = plot_list # 返回所有图片对象
    ))
}

# 使用示例函数
example_usage <- function() {
    cat("使用示例:\n")
    cat("# 1. 基本使用\n")
    cat("plot_mr_performance(summary_stats)\n\n")

    cat("# 2. 自定义方法名称和保存单独指标图\n")
    cat("method_names <- c(\n")
    cat("  'a' = 'IVW',\n")
    cat("  'b' = 'Weighted Median',\n")
    cat("  'c' = 'MR-Egger',\n")
    cat("  'd' = 'MR-PRESSO',\n")
    cat("  'e' = 'Contamination Mixture',\n")
    cat("  'f' = 'Debiased IVW',\n")
    cat("  'g' = 'Robust IVW',\n")
    cat("  'h' = 'Penalized IVW',\n")
    cat("  'i' = 'MR-Lasso'\n")
    cat(")\n")
    cat("result <- plot_mr_performance(summary_stats, method_names = method_names, save_plot = TRUE)\n\n")

    cat("# 3. 访问单独的图片\n")
    cat("# 单个指标图片:\n")
    cat("power_005_plot <- result$plots$power_005\n")
    cat("coverage_95_plot <- result$plots$coverage_95\n")
    cat("bias_plot <- result$plots$bias\n")
    cat("rmse_plot <- result$plots$rmse\n")
    cat("# 综合图片:\n")
    cat("all_metrics_plot <- result$plots$all_metrics\n")
}


plot_bias <- function(summary_stats,
                      method_names = NULL,
                      colors = NULL,
                      label_size = 2.8, # 方法名字体大小参数
                      highlight_methods = NULL, # 需要突出显示的方法
                      show_legend = FALSE, # 新增：是否显示图例
                      save_plot = FALSE,
                      output_file = "bias_plot.png",
                      width = 10,
                      height = 6) {
    # 加载必需的包
    required_packages <- c("ggplot2", "RColorBrewer")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    library(ggplot2)
    library(RColorBrewer)

    # 获取方法列表 - 按照method_names的顺序排列
    if (is.null(method_names)) {
        methods <- names(summary_stats) # 如果没有指定method_names，使用原始顺序
        method_labels <- setNames(methods, methods)
    } else {
        # 使用method_names的顺序作为基准
        methods <- names(method_names)
        method_labels <- method_names

        # 检查method_names中的所有方法是否都存在于summary_stats中
        missing_methods <- setdiff(methods, names(summary_stats))
        if (length(missing_methods) > 0) {
            warning("以下方法在summary_stats中不存在: ", paste(missing_methods, collapse = ", "))
            # 只保留存在的方法，但保持原有顺序
            methods <- methods[methods %in% names(summary_stats)]
            method_labels <- method_labels[methods]
        }
    }

    n_methods <- length(methods)

    # 设置颜色
    if (is.null(colors)) {
        if (n_methods <= 9) {
            colors <- RColorBrewer::brewer.pal(max(3, n_methods), "Set1")[1:n_methods]
        } else {
            colors <- rainbow(n_methods)
        }
    } else {
        if (length(colors) < n_methods) {
            warning("颜色数量不足，将重复使用颜色")
            colors <- rep(colors, length.out = n_methods)
        }
    }

    # 处理突出显示的方法
    if (!is.null(highlight_methods)) {
        # 检查highlight_methods是否有效
        invalid_methods <- setdiff(highlight_methods, methods)
        if (length(invalid_methods) > 0) {
            warning("以下突出显示的方法不存在: ", paste(invalid_methods, collapse = ", "))
            highlight_methods <- intersect(highlight_methods, methods)
        }

        if (length(highlight_methods) > 0) {
            # 创建新的颜色方案：突出显示的方法保持彩色，其他变灰
            highlight_colors <- colors[1:length(highlight_methods)] # 为突出显示的方法分配彩色
            gray_color <- "grey70" # 灰色

            final_colors <- rep(gray_color, n_methods)
            names(final_colors) <- methods

            # 为突出显示的方法分配彩色
            for (i in seq_along(highlight_methods)) {
                method <- highlight_methods[i]
                final_colors[method] <- highlight_colors[i]
            }

            colors <- final_colors
        }
    } else {
        # 如果没有指定突出显示的方法，使用正常颜色方案
        names(colors) <- methods
    }

    # 提取bias数据 - 按照指定的方法顺序
    bias_values <- sapply(methods, function(m) {
        val <- summary_stats[[m]][["bias"]]
        if (is.null(val) || is.na(val)) {
            return(NA)
        }
        return(val)
    })

    # 创建数据框 - 保持方法顺序
    plot_data <- data.frame(
        Method = factor(methods, levels = methods), # 使用factor保持顺序
        Method_Label = method_labels[methods],
        Bias = bias_values,
        stringsAsFactors = FALSE
    )

    # 移除NA值，但保持因子水平顺序
    plot_data <- plot_data[!is.na(plot_data$Bias), ]
    plot_data$Method <- factor(plot_data$Method, levels = methods[methods %in% plot_data$Method])

    if (nrow(plot_data) == 0) {
        stop("没有有效的bias数据可以绘图")
    }

    # 创建bias图
    p <- ggplot(plot_data, aes(x = Method_Label, y = Bias, fill = Method)) +
        geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
        scale_fill_manual(values = colors) +
        scale_x_discrete(limits = plot_data$Method_Label) + # 确保x轴顺序正确
        geom_text(
            aes(label = Method_Label),
            vjust = ifelse(plot_data$Bias >= 0, -0.2, 1.2), # 正值标签在上方，负值标签在下方
            size = label_size, # 使用参数控制字体大小
            fontface = "bold",
            color = "black"
        ) +
        labs(
            title = "Bias Comparison Across MR Methods",
            x = "Method",
            y = "Bias",
            fill = "Method"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 11, face = "bold"),
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
            legend.position = if (show_legend) "bottom" else "none", # 根据参数控制图例显示
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = "grey90")
        )

    # 只有在显示图例时才添加图例设置
    if (show_legend) {
        p <- p + guides(fill = guide_legend(nrow = 2))
    } else {
        p <- p + guides(fill = "none")
    }

    # 调整y轴范围，考虑正负值
    y_max <- max(abs(plot_data$Bias), na.rm = TRUE) * 1.2
    p <- p + ylim(-y_max, y_max)

    # 显示图片
    print(p)

    # 保存图片
    if (save_plot) {
        ggsave(output_file, p, width = width, height = height, dpi = 300)
        cat("Bias图已保存为:", output_file, "\n")
    }

    # 返回图片对象和数据
    return(list(
        plot = p,
        data = plot_data,
        method_labels = method_labels
    ))
}
# 专门绘制Power指标的函数
plot_power <- function(summary_stats,
                       method_names = NULL,
                       colors = NULL,
                       label_size = 2.8, # 方法名字体大小参数
                       highlight_methods = NULL, # 需要突出显示的方法
                       show_legend = TRUE, # 是否显示图例
                       power_type = "both", # 显示哪种power: "005", "001", "both"
                       plot_type = "power", # 绘图类型 "power" 或 "type1error"
                       acceptable_range = c(0.0365, 0.0635), # 新增：可接受范围，默认为(3.65%, 6.35%)
                       save_plot = FALSE,
                       output_file = "power_plot.png",
                       width = 10,
                       height = 6) {
    # 参数验证
    if (plot_type == "type1error" && !is.null(acceptable_range)) {
        if (!is.numeric(acceptable_range) || length(acceptable_range) != 2) {
            stop("acceptable_range 必须是包含两个数值的向量 c(下限, 上限)")
        }
        if (acceptable_range[1] >= acceptable_range[2]) {
            stop("acceptable_range 的下限必须小于上限")
        }
        if (any(acceptable_range < 0) || any(acceptable_range > 1)) {
            stop("acceptable_range 的值必须在 0 和 1 之间")
        }
    }

    # 加载必需的包
    required_packages <- c("ggplot2", "RColorBrewer")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    library(ggplot2)
    library(RColorBrewer)

    # 获取方法列表 - 按照method_names的顺序排列
    if (is.null(method_names)) {
        methods <- names(summary_stats) # 如果没有指定method_names，使用原始顺序
        method_labels <- setNames(methods, methods)
    } else {
        # 使用method_names的顺序作为基准
        methods <- names(method_names)
        method_labels <- method_names

        # 检查method_names中的所有方法是否都存在于summary_stats中
        missing_methods <- setdiff(methods, names(summary_stats))
        if (length(missing_methods) > 0) {
            warning("以下方法在summary_stats中不存在: ", paste(missing_methods, collapse = ", "))
            # 只保留存在的方法，但保持原有顺序
            methods <- methods[methods %in% names(summary_stats)]
            method_labels <- method_labels[methods]
        }
    }

    n_methods <- length(methods)

    # 设置颜色
    if (is.null(colors)) {
        if (n_methods <= 9) {
            colors <- RColorBrewer::brewer.pal(max(3, n_methods), "Set1")[1:n_methods]
        } else {
            colors <- rainbow(n_methods)
        }
    } else {
        if (length(colors) < n_methods) {
            warning("颜色数量不足，将重复使用颜色")
            colors <- rep(colors, length.out = n_methods)
        }
    }

    # 处理突出显示的方法
    if (!is.null(highlight_methods)) {
        # 检查highlight_methods是否有效
        invalid_methods <- setdiff(highlight_methods, methods)
        if (length(invalid_methods) > 0) {
            warning("以下突出显示的方法不存在: ", paste(invalid_methods, collapse = ", "))
            highlight_methods <- intersect(highlight_methods, methods)
        }

        if (length(highlight_methods) > 0) {
            # 创建新的颜色方案：突出显示的方法保持彩色，其他变灰
            highlight_colors <- colors[1:length(highlight_methods)] # 为突出显示的方法分配彩色
            gray_color <- "grey70" # 灰色

            final_colors <- rep(gray_color, n_methods)
            names(final_colors) <- methods

            # 为突出显示的方法分配彩色
            for (i in seq_along(highlight_methods)) {
                method <- highlight_methods[i]
                final_colors[method] <- highlight_colors[i]
            }

            colors <- final_colors
        }
    } else {
        # 如果没有指定突出显示的方法，使用正常颜色方案
        names(colors) <- methods
    }

    # 提取数据，根据plot_type决定提取power还是type I error
    metrics <- c()
    labels <- c()
    y_axis_label <- ""
    plot_title_base <- ""

    if (plot_type == "power") {
        # Power数据
        if (power_type %in% c("005", "both")) {
            metrics <- c(metrics, "power_005")
            labels <- c(labels, "Power (α=0.05)")
        }

        if (power_type %in% c("001", "both")) {
            metrics <- c(metrics, "power_001")
            labels <- c(labels, "Power (α=0.01)")
        }

        y_axis_label <- "Power"
        plot_title_base <- "Power Comparison Across MR Methods"
    } else if (plot_type == "type1error") {
        # 一类错误数据 - 在beta_exp_to_out=0时，power实际上就是一类错误率
        if (power_type %in% c("005", "both")) {
            metrics <- c(metrics, "power_005")
            labels <- c(labels, "Type I Error (α=0.05)")
        }

        if (power_type %in% c("001", "both")) {
            metrics <- c(metrics, "power_001")
            labels <- c(labels, "Type I Error (α=0.01)")
        }

        y_axis_label <- "Type I Error Rate"
        plot_title_base <- "Type I Error Rate Comparison Across MR Methods"
    } else {
        stop("plot_type 必须是 'power' 或 'type1error'")
    }

    if (length(metrics) == 0) {
        stop("power_type 必须是 '005', '001', 或 'both'")
    }

    # 创建数据框 - 按照指定的方法顺序
    plot_data_list <- list()

    for (i in seq_along(metrics)) {
        metric <- metrics[i]
        label <- labels[i]

        values <- sapply(methods, function(m) {
            val <- summary_stats[[m]][[metric]]
            if (is.null(val) || is.na(val)) {
                return(NA)
            }
            return(val)
        })

        temp_data <- data.frame(
            Method = factor(methods, levels = methods), # 使用factor保持顺序
            Method_Label = method_labels[methods],
            Value = values,
            Metric_Type = label,
            stringsAsFactors = FALSE
        )

        plot_data_list[[i]] <- temp_data
    }

    # 合并数据
    plot_data <- do.call(rbind, plot_data_list)

    # 移除NA值，但保持因子水平顺序
    plot_data <- plot_data[!is.na(plot_data$Value), ]
    plot_data$Method <- factor(plot_data$Method, levels = methods[methods %in% unique(plot_data$Method)])

    if (nrow(plot_data) == 0) {
        stop(paste("没有有效的", y_axis_label, "数据可以绘图"))
    }

    # 创建图表
    if (power_type == "both") {
        # 如果显示两种类型，使用facet
        p <- ggplot(plot_data, aes(x = Method_Label, y = Value, fill = Method)) +
            geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
            scale_fill_manual(values = colors) +
            scale_x_discrete(limits = unique(plot_data$Method_Label[order(as.numeric(plot_data$Method))])) + # 确保x轴顺序正确
            geom_text(
                aes(label = Method_Label),
                vjust = -0.2,
                size = label_size,
                fontface = "bold",
                color = "black"
            ) +
            facet_wrap(~Metric_Type, scales = "free_y") +
            labs(
                title = plot_title_base,
                x = "Method",
                y = y_axis_label,
                fill = "Method"
            ) +
            ylim(0, 1.2) # 值在0-1之间，留出标签空间
    } else {
        # 如果只显示一种类型，不使用facet
        p <- ggplot(plot_data, aes(x = Method_Label, y = Value, fill = Method)) +
            geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
            scale_fill_manual(values = colors) +
            scale_x_discrete(limits = plot_data$Method_Label[order(as.numeric(plot_data$Method))]) + # 确保x轴顺序正确
            geom_text(
                aes(label = Method_Label),
                vjust = -0.2,
                size = label_size,
                fontface = "bold",
                color = "black"
            ) +
            labs(
                title = paste(plot_title_base, "-", unique(plot_data$Metric_Type)),
                x = "Method",
                y = y_axis_label,
                fill = "Method"
            ) +
            ylim(0, 1.2) # 值在0-1之间，留出标签空间
    }

    # 添加参考线和区域（根据plot_type）
    if (plot_type == "type1error" && !is.null(acceptable_range)) {
        # 对于一类错误，添加可接受范围的绿色区域
        p <- p + annotate("rect",
            xmin = -Inf, xmax = Inf,
            ymin = acceptable_range[1], ymax = acceptable_range[2],
            fill = "green", alpha = 0.15
        ) # 低饱和度绿色

        # 添加可接受范围的边界线（深绿色虚线）
        p <- p + geom_hline(
            yintercept = acceptable_range[1],
            linetype = "dashed", color = "darkgreen", size = 1, alpha = 0.8
        ) # 下边界
        p <- p + geom_hline(
            yintercept = acceptable_range[2],
            linetype = "dashed", color = "darkgreen", size = 1, alpha = 0.8
        ) # 上边界
    }

    # 添加主题
    p <- p + theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 11, face = "bold"),
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 10, face = "bold"),
            legend.position = if (show_legend) "bottom" else "none",
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = "grey90"),
            strip.background = element_rect(fill = "grey95", color = "grey90")
        )

    # 只有在显示图例时才添加图例设置
    if (show_legend) {
        p <- p + guides(fill = guide_legend(nrow = 2))
    } else {
        p <- p + guides(fill = "none")
    }

    # 显示图片
    print(p)

    # 保存图片
    if (save_plot) {
        ggsave(output_file, p, width = width, height = height, dpi = 300)
        cat(paste(y_axis_label, "图已保存为:", output_file, "\n"))
    }

    # 返回图片对象和数据
    return(list(
        plot = p,
        data = plot_data,
        method_labels = method_labels
    ))
}
# 专门绘制Coverage指标的函数
plot_coverage <- function(summary_stats,
                          method_names = NULL,
                          colors = NULL,
                          label_size = 2.8, # 方法名字体大小参数
                          highlight_methods = NULL, # 需要突出显示的方法
                          show_legend = TRUE, # 是否显示图例
                          coverage_type = "both", # 显示哪种coverage: "95", "99", "both"
                          acceptable_range = c(0.9253, 0.9547), # 可接受范围，默认为95%CI的可接受范围
                          save_plot = FALSE,
                          output_file = "coverage_plot.png",
                          width = 10,
                          height = 6) {
    # 参数验证
    if (!is.null(acceptable_range)) {
        if (!is.numeric(acceptable_range) || length(acceptable_range) != 2) {
            stop("acceptable_range 必须是包含两个数值的向量 c(下限, 上限)")
        }
        if (acceptable_range[1] >= acceptable_range[2]) {
            stop("acceptable_range 的下限必须小于上限")
        }
        if (any(acceptable_range < 0) || any(acceptable_range > 1)) {
            stop("acceptable_range 的值必须在 0 和 1 之间")
        }
    }

    # 加载必需的包
    required_packages <- c("ggplot2", "RColorBrewer")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    library(ggplot2)
    library(RColorBrewer)

    # 获取方法列表 - 按照method_names的顺序排列
    if (is.null(method_names)) {
        methods <- names(summary_stats) # 如果没有指定method_names，使用原始顺序
        method_labels <- setNames(methods, methods)
    } else {
        # 使用method_names的顺序作为基准
        methods <- names(method_names)
        method_labels <- method_names

        # 检查method_names中的所有方法是否都存在于summary_stats中
        missing_methods <- setdiff(methods, names(summary_stats))
        if (length(missing_methods) > 0) {
            warning("以下方法在summary_stats中不存在: ", paste(missing_methods, collapse = ", "))
            # 只保留存在的方法，但保持原有顺序
            methods <- methods[methods %in% names(summary_stats)]
            method_labels <- method_labels[methods]
        }
    }

    n_methods <- length(methods)

    # 设置颜色
    if (is.null(colors)) {
        if (n_methods <= 9) {
            colors <- RColorBrewer::brewer.pal(max(3, n_methods), "Set1")[1:n_methods]
        } else {
            colors <- rainbow(n_methods)
        }
    } else {
        if (length(colors) < n_methods) {
            warning("颜色数量不足，将重复使用颜色")
            colors <- rep(colors, length.out = n_methods)
        }
    }

    # 处理突出显示的方法
    if (!is.null(highlight_methods)) {
        # 检查highlight_methods是否有效
        invalid_methods <- setdiff(highlight_methods, methods)
        if (length(invalid_methods) > 0) {
            warning("以下突出显示的方法不存在: ", paste(invalid_methods, collapse = ", "))
            highlight_methods <- intersect(highlight_methods, methods)
        }

        if (length(highlight_methods) > 0) {
            # 创建新的颜色方案：突出显示的方法保持彩色，其他变灰
            highlight_colors <- colors[1:length(highlight_methods)] # 为突出显示的方法分配彩色
            gray_color <- "grey70" # 灰色

            final_colors <- rep(gray_color, n_methods)
            names(final_colors) <- methods

            # 为突出显示的方法分配彩色
            for (i in seq_along(highlight_methods)) {
                method <- highlight_methods[i]
                final_colors[method] <- highlight_colors[i]
            }

            colors <- final_colors
        }
    } else {
        # 如果没有指定突出显示的方法，使用正常颜色方案
        names(colors) <- methods
    }

    # 提取coverage数据
    metrics <- c()
    labels <- c()

    if (coverage_type %in% c("95", "both")) {
        metrics <- c(metrics, "coverage_95")
        labels <- c(labels, "Coverage (95% CI)")
    }

    if (coverage_type %in% c("99", "both")) {
        metrics <- c(metrics, "coverage_99")
        labels <- c(labels, "Coverage (99% CI)")
    }

    if (length(metrics) == 0) {
        stop("coverage_type 必须是 '95', '99', 或 'both'")
    }

    # 创建数据框 - 按照指定的方法顺序
    plot_data_list <- list()

    for (i in seq_along(metrics)) {
        metric <- metrics[i]
        label <- labels[i]

        values <- sapply(methods, function(m) {
            val <- summary_stats[[m]][[metric]]
            if (is.null(val) || is.na(val)) {
                return(NA)
            }
            return(val)
        })

        temp_data <- data.frame(
            Method = factor(methods, levels = methods), # 使用factor保持顺序
            Method_Label = method_labels[methods],
            Coverage = values,
            Coverage_Type = label,
            stringsAsFactors = FALSE
        )

        plot_data_list[[i]] <- temp_data
    }

    # 合并数据
    plot_data <- do.call(rbind, plot_data_list)

    # 移除NA值，但保持因子水平顺序
    plot_data <- plot_data[!is.na(plot_data$Coverage), ]
    plot_data$Method <- factor(plot_data$Method, levels = methods[methods %in% unique(plot_data$Method)])

    if (nrow(plot_data) == 0) {
        stop("没有有效的coverage数据可以绘图")
    }

    # 创建coverage图
    if (coverage_type == "both") {
        # 如果显示两种coverage，使用facet
        p <- ggplot(plot_data, aes(x = Method_Label, y = Coverage, fill = Method)) +
            geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
            scale_fill_manual(values = colors) +
            scale_x_discrete(limits = unique(plot_data$Method_Label[order(as.numeric(plot_data$Method))])) + # 确保x轴顺序正确
            geom_text(
                aes(label = Method_Label),
                vjust = -0.2,
                size = label_size,
                fontface = "bold",
                color = "black"
            ) +
            facet_wrap(~Coverage_Type, scales = "free_y") +
            labs(
                title = "Coverage Comparison Across MR Methods",
                x = "Method",
                y = "Coverage",
                fill = "Method"
            )
    } else {
        # 如果只显示一种coverage，不使用facet
        p <- ggplot(plot_data, aes(x = Method_Label, y = Coverage, fill = Method)) +
            geom_col(color = "white", size = 0.5, alpha = 0.6, width = 0.9) +
            scale_fill_manual(values = colors) +
            scale_x_discrete(limits = plot_data$Method_Label[order(as.numeric(plot_data$Method))]) + # 确保x轴顺序正确
            geom_text(
                aes(label = Method_Label),
                vjust = -0.2,
                size = label_size,
                fontface = "bold",
                color = "black"
            ) +
            labs(
                title = paste("Coverage Comparison Across MR Methods -", unique(plot_data$Coverage_Type)),
                x = "Method",
                y = "Coverage",
                fill = "Method"
            )
    }

    # 智能设置y轴范围
    min_coverage <- min(plot_data$Coverage, na.rm = TRUE)
    max_coverage <- max(plot_data$Coverage, na.rm = TRUE)

    # 如果有可接受范围，也考虑进去
    if (!is.null(acceptable_range)) {
        min_coverage <- min(min_coverage, acceptable_range[1])
        max_coverage <- max(max_coverage, acceptable_range[2])
    }

    # 设置合理的y轴范围（确保柱子和标签都能显示）
    y_range <- max_coverage - min_coverage
    y_buffer_bottom <- y_range * 0.05 # 底部5%缓冲
    y_buffer_top <- y_range * 0.25 # 顶部25%缓冲（为标签留更多空间）

    y_min <- max(0, min_coverage - y_buffer_bottom) # 不小于0
    y_max <- min(1, max_coverage + y_buffer_top) # 不大于1

    # 如果数据范围太小，设置最小显示范围
    if (y_max - y_min < 0.1) {
        center <- (y_min + y_max) / 2
        y_min <- max(0, center - 0.05)
        y_max <- min(1, center + 0.05)
    }

    p <- p + coord_cartesian(ylim = c(y_min, y_max)) # 使用coord_cartesian而不是ylim

    # 添加可接受范围的绿色区域和边界线
    if (!is.null(acceptable_range)) {
        # 添加可接受范围的绿色区域
        p <- p + annotate("rect",
            xmin = -Inf, xmax = Inf,
            ymin = acceptable_range[1], ymax = acceptable_range[2],
            fill = "green", alpha = 0.15
        ) # 低饱和度绿色

        # 添加可接受范围的边界线（深绿色虚线）
        p <- p + geom_hline(
            yintercept = acceptable_range[1],
            linetype = "dashed", color = "darkgreen", size = 1, alpha = 0.8
        ) # 下边界
        p <- p + geom_hline(
            yintercept = acceptable_range[2],
            linetype = "dashed", color = "darkgreen", size = 1, alpha = 0.8
        ) # 上边界
    }

    # 添加主题
    p <- p + theme_minimal() +
        theme(
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.text.y = element_text(size = 9),
            axis.title = element_text(size = 11, face = "bold"),
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
            strip.text = element_text(size = 10, face = "bold"),
            legend.position = if (show_legend) "bottom" else "none",
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 9),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "white", color = "grey90"),
            strip.background = element_rect(fill = "grey95", color = "grey90")
        )

    # 只有在显示图例时才添加图例设置
    if (show_legend) {
        p <- p + guides(fill = guide_legend(nrow = 2))
    } else {
        p <- p + guides(fill = "none")
    }

    # 显示图片
    print(p)

    # 保存图片
    if (save_plot) {
        ggsave(output_file, p, width = width, height = height, dpi = 300)
        cat("Coverage图已保存为:", output_file, "\n")
    }

    # 返回图片对象和数据
    return(list(
        plot = p,
        data = plot_data,
        method_labels = method_labels
    ))
}

plot_estimates_boxplot <- function(detailed_results,
                                   method_names = NULL,
                                   colors = NULL,
                                   highlight_methods = NULL, # 需要突出显示的方法
                                   show_legend = TRUE, # 是否显示图例
                                   true_value = 0, # 真实值（默认为0）
                                   add_true_line = TRUE, # 是否添加真实值参考线
                                   point_alpha = 0.6, # 散点透明度
                                   point_size = 1.5, # 散点大小
                                   jitter_width = 0.2, # 散点抖动宽度
                                   save_plot = FALSE,
                                   output_file = "estimates_boxplot.png",
                                   width = 12,
                                   height = 8) {
    # 加载必需的包
    required_packages <- c("ggplot2", "dplyr", "tidyr", "RColorBrewer")

    for (pkg in required_packages) {
        if (!requireNamespace(pkg, quietly = TRUE)) {
            stop(paste("请安装必需的包:", pkg, "\n可以使用: install.packages('", pkg, "')", sep = ""))
        }
    }

    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(RColorBrewer)

    # 检查输入数据
    if (is.null(detailed_results) || nrow(detailed_results) == 0) {
        stop("detailed_results 为空或无数据")
    }

    # 识别方法列（以_theta_point结尾的列）
    theta_cols <- grep("_theta_point$", names(detailed_results), value = TRUE)
    if (length(theta_cols) == 0) {
        stop("在 detailed_results 中未找到 _theta_point 列")
    }

    # 提取方法标识符（去掉_theta_point后缀）
    methods <- gsub("_theta_point$", "", theta_cols)
    n_methods <- length(methods)

    # 设置方法名称
    if (is.null(method_names)) {
        method_labels <- setNames(methods, methods) # 直接使用a, b, c...
    } else {
        method_labels <- method_names
        # 确保所有方法都有标签
        missing_methods <- setdiff(methods, names(method_labels))
        if (length(missing_methods) > 0) {
            warning("以下方法缺少标签，将使用默认标签: ", paste(missing_methods, collapse = ", "))
            for (m in missing_methods) {
                method_labels[m] <- m
            }
        }
    }

    # 设置颜色
    if (is.null(colors)) {
        if (n_methods <= 9) {
            colors <- RColorBrewer::brewer.pal(max(3, n_methods), "Set1")[1:n_methods]
        } else {
            colors <- rainbow(n_methods)
        }
    } else {
        if (length(colors) < n_methods) {
            warning("颜色数量不足，将重复使用颜色")
            colors <- rep(colors, length.out = n_methods)
        }
    }

    # 处理突出显示的方法
    if (!is.null(highlight_methods)) {
        # 检查highlight_methods是否有效
        invalid_methods <- setdiff(highlight_methods, methods)
        if (length(invalid_methods) > 0) {
            warning("以下突出显示的方法不存在: ", paste(invalid_methods, collapse = ", "))
            highlight_methods <- intersect(highlight_methods, methods)
        }

        if (length(highlight_methods) > 0) {
            # 创建新的颜色方案：突出显示的方法保持彩色，其他变灰
            highlight_colors <- colors[1:length(highlight_methods)]
            gray_color <- "grey70"

            final_colors <- rep(gray_color, n_methods)
            names(final_colors) <- methods

            # 为突出显示的方法分配彩色
            for (i in seq_along(highlight_methods)) {
                method <- highlight_methods[i]
                final_colors[method] <- highlight_colors[i]
            }

            colors <- final_colors
        }
    } else {
        names(colors) <- methods
    }

    # 转换数据为长格式
    # 先选择需要的列
    selected_cols <- c("simulation_id", theta_cols)
    plot_data_wide <- detailed_results[, selected_cols, drop = FALSE]

    # 转换为长格式
    plot_data <- plot_data_wide %>%
        gather(key = "Method_Raw", value = "Estimate", -simulation_id) %>%
        mutate(
            Method = gsub("_theta_point$", "", Method_Raw),
            Method_Label = method_labels[Method],
            Method_Factor = factor(Method_Label, levels = method_labels[methods]) # 保持顺序
        ) %>%
        filter(!is.na(Estimate)) # 移除NA值

    if (nrow(plot_data) == 0) {
        stop("没有有效的估计值数据")
    }

    # 创建简洁的箱线图+散点图
    p <- ggplot(plot_data, aes(x = Method_Factor, y = Estimate, fill = Method_Label, color = Method_Label)) +

        # 1. 箱线图
        geom_boxplot(alpha = 0.7, width = 0.6, outlier.shape = NA) + # 不显示异常值点，用散点代替

        # 2. 散点图 - 覆盖在箱线图上
        geom_jitter(
            alpha = point_alpha, size = point_size,
            width = jitter_width, height = 0,
            stroke = 0
        ) +

        # 使用方法名称作为颜色标识
        scale_fill_manual(values = setNames(colors, method_labels[methods]), name = "Method") +
        scale_color_manual(values = setNames(colors, method_labels[methods]), name = "Method") +
        labs(
            title = "Distribution of Point Estimates Across Simulations",
            subtitle = paste("Based on", length(unique(plot_data$simulation_id)), "simulations"),
            x = "Method",
            y = "Point Estimate"
        ) +
        theme_minimal() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12, face = "bold"),
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
            legend.position = if (show_legend) "bottom" else "none",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 10),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.background = element_rect(fill = "white", color = "grey90")
        )

    # 添加真实值参考线
    if (add_true_line) {
        p <- p + geom_hline(
            yintercept = true_value,
            linetype = "dashed", color = "darkred",
            size = 1.2, alpha = 0.8
        )
    }

    # 图例设置
    if (show_legend) {
        p <- p + guides(
            fill = guide_legend(nrow = 2, override.aes = list(alpha = 0.8)),
            color = "none" # 只显示fill图例
        )
    } else {
        p <- p + guides(fill = "none", color = "none")
    }

    # 显示图片
    print(p)

    # 保存图片
    if (save_plot) {
        ggsave(output_file, p, width = width, height = height, dpi = 300)
        cat("箱线图已保存为:", output_file, "\n")
    }

    # 计算汇总统计
    summary_stats <- plot_data %>%
        group_by(Method, Method_Label) %>%
        summarise(
            n_sims = n(),
            mean_estimate = mean(Estimate, na.rm = TRUE),
            median_estimate = median(Estimate, na.rm = TRUE),
            sd_estimate = sd(Estimate, na.rm = TRUE),
            q25 = quantile(Estimate, 0.25, na.rm = TRUE),
            q75 = quantile(Estimate, 0.75, na.rm = TRUE),
            min_estimate = min(Estimate, na.rm = TRUE),
            max_estimate = max(Estimate, na.rm = TRUE),
            bias = mean_estimate - true_value,
            .groups = "drop"
        )

    # 返回结果
    return(list(
        plot = p,
        plot_data = plot_data,
        summary_stats = summary_stats,
        method_labels = method_labels
    ))
}

# %% 画图区

method_names <- c(
    "a" = "cML_F",
    "b" = "cML_family",
    "c" = "IVW_F",
    "d" = "MR-Egger_F",
    "e" = "MR-lasso_F",
    "f" = "cML_G",
    "g" = "IVW_G",
    "h" = "MR-Egger_G",
    "i" = "MR-lasso_G"
)
highlight_methods <- c("f", "b", "i")
bias_plot <- plot_bias(b,
    label_size = 3, highlight_methods = highlight_methods,
    method_names = method_names
)
bias_plot
ggsave(a, "MR简单模拟结果/图包/")
power_plot <- plot_power(b,
    label_size = 3, method_names = method_names,
    power_type = "005", plot_type = "type1error", show_legend = FALSE
)

ggsave("MR简单模拟结果/图包/都是好的SNPs真值为0.02.svg",
    plot = power_plot$plot,
    width = 8, height = 6,
    units = "in",
    dpi = 300
)
coverage_plot <- plot_coverage(b,
    label_size = 3, highlight_methods = highlight_methods,
    coverage_type = "95", show_legend = FALSE, method_names = method_names,
    acceptable_range = c(0.9365, 0.9635)
)
coverage_plot
a$simulation_id <- 1:1000
d <- plot_estimates_boxplot(a,
    method_names = method_names,
    true_value = 0.1
)
a
ggsave("MR简单模拟结果/图包/都是好的SNPs真值为0.02（箱线图）.svg",
    plot = d$plot,
    width = 8, height = 6,
    units = "in",
    dpi = 300
)

