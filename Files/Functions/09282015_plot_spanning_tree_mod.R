plot_spanning_tree_mod<-function (cds, x = 1, y = 2, color_by = "State", show_tree = TRUE, 
    show_backbone = TRUE, backbone_color = "black", markers = NULL, 
    show_cell_names = FALSE, cell_name_size = 1,color_custom) 
{
    require(grid)
    lib_info_with_pseudo <- pData(cds)
    S_matrix <- reducedDimS(cds)
    if (is.null(S_matrix)) {
        stop("You must first call reduceDimension() before using this function")
    }
    ica_space_df <- data.frame(t(S_matrix[c(x, y), ]))
    colnames(ica_space_df) <- c("ICA_dim_1", "ICA_dim_2")
    ica_space_df$sample_name <- row.names(ica_space_df)
    ica_space_with_state_df <- merge(ica_space_df, lib_info_with_pseudo, 
        by.x = "sample_name", by.y = "row.names")
    dp_mst <- minSpanningTree(cds)
    if (is.null(dp_mst)) {
        stop("You must first call orderCells() before using this function")
    }
    edge_list <- as.data.frame(get.edgelist(dp_mst))
    colnames(edge_list) <- c("source", "target")
    edge_df <- merge(ica_space_with_state_df, edge_list, by.x = "sample_name", 
        by.y = "source", all = TRUE)
    edge_df <- rename(edge_df, c(ICA_dim_1 = "source_ICA_dim_1", 
        ICA_dim_2 = "source_ICA_dim_2"))
    edge_df <- merge(edge_df, ica_space_with_state_df[, c("sample_name", 
        "ICA_dim_1", "ICA_dim_2")], by.x = "target", by.y = "sample_name", 
        all = TRUE)
    edge_df <- rename(edge_df, c(ICA_dim_1 = "target_ICA_dim_1", 
        ICA_dim_2 = "target_ICA_dim_2"))
    diam <- as.data.frame(as.vector(V(dp_mst)[get.diameter(dp_mst, 
        weights = NA)]$name))
    colnames(diam) <- c("sample_name")
    diam <- arrange(merge(ica_space_with_state_df, diam, by.x = "sample_name", 
        by.y = "sample_name"), Pseudotime)
    markers_exprs <- NULL
    if (is.null(markers) == FALSE) {
        markers_fData <- subset(fData(cds), gene_short_name %in% 
            markers)
        if (nrow(markers_fData) >= 1) {
            markers_exprs <- melt(exprs(cds[row.names(markers_fData), 
                ]))
            markers_exprs <- merge(markers_exprs, markers_fData, 
                by.x = "Var1", by.y = "row.names")
            markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
            markers_exprs$feature_label[is.na(markers_exprs$feature_label)] <- markers_exprs$Var1
        }
    }
    if (is.null(markers_exprs) == FALSE && nrow(markers_exprs) > 
        0) {
        edge_df <- merge(edge_df, markers_exprs, by.x = "sample_name", 
            by.y = "Var2")
        g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
            y = source_ICA_dim_2, size = log10(value + 0.1))) + 
            facet_wrap(~feature_label)
    }
    else {
        g <- ggplot(data = edge_df, aes(x = source_ICA_dim_1, 
            y = source_ICA_dim_2))
    }
    if (show_tree) {
        g <- g + geom_segment(aes_string(xend = "target_ICA_dim_1", 
            yend = "target_ICA_dim_2", color = color_by), size = 0.3, 
            linetype = "solid", na.rm = TRUE)
    }
    g <- g + geom_point(aes_string(color = color_by), na.rm = TRUE,size=4)
    if (show_backbone) {
#         g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2), 
#             color = I(backbone_color), size = 0.75, data = diam, 
#             na.rm = TRUE) + geom_point(aes_string(x = "ICA_dim_1", 
#             y = "ICA_dim_2", color = color_by),size=4, 
#             data = diam, na.rm = TRUE,alpha=0.5) + scale_color_manual(values=color_custom)
              g <- g + geom_path(aes(x = ICA_dim_1, y = ICA_dim_2), 
                  color = I(backbone_color), size = 1.25, alpha=0.8, data = diam, 
                  na.rm = TRUE) + scale_color_manual(values=color_custom)
    }
    if (show_cell_names) {
        g <- g + geom_text(aes(label = sample_name), size = cell_name_size)
    }
    g <- g + theme(panel.border = element_blank(), axis.line = element_line()) + 
        theme(panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank()) + 
        theme(panel.grid.major.x = element_blank(), panel.grid.major.y = element_blank()) + 
        ylab("Component 1") + xlab("Component 2") + theme(legend.position = "top", 
        legend.key.height = unit(0.35, "in")) + theme(legend.key = element_blank()) + 
        theme(panel.background = element_rect(fill = "white"))
    g
}