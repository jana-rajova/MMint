# function takes seurt obj and the deconvolution object. As Seurat can modify the names of the cells, this function makes sure that the deconvolution names become compatible
merge_sample_affix = function(seurat.obj, spacet.obj) {
  if (grepl('_', colnames(seurat.obj)[1])) {
    sample_suffix = t(as.data.frame(strsplit(colnames(seurat.obj), '_')))
    sample_suffix = cbind(sample_suffix, colnames(seurat.obj))
    rownames(sample_suffix) = sample_suffix[, 1]
    sample_suffix = sample_suffix[spacet.obj@input$counts@Dimnames[[2]],]
    spacet.obj@input$counts@Dimnames[[2]] = sample_suffix[, ncol(sample_suffix)]
  }
  return(spacet.obj)
}

# add deconvolution data from files in specified folder to the seurat object
add.spacet = function(seurat.obj,
                      st.folder = '/home/rstudio/data/stdata/',
                      folder = '/home/rstudio/mod/results/spacet/') {
  samples = c(
    'S1_C1',
    'S1_D1',
    'S2_A1',
    'S2_B1',
    'S2_C1',
    'S3_A1',
    'S3_B1',
    'S3_D1',
    'S4_A1' ,
    'S4_B1' ,
    'S4_C1',
    'S4_D1'
  )
  deconv.folder =  folder
  interface = c()
  init = TRUE
  
  for (s in unique(seurat.obj$sample_id)) {
    print(s) #processing sample
    affix = which(s == samples)
    
    SpaCET_obj = readRDS(paste0(deconv.folder, s, '_spacet.rds'))
    
    if (!(grepl('[ACTG]+-1', SpaCET_obj@input$counts@Dimnames[[2]][1]))) {
      SpaCET_obj = translate.coords.to.barcodes(SpaCET_obj, st.folder, s)
    }
    
    # # modification to the naming convention
    # # SpaCET_obj@input$counts@Dimnames[[2]] = paste0(s, '_', SpaCET_obj@input$counts@Dimnames[[2]])
    # # SpaCET_obj@input$counts@Dimnames[[2]] = gsub('$', paste0('_', affix),  SpaCET_obj@input$counts@Dimnames[[2]])
    # # idx_features = paste0(seurat.obj[, seurat.obj$sample_id == s]@assays$SCT@data@Dimnames[[2]], '_', affix)
    # idx_features = seurat.obj[, seurat.obj$sample_id == s]@assays$SCT@data@Dimnames[[2]]
    #
    # filtered.features.names = SpaCET_obj@input$counts[, idx_features]@Dimnames[[2]]
    #
    SpaCET_obj = merge_sample_affix(seurat.obj = seurat.obj[, seurat.obj$sample_id == s],
                                    spacet.obj = SpaCET_obj)
    
    # deconvolution addition
    all.celltypes = rownames(SpaCET_obj@results$deconvolution$propMat)
    main.celltypes = c('Malignant',
                       names(SpaCET_obj@results$deconvolution$Ref$lineageTree))
    sub.celltypes = all.celltypes[!all.celltypes %in% main.celltypes] # all sublineages besides the Malignant substates
    
    colnames(SpaCET_obj@results$deconvolution$propMat) = SpaCET_obj@input$counts@Dimnames[[2]]
    if (init == TRUE) {
      init = FALSE
      propmat = SpaCET_obj@results$deconvolution$propMat[all.celltypes,]
    } else{
      propmat = cbind(propmat, SpaCET_obj@results$deconvolution$propMat[all.celltypes,])
    }
    
    # add interface
    if (!('interface' %in% names(SpaCET_obj@results$CCI))) {
      SpaCET_obj <- SpaCET.identify.interface(SpaCET_obj)
    }
    
    names(SpaCET_obj@results$CCI$interface) = SpaCET_obj@input$counts@Dimnames[[2]]
    interface = c(interface,
                  SpaCET_obj@results$CCI$interface)
  }
  
  seurat.obj@tools$SpaCET.main = main.celltypes
  seurat.obj@tools$SpaCET.sub = sub.celltypes
  
  seurat.obj$interface = interface
  seurat.obj@reductions$SpaCET = propmat
  return(seurat.obj)
}

translate.coords.to.barcodes = function(SpaCET_obj, st.folder, sample) {
  spacet.translation = read.table(paste0(st.folder, sample, '/spatial/tissue_positions_list.csv'),
                                  sep = ',')
  spacet.translation$coord = paste0(spacet.translation$V3, 'x', spacet.translation$V4)
  rownames(spacet.translation) = spacet.translation$coord
  
  # spacet.translation = subset(spacet.translation, subset = spacet.translation$V2 == 1)
  spacet.translation = spacet.translation[SpaCET_obj@input$counts@Dimnames[[2]],]
  
  SpaCET_obj@results$dimnames.original = SpaCET_obj@input$counts@Dimnames[[2]]
  SpaCET_obj@input$counts@Dimnames[[2]] = spacet.translation$V1
  return(SpaCET_obj)
}

plot_stats = function(se.object,
                      sample,
                      output.folder) {
  plot1 = VlnPlot(se.object,
                  features = c("nCount_Spatial"),
                  pt.size = 0.1) + NoLegend()
  
  plot2 = SpatialFeaturePlot(se.object,
                             features = c("nCount_Spatial"),
                             combine = TRUE)  + theme(legend.position = "right")
  
  plot3 = VlnPlot(se.object,
                  features = c('percent.mito'),
                  pt.size = 0.1) + NoLegend()
  
  plot4 = SpatialFeaturePlot(se.object,
                             features = c('percent.mito'),
                             combine = TRUE) + theme(legend.position = "right")
  
  plot5 = VlnPlot(se.object,
                  features = c('percent.ribo'),
                  pt.size = 0.1) + NoLegend()
  
  plot6 = SpatialFeaturePlot(se.object,
                             features = c('percent.ribo'),
                             combine = TRUE) + theme(legend.position = "right")
  
  # wrap_plots(plot1, plot2, plot3, plot4, plot5, plot6, ncol = 2, byrow = TRUE, widths = c(1,5))
  
  pdf(
    file = paste0(output.folder, sample, '_nCounts_mito_ribo.pdf'),
    height = 15,
    width = 7.5
  )
  print(wrap_plots(
    plot1,
    plot2,
    plot3,
    plot4,
    plot5,
    plot6,
    ncol = 2,
    byrow = TRUE,
    widths = c(1, 2)
  ))
  dev.off()
}

# load seurat objet with predefined samples
load.samples = function(visium.folder = '/home/rstudio/data/stdata/',
                        output.folder.data =  '/home/rstudio/mod/data/transcriptome_analysis/',
                        output.folder.plots =  '/home/rstudio/mod/results/transcriptome_analysis/plot/',
                        samples = c(
                          'S1_C1',
                          'S1_D1',
                          'S2_A1',
                          'S2_B1',
                          'S2_C1',
                          'S3_A1' ,
                          'S3_B1' ,
                          'S3_D1' ,
                          'S4_A1' ,
                          'S4_B1' ,
                          'S4_C1',
                          'S4_D1'
                        ),
                        save = TRUE) {
  dir.create(output.folder.data, showWarnings = FALSE)
  dir.create(output.folder.plots, showWarnings = FALSE)
  
  variable_features = c()
  
  slide.section = data.frame(
    'slides' = c(117, 117, 121, 121, 121, 117, 117, 121, 110, 110, 110, 110),
    row.names = c(
      'S1_C1',
      'S1_D1',
      'S2_A1',
      'S2_B1',
      'S2_C1',
      'S3_A1',
      'S3_B1',
      'S3_D1',
      'S4_A1',
      'S4_B1',
      'S4_C1',
      'S4_D1'
    )
  )
  gene_threshold = 5
  
  if (!1 %in% grepl('pre_processed_slides.rds',
                    list.files(output.folder.data))) {
    print('Processing slides anew')
    st.list = list()
    
    for (i in 1:length(samples)) {
      # for (i in 1:5){
      
      print(samples[i])
      st.list[[i]] = Load10X_Spatial(data.dir = paste0(visium.folder, samples[i]), 
                                     slice = samples[i])
      
      # calculate mt and rb
      st.list[[i]]$percent.mito = colSums(st.list[[i]][rownames(st.list[[i]])[grepl('^mt-', rownames(st.list[[i]]), ignore.case = TRUE)],][['Spatial']]$counts) /
        colSums(st.list[[i]][['Spatial']]$counts)
      st.list[[i]]$percent.ribo = colSums(st.list[[i]][rownames(st.list[[i]])[grepl('^rps|^rpl', rownames(st.list[[i]]), ignore.case = TRUE)],][['Spatial']]$counts) /
        colSums(st.list[[i]][['Spatial']]$counts)
      
      # calculate per gene
      genes_QC = t(as.matrix(st.list[[i]]@assays$Spatial@layers$counts))
      colnames(genes_QC) = rownames(st.list[[i]])
      genes_QC = colSums(genes_QC>0) > gene_threshold
      genes_QC = names(genes_QC[genes_QC==TRUE])
      
      # calculate per feature
      st.list[[i]] = st.list[[i]][, st.list[[i]]$nCount_Spatial > 100]
      st.list[[i]] = st.list[[i]][, st.list[[i]]$nFeature_Spatial > 50]
      st.list[[i]] = st.list[[i]][genes_QC,]
      
      
      # calculate mt and rb
      
      # add slide and sample_id
      st.list[[i]]$slide = slide.section[samples[[i]],]
      st.list[[i]]$sample_id = samples[i]
      if (remove.mt.rb == TRUE) {
        st.list[[i]] = st.list[[i]][rownames(st.list[[i]])[!grepl('^mt|^rps|^rpl', rownames(st.list[[i]]), ignore.case = TRUE)],]
      }
      
      plot_stats(st.list[[i]],
                 sample = samples[i],
                 output.folder = output.folder.plots)
      
      print(SpatialFeaturePlot(
        st.list[[i]],
        features = c("MET", "MALAT1", "MLANA", "PMEL"),
        slot = 'counts'
      ))
      
      st.list[[i]] = SCTransform(
        st.list[[i]],
        assay = "Spatial",
        verbose = FALSE,
        vars.to.regress = c('nCount_Spatial'),
        return.only.var.genes = FALSE,
        n_genes = 5000,
        method = "glmGamPoi"
      )
      
      
      st.list[[i]] = NormalizeData(st.list[[i]])
      st.list[[i]] = DimProcess(st.list[[i]])
      st.list[[i]]$orig.ident = samples[i]
      st.list[[i]]$source = gsub('_[A-D][1-2]', '' , st.list[[i]]$orig.ident)
      print(SpatialFeaturePlot(
        st.list[[i]],
        features = c("MET", "MALAT1", "MLANA", "PMEL"),
        slot = 'counts'
      ))
      variable_features = unique(c(st.list[[i]]@assays$SCT@var.features, variable_features))
    }
    
    # saveRDS(st.list, file = paste0(output.folder, 'data/', paste0(sources, collapse ='-'), '/samples_individuals.Rds'))
    
    if (save == TRUE) {
      saveRDS(st.list,
              file = paste0(
                output.folder.data,
                paste(sources, collapse = "-"),
                '_pre_processed_slides.rds'
              ))
      
    }
    
    for (st in st.list) {
      variable_features = unique(c(st@assays$SCT@var.features,
                                   variable_features))
    }
    
    if (length(samples) > 1) {
      se = merge(st.list[[1]], st.list[2:length(samples)])
      se = JoinLayers(se, assay = 'Spatial')
      VariableFeatures(se) = variable_features
      
    } else{
      se = st.list[[1]]
    }
    
    DefaultAssay(se) = "SCT"
    
    se = DimProcess(se)
    
    se$source = gsub('_[A-D][1-2]', '' , se$orig.ident)
    VariableFeatures(se) = variable_features
    
    saveRDS(se,
            file = paste0(
              output.folder.data,
              paste(sources, collapse = "-"),
              '_pre_processed_slides_merged.rds'
            ))
  } else {
    st.list = readRDS(paste0(
      output.folder.data,
      paste(sources, collapse = "-"),
      '_pre_processed_slides.rds'
    ))
    se = readRDS(paste0(
      output.folder.data,
      paste(sources, collapse = "-"),
      '_pre_processed_slides_merged.rds'
    ))
  }
  
  SpatialDimPlot(se,
                 label = TRUE,
                 label.size = 9,
                 ncol = 3)
  gc()
  return(list(se, st.list))
}

DimProcess = function(seurat.obj,
                      var.exp = 0.05,
                      pca_npcs = 100) {
  seurat.obj = RunPCA(seurat.obj,
                      assay = 'SCT',
                      verbose = FALSE,
                      npcs = pca_npcs)
  pca.max = max(which((
    seurat.obj@reductions$pca@stdev[1:99] - seurat.obj@reductions$pca@stdev[2:100]
  ) > var.exp
  ))
  print(pca.max)
  
  seurat.obj = FindNeighbors(seurat.obj,
                             dims = 1:pca.max,
                             prune.SNN = 1 / 15)
  seurat.obj = FindClusters(seurat.obj,
                            resolution = res)
  
  print(seurat.obj)
  seurat.obj = RunUMAP(seurat.obj,
                       dims = 1:pca.max)
  
  return(seurat.obj)
}


plot_qc = function(seurat.obj,
                   output.folder) {
  qc.list = list()
  # for (source in c(sources, list(sources))){ # use if you want an error unless you're running all samples/sources at least
  present_sources = unique(seurat.obj$source)
  
  if (length(present_sources) > 1) {
    present_sources = c(present_sources, list(present_sources))
  }
  for (source in present_sources) {
    plot.name = paste0(source, collapse = '-')
    print(paste('Processing:', plot.name))
    if (plot.name == sources[1]) {
      yax = c(
        "Unique genes per spot",
        "Total counts per spots",
        "Total counts per gene (log10 scale)",
        "Total spots per gene"
      )
    } else{
      yax = c("", "", "", "")
    }
    
    p1 <- ggplot() +
      geom_histogram(
        data = seurat.obj[, seurat.obj$source %in% source][[]],
        aes(nFeature_Spatial),
        fill = "aquamarine",
        alpha = 0.7,
        bins = 50
      ) +
      xlim(min(seurat.obj$nFeature_Spatial),
           max(seurat.obj$nFeature_Spatial)) +
      ylab(yax[1]) +
      ggtitle(plot.name)
    
    p2 <- ggplot() +
      geom_histogram(
        data = seurat.obj[, seurat.obj$source %in% source][[]],
        aes(nCount_Spatial),
        fill = "skyblue",
        alpha = 0.7,
        bins = 50
      ) +
      xlim(min(seurat.obj$nCount_Spatial),
           max(seurat.obj$nCount_Spatial)) +
      ylab(yax[2])
    
    gene_attr_all <-
      data.frame(
        nUMI = Matrix::rowSums(seurat.obj@assays$Spatial@layers$counts),
        nSpots = Matrix::rowSums(seurat.obj@assays$Spatial@layers$counts > 0)
      )
    
    gene_attr <-
      data.frame(
        nUMI = Matrix::rowSums(seurat.obj[, seurat.obj$source %in% source]@assays$Spatial@layers$counts),
        nSpots = Matrix::rowSums(seurat.obj[, seurat.obj$source %in% source]@assays$Spatial@layers$counts > 0)
      )
    p3 <- ggplot() +
      geom_histogram(
        data = gene_attr,
        aes(nUMI),
        fill = "cornflowerblue",
        alpha = 0.7,
        bins = 50
      ) +
      xlim(log10(min(gene_attr_all$nUMI)), log10(max(gene_attr_all$nUMI))) +
      scale_x_log10() +
      ylab(yax[3])
    
    p4 <- ggplot() +
      geom_histogram(
        data = gene_attr,
        aes(nSpots),
        fill = "navy",
        alpha = 0.7,
        bins = 50
      ) +
      xlim(min(gene_attr_all$nSpots), max(gene_attr_all$nSpots)) +
      ylab(yax[4])
    
    qc.list[[plot.name]] = grid.arrange(p1, p2, p3, p4, ncol = 1)
  }
  
  pdf(
    file = paste0(output.folder, paste0(unique(
      seurat.obj$source
    ), collapse = '-'), '_QC.pdf'),
    height = 10,
    width = 15
  )
  marrangeGrob(qc.list,
               ncol = 5,
               nrow = 1,
               top = 'QC per source')
  dev.off()
}

