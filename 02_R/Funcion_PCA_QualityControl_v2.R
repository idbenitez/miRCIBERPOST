
# Principal Component Analysis (PCA) para control de calidad

require(vegan)  # vegan::betadisper -> Multivariate homogeneity of groups dispersions
#BiocManager::install("limma")
library(limma)
library(pROC)
library(dplyr)


pca_QC <- function(datos, patient_id_name, group_name, noms_proteines){
  
  # Parametros de entrada:
  # datos: base de datos
  # patient_id_name: string con el nombre de la columna que contiene el id del paciente
  # group_name: string con el nombre de la columna que contiene el grupo
  # noms_proteines: vector con el nombre de las proteinas a estudiar
  #
  
  # Define inital parameters
  do = "analyze"
  method = "euclidean"
  type = "median"
  coef = 1.5
  labels = FALSE
  
  
  # Data
  to_outliers <- datos[,noms_proteines]

  # Body of the funcion
  dd <- dist(to_outliers, method = method)
  distances <- vegan::betadisper(dd, datos[,group_name], type = type, bias.adjust = FALSE,
                                 sqrt.dist = FALSE, add = FALSE)
  detect_outliers <- data.frame(distances = distances$distances,
                                group = distances$group, sample = datos[,patient_id_name])
  limit <- data.frame(aggregate(detect_outliers$distances,
                                list(detect_outliers$group), function(x) {
                                  quantile(x, 0.75) + coef * IQR(x)
                                }))
  colnames(limit)[1] <- "group"
  detect_outliers <- merge(detect_outliers, limit, by = "group")
  detect_outliers <- detect_outliers %>% mutate(out = as.factor(ifelse(distances > x, 1, 0)))
  final_outliers <- detect_outliers %>% filter(out == 1) %>%
    dplyr::select(sample, group, distances, x) %>% 
    rename(distance_to_centroid = distances, limit_distance = x)
  
  
  if (do == "analyze") {
    vectors <- data.frame(distances$vectors)
    centroids <- data.frame(distances$centroids)
    total_outliers <- vectors %>% tibble::rownames_to_column("sample") %>%
      mutate(Group = datos[,group_name])
    find_hull <- function(x) {
      x[chull(x$PCoA1, x$PCoA2), ]
    }
    hulls <- total_outliers %>% dplyr::select(-sample) %>% group_by(Group) %>%
      dplyr::do(find_hull(.))
    polygon_plot <- ggplot(total_outliers, aes(x = PCoA1,
                                               y = PCoA2)) + geom_polygon(data = hulls, alpha = 0.5,
                                                                          aes(fill = Group)) + {
                                                                            if (!labels)
                                                                              geom_point(aes(shape = Group), size = 3, alpha = 0.7)
                                                                          } + geom_label(data = centroids, aes(x = PCoA1, y = PCoA2,
                                                                                                               color = rownames(centroids), label = rownames(centroids)),
                                                                                         show.legend = FALSE) + {
                                                                                           if (labels)
                                                                                             geom_text(aes(label = sample))
                                                                                         } + theme_bw()
    distance_boxplot <- ggplot(detect_outliers, aes(group,
                                                    distances, fill = group)) + geom_boxplot(coef = coef,
                                                                                              alpha = 0.8) + ylab("Distance to group centroid") +
      xlab("") + {
        if (labels)
          ggrepel::geom_label_repel(data = detect_outliers[detect_outliers$out ==
                                                             1, ], aes(label = sample), na.rm = TRUE, size = 4,
                                    show.legend = FALSE)
      } + theme_bw() + theme(axis.text.x = element_text(angle = 45,
                                                        hjust = 1))
    # return(list(polygon_plot = polygon_plot, distance_boxplot = distance_boxplot,
    #     outliers = final_outliers))
  }
  
  # Results to return
  list(polygon_plot, 
       final_outliers,
       distance_boxplot)

}


## EXPRESION DIFERENCIAL CRUDA/AJUSTADA
## VARIABLE GRUPO 2 categorias

expresion_diferencial <- function(datos, var_grupo, vars_ajuste, noms_proteines,names_id,  base = NULL) {
  
  ## Table: Models have been adjusted for age, sex and BMI.
  datos <- datos[complete.cases(datos[,c( var_grupo, vars_ajuste,names_id)]),]
  
  # "unadjusted"
  if(is.null(vars_ajuste)){
    design <- model.matrix(as.formula(paste0("~  0 + ", var_grupo)), data = datos)
    colnames(design) <- c("gr1","gr2")
    rownames(design) <- datos[,names_id]
    cont.matrix <- makeContrasts(DIF = gr2-gr1, levels=design) 
  }
  # if "adjusted"
  else{
    design <- model.matrix(as.formula(paste0("~  ", var_grupo, "+ ", paste0(vars_ajuste, collapse="+"))), data = datos)
    colnames(design)[1:2] <- c("INTERCEPT", "gr2")
    rownames(design) <- datos[,names_id]
    cont.matrix <- makeContrasts(gr2, levels=design) 
    }
  
  aux <- t(datos[, noms_proteines])
  rownames(aux) <- noms_proteines
  lmF <- lmFit(aux, design)
  
  ## definir los contrastes ##
  cont.fit <- contrasts.fit(lmF, cont.matrix)
  eBa <- eBayes(cont.fit)
  
  # Obtencion de la tabla de resultados con los 10^-ddCt y los p-values ##
  limma <- topTable(eBa,coef = 1, adjust.method = "fdr", number=length(noms_proteines), sort.by="p")
  # logFC: estimate of the log10-fold-change corresponding to the effect or contrast
  #head(limma) # logFC coincide con el slope que reporta Olink en el informe
  #head(limma[order(limma$P.Value, decreasing=FALSE),])
  #limma[row.names(limma)=="Infl_IL_20",]

  ddCt <- data.frame(rownames(limma),base^(limma$logFC),limma$P.Value,limma$adj.P.Val)
  colnames(ddCt) <- c("Names","FC","p.value","FDR")
  
  pval_FDR0.20 <- ifelse(length(ddCt$p.value[ddCt$FDR<0.20])==0, NA, max(ddCt$p.value[ddCt$FDR<0.20]))
   
  aux <- data.frame(datos %>% dplyr::select(noms_proteines))
  res <- apply(aux,2,function(x){auc(datos[,var_grupo], x)})
  res <- ifelse(res<0.5,1-res,res)
  res <- data.frame(Names =names(res),AUC = res)
  ddCt <- merge(ddCt,res,by = "Names")
  #ddCt$Gene_ID <- gsub("\\.","-",ddCt$Gene_ID)
  #Reordenaci?n de la tabla de resultados ddCt
  ddCt <- ddCt[order(ddCt$p.value, decreasing=FALSE),]
  
  #DT::datatable(ddCt, extensions = "Buttons", options = list(pageLength = 20, dom = "Bftip", buttons = c("copy", "csv", "excel", "pdf")), caption = "Differentially expressed") %>% DT::formatRound(names(ddCt)[-1], 3)
  #DE <- ddCt[which((ddCt$FC>1.25 | ddCt$FC<0.8) & ddCt$p.value <0.05),"Names"]
  
  
  ### Volcano plot (adjusted Filtrado dicotomizado)
  
  # Puntos de corte:
  #  -ln(pvalue): (-1)*log(0.05)
  #             : (-1)*log(max(ddCt$p.value[ddCt$FDR<0.20]))
  #  log10(Fold Change): (-1)*log10(1.15)
  #                   : log10(1.15)
  # x = log10(FC) = log10(ddCt$FC)
  # y = (-1)log(pval) = (-1)log(ddCt$p.value)
  
  
  
  cols <- c("Downregulated" = "red", "Nonsignificant" = "black", "Upregulated" = "blue")
  ss <- data.frame(x=eBa$coefficients[,1], y=c(-1)*log(eBa$p.value[,1]),names = rownames(eBa$lods))
  ss$color <- ifelse(ss$y > c(-1)*log(0.05) & ss$x < c(- 0.06069784),"Downregulated",
                     ifelse(ss$y > c(-1)*log(0.05) & ss$x >  0.06069784 ,"Upnregulated","Nonsignificant"))
  
  ss$names[ss$y < (-1)*log(0.05) & ss$x >= -(1)*log10(1.15)] <- NA
  ss$names[ss$y < (-1)*log(0.05) & ss$x < log10(1.15)] <- NA
  
  #cairo_pdf("Olink post-covid/09_results/Volcano.pdf",width = 12,height = 10)
  bp <- ggplot(ss, aes(x = x, y = y, labels=names,color = color)) + geom_point(size = 2.5, alpha = 1, na.rm = T)  +
    scale_colour_manual(values = cols) +
    ggtitle(label = "", subtitle = "") +
    geom_point(size = 2.5, alpha = 1, na.rm = T) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top",
          legend.title = element_blank(),
          axis.text=element_text(size=23),
          axis.title=element_text(size=23,face="bold"),
          legend.text = element_text(size=23)) +
    xlab(expression(log10("Fold Change"))) +
    ylab(expression(-ln("p value"))) +
    geom_hline(yintercept = (-1)*log(0.05), colour="blue", linetype="dashed") +
    geom_vline(xintercept = (-1)*log10(1.15), colour="#990000", linetype="dashed") +
    geom_vline(xintercept = log10(1.15), colour="#990000", linetype="dashed") +
    scale_y_continuous(trans = "log1p") + ggrepel::geom_text_repel(aes(x = x, y = y, label=names))
  
  
  
  #+ xlim(-0.1,0.1)
    # + ggrepel::geom_text_repel(aes(x = x, y = y, label=names)) 
  #if(!is.na(pval_FDR0.20)) bp  <- bp + geom_hline(yintercept = (-1)*log(pval_FDR0.20), colour="darkgreen", linetype="dashed")

  #bp
  #dev.off()
  list(limma %>% mutate_if(is.numeric, round, digits=5),
       ddCt  %>% mutate_if(is.numeric, round, digits=3), 
       bp,
       eBa)
}









































cors <- function(df) { 
  # turn all three matrices (r, n, and P into a data frame)
  M <- Hmisc::rcorr(as.matrix(df))
  # return the three data frames in a list return(Mdf)
  Mdf <- purrr::map(M, ~data.frame(.x))
}
#cors(sel_patients_olink[,c(vars.sleep, vars.prot)]) %>% first() %>% head() %>% kable()

formatted_cors <- function(df){
  res <- cors(df) %>%
    purrr::map(~tibble::rownames_to_column(.x, var="measure1")) %>%
    purrr::map(~tidyr::pivot_longer(.x, -measure1, "measure2")) %>% 
    bind_rows(.id = "id") %>%
    tidyr::pivot_wider(names_from = id, values_from = value) %>%
    mutate(sig_p = ifelse(P < .05, T, F), p_if_sig = ifelse(P <.05, P, NA), 
           r_if_sig = ifelse(P <.05, r, NA), pval=ifelse(P<0.001,"***", ifelse(P<0.01,"**",ifelse(P<0.05,"*",NA))))
  res %>% filter(measure1 %in% vars.sleep) %>% filter(!measure2 %in% vars.sleep)
}
#formatted_cors(sel_patients_olink[,c(vars.sleep, vars.prot)]) %>% head() %>% kable()
#formatted_cors(sel_patients_olink[,c(vars.sleep, vars.prot)]) %>% filter(P < 0.001)


# Correlaciones CRUDAS. Variables sleep y proteinas (p.value < 0.05)