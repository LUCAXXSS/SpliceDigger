
#' Title
#'
#' @param SD_object
#' @param use_data
#' @param tSNE_perplexity
#' @param kmeans_centers
#' @param seed
#'
#' @return
#' @export
#'
#' @examples
SpliceDigger.calculate_tSNE <- function(SD_object=A,use_data="iso_exp_pct_filtered",
                                        tSNE_perplexity=8,kmeans_centers = 3,seed="42"){

  exp_to_use <- SD_object[[use_data]]



  exp_to_use <- scale(exp_to_use)%>%as.data.frame()

  exp_to_use <-exp_to_use%>%t()%>%as.data.frame()

  set.seed(seed)

  safevalue <- (length(row.names(exp_to_use))-1)/3
  tSNE_perplexity_safevalue <- round(safevalue,digits = 0)
  if (tSNE_perplexity>tSNE_perplexity_safevalue) {
    stop(paste0("tSNE_perplexity must be smaller than ", tSNE_perplexity_safevalue))

  }


  tsne_out <- Rtsne::Rtsne(exp_to_use,perplexity=tSNE_perplexity) # 运行tSNE

  tsne_out$Y <-as.data.frame(tsne_out$Y)
  SD_object$tsne_res[[use_data]] <- tsne_out


  kmeans_result <- stats::kmeans(tsne_out$Y, centers = kmeans_centers)
  clusters <- kmeans_result$cluster
  clusters1 <- paste0("cluster_",clusters)

  anoo_kmeans <- data.frame(clusters1,row.names =row.names(exp_to_use) )
  colnames(anoo_kmeans) <- paste0(use_data,"_kmeans_cluster")

  SD_object <- SpliceDigger.add_annotation(SD_object,anoo_kmeans)
  SD_object
}


#' Title SpliceDigger.plot_tSNE
#'
#' @param SD_object
#' @param plot_dataset
#' @param anno_to_show
#'
#' @return
#' @export
#'
#' @examples
#'


SpliceDigger.plot_tSNE <- function(SD_object=A,plot_dataset="iso_exp_pct_filtered",
                                        anno_to_show="sort"){

   tsne_out<-SD_object$tsne_res[[plot_dataset]]

   annos_2show <- SD_object$anno[,anno_to_show]

  p <- ggplot(data=tsne_out$Y)+
    geom_point(aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2],color=annos_2show),size=3)+theme_classic()

 p
}
