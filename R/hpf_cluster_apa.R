#' Cluster the two-dimensional coordinate files (RC and SA) to find the APA clusters
#'
#' @param RC_prebed RC file of one strand generated from WriteEndlist()
#' @param SA_prebed SA file of one strand generated from WriteEndlist()
#' @param out_dir The directory where the output files are stored
#' @param n_endloc A approximate number of samples used in the first clustering, changes when processing (default and recommended 100000)
#' @param blocksize A approximate gap length and size of clusters, changes with the situation during clustering (default and recommended 100)
#' @param lessmerge_n The minimum size of APA clusters (default 24)
#' @param SA_weight The weight increase of SA file (default 10)
#'
#' @return NULL. writes output to APA_clusters.out in the out_dir
#' @export
#'
#' @import data.table kpeaks Ckmeans.1d.dp
#' @rawNamespace import(GenomicAlignments, except = c(first, second, last))
#'
#' @examples
#' ClusterAPA("output/negative_strand/sa.prebed", "output/negative_strand/sa.prebed", "output/negative_strand/")
#' ClusterAPA("output/positive_strand/sa.prebed", "output/positive_strand/sa.prebed", "output/positive_strand/")
ClusterAPA <- function(RC_prebed, SA_prebed, out_dir, n_endloc = 100000, blocksize = 100, lessmerge_n = 24, SA_weight = 10){
  # start, need a cup of tea


  # read the files and merge
  RC <- fread(RC_prebed, header = F, drop = 4)
  SA <- fread(SA_prebed, header = F, drop = 4)
  SA <- SA[, V3 := V3 * SA_weight]
  CLH <- rbindlist(list(RC, SA))
  all_clu_bed <- list()
  rm(SA)
  rm(RC)
  gc()


  # loop to find APA clusters for every chr
  for (chrno in unique(CLH$V1)) {


    # prepare deduplicated data.table
    print(paste0("prebed of ",chrno," is loading"))
    CLH_chr <- CLH[V1 == chrno]
    setkey(CLH_chr, V2)
    setorder(CLH_chr, V2)
    CLH_chr <- CLH_chr[, sum(V3), by = V2]


    # find the nearest gap before yield location, and separate prebed
    start <- c()
    end <- c()
    n_tpcall <- ceiling(nrow(CLH_chr)/n_endloc)
    l_gap <- which(CLH_chr[2 : nrow(CLH_chr), V2] - CLH_chr[1 : (nrow(CLH_chr) - 1), V2]  >= blocksize)
    end_tpcall <- sapply((1 : n_tpcall), function(x) l_gap[max(which((x - l_gap/n_endloc) >= 0),1)])
    end_tpcall[n_tpcall] <- nrow(CLH_chr)  # make the last position of the chr included
    start_tpcall <- c(1, end_tpcall + 1)
    print(paste0("cluster the nodes of ",chrno," roughly"))


    # cluster to find possible rough clusters for every separated prebed
    for (i in 1 : n_tpcall) {

      # load node list
      clu_dt <- CLH_chr[start_tpcall[i] : end_tpcall[i]]
      setkey(CLH_chr, V2)
      clu_dt <- clu_dt[complete.cases(clu_dt[, 'V2'])] # avoid na

      # pre eval the data to get the K
      if (nrow(CLH_chr) <= 100) {
        tmpk <- ceiling(nrow(CLH_chr)/10)
      }else{
        histsmooth <- suppressWarnings(genpolygon(rep(clu_dt$V2, clu_dt$V1),   # suppressWarnings!!!!!!
                                                  binrule = "usr",
                                                  nbins = (clu_dt[.N, V2]-clu_dt[1, V2])/blocksize,
                                                  disp = FALSE))
        tmp <- findpolypeaks(histsmooth$mids, histsmooth$freqs, tcmethod="min")
        tmpk <- tmp$np
      }

      # cluster
      ckmean_1 <- suppressWarnings(Ckmeans.1d.dp(as.numeric(na.omit(clu_dt$V2)),    # suppressWarnings!!!!!!
                                                 k = tmpk, as.numeric(na.omit(clu_dt$V3))))
      clu_dt[, group := ckmean_1$cluster]

      # get the start & end of rough clusters
      start <- c(start,clu_dt[, min(V2), by = group]$V1)
      end <- c(end,clu_dt[, max(V2), by = group]$V1)
    }


    # get the data.table of rough clusters
    clu_bed <- data.table(Chr = chrno, clu_startposi = start, clu_endposi = end)


    # cluster to separate the rough clusters into APA clusters
    print(paste0("cluster the nodes of ",chrno," to get APA clusters"))
    start <- c()
    end <- c()
    for (i in 1 : nrow(clu_bed)) {

      # load the nodes of the rough clusters
      clu_dt <- CLH_chr[V2 %between% c(clu_bed[i, clu_startposi], clu_bed[i, clu_endposi])]

      # detect dump region, get the length of no-dump region to get the range of k
      ldump <- which((clu_dt[2 : .N]$V2 - clu_dt[1 : .N-1]$V2) > blocksize )
      ndump <- length(ldump)
      dedumplen <-  (clu_dt[.N, V2] - clu_dt[1, V2]) - sum(clu_dt[ldump+1, V2] - clu_dt[ldump, V2])
      kmax <- ceiling(dedumplen/blocksize)
      if ((clu_dt[.N, V2] - clu_dt[1, V2]) <= blocksize) {

        # too small to separate, rough clusters remain the same
        start <- c(start, clu_bed[i, clu_startposi])
        end <- c(end, clu_bed[i, clu_endposi])
      }else {

        # eval the range of k, to make long block separated and long dump ignored
        kdump <- ifelse(dedumplen > blocksize, (ndump+2), (ndump+1))
        k_range <- c(kdump : (ndump+kmax+1))

        # cluster
        ckmean_1 <- suppressWarnings(Ckmeans.1d.dp(as.numeric(na.omit(clu_dt$V2)),    # suppressWarnings!!!!!!
                                                   k = k_range,
                                                   as.numeric(na.omit(clu_dt$V3))))

        # get the start & end of APA clusters
        clu_dt[, group := ckmean_1$cluster]
        start <- c(start, clu_dt[, min(V2), by = group]$V1)
        end <- c(end, clu_dt[, max(V2), by = group]$V1)
      }
    }


    # get the data.table of APA clusters
    s_clu_bed <- data.table(Chr = chrno, clu_startposi = start, clu_endposi = end)


    # adjust the clusters size, >= 24 by defult
    lessmerge_l <- which(s_clu_bed[, 3] - s_clu_bed[, 2] < lessmerge_n)
    s_clu_bed[lessmerge_l, ]$clu_startposi <- floor((s_clu_bed[lessmerge_l, ]$clu_startposi +
                                                       s_clu_bed[lessmerge_l,]$clu_endposi)/2 - floor(lessmerge_n/2))
    s_clu_bed[lessmerge_l, ]$clu_endposi <- floor((s_clu_bed[lessmerge_l, ]$clu_startposi +
                                                     s_clu_bed[lessmerge_l,]$clu_endposi)/2 + floor(lessmerge_n/2))
    OLtmp <- which(s_clu_bed[2 : .N,clu_startposi] -   # check the overlap after increase in size
                     s_clu_bed[1 : .N-1,clu_endposi] < 0)
    s_clu_bed[OLtmp, clu_endposi :=   # adjust overlap
                (round((s_clu_bed[OLtmp, clu_endposi] + s_clu_bed[OLtmp+1, clu_startposi])/2))]
    s_clu_bed[OLtmp+1, clu_startposi := (s_clu_bed[OLtmp,clu_endposi])]


    # get the data.table of adjusted APA clusters
    all_clu_bed[[chrno]] <- s_clu_bed
  }


  # write down
  all_clu_bed_DT <- rbindlist(all_clu_bed)
  fwrite(all_clu_bed_DT,paste(out_dir, "APA_clusters.out", sep = ""),
         quote = FALSE, sep = "\t", col.names = FALSE)

}




