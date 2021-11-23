library(SCOPE)
library(BSgenome.Mmusculus.UCSC.mm10)
bamFile <- list.files('/nfs02data1/home/xuwanxing/test/cnv/SCOPE/dedup', pattern = '*.dedup.bam$')
bamdir <- file.path('/nfs02data1/home/xuwanxing/test/cnv/SCOPE/dedup', bamFile)
sampname_raw <- sapply(strsplit(bamFile, ".", fixed = TRUE), "[", 1)
bambedObj <- get_bam_bed(bamdir = bamdir, sampname = sampname_raw, 
                        hgref = "mm10")
bamdir <- bambedObj$bamdir
sampname_raw <- bambedObj$sampname
ref_raw <- bambedObj$ref

mapp <- get_mapp(ref_raw, hgref = "mm10")
head(mapp)
gc <- get_gc(ref_raw, hgref = "mm10")
values(ref_raw) <- cbind(values(ref_raw), DataFrame(gc, mapp))
ref_raw

coverageObj <- get_coverage_scDNA(bambedObj, mapqthres = 40,     》》》用时长
                                seq = 'paired-end', hgref = "mm10")
Y_raw <- coverageObj$Y

QCmetric_raw <- get_samp_QC(bambedObj)   》》》 用时长
qcObj <- perform_qc(Y_raw = Y_raw, 
    sampname_raw = sampname_raw, ref_raw = ref_raw, 
    QCmetric_raw = QCmetric_raw)
Y <- qcObj$Y
sampname <- qcObj$sampname
ref <- qcObj$ref
QCmetric <- qcObj$QCmetric

Gini <- get_gini(Y_sim)
normObj.sim <- normalize_codex2_ns_noK(Y_qc = Y_sim,
                                        gc_qc = ref_sim$gc,
                                        norm_index = which(Gini<=0.12))
Yhat.noK.sim <- normObj.sim$Yhat
beta.hat.noK.sim <- normObj.sim$beta.hat
fGC.hat.noK.sim <- normObj.sim$fGC.hat
N.sim <- normObj.sim$N

ploidy.sim <- initialize_ploidy(Y = Y_sim, Yhat = Yhat.noK.sim, ref = ref_sim)
normObj.scope.sim <- normalize_scope_foreach(Y_qc = Y_sim, gc_qc = ref_sim$gc,
    K = 1, ploidyInt = ploidy.sim,
    norm_index = which(Gini<=0.12), T = 1:5,
    beta0 = beta.hat.noK.sim, nCores = 2)

Yhat.sim <- normObj.scope.sim$Yhat[[which.max(normObj.scope.sim$BIC)]]
fGC.hat.sim <- normObj.scope.sim$fGC.hat[[which.max(normObj.scope.sim$BIC)]]

clones <- c("normal", "tumor1", "normal", "tumor1", "tumor1")
ploidy.sim.group <- initialize_ploidy_group(Y = Y_sim, Yhat = Yhat.noK.sim,
                                ref = ref_sim, groups = clones)
ploidy.sim.group

normObj.scope.sim.group <- normalize_scope_group(Y_qc = Y_sim,
                                    gc_qc = ref_sim$gc,
                                    K = 1, ploidyInt = ploidy.sim.group,
                                    norm_index = which(clones=="normal"),
                                    groups = clones,
                                    T = 1:5,
                                    beta0 = beta.hat.noK.sim)
Yhat.sim.group <- normObj.scope.sim.group$Yhat[[which.max(
                                    normObj.scope.sim.group$BIC)]]
fGC.hat.sim.group <- normObj.scope.sim.group$fGC.hat[[which.max(
                                    normObj.scope.sim.group$BIC)]]

plot_EM_fit(Y_qc = Y_sim, gc_qc = ref_sim$gc, norm_index = which(Gini<=0.12), 
            T = 1:5,
            ploidyInt = ploidy.sim, beta0 = beta.hat.noK.sim,
            filename = "plot_EM_fit_demo.pdf")

chrs <- unique(as.character(seqnames(ref_sim)))
segment_cs <- vector('list',length = length(chrs))
names(segment_cs) <- chrs
for (chri in chrs) {
    message('\n', chri, '\n')
    segment_cs[[chri]] <- segment_CBScs(Y = Y_sim,
                                    Yhat = Yhat.sim,
                                    sampname = colnames(Y_sim),
                                    ref = ref_sim,
                                    chr = chri,
                                    mode = "integer", max.ns = 1)
}
iCN_sim <- do.call(rbind, lapply(segment_cs, function(z){z[["iCN"]]}))

plot_iCN(iCNmat = iCN_sim, ref = ref_sim, Gini = Gini, 
        filename = "plot_iCN_demo")

