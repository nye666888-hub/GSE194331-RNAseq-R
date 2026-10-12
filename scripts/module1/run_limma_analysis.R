#!/usr/bin/env Rscript
# User-directed exploratory limma inference on the supplied quantile-normalized matrix.
# The uploaded source units are unknown; log2(x+1) is an explicit analytic assumption,
# not an inferred or authenticated upstream transform. No samples/genes are silently removed.
options(stringsAsFactors=FALSE,warn=1)
set.seed(194331)
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(pheatmap);library(limma)})
root='/workspace/GSE194331-RNAseq-R';setwd(root)
for(d in c('results/module1_deg','results/module1_severity','results/module1_downstream','reports/module1'))dir.create(d,recursive=TRUE,showWarnings=FALSE)
zip='GSE示例_归一后.zip';ent=unzip(zip,list=TRUE)$Name
tab=read.csv(unz(zip,ent,encoding='UTF-8'),check.names=FALSE)
gene=as.character(tab[[1]]);sid=names(tab)[-1]
expr=as.matrix(tab[,-1,drop=FALSE]);storage.mode(expr)='double';rownames(expr)=gene
meta=fread('轻中重分组.csv',encoding='UTF-8')
stopifnot(identical(sid,meta$id),!anyDuplicated(gene),!anyDuplicated(sid),all(is.finite(expr)),all(expr>=0),nrow(meta)==119)
meta$group=factor(ifelse(meta$group=='H','HC',meta$group),levels=c('HC','MAP','MSAP','SAP'))
stopifnot(identical(as.integer(table(meta$group)),c(32L,57L,20L,10L)))
# The preceding file-to-file analysis independently recovered mean>=0.1 filtering and
# quantile normalization. Do not repeat either operation on this already-normalized matrix.
y=log2(expr+1)
fwrite(data.frame(gene_id=gene,as.data.frame(y,check.names=FALSE),check.names=FALSE),
       'results/module1_downstream/exploratory_log2_quantile_expression_plus1.tsv.gz',sep='\t')
write.csv(data.frame(parameter=c('input','input_scale','preprocessing','model','variance_moderation','design_covariates','seed',
    'descriptive_DEG_filter','null_hypothesis','status'),
    value=c(paste0(zip,'::',ent),'unlogged quantile-normalized, source units unknown',
    'existing filter/quantile normalization preserved; exploratory log2(expression+1); no new gene filtering',
    'limma linear model, contrast-specific binary group-only designs',
    'eBayes(robust=TRUE, trend=FALSE); sensitivity eBayes(robust=TRUE, trend=TRUE)',
    'group only: age/sex/batch unavailable','194331','BH padj<0.05 AND |estimated log2FC|>1; post-test descriptive filter only',
    'two-sided H0 log2FC=0','exploratory because source units and upstream origin remain unknown')),
    'results/module1_deg/limma_parameters.csv',row.names=FALSE)
writeLines(capture.output(sessionInfo()),'reports/module1/sessionInfo_limma.txt')
write.csv(data.frame(package=c('R','limma','data.table','ggplot2','pheatmap','statmod'),
    version=c(as.character(getRversion()),as.character(packageVersion('limma')),as.character(packageVersion('data.table')),
    as.character(packageVersion('ggplot2')),as.character(packageVersion('pheatmap')),as.character(packageVersion('statmod')))),
    'results/module1_deg/limma_package_versions.csv',row.names=FALSE)
contrasts_plan=list(
 A_AP_vs_HC=list(case='AP',ref='HC',groups=c('HC','MAP','MSAP','SAP'),role='primary'),
 B_MSAP_vs_MAP=list(case='MSAP',ref='MAP',groups=c('MAP','MSAP'),role='primary'),
 C_SAP_vs_MAP=list(case='SAP',ref='MAP',groups=c('MAP','SAP'),role='primary'),
 D_SAP_vs_MSAP=list(case='SAP',ref='MSAP',groups=c('MSAP','SAP'),role='exploratory'),
 E_MSAP_SAP_vs_MAP=list(case='severe',ref='MAP',groups=c('MAP','MSAP','SAP'),role='severity_sensitivity'))
colors=c(HC='#0072B2',MAP='#009E73',MSAP='#E69F00',SAP='#D55E00',AP='#CC79A7',severe='#CC79A7')
save_plot=function(p,path,w=7,h=5){ggsave(paste0(path,'.pdf'),p,width=w,height=h);ggsave(paste0(path,'.png'),p,width=w,height=h,dpi=300)}
summaries=list(); fits_primary=list(); all_results=list()
for(nm in names(contrasts_plan)){
    pl=contrasts_plan[[nm]]; sel=meta$group %in% pl$groups; ix=which(sel)
    g=as.character(meta$group[sel]);is_case=if(nm=='A_AP_vs_HC')g!='HC' else if(nm=='E_MSAP_SAP_vs_MAP')g!='MAP' else g==pl$case
    binary=factor(ifelse(is_case,pl$case,pl$ref),levels=c(pl$ref,pl$case))
    design=model.matrix(~0+binary);colnames(design)=levels(binary)
    stopifnot(qr(design)$rank==2L,nrow(design)==length(ix))
    contrast=makeContrasts(contrasts=paste0(pl$case,'-',pl$ref),levels=design)
    fit0=lmFit(y[,ix,drop=FALSE],design)
    fitc=contrasts.fit(fit0,contrast)
    fit=eBayes(fitc,robust=TRUE,trend=FALSE);fit_trend=eBayes(fitc,robust=TRUE,trend=TRUE)
    tt=topTable(fit,coef=1,number=Inf,sort.by='none',adjust.method='BH')
    ts=topTable(fit_trend,coef=1,number=Inf,sort.by='none',adjust.method='BH')
    stopifnot(nrow(tt)==nrow(y),nrow(ts)==nrow(y))
    se=fit$stdev.unscaled[,1]*sqrt(fit$s2.post);ci=qt(.975,df=fit$df.total)*se
    res=data.frame(gene_id=gene,gene_symbol=gene,gene_symbol_status='source gene label; annotation version unverified',
      baseMean=tt$AveExpr,AveExpr=tt$AveExpr,log2FoldChange=tt$logFC,standard_error=se,
      moderated_t=tt$t,pvalue=tt$P.Value,padj=tt$adj.P.Val,CI95_low=tt$logFC-ci,CI95_high=tt$logFC+ci,
      B=tt$B,trend_moderated_t=ts$t,trend_pvalue=ts$P.Value,trend_padj=ts$adj.P.Val,
      contrast=paste0(pl$case,' vs ',pl$ref),n_case=sum(is_case),n_reference=sum(!is_case),
      model_scale='log2(existing QN value + 1); exploratory assumption',check.names=FALSE)
    out=file.path('results/module1_deg',paste0(nm,'.csv.gz'))
    fwrite(res,out)
    all_results[[nm]]=res;fits_primary[[nm]]=fit
    ranked=res[is.finite(res$moderated_t),c('gene_id','gene_symbol','moderated_t')]
    ranked=ranked[order(-ranked$moderated_t,ranked$gene_id,method='radix'),]
    write.table(ranked,file.path('results/module1_downstream',paste0(nm,'_ranked_all_genes.tsv')),
                sep='\t',row.names=FALSE,quote=FALSE)
    sig=!is.na(res$padj)&res$padj<.05&abs(res$log2FoldChange)>1
    up=sum(sig&res$log2FoldChange>1);down=sum(sig&res$log2FoldChange< -1)
    summaries[[nm]]=data.frame(comparison=nm,contrast=paste0(pl$case,' vs ',pl$ref),role=pl$role,
        n_case=sum(is_case),n_reference=sum(!is_case),genes_total=nrow(res),genes_tested=sum(is.finite(res$pvalue)),
        genes_padj_lt_0_05=sum(res$padj<.05,na.rm=TRUE),descriptive_DEG_up=up,descriptive_DEG_down=down,
        DEG_rule='BH padj<0.05 and abs(log2FC)>1',logFC_null_test='eBayes H0=0; no treat/lfcThreshold test')
    stab=data.frame(gene_id=gene,standard_logFC=tt$logFC,trend_logFC=ts$logFC,
      standard_padj=tt$adj.P.Val,trend_padj=ts$adj.P.Val,
      standard_sig=!is.na(tt$adj.P.Val)&tt$adj.P.Val<.05&abs(tt$logFC)>1,
      trend_sig=!is.na(ts$adj.P.Val)&ts$adj.P.Val<.05&abs(ts$logFC)>1)
    write.csv(data.frame(metric=c('logFC_pearson','logFC_spearman','same_direction_fraction',
        'standard_DEG_count','trend_DEG_count','DEG_intersection','DEG_union'),
        value=c(cor(tt$logFC,ts$logFC),cor(tt$logFC,ts$logFC,method='spearman'),
            mean(sign(tt$logFC)==sign(ts$logFC)),sum(stab$standard_sig),sum(stab$trend_sig),
            sum(stab$standard_sig&stab$trend_sig),sum(stab$standard_sig|stab$trend_sig))),
        file.path('results/module1_deg',paste0(nm,'_method_stability.csv')),row.names=FALSE)
    write.csv(stab,file.path('results/module1_deg',paste0(nm,'_method_stability_by_gene.csv.gz')),row.names=FALSE)
    volcano=data.frame(logFC=res$log2FoldChange,minusLog10P=-log10(pmax(res$padj,.Machine$double.xmin)),
                       padj=res$padj,DEG=sig)
    save_plot(ggplot(volcano,aes(logFC,minusLog10P,color=DEG))+geom_point(size=.7,alpha=.55)+
        scale_color_manual(values=c('FALSE'='grey65','TRUE'='#D55E00'))+theme_bw()+
        labs(title=paste(nm,pl$case,'vs',pl$ref),subtitle='Exploratory limma on log2(QN expression + 1)',
             x='log2 fold change on transformed QN values',y='-log10(BH adjusted P)',color='padj<0.05 & |LFC|>1'),
        file.path('results/module1_deg',paste0(nm,'_volcano')),7,5)
    top=order(ifelse(is.na(res$padj),Inf,res$padj),-abs(res$log2FoldChange),gene)[seq_len(min(50,nrow(res)))]
    mat=y[top,ix,drop=FALSE]
    ann=data.frame(group=meta$group[sel],row.names=meta$id[sel])
    pdf(file.path('results/module1_deg',paste0(nm,'_top50_heatmap.pdf')),width=11,height=10)
    pheatmap(mat,scale='row',annotation_col=ann,annotation_colors=list(group=colors),
             show_colnames=FALSE,fontsize_row=5,main=paste(nm,'top 50 by exploratory limma p-value'));dev.off()
    png(file.path('results/module1_deg',paste0(nm,'_top50_heatmap.png')),width=3300,height=3000,res=300)
    pheatmap(mat,scale='row',annotation_col=ann,annotation_colors=list(group=colors),
             show_colnames=FALSE,fontsize_row=5,main=paste(nm,'top 50 by exploratory limma p-value'));dev.off()
}
summary=do.call(rbind,summaries)
write.csv(summary,'results/module1_deg/limma_DEG_summary.csv',row.names=FALSE)
meta_table=as.data.frame(meta)
meta_table$analysis_group=ifelse(meta_table$group=='H','HC',as.character(meta_table$group))
trend_sample=meta_table$analysis_group %in% c('MAP','MSAP','SAP')
apmeta=meta_table[trend_sample,,drop=FALSE];apidx=which(trend_sample)
apmeta$severity=factor(apmeta$analysis_group,levels=c('MAP','MSAP','SAP'))
apmeta$severity_score=as.integer(apmeta$severity)
stopifnot(identical(as.integer(table(apmeta$severity)),c(57L,20L,10L)))
# Primary ordinal score association: per one-category step, with no linear-dose causal claim.
design_trend=model.matrix(~severity_score,apmeta);stopifnot(qr(design_trend)$rank==2L)
ftrend=eBayes(lmFit(y[,apidx,drop=FALSE],design_trend),robust=TRUE,trend=FALSE)
ttrend=topTable(ftrend,coef='severity_score',number=Inf,sort.by='none',adjust.method='BH')
se=ftrend$stdev.unscaled[,2]*sqrt(ftrend$s2.post);ci=qt(.975,df=ftrend$df.total)*se
trend_res=data.frame(gene_id=gene,gene_symbol=gene,gene_symbol_status='source gene label; annotation version unverified',
  AveExpr=ttrend$AveExpr,log2FC_per_score_step=ttrend$logFC,standard_error=se,moderated_t=ttrend$t,
  pvalue=ttrend$P.Value,padj=ttrend$adj.P.Val,CI95_low=ttrend$logFC-ci,CI95_high=ttrend$logFC+ci,
  score_mapping='MAP=1, MSAP=2, SAP=3',test='moderated limma t for ordinal score slope; H0 slope=0',
  interpretation='association per score step; not a biological linear-dose proof')
fwrite(trend_res,'results/module1_severity/ordinal_score_trend_all_genes.csv.gz')
trend_rank=trend_res[is.finite(trend_res$moderated_t),c('gene_id','gene_symbol','moderated_t')]
trend_rank=trend_rank[order(-trend_rank$moderated_t,trend_rank$gene_id,method='radix'),]
write.table(trend_rank,'results/module1_downstream/AP_ordinal_trend_ranked_all_genes.tsv',sep='\t',row.names=FALSE,quote=FALSE)
# Categorical AP model: global 2-df moderated F test, sensitivity to linear-score assumption.
design_cat=model.matrix(~severity,apmeta);stopifnot(qr(design_cat)$rank==3L)
fcat=eBayes(lmFit(y[,apidx,drop=FALSE],design_cat),robust=TRUE,trend=FALSE)
cat_tt=topTable(fcat,coef=2:3,number=Inf,sort.by='none',adjust.method='BH')
cat_se=data.frame(gene_id=gene,gene_symbol=gene,AveExpr=cat_tt$AveExpr,F=cat_tt$F,
                  pvalue=cat_tt$P.Value,padj=cat_tt$adj.P.Val,n_MAP=57,n_MSAP=20,n_SAP=10,
                  test='global 2-df moderated F test: category model ~ severity')
fwrite(cat_se,'results/module1_severity/category_global_F_sensitivity_all_genes.csv.gz')
write.csv(data.frame(model=c('ordinal_score','AP_category_global_F'),n=c(87,87),rank=c(qr(design_trend)$rank,qr(design_cat)$rank),
  residual_df=c(nrow(design_trend)-qr(design_trend)$rank,nrow(design_cat)-qr(design_cat)$rank),
  significant_BH_lt_0_05=c(sum(trend_res$padj<.05,na.rm=TRUE),sum(cat_se$padj<.05,na.rm=TRUE))),
  'results/module1_severity/severity_model_summary.csv',row.names=FALSE)
writeLines(c('# Severity model notes',
    '- AP-only samples: MAP 57, MSAP 20, SAP 10. Input values are log2(existing QN expression + 1), explicitly exploratory due unknown source units.',
    '- Primary trend design: ~ severity_score with MAP=1, MSAP=2, SAP=3; limma robust empirical Bayes moderated t on per-step slope.',
    '- This tests association with an ordinal score under a constant per-step log-expression slope; it does not establish a linear biological dose response.',
    '- Sensitivity model: ~ severity categorical with a global 2-df moderated F test. Comparisons B, C, D provide pairwise AP category contrasts.',
    '- All 87 AP samples retained. SAP has only 10 samples; effect sizes/95% intervals accompany the trend and pairwise tables; nonsignificance is not absence.',
    '- No clinical covariates or batches were available.'),'results/module1_severity/model_notes.md')
# Stable gene sets and result summaries.
write.csv(summary,'results/module1_deg/DEG_direction_and_effect_summary.csv',row.names=FALSE)
cat('EXPLORATORY_LIMMA_COMPLETE; log2(QN+1), five contrasts; 19,645 genes; severity ordinal plus categorical sensitivity.\n')
print(summary,row.names=FALSE)
cat('Severity ordinal BH<0.05:',sum(trend_res$padj<.05,na.rm=TRUE),'; category global BH<0.05:',sum(cat_se$padj<.05,na.rm=TRUE),'\n')
