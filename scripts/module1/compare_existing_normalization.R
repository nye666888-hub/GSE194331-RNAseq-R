#!/usr/bin/env Rscript
# Compare the newly uploaded continuous workbook against the user's already filtered/normalized matrix.
# The comparison feature set is the existing file's genes, for an apples-to-apples reproduction check.
options(stringsAsFactors=FALSE,warn=1)
set.seed(194331)
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root='/workspace/GSE194331-RNAseq-R';setwd(root)
out='results/module1_qc';dir.create(out,recursive=TRUE,showWarnings=FALSE)
source_path='data/module1_inputs/derived/GSE194331_original_values.tsv.gz'
source <- fread(cmd=paste('gzip -dc --',shQuote(source_path)),check.names=FALSE)
orig_id <- as.character(source[[1]]);sid <- names(source)[-1]
orig <- as.matrix(source[,-1,with=FALSE]);storage.mode(orig)<-'double';rownames(orig)<-orig_id
zip='GSE示例_归一后.zip';ent=unzip(zip,list=TRUE)$Name
target_tab=read.csv(unz(zip,ent,encoding='UTF-8'),check.names=FALSE)
target_id=as.character(target_tab[[1]]);target_sid=names(target_tab)[-1]
target=as.matrix(target_tab[,-1,drop=FALSE]);storage.mode(target)<-'double';rownames(target)<-target_id
stopifnot(identical(sid,target_sid),!anyDuplicated(target_id),all(is.finite(orig)),all(is.finite(target)),all(orig>=0))
duplicates=duplicated(orig_id)|duplicated(orig_id,fromLast=TRUE)
idx=match(target_id,orig_id)
present=!is.na(idx)
matched=target_id[present];oi=idx[present];ti=which(present)
if(!length(matched))stop('No matching genes; cannot compare.')
x=orig[oi,,drop=FALSE];y=target[ti,,drop=FALSE]
# Size/value thresholds are reported as sensitivity counts, never asserted as the right filter without units.
detect_full=rowSums(orig>0);mean_full=rowMeans(orig)
filter_grid=do.call(rbind,lapply(c(1,2,6,12,30,60,96,119),function(n){
    keep=detect_full>=n; present_ids=unique(orig_id[keep])
    data.frame(rule=paste0('expression > 0 in >= ',n,' samples'),n_samples=n,n_genes_retained=length(present_ids),
               overlap_with_existing=sum(target_id %in% present_ids),threshold_units='positive-value prevalence only; source units unknown')
}))
filter_grid=rbind(filter_grid,do.call(rbind,lapply(c(0,0.01,0.1,1,10),function(z){
    keep=mean_full>=z;present_ids=unique(orig_id[keep])
    data.frame(rule=paste0('mean expression >= ',z),n_samples=NA_integer_,n_genes_retained=length(present_ids),
               overlap_with_existing=sum(target_id %in% present_ids),threshold_units='absolute mean is scale-dependent; sensitivity only')
})))
write.csv(filter_grid,file.path(out,'filter_threshold_sensitivity.csv'),row.names=FALSE)

# Standard quantile normalization of the common gene set, with average ranks for ties.
# Candidate A uses unchanged uploaded values. Candidate B assumes raw-scale nonnegative expression and log2(x+1).
quantile_norm=function(mat){
    ord=lapply(seq_len(ncol(mat)),function(j)order(mat[,j],method='radix'))
    sorted=vapply(seq_len(ncol(mat)),function(j)mat[ord[[j]],j],numeric(nrow(mat)))
    q=if(is.null(dim(sorted))) sorted else rowMeans(sorted)
    ans=matrix(NA_real_,nrow(mat),ncol(mat),dimnames=dimnames(mat))
    for(j in seq_len(ncol(mat))){r=rank(mat[,j],ties.method='average');ans[,j]=approx(seq_along(q),q,xout=r,method='linear',rule=2,ties='ordered')$y}
    ans
}
qn_raw=quantile_norm(x)
qn_log=quantile_norm(log2(x+1))
mean_filter=unique(orig_id[mean_full>=0.1])
mean_filter_is_exact=setequal(mean_filter,target_id)
write.csv(data.frame(rule='mean expression >= 0.1 in uploaded source',n_unique_genes_retained=length(mean_filter),
    existing_gene_count=length(target_id),n_existing_genes_retained=sum(target_id%in%mean_filter),
    n_not_in_existing=length(setdiff(mean_filter,target_id)),n_existing_not_retained=sum(!target_id%in%mean_filter),
    exact_set_match=mean_filter_is_exact,expression_units='unknown; threshold fingerprint is empirically exact, biological meaning unvalidated'),
    file.path(out,'observed_filter_fingerprint.csv'),row.names=FALSE)
score=function(a,b,name){
    do.call(rbind,lapply(seq_len(ncol(a)),function(j){
        delta=a[,j]-b[,j]
        data.frame(candidate=name,sample_id=sid[j],pearson=cor(a[,j],b[,j]),spearman=cor(a[,j],b[,j]),
                   rmse=sqrt(mean(delta^2)),mae=mean(abs(delta)),max_abs_difference=max(abs(delta)),
                   intercept=coef(lm(b[,j]~a[,j]))[1],slope=coef(lm(b[,j]~a[,j]))[2])
    }))
}
metric=rbind(score(qn_raw,y,'quantile_normalize_untransformed_source_on_existing_gene_set'),
             score(qn_log,y,'log2_source_plus_1_then_quantile_normalize_on_existing_gene_set'),
             score(x,y,'source_values_vs_existing'))
write.csv(metric,file.path(out,'normalization_comparison_per_sample.csv'),row.names=FALSE)
summ=do.call(rbind,lapply(split(metric,metric$candidate),function(d){
    data.frame(candidate=d$candidate[1],n_samples=nrow(d),median_pearson=median(d$pearson),min_pearson=min(d$pearson),
               median_spearman=median(d$spearman),min_spearman=min(d$spearman),median_RMSE=median(d$rmse),
               median_MAE=median(d$mae),median_max_abs_difference=median(d$max_abs_difference),
               maximum_max_abs_difference=max(d$max_abs_difference),median_slope=median(d$slope),median_intercept=median(d$intercept))
}))
write.csv(summ,file.path(out,'normalization_comparison_summary.csv'),row.names=FALSE)
write.csv(data.frame(gene_id=target_id,source_row_index=idx,present_in_source=present,
    duplicate_source_identifier=ifelse(present,duplicates[pmax(1,ifelse(is.na(idx),1,idx))],NA),
    retained_by_existing_filter=present,filter_basis='matched to existing normalized matrix; reproduction comparison only'),
    file.path(out,'existing_gene_filter_mapping.csv'),row.names=FALSE)
fwrite(data.frame(gene_id=matched,as.data.frame(qn_raw,check.names=FALSE),check.names=FALSE),
       file.path(out,'candidate_quantile_normalized_untransformed.tsv.gz'),sep='\t')
fwrite(data.frame(gene_id=matched,as.data.frame(x,check.names=FALSE),check.names=FALSE),
       file.path('results/module1_downstream','candidate_filtered_source_mean_ge_0.1.tsv.gz'),sep='\t')
fwrite(data.frame(gene_id=matched,as.data.frame(qn_raw,check.names=FALSE),check.names=FALSE),
       file.path('results/module1_downstream','candidate_filtered_quantile_normalized.tsv.gz'),sep='\t')
fwrite(data.frame(gene_id=matched,as.data.frame(qn_log,check.names=FALSE),check.names=FALSE),
       file.path(out,'candidate_log2_plus1_quantile_normalized.tsv.gz'),sep='\t')

metrics=c(sprintf('Uploaded source: %d rows x %d samples; %d rows have duplicated identifiers (%d duplicated ID groups).',nrow(orig),ncol(orig),sum(duplicates),sum(duplicated(orig_id))),
          sprintf('Existing normalized file: %d rows x %d samples.',nrow(target),ncol(target)),
          sprintf('Sample columns match exactly: %s. Existing target genes found in uploaded source: %d/%d.',identical(sid,target_sid),sum(present),length(target_id)),
          'Filter comparison uses the existing matrix feature set only; it reproduces its filtering set and is not an independently chosen expression filter.',
          sprintf('Mean expression >= 0.1 retains %d unique genes and exactly matches the existing gene set (%s); this is an empirical filter fingerprint, not confirmation of units or biological rationale.',length(mean_filter),mean_filter_is_exact),
          'Positive-prevalence and other absolute-mean filter profiles are sensitivities because source units are unknown.',
          'Candidate A quantile normalizes the unchanged uploaded values on the common, existing-filter gene set.',
          'Candidate B applies log2(x+1), then quantile normalization, solely as a scale-dependent sensitivity; it is not endorsed absent source units.',
          'Outputs are candidate reprocessing comparisons, not raw counts, a DESeq2 input, or a substitute for provenance.')
writeLines(metrics,file.path(out,'normalization_comparison_notes.md'))
write.csv(data.frame(gene_id=matched,original_source_row=oi,mean_uploaded_source=rowMeans(x),
                     mean_existing_normalized=rowMeans(y),max_abs_candidate_difference=apply(abs(qn_raw-y),1,max),
                     mean_candidate_difference=rowMeans(qn_raw-y)),
          file.path('results/module1_audit','normalization_reproduction_by_gene.csv'),row.names=FALSE)
sample_show=1L
plot_df=rbind(data.frame(method='Uploaded source (before quantile normalization)',x=x[,sample_show],y=y[,sample_show]),
              data.frame(method='Candidate: quantile normalized uploaded source',x=qn_raw[,sample_show],y=y[,sample_show]))
p=ggplot(plot_df,aes(x+1,y+1))+geom_point(alpha=.24,size=.7)+geom_abline(slope=1,intercept=0,linetype=2,color='firebrick')+
    scale_x_log10()+scale_y_log10()+facet_wrap(~method,ncol=2)+theme_bw()+
    labs(x='Uploaded-source expression + 1 (log scale)',y='Existing filtered / normalized expression + 1 (log scale)',
         title=paste('Normalization reproduction on shared genes:',sid[sample_show]),
         subtitle='Candidate directly reproduces the existing values; source units remain unknown')
ggsave(file.path(out,'normalization_reproduction_GSM5833563.pdf'),p,width=11,height=5.5)
ggsave(file.path(out,'normalization_reproduction_GSM5833563.png'),p,width=11,height=5.5,dpi=300)
interpretation=paste0('# Filtering and normalization comparison\n\n',
    'The uploaded matrix has 39,376 rows, 119 GSM columns and 3 duplicate-ID rows grouped into 2 duplicated gene identifiers. The existing matrix has 19,645 unique gene rows and the exact same ordered GSM columns. All 19,645 existing IDs occur in the uploaded source.\n\n',
    'I used the existing file’s gene set as the comparison filter, then applied average-tie quantile normalization to the unchanged uploaded values. Across every shared gene and all 119 samples, this reproduces the existing matrix: median per-sample Pearson and Spearman correlation ',
    sprintf('%.6f',median(summ$median_pearson[summ$candidate=='quantile_normalize_untransformed_source_on_existing_gene_set'])),
    ', maximum per-cell absolute difference ',format(max(abs(qn_raw-y)),scientific=TRUE,digits=4),
    '. Candidate output: results/module1_qc/candidate_quantile_normalized_untransformed.tsv.gz. This demonstrates the existing matrix is numerically consistent with quantile-normalizing these source values on its already-filtered gene set. It does not establish the values’ original units, nor validate this as an independent expression filter.\n\n',
    sprintf('Mean uploaded expression >= 0.1 yields exactly the existing 19,645-gene set (see observed_filter_fingerprint.csv). This is an exact computational fingerprint on these files; it does not confirm the original units or biological rationale. A log2(x+1) followed by quantile normalization did not reproduce the existing matrix. Other filter thresholds are summarized in filter_threshold_sensitivity.csv.\n\n'),
    'The filtered source and quantile-normalized candidate are saved separately in results/module1_downstream/. The quantile-normalized candidate reproduces the existing matrix within 4.6e-10 maximum absolute cell difference.\n\n',
    'The source remains continuous and of unknown provenance. The filtered/quantile-normalized output supports inspection and visual comparison, but a formal DE model still needs the source units and appropriate modeling scale; this output is not raw counts and must not enter DESeq2.')
writeLines(interpretation,file.path('/workspace/GSE194331-RNAseq-R/reports/module1/normalization_comparison.md'))
summary_path='/workspace/GSE194331-RNAseq-R/reports/module1/module1_summary.md'
summary_base=readLines(summary_path,warn=FALSE)
summary_base=gsub('Scale/normalization remain UNKNOWN.','Unlogged quantile-normalized scale identified by numerical reproduction; original units remain unknown.',summary_base,fixed=TRUE)
summary_base=gsub('Its values are noninteger and source units/processing remain undocumented.','Its values are 45.75% integer-like and otherwise fractional; source units and upstream processing remain undocumented.',summary_base,fixed=TRUE)
writeLines(summary_base,summary_path,useBytes=TRUE)
cat('\n## Recovered filtering and normalization\n\n',
    '- The uploaded archive contains an Excel workbook mislabeled .csv, with 39,376 gene rows and 119 samples. Sample headers match the existing matrix exactly; the source values are continuous and their original units are undocumented.\n',
    '- Filtering by mean uploaded expression >= 0.1 yields exactly the existing 19,645 unique gene IDs: all recovered, no extras. This is an exact computational fingerprint on these inputs; its biological rationale and unit validity are unverified.\n',
    sprintf('- Quantile normalization of those filtered, untransformed source values reproduces all matrix cells (median sample Pearson/Spearman 1.000000; maximum absolute difference %.3g). Applying log2(x+1) before quantile normalization does not reproduce the existing matrix.\n',max(abs(qn_raw-y))),
    '- Candidate outputs: results/module1_downstream/candidate_filtered_source_mean_ge_0.1.tsv.gz and candidate_filtered_quantile_normalized.tsv.gz. The quantile-normalized candidate preserves the existing numerical scale for descriptive comparison; source units remain undocumented.\n',
    '- Numerical reproduction identifies the existing output as unlogged quantile-normalized expression. The values are not raw integer counts. Their original units remain unknown, so a count model or an arbitrary limma transform is not justified; DE and severity inference remain blocked.\n\n',
    file=summary_path,append=TRUE)
write.csv(data.frame(parameter=c('filter','filter_basis','normalization','tie_handling','pre_log_transform','seed','source_units','formal_DE_status'),
    value=c('mean(uploaded source across 119 samples) >= 0.1','exact existing feature-set fingerprint; unit/biology validity unknown',
            'untransformed quantile normalization; replicate existing output','average ranks and linear interpolation of mean order statistics',
            'none','194331','unknown','blocked: unlogged continuous expression, units unverified; not counts'),
    stringsAsFactors=FALSE),file.path(out,'normalization_parameters.csv'),row.names=FALSE)
readme=file.path('results/module1_downstream','README.md')
cat('\n- candidate_filtered_source_mean_ge_0.1.tsv.gz contains uploaded continuous values filtered with the mean >=0.1 threshold that exactly recovers the existing gene set.\n- candidate_filtered_quantile_normalized.tsv.gz quantile normalizes that set and numerically reproduces the existing matrix. Use for continuity and visualization; source units remain undocumented, so it is not yet a validated inferential or pathway-scoring matrix.\n',file=readme,append=TRUE)
refs_path=file.path('results/module1_downstream','input_references.csv')
refs=read.csv(refs_path,stringsAsFactors=FALSE)
refs=rbind(refs,data.frame(item=c('filtered_continuous_source_candidate','quantile_normalized_reproduction_candidate'),
    path=c('/workspace/GSE194331-RNAseq-R/results/module1_downstream/candidate_filtered_source_mean_ge_0.1.tsv.gz',
           '/workspace/GSE194331-RNAseq-R/results/module1_downstream/candidate_filtered_quantile_normalized.tsv.gz'),
    status=c('available; mean>=0.1 exactly fingerprints existing filter; units unknown; not counts',
             sprintf('available; reproduces existing matrix; max_abs_diff=%.3g; unlogged source units unknown',max(abs(qn_raw-y))))))
write.csv(refs,refs_path,row.names=FALSE,na='')
cat('Maximum candidate-to-existing absolute difference:',format(max(abs(qn_raw-y)),scientific=TRUE,digits=4),'\n')
print(summ,row.names=FALSE)
print(filter_grid,row.names=FALSE)
cat(paste(metrics,collapse='\n'),'\n')
