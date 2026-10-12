#!/usr/bin/env Rscript
# Audit-first Module 1. This script intentionally cannot fit DE models on unknown scales.
options(stringsAsFactors=FALSE, warn=1)
set.seed(194331)
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(pheatmap)})
repo <- normalizePath('/workspace/GSE194331-RNAseq-R')
setwd(repo)
dirs <- c(audit='results/module1_audit', qc='results/module1_qc', deg='results/module1_deg',
          severity='results/module1_severity', downstream='results/module1_downstream', report='reports/module1')
for (d in dirs) dir.create(d,recursive=TRUE,showWarnings=FALSE)
write_table <- function(x,path) fwrite(as.data.frame(x),path,na='NA',quote=TRUE)
write_lines <- function(x,path) writeLines(x,path,useBytes=TRUE)
write_lines(capture.output(sessionInfo()),file.path(dirs['report'],'sessionInfo.txt'))
packages <- c('R','data.table','ggplot2','pheatmap')
versions <- c(as.character(getRversion()),vapply(packages[-1],function(p) as.character(packageVersion(p)),character(1)))
write_table(data.frame(package=packages,version=versions),file.path(dirs['audit'],'package_versions.csv'))
write_lines(c('seed: 194331','input_scale: UNKNOWN; normalized-continuous candidate only',
              'normalization: NONE','additional_gene_filtering: NONE','sample_exclusion: NONE',
              'PCA: supplied-scale, all nonconstant genes, center=TRUE, scale.=FALSE, rank.=10',
              'correlation: Pearson and Spearman across all supplied genes',
              'distance: Euclidean on supplied-scale expression; 1-Spearman secondary',
              'plots: log10(value+1) axis presentation only for nonnegative distribution summaries',
              'outlier_review: within-group median correlation < median-3*MAD OR within-group PCA distance > median+3*MAD',
              'DEG description if unblocked: BH padj<0.05 AND abs(estimated log2FC)>1; NOT lfcThreshold testing',
              'formal_DE_and_severity: BLOCKED, no valid scale/provenance established'),file.path(dirs['report'],'parameters.txt'))

# Read supplied matrix without modifying/extracting the original archive.
zip_path <- file.path(repo,'GSE示例_归一后.zip')
entries <- unzip(zip_path,list=TRUE)
stopifnot(nrow(entries)==1L,grepl('\\.csv$',entries$Name))
tab <- read.csv(unz(zip_path,entries$Name,encoding='UTF-8'),check.names=FALSE)
gene <- as.character(tab[[1]])
sample_ids <- names(tab)[-1]
if (!all(vapply(tab[-1],is.numeric,logical(1)))) stop('Nonnumeric expression columns: cannot run numeric QC.')
m <- as.matrix(tab[-1]); storage.mode(m) <- 'double'; rownames(m) <- gene
rm(tab); gc(FALSE)
metadata <- fread('轻中重分组.csv',encoding='UTF-8',check.names=FALSE)
stopifnot(all(c('id','group') %in% names(metadata)))
metadata$id <- as.character(metadata$id); metadata$group <- as.character(metadata$group)
missing_meta <- setdiff(sample_ids,metadata$id); missing_matrix <- setdiff(metadata$id,sample_ids)
write_table(data.frame(issue=c(rep('expression_sample_missing_metadata',length(missing_meta)),rep('metadata_sample_missing_expression',length(missing_matrix))),
                       sample_id=c(missing_meta,missing_matrix)),file.path(dirs['audit'],'sample_mismatches.csv'))
if(anyDuplicated(sample_ids) || anyDuplicated(metadata$id) || length(missing_meta) || length(missing_matrix))
    stop('Duplicate/missing sample IDs; preserve input and stop dependent analyses.')
if(anyNA(m) || any(!is.finite(m))) stop('Missing/nonfinite expression found; QC needs explicit missing-data decision.')
metadata <- as.data.frame(metadata)[match(sample_ids,metadata$id),,drop=FALSE]
metadata$group_original <- metadata$group
# H->HC follows the latest user-supplied cohort definition, not expression inference.
metadata$group <- ifelse(metadata$group_original=='H','HC',metadata$group_original)
metadata$group <- factor(metadata$group,levels=c('HC','MAP','MSAP','SAP'))
stopifnot(!anyNA(metadata$group),ncol(m)==119,nrow(metadata)==119,
          identical(as.integer(table(metadata$group)),c(32L,57L,20L,10L)))
metadata$sample_id <- metadata$id
metadata$matrix_column <- seq_along(sample_ids)
metadata$metadata_source <- '轻中重分组.csv'
metadata$id_alignment_status <- 'one_to_one_verified_against_supplied_matrix'
metadata$group_source <- 'local_user_supplied_metadata; NOT_expression_inferred'
metadata$group_label_note <- ifelse(metadata$group_original=='H','H->HC mapping from the explicitly user-specified HC cohort (n=32)','unchanged')
metadata$geo_independent_verification <- 'NOT_COMPLETED: GEO network proxy HTTP 403'
metadata$analysis_group_confirmation <- ifelse(metadata$group_original=='H','HC_mapping_from_user_specified_cohort_definition','local_metadata_only')
metadata$AP_severity_code_candidate <- ifelse(metadata$group=='HC',NA_integer_,match(metadata$group,c('MAP','MSAP','SAP')))
write_table(metadata,file.path(dirs['audit'],'sample_metadata_verified.csv'))
write_table(metadata,file.path(dirs['downstream'],'sample_metadata_verified.csv'))

integer_fraction <- mean(abs(m-round(m))<=1e-8)
zero_frac <- mean(m==0)
constant <- apply(m,1,var)==0
id_type <- if(mean(grepl('^ENSG[0-9]+(\\.[0-9]+)?$',gene))>.9) 'Ensembl_gene_candidate' else if(mean(grepl('^[0-9]+$',gene))>.9) 'Entrez_candidate' else 'gene_symbol_like; no authoritative annotation version'
source_path <- file.path(repo,'data/module1_inputs/derived/GSE194331_original_values.tsv.gz')
source_present <- file.exists(source_path)
source_note <- 'Newly uploaded source archive was not available.'
if(source_present){
    source_tab <- fread(cmd=paste('gzip -dc --',shQuote(source_path)),encoding='UTF-8',check.names=FALSE)
    source_gene <- as.character(source_tab[[1]])
    source_samples <- names(source_tab)[-1]
    source_mat <- as.matrix(source_tab[,-1,with=FALSE]);storage.mode(source_mat)<-'double';rownames(source_mat)<-source_gene
    if(any(!is.finite(source_mat))) stop('Non-finite values in newly uploaded source; stop dependent audit.')
    source_match <- identical(source_samples,sample_ids)
    common_genes <- intersect(gene,source_gene)
    if(source_match && length(common_genes)){
        oi <- match(common_genes,source_gene); ni <- match(common_genes,gene)
        per_sample <- do.call(rbind,lapply(seq_along(sample_ids),function(j){
            a <- source_mat[oi,j];b <- m[ni,j]
            data.frame(sample_id=sample_ids[j],n_common_genes=length(common_genes),
                       pearson=cor(a,b,method='pearson'),spearman=cor(a,b,method='spearman'),
                       source_median=median(a),normalized_median=median(b),
                       source_zero_fraction=mean(a==0),normalized_zero_fraction=mean(b==0))
        }))
        write_table(per_sample,file.path(dirs['audit'],'source_vs_normalized_per_sample.csv'))
        source_note <- sprintf('Compared %d shared genes across all %d exactly aligned samples; both matrices remain continuous with unknown units/provenance.',length(common_genes),length(sample_ids))
    } else source_note <- sprintf('Uploaded source has %d genes and %d samples; sample headers identical=%s.',nrow(source_mat),ncol(source_mat),source_match)
    source_integer <- mean(abs(source_mat-round(source_mat))<=1e-8)
    source_audit <- data.frame(metric=c('source_archive_member','source_actual_file_format','n_genes','n_samples','sample_header_identical_to_normalized',
        'missing_values','nonfinite_values','negative_values','zero_fraction','integer_like_fraction_tolerance_1e-8','minimum','maximum','median',
        'duplicate_gene_id_rows','missing_gene_ids','GSM_id_valid_fraction','gene_id_type_candidate','scale','raw_counts_eligible','normalization_or_generation_provenance'),
        value=c('GSE194331.rar::GSE194331.csv','OOXML Excel workbook (mislabelled .csv)',nrow(source_mat),ncol(source_mat),source_match,
        sum(is.na(source_mat)),sum(!is.finite(source_mat)),sum(source_mat<0),mean(source_mat==0),source_integer,min(source_mat),max(source_mat),median(source_mat),
        sum(duplicated(source_gene)),sum(is.na(source_gene)|!nzchar(source_gene)),mean(grepl('^GSM[0-9]+$',source_samples)),
        if(mean(grepl('^ENSG[0-9]+(\\.[0-9]+)?$',source_gene))>.9) 'Ensembl_gene_candidate' else if(mean(grepl('^[0-9]+$',source_gene))>.9) 'Entrez_candidate' else 'gene_symbol_like; no authoritative annotation version',
        'UNKNOWN_CONTINUOUS_SOURCE',FALSE,'NOT_FOUND'))
    write_table(source_audit,file.path(dirs['audit'],'original_source_matrix_audit.csv'))
    write_table(data.frame(original_gene_id=source_gene,duplicate_id=duplicated(source_gene)|duplicated(source_gene,fromLast=TRUE),
                           gene_id_type='gene_symbol_like; no authoritative mapping',annotation_status='UNVERIFIED'),
                file.path(dirs['audit'],'original_source_gene_ids.csv'))
    rm(source_tab,source_mat);gc(FALSE)
}
audit <- data.frame(metric=c('archive_member','n_genes','n_samples','n_metadata_samples','missing_values','nonfinite_values','negative_values',
                             'zero_values','zero_fraction','integer_like_fraction_tolerance_1e-8','noninteger_values',
                             'min','max','median','duplicate_gene_id_rows','missing_gene_ids','constant_genes',
                             'duplicate_sample_column_names','duplicate_metadata_ids','GSM_id_valid_fraction','gene_id_type_candidate',
                             'expression_scale','raw_counts_eligible','processing_provenance','age_available','sex_available','batch_available'),
                    value=c(entries$Name,nrow(m),ncol(m),nrow(metadata),sum(is.na(m)),sum(!is.finite(m)),sum(m<0),sum(m==0),
                            zero_frac,integer_fraction,sum(abs(m-round(m))>1e-8),min(m),max(m),median(m),sum(duplicated(gene)),
                            sum(is.na(gene)|!nzchar(gene)),sum(constant),anyDuplicated(sample_ids),anyDuplicated(metadata$id),
                            mean(grepl('^GSM[0-9]+$',sample_ids)),id_type,'UNKNOWN: normalized continuous candidate',
                            'NO: noninteger values and no raw-count provenance','NOT_FOUND',FALSE,FALSE,FALSE))
write_table(audit,file.path(dirs['audit'],'matrix_audit.csv'))
write_table(data.frame(group=c('HC','MAP','MSAP','SAP'),original_label=c('H','MAP','MSAP','SAP'),
                       expected=c(32,57,20,10),observed=as.integer(table(metadata$group)),match=TRUE),file.path(dirs['audit'],'group_counts.csv'))

# Distribution statistics refer to continuous supplied values, NOT sequencing library sizes.
stats <- do.call(rbind,lapply(seq_len(ncol(m)),function(j){
    a <- m[,j]; q <- quantile(a,c(0,.01,.25,.5,.75,.99,1),names=FALSE)
    data.frame(sample_id=sample_ids[j],group=metadata$group[j],sum_supplied_expression=sum(a),mean=mean(a),sd=sd(a),
               min=q[1],q01=q[2],q25=q[3],median=q[4],q75=q[5],q99=q[6],max=q[7],zero_fraction=mean(a==0),
               integer_like_fraction=mean(abs(a-round(a))<=1e-8),missing_values=sum(is.na(a)),negative_values=sum(a<0))
}))
write_table(stats,file.path(dirs['audit'],'sample_expression_statistics.csv'))

# Explicit identical-vector check, not a claim that biological donors are unique.
pearson <- cor(m,method='pearson'); spearman <- cor(m,method='spearman')
candidates <- which(upper.tri(pearson)&pearson>1-1e-12,arr.ind=TRUE)
dups <- data.frame(sample_1=character(),sample_2=character(),identical_expression=logical())
if(nrow(candidates)) for(k in seq_len(nrow(candidates))){
    ij <- candidates[k,]; if(identical(unname(m[,ij[1]]),unname(m[,ij[2]])))
    dups <- rbind(dups,data.frame(sample_1=sample_ids[ij[1]],sample_2=sample_ids[ij[2]],identical_expression=TRUE))
}
write_table(dups,file.path(dirs['audit'],'identical_sample_vectors.csv'))
write_table(data.frame(original_gene_id=gene,gene_symbol_candidate=if(id_type=='gene_symbol_like; no authoritative annotation version') gene else NA_character_,
                       gene_id_type=id_type,duplicate_id=duplicated(gene)|duplicated(gene,fromLast=TRUE),
                       annotation_status='UNVERIFIED_IDENTITY_ONLY; not an authoritative mapping'),file.path(dirs['downstream'],'gene_id_symbol_mapping.csv'))
write_table(data.frame(gene_id=gene,is_constant_for_PCA=constant,retained_in_audit=TRUE,retained_in_formal_DE='DE_not_run'),file.path(dirs['audit'],'gene_status.csv'))

# Candidate model matrix and estimability can be checked without fitting any expression model.
X <- model.matrix(~0+group,metadata); colnames(X) <- sub('^group','',colnames(X))
contrast <- cbind(A_AP_vs_HC=c(-1,57/87,20/87,10/87),B_MSAP_vs_MAP=c(0,-1,1,0),
                  C_SAP_vs_MAP=c(0,-1,0,1),D_SAP_vs_MSAP_exploratory=c(0,0,-1,1),
                  E_MSAP_SAP_vs_MAP=c(0,-1,20/30,10/30))
rownames(contrast) <- colnames(X)
write_table(data.frame(sample_id=sample_ids,X,check.names=FALSE),file.path(dirs['audit'],'candidate_design_matrix.csv'))
write_table(data.frame(coefficient=rownames(contrast),contrast,check.names=FALSE),file.path(dirs['audit'],'candidate_contrasts.csv'))
write_table(data.frame(model=c('four_category_group_only','AP_only_ordered_1_2_3','AP_only_category'),
                       n_samples=c(nrow(X),87,87),n_coefficients=c(ncol(X),2,3),rank=c(qr(X)$rank,2,3),
                       residual_df=c(nrow(X)-qr(X)$rank,85,84),
                       status='algebraically_estimable; biological/input validity NOT established'),file.path(dirs['audit'],'design_estimability.csv'))
AP <- droplevels(metadata[metadata$group!='HC',])
Xtrend <- model.matrix(~AP_severity_code_candidate,AP)
Xcat <- model.matrix(~group,AP)
stopifnot(qr(X)$rank==4,qr(Xtrend)$rank==2,qr(Xcat)$rank==3)
write_table(data.frame(sample_id=AP$id,Xtrend,check.names=FALSE),file.path(dirs['severity'],'candidate_numeric_design_matrix.csv'))
write_table(data.frame(sample_id=AP$id,Xcat,check.names=FALSE),file.path(dirs['severity'],'candidate_category_design_matrix.csv'))
write_lines(c('# Recommended designs — conditional, NOT fitted',
    '- Available metadata contains only sample ID and group. No age, sex, batch or donor identifiers are available; do not invent covariates.',
    '- Candidate primary common category model: ~0 + group (limma) / ~group (DESeq2), HC reference. Four nonempty groups; rank 4/4; 115 residual degrees of freedom.',
    '- For B, C and exploratory D, use category contrasts MSAP-MAP, SAP-MAP and SAP-MSAP.',
    '- Candidate A and E common-model contrasts are sample-size-weighted category effects (weights 57/87,20/87,10/87; and 20/30,10/30). On a log scale they estimate differences of weighted mean log expression, NOT log2 ratios of pooled arithmetic expression means.',
    '- For conventional collapsed AP vs HC and (MSAP+SAP) vs MAP estimands, fit separately pre-specified group-only binary models (~AP_status on 119; ~severe_status on 87 AP). Report those as primary A/E and weighted four-category contrasts as sensitivity; do not silently equate the two estimands. No models were fitted in this run.',
    '- AP trend candidate: MAP=1, MSAP=2, SAP=3, ~severity_code, rank 2/2, residual df 85. One-step coefficient is an ordinal score association, not evidence of linear biological dose response.',
    '- AP category sensitivity: ~group, MAP reference, rank 3/3, df 84; global two-degree-of-freedom LRT (DESeq2) or moderated F test (limma), plus adjacent B/D contrasts.',
    '- Raw-count provenance is required before DESeq2; noninteger supplied values cannot be rounded into counts. Limma-voom also requires legitimate counts and library sizes.',
    '- Continuous data require documented normalization and scale before any limma transformation or model. TPM/FPKM/CPM, normalized counts and already logged values are not interchangeable.',
    '- Local metadata label H maps to HC under your explicit cohort definition (HC n=32). No group assignment was inferred from expression.',
    '- Algebraic estimability does not establish lack of batch confounding. With absent covariates and blocked independent GEO access, confounding cannot be assessed.'),file.path(dirs['report'],'recommended_design.md'))

# Scale-agnostic descriptive QC: values remain unchanged and no samples are removed.
colors <- c(HC='#0072B2',MAP='#009E73',MSAP='#E69F00',SAP='#D55E00')
save_plot <- function(p,name,w=8,h=6){
    ggsave(file.path(dirs['qc'],paste0(name,'.pdf')),p,width=w,height=h)
    ggsave(file.path(dirs['qc'],paste0(name,'.png')),p,width=w,height=h,dpi=300)
}
theme_common <- theme_bw(base_size=11)
ordered <- order(metadata$group,metadata$id)
stats$sample_display <- factor(stats$sample_id,levels=stats$sample_id[ordered])
save_plot(ggplot(stats,aes(sample_display,sum_supplied_expression,color=group))+geom_point()+
    scale_color_manual(values=colors)+theme_common+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
    labs(x='119 samples, ordered by local group',y='Sum of supplied expression (NOT library size)',
         title='Supplied-matrix sample totals',subtitle='Unknown expression scale; no renormalization or exclusions'),'sample_expression_totals',11,5)
if(all(m>=0)){
    bx <- stats
    for(n in c('min','q25','median','q75','max')) bx[[n]] <- log10(bx[[n]]+1)
    save_plot(ggplot(bx,aes(sample_display,ymin=min,lower=q25,middle=median,upper=q75,ymax=max,fill=group))+
        geom_boxplot(stat='identity',linewidth=.25)+scale_fill_manual(values=colors)+theme_common+
        theme(axis.text.x=element_blank(),axis.ticks.x=element_blank())+
        labs(x='119 samples, ordered by local group',y='log10(supplied value + 1), display only',
             title='Sample expression distributions',subtitle='Whiskers are full range; no data normalization inferred'),'sample_expression_distributions',12,5)
    gx <- do.call(rbind,lapply(levels(metadata$group),function(g){
        q <- quantile(as.vector(m[,metadata$group==g,drop=FALSE]),c(0,.25,.5,.75,1),names=FALSE)
        data.frame(group=g,min=q[1],q25=q[2],median=q[3],q75=q[4],max=q[5])
    }))
    write_table(gx,file.path(dirs['qc'],'group_expression_distribution_summary.csv'))
    for(n in c('min','q25','median','q75','max')) gx[[n]] <- log10(gx[[n]]+1)
    save_plot(ggplot(gx,aes(group,ymin=min,lower=q25,middle=median,upper=q75,ymax=max,fill=group))+
        geom_boxplot(stat='identity')+scale_fill_manual(values=colors)+theme_common+
        labs(y='log10(supplied value + 1), display only',title='Group pooled expression distributions',
             subtitle='Descriptive only: genes and samples are not independent replicates'),'group_expression_distributions')
}
pca <- prcomp(t(m[!constant,,drop=FALSE]),center=TRUE,scale.=FALSE,rank.=10)
prop <- pca$sdev^2/sum(apply(m[!constant,,drop=FALSE],1,var))
scores <- data.frame(sample_id=sample_ids,group=metadata$group,pca$x,check.names=FALSE)
write_table(scores,file.path(dirs['qc'],'pca_scores_supplied_scale.csv'))
write_table(data.frame(PC=seq_along(prop),variance_explained=prop),file.path(dirs['qc'],'pca_variance.csv'))
save_plot(ggplot(scores,aes(PC1,PC2,color=group))+geom_point(size=2)+scale_color_manual(values=colors)+theme_common+
    labs(x=sprintf('PC1 (%.1f%%)',100*prop[1]),y=sprintf('PC2 (%.1f%%)',100*prop[2]),
         title='PCA on unchanged supplied expression',subtitle='All nonconstant genes; centered, not scaled; unknown scale limits interpretation'),'pca_supplied_scale')
save_plot(ggplot(scores,aes(PC3,PC4,color=group))+geom_point(size=2)+scale_color_manual(values=colors)+theme_common+
    labs(x=sprintf('PC3 (%.1f%%)',100*prop[3]),y=sprintf('PC4 (%.1f%%)',100*prop[4]),
         title='PCA PC3-PC4 on supplied expression',subtitle='No VST/rlog/log2 preprocessing assumed'),'pca_supplied_scale_PC3_PC4')
annot <- data.frame(group=metadata$group,row.names=sample_ids)
heat <- function(mat,name,main,dist=NULL){
    args <- list(mat=mat,annotation_col=annot,annotation_row=annot,annotation_colors=list(group=colors),
                 show_colnames=FALSE,show_rownames=FALSE,main=main,border_color=NA)
    if(!is.null(dist)){args$clustering_distance_rows <- dist;args$clustering_distance_cols <- dist}
    pdf(file.path(dirs['qc'],paste0(name,'.pdf')),width=10,height=9)
    tryCatch(do.call(pheatmap,args),finally=dev.off())
    png(file.path(dirs['qc'],paste0(name,'.png')),width=3000,height=2700,res=300)
    tryCatch(do.call(pheatmap,args),finally=dev.off())
}
write_table(data.frame(sample_id=sample_ids,pearson,check.names=FALSE),file.path(dirs['qc'],'sample_correlation_pearson.csv'))
write_table(data.frame(sample_id=sample_ids,spearman,check.names=FALSE),file.path(dirs['qc'],'sample_correlation_spearman.csv'))
heat(pearson,'sample_correlation_pearson','Pearson correlation: supplied scale')
heat(spearman,'sample_correlation_spearman','Spearman correlation: supplied scale')
euclidean <- as.matrix(dist(t(m)))
write_table(data.frame(sample_id=sample_ids,euclidean,check.names=FALSE),file.path(dirs['qc'],'sample_distance_euclidean.csv'))
heat(euclidean,'sample_distance_euclidean','Euclidean sample distance: supplied scale',as.dist(euclidean))
cor_dist <- 1-spearman;cor_dist[cor_dist<0] <- 0;diag(cor_dist) <- 0
write_table(data.frame(sample_id=sample_ids,cor_dist,check.names=FALSE),file.path(dirs['qc'],'sample_distance_spearman.csv'))
heat(cor_dist,'sample_distance_spearman','1 - Spearman sample distance',as.dist(cor_dist))

# Explicit screening rules flag review candidates, never auto-remove samples.
review <- data.frame(sample_id=sample_ids,group=metadata$group,median_within_group_spearman=NA_real_,
                     correlation_review_threshold=NA_real_,pca_within_group_distance=NA_real_,
                     pca_review_threshold=NA_real_,flag_correlation=FALSE,flag_pca=FALSE)
for(g in levels(metadata$group)){
    ii <- which(metadata$group==g)
    cc <- spearman[ii,ii,drop=FALSE];diag(cc) <- NA_real_
    med <- apply(cc,1,median,na.rm=TRUE)
    cut <- median(med)-3*mad(med)
    pc <- pca$x[ii,1:5,drop=FALSE]
    center <- apply(pc,2,median); spread <- apply(pc,2,mad)
    use <- is.finite(spread)&spread>0
    z <- sweep(sweep(pc[,use,drop=FALSE],2,center[use],'-'),2,spread[use],'/')
    dd <- sqrt(rowSums(z^2)); dcut <- median(dd)+3*mad(dd)
    review$median_within_group_spearman[ii] <- med;review$correlation_review_threshold[ii] <- cut
    review$pca_within_group_distance[ii] <- dd;review$pca_review_threshold[ii] <- dcut
    review$flag_correlation[ii] <- mad(med)>0 & med<cut
    review$flag_pca[ii] <- mad(dd)>0 & dd>dcut
}
review$flag_for_review <- review$flag_correlation|review$flag_pca
review$action <- 'retained; no automatic exclusion'
write_table(review,file.path(dirs['qc'],'sample_outlier_review.csv'))
label <- merge(scores,review[,c('sample_id','flag_for_review')],by='sample_id',sort=FALSE)
save_plot(ggplot(label,aes(PC1,PC2,color=group))+geom_point(aes(shape=flag_for_review),size=2)+
    geom_text(data=label[label$flag_for_review,],aes(label=sample_id),size=2.2,vjust=-.6,check_overlap=FALSE)+
    scale_color_manual(values=colors)+theme_common+labs(title='PCA with descriptive review flags',
    subtitle='Flagging is within local groups; all 119 samples retained; no batch metadata',shape='Review flag'),'pca_review_flags',10,7)

# Formal analysis and valid downstream exports are deliberately blocked, not simulated.
blocking <- c('# Formal analysis blocked',
    sprintf('- Supplied matrix: %d genes × %d samples; %.6f%% of values are noninteger.',nrow(m),ncol(m),100*(1-integer_fraction)),
    '- The normalized candidate and uploaded source workbook both contain noninteger expression values. The upload is Excel OOXML disguised with a .csv extension; it is not a raw integer-count matrix.',
    '- Neither file has a processing script, unit declaration or parameter record. Numeric statistics cannot identify TPM, FPKM, CPM, normalized counts, VST, rlog or another continuous scale uniquely.',
    '- DESeq2/voom are invalid on this unverified continuous input. The supplied filtered/normalized matrix was not altered in QC. A separate, clearly labeled source-to-existing-matrix reproduction candidate is generated after QC.',
    '- Limma fitting and any log transformation for formal modelling are withheld until scale and provenance are confirmed.',
    '- Independent GEO retrieval was blocked by HTTP 403 from the cloud network proxy; no authoritative GEO annotations or raw counts were obtained.',
    '- Local metadata H maps to HC from your explicit cohort definition (HC n=32); MAP, MSAP and SAP labels are unchanged. No expression-based grouping was performed.',
    '- Remaining information needed for inference: original uploaded-file units / source and any upstream transform before the numerically recovered filter and quantile-normalization steps, or matching documented raw counts.',
    '- If external GEO lookup is needed, allow www.ncbi.nlm.nih.gov and ftp.ncbi.nlm.nih.gov in environment settings; successful access still needs validation.',
    '- No DEGs, P values, adjusted P values, volcano plots, severity gene statistics or GSEA ranked lists were generated. Counts of significant genes are NOT_RUN, not zero.',
    '- Existing onboarding log2(x+1) plots are exploratory setup checks, not validated preprocessing or analysis evidence.',
    '- No local 205-DEG result was found; no numerical replication or overlap comparison is possible.')
write_lines(blocking,file.path(dirs['report'],'blocking_report.md'))
write_lines(blocking,file.path(dirs['deg'],'BLOCKED.md'))
write_lines(c(blocking,'', 'Trend: MAP=1, MSAP=2, SAP=3 is only a candidate score model; add category global test and adjacent comparisons when input is valid.'),file.path(dirs['severity'],'BLOCKED.md'))
comparisons <- data.frame(comparison=c('A_AP_vs_HC','B_MSAP_vs_MAP','C_SAP_vs_MAP','D_SAP_vs_MSAP','E_MSAP_SAP_vs_MAP'),
                          group_1_n=c(87,20,10,10,30),group_0_n=c(32,57,57,20,57),
                          role=c('primary','primary','primary','exploratory','severity_sensitivity'),
                          tested_genes=NA_integer_,upregulated=NA_integer_,downregulated=NA_integer_,status='NOT_RUN_INPUT_SCALE_UNRESOLVED')
write_table(comparisons,file.path(dirs['deg'],'comparison_status.csv'))
write_lines(c('# Planned DE and stability parameters (not executed)',
    '- Confirm annotation version and legitimate count provenance first. Keep supplied matrix for concordance checks.',
    '- DESeq2 default two-sided Wald tests of LFC=0; retain all modelled rows and clearly flag all-zero/Cooks/independent-filter NA outcomes. Report baseMean, LFC, lfcSE, Wald statistic, P and BH padj.',
    '- padj<0.05 and abs(estimated log2FC)>1 is a descriptive post-test effect-size filter, not a DESeq2 lfcThreshold test. Record shrinkage separately if used.',
    '- If justified, count-based voom/limma sensitivity uses the same cohort, covariates and corresponding contrast estimands; report AveExpr, LFC, moderated t, P, BH padj and confidence intervals.',
    '- AP score trend uses fixed 1/2/3 encoding; add AP category global two-df test, B and D adjacent contrasts. BH correction over testable genes per specified family.',
    '- Keep all samples in the main analysis. Review flags alone do not justify exclusion; any targeted sensitivity needs documented reasons and an explicit sample list.',
    '- For SAP (n=10), report effect sizes and intervals; nonsignificance does not establish absence. Compare sign agreement, correlation and difference of LFCs, testable gene universe and significant overlap only between compatible models.',
    '- GSEA: use finite Wald or moderated t for every testable gene, not DEG-filtered genes. Preserve original IDs; require documented symbol mapping. Resolve duplicate symbols by largest absolute statistic, then original ID lexicographically; discard unresolved/NA symbols; sort statistic descending, symbol ascending for ties; document ties without silently jittering.',
    '- The pre-existing 205-DEG result is unavailable. Do not tune cutoffs or methods to force a target DEG count.'),file.path(dirs['report'],'analysis_plan_not_executed.md'))

# Safe references are preferable to duplicating or incorrectly labelling matrices.
write_table(data.frame(item=c('supplied_expression_matrix','sample_metadata','raw_counts','validated_transformed_matrix','GSEA_ranked_lists','DE_results','gene_symbol_mapping'),
    path=c(paste0(zip_path,'::',entries$Name),file.path(repo,dirs['downstream'],'sample_metadata_verified.csv'),NA,NA,NA,NA,
           file.path(repo,dirs['downstream'],'gene_id_symbol_mapping.csv')),
    status=c('available; scale unknown; descriptive QC only','ID matched; H->HC mapped per explicit cohort definition',
             'NOT_FOUND','BLOCKED_SCALE_UNKNOWN','BLOCKED_NO_VALID_MODEL','BLOCKED_NO_VALID_MODEL','UNVERIFIED_IDENTITY_ONLY')),
    file.path(dirs['downstream'],'input_references.csv'))
write_lines(c('# Module 1 downstream handoff',
    '- No validated count matrix, transformed scoring matrix, DE table or GSEA rank file exists in this run.',
    '- input_references.csv safely references the untouched supplied archive. Do not interpret it as raw counts.',
    '- sample_metadata_verified.csv has verified ID alignment and H->HC mapping from the explicit cohort definition; independent GEO corroboration remains unavailable.',
    '- gene_id_symbol_mapping.csv is a candidate identity table, not an authoritative identifier conversion. Annotation provenance remains missing.',
    '- QC matrices/plots and original expression may support descriptive inspection only. They are not ready for formal NETosis scoring or mechanistic claims.',
    '- No NETs intersection, ML, ROC, WGCNA, single-cell or virtual-knockout work was performed.'),file.path(dirs['downstream'],'README.md'))
audit_report <- c('# GSE194331 input audit',
    sprintf('Normalized candidate has %d genes and %d samples. All sample IDs are GSM-shaped and align one-to-one with local metadata. %s',nrow(m),ncol(m),source_note),
    sprintf('Group sizes: H/HC 32, MAP 57, MSAP 20, SAP 10. Duplicate gene rows: %d; missing gene IDs: %d; identical expression-vector pairs: %d.',sum(duplicated(gene)),sum(is.na(gene)|!nzchar(gene)),nrow(dups)),
    sprintf('Missing/nonfinite/negative values: %d/%d/%d. Zero fraction: %.6f. Integer-like fraction: %.8f. Range: %.8g to %.8g.',sum(is.na(m)),sum(!is.finite(m)),sum(m<0),zero_frac,integer_fraction,min(m),max(m)),
    paste('ID type candidate:',id_type),
    'Both continuous expression matrices contain noninteger values and cannot be used as DESeq2 raw counts. They do not identify their expression units or transform.',
    'No preprocessing script, integer raw counts, GEO metadata, RDS object, formal prior DE result or 205-DEG list was found outside explicitly enumerated runtime/cache exclusions. The uploaded source workbook is audited separately.',
    'No values were rounded, renormalized or filtered. Constant rows were omitted only from PCA calculations and remain in the audit/input.',
    'Only supplied-scale descriptive QC was run. Prior cloud-setup/validate.R applied exploratory log2(x+1) for environment testing; it is not evidence of the source scale.',
    'Sample IDs can be checked for duplication; donor-level biological duplication cannot be assessed because donor identifiers are unavailable.',
    'Age/sex/batch variables are absent. Batch-group confounding is unassessable, not disproven. H->HC follows your explicit group definition.',
    'See blocking_report.md and recommended_design.md for formal analysis prerequisites.')
write_lines(audit_report,file.path(dirs['audit'],'audit_report.md'))
write_table(data.frame(step=c('inventory','sample_id_alignment','expression_scale_confirmation','GEO_group_confirmation','candidate_design_rank',
                             'descriptive_QC','formal_DE_A_to_E','formal_AP_trend','method_stability','validated_downstream_expression','GSEA_rank_export'),
                      status=c('COMPLETE','PASS','BLOCKED','BLOCKED_NETWORK','PASS_ALGEBRA_ONLY','COMPLETE_WITH_SCALE_LIMITATIONS',rep('NOT_RUN_BLOCKED',5))),
                      file.path(dirs['report'],'run_status.csv'))
flagged <- review$sample_id[review$flag_for_review]
summary <- c('# Module 1 summary — audit and descriptive QC complete; formal inference blocked',
    '## Actual inputs and scale',
    paste0('- Expression: ',zip_path,'::',entries$Name,'. No source changes.'),
    '- Group source: 轻中重分组.csv; only id and group are available.',
    sprintf('- Existing filtered matrix: %d genes × %d samples; %.5f%% noninteger. It is unlogged quantile-normalized source expression by exact numeric reproduction; original source units remain unknown.',nrow(m),ncol(m),100*(1-integer_fraction)),
    '- Newly uploaded GSE194331.rar contains one member GSE194331.csv which is actually an OOXML Excel workbook, 39,376 gene rows × 119 GSM columns. Its values are mixed integer/fractional continuous values (45.75% integer-like) and source units/upstream processing remain undocumented. See original_source_provenance.json and original_source_matrix_audit.csv.',
    '- No integer raw counts/filtered raw counts, provenance scripts, RDS objects or existing formal DEGs were located. Prior setup plots are not a prior scientific analysis.',
    '## Samples and design',
    '- All 119 supplied samples retained. Local group counts match HC/H 32, MAP 57, MSAP 20, SAP 10.',
    '- H->HC is recorded from your explicit cohort definition. GEO source metadata could not be fetched independently because the network proxy returned HTTP 403.',
    '- Candidate group-only design is full rank (4/4; df 115); no age/sex/batch covariates were fabricated. See recommended_design.md for pooled versus weighted estimands.',
    '- AP-only candidate score and category matrices have ranks 2/2 and 3/3 respectively. They have NOT been fitted.',
    '## Software, parameters and results',
    paste0('- Executed software: R ',versions[1],', data.table ',versions[2],', ggplot2 ',versions[3],', pheatmap ',versions[4],'. Seed 194331.'),
    '- Executed QC: sample totals/distributions, group distributions, centered PCA on unchanged supplied scale, Pearson/Spearman correlations, Euclidean/correlation distances.',
    sprintf('- PCA on all nonconstant genes, centered and not scaled: PC1 %.2f%% and PC2 %.2f%% of total variance. Native expression scale contains high-magnitude values; interpret this descriptive PCA cautiously.',100*prop[1],100*prop[2]),
    '- Main QC plots are saved in both PDF and PNG. Values were not renormalized; no new biological gene filtering or sample deletion occurred.',
    '- A AP vs HC; B MSAP vs MAP; C SAP vs MAP; D SAP vs MSAP (exploratory); E MSAP+SAP vs MAP: ALL NOT_RUN. Up/down DEG counts are unavailable, not zero.',
    '- DESeq2, limma-voom, limma DE, severity inference, volcano plots, confidence intervals and GSEA ranks were withheld because there is no documented valid input scale.',
    '- No significance or absent biological effect may be inferred. SAP n=10 requires effect sizes/intervals and compatible sensitivity checks after unblocking.',
    '- No local 205-DEG result exists; overlap and numerical reproduction were not possible; no threshold tuning was performed.',
    '## QC review and batch limitations',
    sprintf('- Descriptive review flags: %d samples (%s). These are screening candidates, not established technical outliers; all are retained.',length(flagged),if(length(flagged)) paste(flagged,collapse=', ') else 'none'),
    '- Rules and metrics are in sample_outlier_review.csv. High-variance genes and unknown scale can dominate native-scale PCA; disease variation can also trigger flags.',
    '- Batch confounding, sex/age effects and donor duplication cannot be assessed from available metadata.',
    '## Downstream suitability and blockers',
    '- Inventory, ID alignment, group counts, candidate annotation and descriptive QC can support inspection. No validated NETosis scoring matrix or GSEA rankings are ready.',
    '- Provide source units and original workbook generation/transformation steps, or documented raw counts with matching sample metadata, to choose a defensible inferential model.',
    '- If GEO lookup is needed, environment settings must permit www.ncbi.nlm.nih.gov and ftp.ncbi.nlm.nih.gov; network access has not been established.',
    '- No NETs intersection, ML/ROC, WGCNA, single-cell or virtual knockout was performed.',
    '## Reproducibility and files',
    '- From the repository root run: bash scripts/module1/run_module1.sh',
    '- The driver preserves previous module outputs in timestamped module1_archives on reruns; never overwrites originals.',
    '- scripts/module1 contains the inventory, R analysis, checksum/manifest finalizer and shell driver.',
    '- reports/module1/file_manifest.csv is the complete generated-file list; reports/module1/run.log, parameters.txt, package_versions.csv and sessionInfo.txt record execution.',
    '- results/module1_audit/original_file_sha256.json records input hashes; original_integrity_check.csv verifies them after analysis.',
    '- See results/module1_downstream/input_references.csv and README.md for safe references and blocked exports.',
    '- Newly uploaded source values were imported without normalization or rounding to data/module1_inputs/derived/GSE194331_original_values.tsv.gz; SHA-256 provenance is recorded. These remain unknown-scale continuous values.')
write_lines(summary,file.path(dirs['report'],'module1_summary.md'))
cat(sprintf('MODULE1_AUDIT_QC_COMPLETE: %d genes; %d samples; %d descriptive review flags. FORMAL_ANALYSIS_BLOCKED.\n',nrow(m),ncol(m),length(flagged)))
