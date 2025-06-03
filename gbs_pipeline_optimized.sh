#!/usr/bin/env bash

# Optimized generalized GBS population genomics pipeline
# Provides modular functions, dependency checks and cleaner logging.

set -euo pipefail
IFS=$'\n\t'

#==================== CONFIGURATION ==========================================
RAW_VCF="Norwegian_Only.vcf.gz"
SAMPLE_LIST="Norwegian_SampleIDs.txt"
POP_FILE="Silene_Norwegian_Pop366.txt"
METADATA="sample_metadata.txt"
ENV_FILE="sample_environment.txt"
ADMIX_PLOT_SCRIPT="scripts/plotADMIXTURE_updated3.R"

RESULTS_DIR="results"
QC_DIR="${RESULTS_DIR}/qc"
FILTERED_DIR="${RESULTS_DIR}/filtered"
PLINK_DIR="${RESULTS_DIR}/plink"
PCA_DIR="${RESULTS_DIR}/pca"
ADMIX_DIR="${RESULTS_DIR}/Admixture_results"
ADMIX_LOG_DIR="${ADMIX_DIR}/logs"
ADMIX_CV_DIR="${ADMIX_DIR}/cv"
ADMIX_BARPLOT_DIR="${ADMIX_DIR}/barplots"
DIVERSITY_DIR="${RESULTS_DIR}/diversity"
FST_DIR="${RESULTS_DIR}/fst"
IBD_DIR="${RESULTS_DIR}/ibd"

MAF=0.05
MAX_MISSING=0.9
GENO=0.1
MIND=0.15
N_PCS=10
K_MIN=2
K_MAX=10
THREADS=4

VERSION_LOG="${QC_DIR}/software_versions.txt"
PIPELINE_LOG="${QC_DIR}/pipeline_run.log"

#==================== HELPER FUNCTIONS ======================================
log(){ echo -e "$1" | tee -a "${PIPELINE_LOG}"; }

check_dep(){
  for c in "$@"; do
    command -v "$c" >/dev/null 2>&1 || { echo "Required command $c not found" >&2; exit 1; }
  done
}

run_if_missing(){
  local target="$1"; shift
  if [[ ! -e "$target" ]]; then
    "$@"
  else
    log "✔ Skipping $(basename "$target") (already exists)"
  fi
}

#==================== SETUP =================================================
mkdir -p "${QC_DIR}" "${FILTERED_DIR}" "${PLINK_DIR}" "${PCA_DIR}" \
  "${ADMIX_LOG_DIR}" "${ADMIX_CV_DIR}" "${ADMIX_BARPLOT_DIR}" \
  "${DIVERSITY_DIR}" "${FST_DIR}" "${IBD_DIR}"

START_TIME=$(date '+%F %T')
log "Pipeline started at ${START_TIME}"

check_dep vcftools bcftools plink admixture Rscript

#==================== STEP 0: SOFTWARE VERSIONS ==============================
run_if_missing "${VERSION_LOG}" bash -c '
  {
    echo "=== Software Versions ===";
    vcftools --version 2>&1 | head -n1;
    bcftools --version | head -n1;
    plink --version | head -n1;
    admixture --version 2>&1 | head -n1;
    Rscript --version;
  } > "${VERSION_LOG}"
'

#==================== STEP 1: INITIAL QC ====================================
run_if_missing "${QC_DIR}/total_snps.txt" bash -c "zgrep -v '^#' '${RAW_VCF}' | wc -l > '${QC_DIR}/total_snps.txt'"
run_if_missing "${QC_DIR}/total_samples.txt" bash -c "zgrep '^#CHROM' '${RAW_VCF}' | awk '{print NF-9}' > '${QC_DIR}/total_samples.txt'"
run_if_missing "${QC_DIR}/sample_list_from_vcf.txt" bash -c "zgrep '^#CHROM' '${RAW_VCF}' | cut -f10- | tr '\t' '\n' > '${QC_DIR}/sample_list_from_vcf.txt'"

#==================== STEP 2: FILTER RAW VCF =================================
FILTERED_VCF_PREFIX="${FILTERED_DIR}/filtered_MAF${MAF}_miss${MAX_MISSING}"
FILTERED_VCF_GZ="${FILTERED_VCF_PREFIX}.vcf.gz"
run_if_missing "${FILTERED_VCF_GZ}" bash -c '
  vcftools --gzvcf "${RAW_VCF}" --maf "${MAF}" --max-missing "${MAX_MISSING}" \
    --recode --recode-INFO-all --out "${FILTERED_VCF_PREFIX}" && \
  mv "${FILTERED_VCF_PREFIX}.recode.vcf" "${FILTERED_VCF_PREFIX}.vcf" && \
  bgzip -f "${FILTERED_VCF_PREFIX}.vcf" && \
  tabix -p vcf "${FILTERED_VCF_GZ}"
'

#==================== STEP 3: CONVERT TO PLINK ===============================
PLINK_RAW_PREFIX="${PLINK_DIR}/raw"
run_if_missing "${PLINK_RAW_PREFIX}.bed" plink --vcf "${FILTERED_VCF_GZ}" --double-id --make-bed --out "${PLINK_RAW_PREFIX}" --threads "${THREADS}"

#==================== STEP 4: PLINK GENO/MIND FILTER ========================
PLINK_FILTERED_PREFIX="${PLINK_DIR}/filtered"
run_if_missing "${PLINK_FILTERED_PREFIX}.bed" plink --bfile "${PLINK_RAW_PREFIX}" --geno "${GENO}" --mind "${MIND}" --make-bed --out "${PLINK_FILTERED_PREFIX}" --threads "${THREADS}"

#==================== STEP 5: PCA ===========================================
PCA_PREFIX="${PCA_DIR}/pca"
run_if_missing "${PCA_PREFIX}.eigenvec" plink --bfile "${PLINK_FILTERED_PREFIX}" --pca "${N_PCS}" --out "${PCA_PREFIX}" --threads "${THREADS}"

cat <<'EOS' > "${PCA_DIR}/plot_pca.R"
# Simple PCA plot
library(data.table); library(ggplot2); library(dplyr)
vals <- fread('pca.eigenval', header=FALSE)
vecs <- fread('pca.eigenvec', header=FALSE)
pop  <- fread('../../Silene_Norwegian_Pop366.txt', header=FALSE)
setnames(pop, c('IID','POP'))
vecs <- vecs[,1:4]; setnames(vecs, c('FID','IID','PC1','PC2'))
merged <- inner_join(vecs, pop, by='IID')
var_exp <- vals$V1/sum(vals$V1)*100
p <- ggplot(merged, aes(PC1, PC2, color=POP)) + geom_point() +
  xlab(paste0('PC1 (', round(var_exp[1],1),'%)')) +
  ylab(paste0('PC2 (', round(var_exp[2],1),'%)'))+
  theme_minimal()
ggsave('PCA_PC1_vs_PC2.png', p, width=6, height=5)
EOS
run_if_missing "${PCA_DIR}/PCA_PC1_vs_PC2.png" bash -c "cd '${PCA_DIR}' && Rscript plot_pca.R"

#==================== STEP 6: ADMIXTURE ======================================
for K in $(seq "${K_MIN}" "${K_MAX}"); do
  LOG="${ADMIX_LOG_DIR}/admix_K${K}.log"
  if [[ ! -e "$LOG" ]]; then
    admixture -j"${THREADS}" --cv "${PLINK_FILTERED_PREFIX}.bed" "$K" 2>&1 | tee "$LOG"
  else
    log "✔ Skipping ADMIXTURE K=$K (log exists)"
  fi
 done

run_if_missing "${ADMIX_CV_DIR}/cv_errors.txt" bash -c '
  echo -e "K\tCV_Error" > "${ADMIX_CV_DIR}/cv_errors.txt"
  for K in $(seq "${K_MIN}" "${K_MAX}"); do
    CV=$(grep "CV error" "${ADMIX_LOG_DIR}/admix_K${K}.log" | awk -F": " "{print $2}")
    echo -e "${K}\t${CV:-NA}" >> "${ADMIX_CV_DIR}/cv_errors.txt"
  done
'

cat <<'EOS' > "${ADMIX_CV_DIR}/plot_cv.R"
library(data.table); library(ggplot2)
cv <- fread('cv_errors.txt')
cv[,Delta:=c(NA, abs(diff(CV_Error)))]
q <- ggplot(cv, aes(K, CV_Error)) + geom_line() + geom_point() + theme_minimal()
qgg <- ggplot(cv[!is.na(Delta)], aes(K, Delta)) + geom_line(color='steelblue')+geom_point(color='steelblue')+ theme_minimal()
ggsave('ADMIXTURE_CV_Error.png', q, width=6, height=4)
ggsave('ADMIXTURE_Delta_CV.png', qgg, width=6, height=4)
EOS
run_if_missing "${ADMIX_CV_DIR}/ADMIXTURE_CV_Error.png" bash -c "cd '${ADMIX_CV_DIR}' && Rscript plot_cv.R"

if [[ -f "${ADMIX_PLOT_SCRIPT}" ]]; then
  POP_LIST=$(cut -f2 "${POP_FILE}" | sort -u | paste -sd ',' -)
  if [[ -z $(ls "${ADMIX_BARPLOT_DIR}"/*.png 2>/dev/null) ]]; then
    (cd "${ADMIX_BARPLOT_DIR}" && Rscript "../../${ADMIX_PLOT_SCRIPT}" -p "../../${PLINK_FILTERED_PREFIX}" -i "../../${POP_FILE}" -k "${K_MAX}" -l "${POP_LIST}")
  fi
fi

#==================== STEP 7: DIVERSITY ======================================
run_if_missing "${DIVERSITY_DIR}/individual_heterozygosity.het" vcftools --gzvcf "${FILTERED_VCF_GZ}" --het --out "${DIVERSITY_DIR}/individual_heterozygosity"

POP_NAMES=( $(cut -f2 "${POP_FILE}" | sort -u) )
for POP in "${POP_NAMES[@]}"; do
  run_if_missing "${DIVERSITY_DIR}/${POP}_pi.sites.pi" bash -c '
    awk -F"\t" -v p="'"${POP}"'" "$2==p{print $1}" "${POP_FILE}" > "${DIVERSITY_DIR}/${POP}_samples.txt" && \
    bcftools view -Ov -S "${DIVERSITY_DIR}/${POP}_samples.txt" "${FILTERED_VCF_GZ}" -o "${DIVERSITY_DIR}/${POP}.vcf" && \
    bgzip -f "${DIVERSITY_DIR}/${POP}.vcf" && tabix -p vcf "${DIVERSITY_DIR}/${POP}.vcf.gz" && \
    vcftools --gzvcf "${DIVERSITY_DIR}/${POP}.vcf.gz" --site-pi --out "${DIVERSITY_DIR}/${POP}_pi"
  '
  rm -f "${DIVERSITY_DIR}/${POP}_samples.txt" "${DIVERSITY_DIR}/${POP}.vcf.gz" "${DIVERSITY_DIR}/${POP}.vcf.gz.tbi"
done

run_if_missing "${DIVERSITY_DIR}/pop_heterozygosity.het" bash -c '
  awk "BEGIN{OFS=\t}{print $1,$1,$2}" "${POP_FILE}" > "${DIVERSITY_DIR}/within_pop.txt" && \
  plink --bfile "${PLINK_FILTERED_PREFIX}" --within "${DIVERSITY_DIR}/within_pop.txt" --het --out "${DIVERSITY_DIR}/pop_heterozygosity" --threads "${THREADS}"
'

cat <<'EOS' > "${DIVERSITY_DIR}/summarize_fis.R"
library(data.table)
het <- fread('pop_heterozygosity.het')
popinfo <- fread('../../Silene_Norwegian_Pop366.txt', header=FALSE)
setnames(popinfo, c('IID','POP'))
het_pop <- merge(het, popinfo, by='IID')
summary_pop <- het_pop[,.(Mean_F=mean(F,na.rm=TRUE),SD_F=sd(F,na.rm=TRUE)), by=POP]
fwrite(summary_pop,'FIS_summary.txt',sep='\t')
EOS
run_if_missing "${DIVERSITY_DIR}/FIS_summary.txt" bash -c "cd '${DIVERSITY_DIR}' && Rscript summarize_fis.R"

#==================== STEP 8: PAIRWISE FST ==================================
if [[ -z $(ls ${FST_DIR}/*.weir.fst 2>/dev/null) ]]; then
  for ((i=0;i<${#POP_NAMES[@]};i++)); do
    for ((j=i+1;j<${#POP_NAMES[@]};j++)); do
      P1="${POP_NAMES[i]}"; P2="${POP_NAMES[j]}";
      OUT="${FST_DIR}/FST_${P1}_vs_${P2}"
      awk -F"\t" -v p="$P1" '$2==p{print $1}' "${POP_FILE}" > "${FST_DIR}/${P1}.samples"
      awk -F"\t" -v p="$P2" '$2==p{print $1}' "${POP_FILE}" > "${FST_DIR}/${P2}.samples"
      vcftools --gzvcf "${FILTERED_VCF_GZ}" --weir-fst-pop "${FST_DIR}/${P1}.samples" --weir-fst-pop "${FST_DIR}/${P2}.samples" --out "$OUT"
      rm -f "${FST_DIR}/${P1}.samples" "${FST_DIR}/${P2}.samples"
    done
  done
fi

#==================== STEP 9: IBD ===========================================
run_if_missing "${IBD_DIR}/GeneticDist_IBS.dist" plink --bfile "${PLINK_FILTERED_PREFIX}" --distance square ibs --out "${IBD_DIR}/GeneticDist_IBS" --threads "${THREADS}"
if [[ -f "${METADATA}" ]]; then
  cat <<'EOS' > "${IBD_DIR}/compute_ibd.R"
library(data.table); library(geosphere); library(vegan)
meta <- fread('../sample_metadata.txt')
fam <- fread('../results/plink/filtered.fam', header=FALSE)
setnames(fam, c('FID','IID','PAT','MAT','SEX','PHEN'))
meta <- meta[IID %in% fam$V2]
meta <- meta[match(fam$V2, meta$IID)]
coords <- as.matrix(meta[,.(lon,lat)])
geo <- distm(coords, fun=distHaversine)/1000
gen <- as.matrix(fread('GeneticDist_IBS.dist', header=FALSE))
gen <- 1 - gen
mantel_res <- mantel(as.dist(gen), as.dist(geo), method='pearson', permutations=9999)
sink('Mantel_IBD_results.txt'); print(mantel_res); sink()
EOS
  run_if_missing "${IBD_DIR}/Mantel_IBD_results.txt" bash -c "cd '${IBD_DIR}' && Rscript compute_ibd.R"
fi

END_TIME=$(date '+%F %T')
log "Pipeline finished at ${END_TIME}"
