require(data.table)

#use all threads
setDTthreads(parallel::detectCores())
#setDTthreads(64)

#set directory as script path
setwd(system("pwd", intern = T))
#setwd("/media/ozcelik2/TO-11/aavc/")

#parse arguments
args <- commandArgs(trailingOnly = TRUE)

#define function to raise errors
printUsageAndExit <- function() {
  cat("Usage: Rscript vcf_parser.R --vcf_dir variant.vcf [--cohort_dir cohort.txt] [--keep_info]\n")
  quit(status = 1)
}

#define arguments
vcf_dir <- NULL
cohort_dir <- NULL
keep_info <- FALSE

#pass arguments
for (i in 1:(length(args) - 1)) {
  if (args[i] %in% c("--vcf_dir", "-v")) {
    vcf_dir <- args[i + 1]
  } else if (args[i] %in% c("--cohort_dir", "-c")) {
    cohort_dir <- args[i + 1]
  } else if (args[i] %in% c("--keep_info", "-k")) {
    keep_original <- TRUE
  } else if (args[i] == "--help" || args[i] == "-h") {
    printUsageAndExit()
  }
}

#check for required arguments
if (is.null(vcf_dir)) {
  cat("Error: Missing required argument --vcf_dir\n")
  printUsageAndExit()
}

#check for invalid arguments
invalid_args <- setdiff(args, c("--vcf_dir", "-v", "--cohort_dir", "-c", "--keep_info", "-k", "--help", "-h", vcf_dir, cohort_dir))
if (length(invalid_args) > 0) {
  cat("Error: Invalid argument(s):", paste(invalid_args, collapse = " "), "\n")
  printUsageAndExit()
}

#read config file
vcf_config <- readLines("vcf_config.txt")

#get field names
grpmax <- sub(".*grpmax=(.*)", "\\1", vcf_config[grep("grpmax=", vcf_config)])
if (grpmax == "") {stop("grpmax is not defined. Please define the corresponding field name in the vcf_config.txt.")}
nhomalt <- sub(".*nhomalt=(.*)", "\\1", vcf_config[grep("nhomalt=", vcf_config)])
if (nhomalt == "") {stop("nhomalt is not defined. Please define the corresponding field name in the vcf_config.txt.")}
spliceai <- sub(".*spliceai=(.*)", "\\1", vcf_config[grep("spliceai=", vcf_config)])
if (spliceai == "") {stop("spliceai is not defined. Please define the corresponding field name in the vcf_config.txt.")}
phylop <- sub(".*phylop=(.*)", "\\1", vcf_config[grep("phylop=", vcf_config)])
if (phylop == "") {stop("phylop is not defined. Please define the corresponding field name in the vcf_config.txt.")}
predictor <- sub(".*predictor=(.*)", "\\1", vcf_config[grep("predictor=", vcf_config)])
if (predictor == "") {stop("predictor is not defined. Please define the corresponding field name in the vcf_config.txt.")}

#collect fields to be extracted
field_names <- list(grpmax = grpmax, nhomalt = nhomalt, predictor = predictor, spliceai = spliceai, phylop = phylop)

#vcf <- fread("/media/ozcelik2/TO-12/gnomAD4_4.1/protein_coding/csq/Obesity/try.txt", header=F, data.table=T, sep="\t")

#read vcf and get unique columns
vcf <- fread(vcf_dir, header=F, data.table=T, sep="\t")
vcf <- unique(vcf, by = c("V1", "V2", "V4", "V5"))

#fix chromosome names
if (grepl("chr",vcf[5,1])) {vcf[, V1 := substr(V1,4,nchar(V1))]}

#extract gene symbol
vcf[, gene := sub(".*;CSQ=.*\\|(LOW|MODERATE|MODIFIER|HIGH)\\|(.+?)\\|ENSG.*","\\2", V8)]
vcf <- vcf[nchar(gene) < 20,]

#extract gene and tcpt identifiers
vcf[, ensg_id:= sub(".*;CSQ=.*(ENSG[0-9]+).*","\\1", V8)]
vcf[, enst_id := sub(".*;CSQ=.*(ENST[0-9]+).*","\\1", V8)]

#extract nucleotide change
vcf[, nt_change := sub(".*;CSQ=.*:c.(.*?)(?::p|\\|).*","\\1", V8)]
vcf[nchar(nt_change) >= 200, nt_change := sub(".*;CSQ=.*c\\.(\\d+[+-]\\d+[ACGT]>[ACGT]).*","\\1", V8)]
vcf[nchar(nt_change) >= 200, nt_change := "."]

#extract aminoacid change
vcf[, aa_change := sub(".*;CSQ=.*:p.(.*?)\\|.*","\\1", V8)]
vcf[nchar(aa_change) >= 200, aa_change := "."]

#extract variant consequence
vcf[, effect := sub(".*;CSQ=.*?\\|(.*?)\\|.*","\\1", V8)]
vcf[effect == "character(0)", effect := "."]

#define extractor function and extract annotations
extract_annotation <- function(column, field) {
  pattern <- paste0(".*;",field,"=(.*?);.*")
  vcf[, (column) := as.numeric(sub(pattern, "\\1", V8))]
  vcf[nchar(get(column)) >= 200, (column) := NA]
}
for (column_name in names(field_names)) {extract_annotation(column_name, field_names[[column_name]])}

#remove redundant columns and rename them
colnames(vcf)[1:8] <- c("chr","pos","id","ref","alt","format","filter","info")
if (!keep_info) {svcf <- cbind(vcf[, c(1, 2, 4, 5)], vcf[, 9:ncol(vcf)])
} else {svcf <- cbind(vcf[, c(1, 2, 4, 5)], vcf[, 8:ncol(vcf)])}

#run in case-control mode
if (!is.null(cohort_dir)) {
  
  #load cohort files
  cohort <- fread(cohort_dir, header=F, data.table = T, sep= "\t", na.strings = "")
  
  #assign indiv ids
  colnames(vcf)[10:ncol(vcf)] <- cohort$V2
  
  #extract genotypes
  vcf[, (10:ncol(vcf)) := lapply(.SD, function(x) substring(x, 1, 3)), .SDcols = 10:ncol(vcf)]
  
  #get indices for case and control
  case_ind <- as.vector(which(cohort$V6 == 1) + as.integer(which(colnames(vcf) == colnames(vcf)[10])) - 1)
  cont_ind <- as.vector(which(cohort$V6 == 0) + as.integer(which(colnames(vcf) == colnames(vcf)[10])) - 1)
  uq_case_fam <- nlevels(factor(cohort[cohort$V6 == 1]$V1))
  uq_cont_fam <- nlevels(factor(cohort[cohort$V6 == 0]$V1))
  
  #subset data by case/control
  case <- vcf[, ..case_ind]
  cont <- vcf[, ..cont_ind]
  
  #find num of hom/het indivs
  case_hom <- rowSums(sapply(case, function(x) substr(x, 1, 3) %in% c("1/1","1|1")))
  case_het <- rowSums(sapply(case, function(x) substr(x, 1, 3) %in% c("0/1","0|1")))
  case_ref <- rowSums(sapply(case, function(x) substr(x, 1, 3) %in% c("0/0","0|0")))
  cont_hom <- rowSums(sapply(cont, function(x) substr(x, 1, 3) %in% c("1/1","1|1")))
  cont_het <- rowSums(sapply(cont, function(x) substr(x, 1, 3) %in% c("0/1","0|1")))
  cont_ref <- rowSums(sapply(cont, function(x) substr(x, 1, 3) %in% c("0/0","0|0")))
  
  #write results onto table
  vcf[, Case := paste(case_hom,case_het,sep="|")]
  vcf[, Control := paste(cont_hom,cont_het,sep="|")]
  
  #get carriers by each group
  group_carriers <- function(row, fam_wise = F) {
    
    cols <- grep("1\\/1|^1\\|1|^0\\/1|^0\\|1", row, value = T)
    #cols <- grep("1\\/1|^0\\/1", row, value = T)
    indivs <- as.vector(names(cols))
    
    case_indiv <- c()
    cont_indiv <- c()
    case_famil <- c()
    cont_famil <- c()
    
    for (i in indivs) {
      
      if (cohort[cohort$V2 == i, "V6"] == 1) {
        
        fam <- cohort$V1[cohort$V2 == i]
        case_famil <- c(case_famil, fam)
        case_indiv <- c(case_indiv, i)
        
        
      }else if (cohort[cohort$V2 == i, "V6"] == 0) {
        
        fam <- cohort$V1[cohort$V2 == i]
        cont_famil <- c(cont_famil, fam)
        cont_indiv <- c(cont_indiv, i)
        
      }
      
    }
    
    case_indiv <- na.omit(case_indiv)
    case_indiv <- unique(case_indiv)
    case_indiv <- paste(case_indiv, collapse=",")
    if (case_indiv == ""){case_indiv = "." }
    
    cont_indiv <- na.omit(cont_indiv)
    cont_indiv <- unique(cont_indiv)
    cont_indiv <- paste(cont_indiv, collapse=",")
    if (cont_indiv == ""){cont_indiv = "." }
    
    case_famil <- na.omit(case_famil)
    case_famil <- unique(case_famil)
    case_famil <- paste(case_famil, collapse=",")
    if (case_famil == ""){case_famil = "." }
    
    cont_famil <- na.omit(cont_famil)
    cont_famil <- unique(cont_famil)
    cont_famil <- paste(cont_famil, collapse=",")
    if (cont_famil == ""){cont_famil = "." }
    
    return(list(case_indiv,cont_indiv,case_famil,cont_famil))
    
  }
  list_homozygotes <- function(row) {
    
    cols <- grep("1\\/1|^1\\|1", row, value = T)
    indivs <- as.vector(names(cols))
    
    hom_indivs <- c()
    
    for (i in indivs) {
      
      hom_indivs <- c(hom_indivs, i)
      
    }
    
    hom_indivs <- na.omit(hom_indivs)
    hom_indivs <- unique(hom_indivs)
    hom_indivs <- paste(hom_indivs, collapse=",")
    if (hom_indivs == "") {hom_indivs = "."}
    
  }
  
  vcf[, AF := (case_hom+case_het+cont_hom+cont_het)/(case_hom+case_het+cont_hom+cont_het+case_ref+cont_ref)]
  
  vcf[AF > 0.1,  c("Case_Indiv","Cont_Indiv","Case_Famil","Cont_Famil","Hom_Indiv") := "..."]
  vcf[AF <= 0.1, c("Case_Indiv","Cont_Indiv","Case_Famil","Cont_Famil") := as.list(group_carriers(.SD)), by = list(1:nrow(vcf[AF <= 0.1,])), .SDcols = 10:(nrow(cohort) + 9)]
  vcf[AF <= 0.1, Hom_Indiv := list_homozygotes(.SD), by = list(1:nrow(vcf[AF <= 0.1,])), .SDcols = 10:(nrow(cohort) + 9)]

  #calculate ratio for unrelated cohort members
  vcf[, CC_Fam := paste(
    ifelse(
      Case_Famil == ".",0,
      ifelse(Case_Famil == "...","?",
             sapply(
               strsplit(as.character(Case_Famil), ","), length))),
    ifelse(
      Cont_Famil == "...","?",
      ifelse(Cont_Famil == ".",0,
             sapply(
               strsplit(as.character(Cont_Famil), ","), length))),sep="_")]
  
  vcf[, UR_Case := sapply(CC_Fam, function(x) as.numeric(strsplit(x, "_")[[1]][1]))]
  vcf[, WT_Case := uq_case_fam - UR_Case]
  vcf[, UR_Cont := sapply(CC_Fam, function(x) as.numeric(strsplit(x, "_")[[1]][2]))]
  vcf[, WT_Cont := uq_cont_fam - UR_Cont]
  
  #calculate odds ratio & statistics
  vcf[, OR := ifelse(UR_Case  == 0 | WT_Case == 0 | UR_Cont == 0 | WT_Cont == 0, round((UR_Case+0.5)*(WT_Cont+0.5)/((WT_Case+0.5)*(UR_Cont+0.5)),2), round((UR_Case)*(WT_Cont)/((WT_Case)*(UR_Cont)),2))]
  vcf[, SE := ifelse(UR_Case  == 0 | WT_Case == 0 | UR_Cont == 0 | WT_Cont == 0, sqrt(1/(UR_Case+0.5) + 1/(WT_Case+0.5) + 1/(UR_Cont+0.5) + 1/(WT_Cont)), sqrt(1/(UR_Case) + 1/(WT_Case) + 1/(UR_Cont) + 1/(WT_Cont)))]
  vcf[, LL := exp(log(OR) - 1.96 * SE)]
  vcf[, UL := exp(log(OR) + 1.96 * SE)]
  
  #calculate PS4 status
  vcf[,PS4:= ifelse(OR>5 & !(1 >= LL & 1 <= UL), "PS4_S", ifelse(Control == "0|0" & Case != "0|0" & CC_Fam != "1_0", "PS4_M","N/A"))]
  vcf[AF > 0.1, PS4 := "N/A"]
  
  #remove stats columns
  vcf[, c("UR_Case","WT_Case", "UR_Cont", "WT_Cont", "SE","LL","UL") := NULL]
  
}

#get output name and write it
o_name <- paste0(strsplit(vcf_dir, "\\.")[[1]][1],".svcf")
fwrite(svcf, o_name, quote=F, sep="\t", row.names = F, col.names = T)
