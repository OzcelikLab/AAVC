require(data.table)
require(R.utils)
#1) Set directories ----

#Working directory
setwd("/path/to/your/directory/")

#pre_annotation_db
database_path <- "/path/to/the/pre_annotation_db/"

#Ensembl VEP directories
vep_directory <- "/path/to/the/ensembl/vep" # "usually ends with /ensembl/ensembl-vep/vep"
vep_cache_directory <- "/path/to/the/ensembl/vep" # "usually ends with /ensembl/ensembl-vep"
reference_genome_directory <- "/path/to/the/reference/genome" #usually a fasta file

#2) Select your input file as Chr-Pos-Ref-Alt format (See example_input.txt) ----
#   Go to section 3 and read specify your .vcf file if the file is already in the VCF format

#Select your input file
input_file <- fread("input.txt", header = F) 

#Convert to the VCF format
input_file[, c("CHROM", "POS", "REF", "ALT") := tstrsplit(V1, "-", fixed = TRUE)]
vcf <- input_file[, .(
  `#CHROM` = CHROM,
  POS = as.integer(POS),
  ID = ".",
  REF,
  ALT,
  QUAL = ".",
  FILTER = ".",
  INFO = "."
)]

fwrite(vcf, "example_vcf.vcf", sep="\t")
rm(input_file, vcf)

#3)Name of the vcf file ----
vcf_file <- "example_vcf.vcf"

#4)Run VEP ----
cmd <- paste0(vep_directory," --assembly GRCh38 --cache  --force_overwrite --input_file ",
              vcf_file, " --output_file aavc_temp_input --species homo_sapiens --dir_cache ",
              vep_cache_directory, " --canonical --hgvs --fasta ",
              reference_genome_directory, " --pick --pick_order mane_select,canonical,tsl,biotype,ccds,rank,length --vcf"
)
system(cmd)
rm(cmd)

#4)Annotate with additional scores

vcf <- fread("aavc_temp_input")
colnames(vcf)[[1]] <- "CHROM"

#Annotate w/ revel
file_path <- paste0(database_path,"/revel.tsv") #please download genome-wide revel predictions from https://www.dbnsfp.org/
revel <- fread(file_path, header = T)
vcf[revel, on = .(CHROM, POS, REF, ALT), revel := revel]
rm(revel)
gc()

#Annotate w/gnomad
file_path <- paste0(database_path,"/gnomad.tsv") #please download full gnomAD database from https://gnomad.broadinstitute.org/data#v4
gnomad <- fread(file_path, header=T)
vcf[gnomad, on = .(CHROM, POS, REF, ALT), `:=` (maf=maf, nhomalt=nhomalt, grpmax=grpmax)]
rm(gnomad)
gc()

#Annotate w/phyloP
file_path <- paste0(database_path, "/phylop.tsv") #please download full phyloP calculations from https://www.dbnsfp.org/ or https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/
phylop <- fread(file_path, header=T)
vcf[phylop, on = .(CHROM, POS, REF, ALT), phylop := phylop]
rm(phylop)
gc()

#Annotate w/spliceai
file_path <- paste0(database_path, "/spliceai.tsv") #please download full SpliceAI calculations from https://basespace.illumina.com/s/otSPW8hnhaZR
spliceai <- fread(file_path, header = T)
vcf[spliceai, on = .(CHROM, POS, REF, ALT), spliceai := spliceai]
rm(spliceai)
gc()

#5)Write the output
vcf[,INFO:=paste0("PR;",INFO,";revel_max=",revel,";maf=",maf,";nhomalt=",nhomalt,";fafmax_faf95_max=",grpmax,";phylop=",phylop,";spliceai_ds_max=",spliceai,";END")]
gc()
vcf <- vcf[, c(1:8)]
fwrite(vcf, "example_AAVC_input.txt", col.names = F, sep="\t")


