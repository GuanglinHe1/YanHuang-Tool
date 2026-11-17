# Check independently calculated results and deliverable consistency.
samples <- read.delim('input/samples.tsv')
for (level in 1:3) {
 tab <- read.delim(sprintf('output/result/2-frequency/freq_matrix_HH%02d.tsv',level),check.names=FALSE)
 expected <- prop.table(table(samples$population,substr(samples$haplogroup_full,1,level)),1)
 actual <- as.matrix(tab[,-1,drop=FALSE]);rownames(actual)<-tab$population
 stopifnot(max(abs(actual-expected[rownames(actual),colnames(actual)])) < 1e-12)
 for(metric in c('TVD','JSDist','NeiD','DA','BrayCurtis','Euclidean','Manhattan')) {
  z <- read.delim(sprintf('output/result/3-distance/pairwise_%s_HH%02d.tsv',metric,level),check.names=FALSE)
  D <- as.matrix(z[,-1,drop=FALSE])
  stopifnot(max(abs(D-t(D)))<1e-12,all(diag(D)==0))
  p <- actual[1,];q <- actual[2,];m <- (p+q)/2
  value <- switch(metric,TVD=sum(abs(p-q))/2,DA=1-sum(sqrt(p*q)),NeiD=-log(sum(p*q)/sqrt(sum(p*p)*sum(q*q))),BrayCurtis=sum(abs(p-q))/sum(p+q),Euclidean=sqrt(sum((p-q)^2)),Manhattan=sum(abs(p-q)),JSDist={p<-p+1e-12;p<-p/sum(p);q<-q+1e-12;q<-q/sum(q);m<-(p+q)/2;sqrt(sum(p*log(p/m)+q*log(q/m))/2)})
  stopifnot(abs(D[1,2]-value)<1e-12)
 }
 stopifnot(file.exists(sprintf('output/result/5-bootstrap/HH%02d/bootstrap/haplogroup_bootstrap_freq_long.tsv',level)))
}
pdfs <- list.files('output/figure',pattern='[.]pdf$',recursive=TRUE,full.names=TRUE)
stopifnot(length(pdfs)>0,all(file.exists(sub('[.]pdf$','.tsv',pdfs))))
for(f in list.files('R',pattern='[.]R$',full.names=TRUE)) parse(f)
cat('PASS: frequencies, seven distance formulas, symmetry, bootstrap levels, PDF sidecars, R syntax\n')

config <- yaml::read_yaml('conf/1-preprocess.yaml')
run <- function(config,module) {
 values <- unlist(config[c('paths','runtime','analysis')],recursive=TRUE,use.names=TRUE)
 names(values) <- sub('^[^.]+[.]','',names(values))
 args <- as.vector(rbind(paste0('--',names(values)),shQuote(as.character(values))))
 suppressWarnings(system2('Rscript',c(module,args),stdout=FALSE,stderr=FALSE))
}
config$paths$output <- 'temp/1-preprocess/edge-result'
config$paths$figure <- 'temp/1-preprocess/edge-figure'
config$runtime$overwrite <- TRUE
z <- read.delim('input/samples.tsv')
for(kind in c('empty','duplicate','missing')) {
 bad <- switch(kind,empty=z[FALSE,],duplicate=rbind(z,z[1,]),missing=z[,-1])
 write.table(bad,'temp/1-preprocess/invalid.tsv',sep='\t',quote=FALSE,row.names=FALSE)
 config$paths$input <- 'temp/1-preprocess/invalid.tsv'
 stopifnot(run(config,'R/1-1-preprocess.R')!=0)
}
config$paths$input <- 'input/samples.tsv'
stopifnot(run(config,'R/1-1-preprocess.R')==0)
config$runtime$overwrite <- FALSE
stopifnot(run(config,'R/1-1-preprocess.R')!=0)
config <- yaml::read_yaml('conf/3-distance.yaml')
config$paths$output <- 'temp/1-preprocess/edge-distance'
config$paths$figure <- 'temp/1-preprocess/edge-distance-figure'
config$runtime$overwrite <- TRUE
config$paths$metadata <- '';config$paths$ystr <- ''
stopifnot(run(config,'R/3-1-distance.R')==0)
config$paths$metadata <- 'input/nonexistent.tsv'
stopifnot(run(config,'R/3-1-distance.R')!=0)
cat('PASS: empty input, duplicate IDs, missing columns, overwrite protection, omitted optional inputs, missing configured metadata\n')
