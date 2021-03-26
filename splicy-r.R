library(data.table)

x <- fread("annotated-counts.tsv", header=T)
x <- x[classification != 'unannotated']
x <- x[classification != 'intergenic']
x.cols <- colnames(x)
x.cols <- x.cols[x.cols != 'extra']
x <- unique(x[, ..x.cols])

all.samples <- x[, unique(sample)]
all.classes <- x[, unique(classification)]

gene.counts <- x[, .(total = sum(coverage)), by=.(gene, classification)]

gene.totals <- gene.counts[, .(total = sum(total)), by=.(gene)]

all.genes <- gene.totals[total >= 100, gene]

x <- x[gene %in% all.genes]

gene.counts <- CJ(gene=all.genes, classification=all.classes)
gene.counts <- merge(gene.counts, x[, .(total = sum(coverage)), by=.(gene, classification)], by=c('gene', 'classification'), all.x=TRUE)
gene.counts[is.na(total), total := 0]
gene.counts[, q := total/max(1, sum(total)), by=.(gene)]

sample.gene.counts <- CJ(sample=all.samples, gene=all.genes, classification=all.classes)
sample.gene.counts <- merge(sample.gene.counts, x[, .(coverage = sum(coverage)), by=.(sample, gene, classification)], by=c('sample', 'gene', 'classification'), all.x=TRUE)
sample.gene.counts[is.na(coverage), coverage := 0]
sample.gene.counts[, p := coverage/max(1,sum(coverage)), by=.(sample, gene)]
sample.gene.counts <- merge(sample.gene.counts, gene.counts, by=c('gene', 'classification'))[, .(sample, gene, classification, coverage, p, kld = ifelse(p > 0, p*log(p/q), 0))]

gene.class.z <- sample.gene.counts[, .(p.avg = mean(p), p.dev = sd(p)), by=.(gene, classification)]
sample.gene.counts <- merge(sample.gene.counts, gene.class.z, by=c('gene', 'classification'))[, .(sample, gene, classification, coverage, p, kld, z = ifelse(p.dev > 0, (p - p.avg)/p.dev, 0))]

sample.gene.kld <- sample.gene.counts[, .(coverage = sum(coverage), abnormal.coverage = sum(coverage * (classification != 'normal')), sum.z2 = sum(z^2), kld = sum(kld)), by=.(sample, gene)]

gene.kld <- sample.gene.kld[kld > 0, .(sx = sum(kld), slx = sum(log(kld)), sxlx = sum(kld * log(kld)), n = length(kld)), by=.(gene)]
gene.kld[, k := n*sx / (n*sxlx - sx*slx)]
gene.kld[, theta := (1/n^2) * (n*sxlx - sx*slx)]

sample.gene.kld <- merge(sample.gene.kld, gene.kld, by='gene')[, .(sample, gene, coverage, abnormal.coverage, sum.z2, pval = pgamma(kld, shape = k, scale = theta, lower.tail=FALSE))]
sample.gene.kld[, sample.pval.rank := rank(pval), by=.(sample)]
sample.gene.kld[, sample.z2.rank := rank(-sum.z2), by=.(sample)]
sample.gene.kld[, sample.rank.prod := sqrt(sample.pval.rank * sample.z2.rank)]

if (TRUE) {
    genc = rtracklayer::import("/Users/tom.conway/data/hg38/gencode.v37.basic.annotation.gtf.gz")
    gx <- unique(data.table(gene = genc$gene_name, type = genc$gene_type))
    sample.gene.kld <- merge(sample.gene.kld, gx, by='gene')
}

top.hits.summary <- data.table()
top.hits.junctions <- data.table()

for (sn in all.samples) {
    ss <- sample.gene.kld[sample==sn & coverage >= 50 & abnormal.coverage >= 10 & sample.rank.prod < 10]
    xx <- x[sample==sn & gene %in% ss$gene]
    top.hits.summary <- rbind(top.hits.summary, ss)
    top.hits.junctions <- rbind(top.hits.junctions, xx)
}
top.hits.summary <- top.hits.summary[, .(sample, gene, coverage, abnormal.coverage, sum.z2, pval, sample.z2.rank, sample.pval.rank, sample.rank.prod, type)][order(sample, sample.rank.prod)]

fwrite(top.hits.summary, "top-hits-summary.tsv", sep="\t", col.names=T)
fwrite(top.hits.junctions, "top-hits-junctions.tsv", sep="\t", col.names=T)
