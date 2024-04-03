library(data.table)
library(tidyverse)

# Read sumstats.
# For some reason fread couldn't read the header, so I read it separately
# and added to the sumstats.
header = read.table("/PATH/TO/DIR/rhythym_sumstat.txt", header = TRUE, nrow = 1)
rhy_sumstat = fread("/PATH/TO/DIR/rhythym_sumstat.txt", skip = 1, header = FALSE)
setnames(rhy_sumstat, colnames(header))

# Calculate chisquares from existing pvalues,
# divide them by the genomic inflation factor,
# recalculate pvalues and add them to the sumstats.
chisq = qchisq(rhy_sumstat$pvalue, 1, lower.tail = FALSE)
newchisq = chisq/1.474
rhy_sumstat$pvalue = pchisq(newchisq, df = 1, lower.tail = FALSE)

# write
fwrite(rhy_sumstat, "/PATH/TO/DIR/rhythm_sumstat_wGCcorr.tab", quote = F, sep = "\t", row.names = F, na = "NA")
