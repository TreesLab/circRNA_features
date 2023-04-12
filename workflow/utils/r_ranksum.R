#! /usr/bin/env Rscript

con <- file('stdin')
on.exit(close(con))

data <- readLines(con)

data_1 <- as.numeric(unlist(strsplit(data[1], ',')))
data_2 <- as.numeric(unlist(strsplit(data[2], ',')))
all_data <- data.frame(x=c(data_1, data_2), y=factor(rep(c(1, 0), c(length(data_1), length(data_2)))))


# wilcox.test
res <- wilcox.test(x=data_1, y=data_2, conf.int=TRUE)

# wilcox_effsize
res_rstatix <- rstatix::wilcox_effsize(all_data, x ~ y, ci=TRUE)

cat(res$p.value, res$conf.int[1], res$conf.int[2], as.numeric(res$estimate), 
    res_rstatix$effsize, res_rstatix$conf.low, res_rstatix$conf.high, '\n')
