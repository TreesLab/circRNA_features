#! /usr/bin/env Rscript

con <- file('stdin')
on.exit(close(con))

data <- readLines(con)

data_1 <- as.numeric(unlist(strsplit(data[1], ',')))
data_2 <- as.numeric(unlist(strsplit(data[2], ',')))


# wilcox.test
res <- wilcox.test(x=data_1, y=data_2, conf.int=TRUE)


# coin::wilcox_test
all_data <- data.frame(x=c(data_1, data_2), y=factor(rep(c(1, 0), c(length(data_1), length(data_2)))))
res_coin <- coin::wilcox_test(x ~ y, data = all_data, distribution = "exact", conf.int = TRUE)

cat(res$p.value, res$conf.int[1], res$conf.int[2], as.numeric(res$estimate), 
    statistic(res_coin), pvalue(res_coin), confint(res_coin), '\n')
