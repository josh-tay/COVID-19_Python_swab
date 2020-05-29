library(xlsx)
p <- read.xlsx("ct_masked.xlsx", 1, as.data.frame=TRUE, header = TRUE)

# create categorical outcomes for each swab
# FLOQSwab
p$COBAS1_FLOQ_A_cat <- NA
p$COBAS1_FLOQ_A_cat[p$COBAS1_FLOQ_A != "Neg"] <- "Pos"
p$COBAS1_FLOQ_A_cat[p$COBAS1_FLOQ_A == "Neg"] <- "Neg"
p$COBAS2_FLOQ_A_cat <- NA
p$COBAS2_FLOQ_A_cat[p$COBAS2_FLOQ_A != "Neg"] <- "Pos"
p$COBAS2_FLOQ_A_cat[p$COBAS2_FLOQ_A == "Neg"] <- "Neg"
# overall result for FLOQSwab
p$Result_FLOQ_final <- NA
p$Result_FLOQ_final[p$COBAS1_FLOQ_A != "Neg"] <- "Pos"
p$Result_FLOQ_final[p$COBAS1_FLOQ_A == "Neg"] <- "Neg"
p$Result_FLOQ_final[p$COBAS1_FLOQ_A == "Neg" & p$COBAS2_FLOQ_A != "Neg" 
                    & p$COBAS2_FLOQ_A_repeat != "Neg"] <- "Pos"

# Python swab
p$COBAS1_POLY_A_cat <- NA
p$COBAS1_POLY_A_cat[p$COBAS1_POLY_A != "Neg"] <- "Pos"
p$COBAS1_POLY_A_cat[p$COBAS1_POLY_A == "Neg"] <- "Neg"
p$COBAS2_POLY_A_cat <- NA
p$COBAS2_POLY_A_cat[p$COBAS2_POLY_A != "Neg"] <- "Pos"
p$COBAS2_POLY_A_cat[p$COBAS2_POLY_A == "Neg"] <- "Neg"
# overall result for Python swab
p$Result_POLY_final <- NA
p$Result_POLY_final[p$COBAS1_POLY_A != "Neg"] <- "Pos"
p$Result_POLY_final[p$COBAS1_POLY_A == "Neg"] <- "Neg"
p$Result_POLY_final[p$COBAS1_POLY_A == "Neg" & p$COBAS2_POLY_A != "Neg" 
                    & p$COBAS2_POLY_A_repeat != "Neg"] <- "Pos"

# Concordance between the swabs
p$Result_concordance_final <- NA
p$Result_concordance_final[p$Result_FLOQ_final == "Pos" & p$Result_POLY_final == "Pos"] <- "True Positive"
p$Result_concordance_final[p$Result_FLOQ_final == "Neg" & p$Result_POLY_final == "Neg"] <- "True Negative"
p$Result_concordance_final[p$Result_FLOQ_final == "Pos" & p$Result_POLY_final == "Neg"] <- "False Negative"


# Cross-tabs --------------------------------------------------------------

library(fmsb)

# cross-tabs cases and control (final)
t <- table(p$Result_FLOQ_final, p$Result_POLY_final)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# cross-tabs cases only (final)
p_cases <- subset(p, p$Group == "case")
t <- table(p_cases$Result_FLOQ_final, p_cases$Result_POLY_final)
Kappa.test(t, y=NULL, conf.level=0.95)
t
# 1st week cases only (final)
p1 <- subset(p_cases, p_cases$Day_of_illness < 8)
t <- table(p1$Result_FLOQ_final, p1$Result_POLY_final)
t
# 2nd week cases only (final)
p2 <- subset(p_cases, p_cases$Day_of_illness > 7)
t <- table(p2$Result_FLOQ_final, p2$Result_POLY_final)
t
Kappa.test(t, y=NULL, conf.level=0.95)

# cross-tabs for ORF1/a target
# cases and controls
t <- table(p$COBAS1_FLOQ_A_cat, p$COBAS1_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# cases only
p_cases <- subset(p, p$Group == "case")
t <- table(p_cases$COBAS1_FLOQ_A_cat, p_cases$COBAS1_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# first week
p1 <- subset(p_cases, p_cases$Day_of_illness < 8)
t <- table(p1$COBAS1_FLOQ_A_cat, p1$COBAS1_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# second week
p2 <- subset(p_cases, p_cases$Day_of_illness > 7)
t <- table(p2$COBAS1_FLOQ_A_cat, p2$COBAS1_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)


# cross-tabs for E-gene target
# cases and controls
t <- table(p$COBAS2_FLOQ_A_cat, p$COBAS2_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# cases only
p_cases <- subset(p, p$Group == "case")
t <- table(p_cases$COBAS2_FLOQ_A_cat, p_cases$COBAS2_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)
# 1st week
p1 <- subset(p_cases, p_cases$Day_of_illness < 8)
t <- table(p1$COBAS2_FLOQ_A_cat, p1$COBAS2_POLY_A_cat)
t
# 2nd week
p2 <- subset(p_cases, p_cases$Day_of_illness > 7)
t <- table(p2$COBAS2_FLOQ_A_cat, p2$COBAS2_POLY_A_cat)
t
Kappa.test(t, y=NULL, conf.level=0.95)


# Correlation plots -------------------------------------------------------

# assign CT-values for negative cases
neg_CT_value_COBAS1 <- 36 
neg_CT_value_COBAS2 <- 41 
p$COBAS1_FLOQ_A_num <- p$COBAS1_FLOQ_A
p$COBAS1_FLOQ_A_num[p$COBAS1_FLOQ_A == "Neg"] <- neg_CT_value_COBAS1
p$COBAS2_FLOQ_A_num <- p$COBAS2_FLOQ_A
p$COBAS2_FLOQ_A_num[p$COBAS2_FLOQ_A == "Neg"] <- neg_CT_value_COBAS2
p$COBAS1_POLY_A_num <- p$COBAS1_POLY_A
p$COBAS1_POLY_A_num[p$COBAS1_POLY_A == "Neg"] <- neg_CT_value_COBAS1
p$COBAS2_POLY_A_num <- p$COBAS2_POLY_A
p$COBAS2_POLY_A_num[p$COBAS2_POLY_A == "Neg"] <- neg_CT_value_COBAS2


# create categorical groups
# positive result on either swab, for either target
p$any_positive <- NA
p$any_positive[p$COBAS1_FLOQ_A_cat == "Pos" | p$COBAS1_POLY_A_cat == "Pos" | 
                 p$COBAS2_FLOQ_A_cat == "Pos" | p$COBAS2_POLY_A_cat == "Pos"] <- "yes"
p$any_positive[!(p$COBAS1_FLOQ_A_cat == "Pos" | p$COBAS1_POLY_A_cat == "Pos" | 
                 p$COBAS2_FLOQ_A_cat == "Pos" | p$COBAS2_POLY_A_cat == "Pos")] <- "no"

# positive result on either swab, for ORF1/a target
p$any_COBAS1_positive[p$COBAS1_FLOQ_A_cat == "Pos" | p$COBAS1_POLY_A_cat == "Pos" ] <- "yes"
p$any_COBAS1_positive[!(p$COBAS1_FLOQ_A_cat == "Pos" | p$COBAS1_POLY_A_cat == "Pos")] <- "no"
# positive result on either swab, for E-gene target
p$any_COBAS2_positive[p$COBAS2_FLOQ_A_cat == "Pos" | p$COBAS2_POLY_A_cat == "Pos" ] <- "yes"
p$any_COBAS2_positive[!(p$COBAS2_FLOQ_A_cat == "Pos" | p$COBAS2_POLY_A_cat == "Pos")] <- "no"

# positive results for both swabs, for ORF1/a target
p$both_COBAS1_positive[p$COBAS1_FLOQ_A_cat == "Pos" & p$COBAS1_POLY_A_cat == "Pos" ] <- "yes"
p$both_COBAS1_positive[!(p$COBAS1_FLOQ_A_cat == "Pos" & p$COBAS1_POLY_A_cat == "Pos")] <- "no"
# positive results for both swabs, for E-gene target 
p$both_COBAS2_positive[p$COBAS2_FLOQ_A_cat == "Pos" & p$COBAS2_POLY_A_cat == "Pos" ] <- "yes"
p$both_COBAS2_positive[!(p$COBAS2_FLOQ_A_cat == "Pos" & p$COBAS2_POLY_A_cat == "Pos")] <- "no"


# Cases positive for ORF1/a on either swab
p_COBAS1_pos <- subset(p, p$any_COBAS1_positive == "yes")
x_temp <- as.numeric(p_COBAS1_pos$COBAS1_FLOQ_A_num)
y_temp <- as.numeric(p_COBAS1_pos$COBAS1_POLY_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "FLOQSwab", ylab = "Python swab",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(16,16)[as.factor(p_COBAS1_pos$COBAS1_POLY_A_cat)],
     col = c("blue", "black")[as.factor(p_COBAS1_pos$COBAS1_POLY_A_cat)],
     xlim = c(16, 37), ylim =c(16,37),
     main = "Ct values for ORF1/a (n = 29)")
legend("bottomright", pch = c(16,16), col=c("black", "blue"), 
     legend = c("ORF1/a positive for both swabs", "ORF1/a positive for FLOQSwab only"), cex=0.7)
abline(coef = c(0,1), col = "grey40", lty = "dashed")
abline(lm(y_temp ~ x_temp), col = "red", lty = "dotted")

# the same plot, but showing the order of swabs
plot(x_temp,y_temp, xlab = "FLOQSwab", ylab = "Python swab", 
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     col = c("green","blue")[as.numeric(p_COBAS1_pos$Odd_or_Even)],
     pch = c(1,16)[as.factor(p_COBAS1_pos$COBAS1_POLY_A_cat)],
     main = "Ct values for ORF1/a (n = 29)", xlim = c(16, 37), ylim =c(16,37))
legend("bottomright", pch = c(16,16), col=c("green", "blue"),
      legend = c("Python swab first", "FLOQSwab first"), cex=1)
abline(coef = c(0,1), col = "grey40", lty = "dashed")
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")

# Cases positive for E-gene on either swab
p_COBAS2_pos <- subset(p, p$any_COBAS2_positive == "yes")
# create new variable to mark out COBAS1 polymer negative cases
p_COBAS2_pos$COBAS2_POLY_A_cat_both <- p_COBAS2_pos$COBAS2_POLY_A_cat
p_COBAS2_pos$COBAS2_POLY_A_cat_both[p_COBAS2_pos$COBAS1_POLY_A_cat == "Neg" & p_COBAS2_pos$COBAS2_POLY_A_cat == "Pos"] <- "ORF1a discordant" 
x_temp <- as.numeric(p_COBAS2_pos$COBAS2_FLOQ_A_num)
y_temp <- as.numeric(p_COBAS2_pos$COBAS2_POLY_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "FLOQSwab", ylab = "Python swab",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(1,16,16)[as.factor(p_COBAS2_pos$COBAS2_POLY_A_cat_both)],
     col = c("red","blue","black")[as.factor(p_COBAS2_pos$COBAS2_POLY_A_cat_both)],
     main = "Ct values for E-gene (n = 32)", xlim = c(16, 42), ylim =c(16,42) )
legend("bottomright", pch = c(16,1), col=c("black", "red"), 
       legend = c("E-gene positive for both swabs", "E-gene positive for FLOQSwab only"), cex=0.7)
abline(coef = c(0,1), col = "grey40", lty = "dashed")
abline(lm(y_temp ~ x_temp), col = "red", lty = "dotted")

# the same plot, but showing the order of swabs
plot(x_temp,y_temp, xlab = "FLOQSwab", ylab = "Python swab",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(1,16,16)[as.factor(p_COBAS2_pos$COBAS2_POLY_A_cat_both)],
     col = c("green","blue")[as.numeric(p_COBAS2_pos$Odd_or_Even)],
     main = "Ct values for E-gene (n = 32)", xlim = c(16, 42), ylim =c(16,42) )
legend("bottomright", pch = c(16,16), col=c("green", "blue"), 
       legend = c("Python swab first", "FLOQSwab first"), cex=1)
abline(coef = c(0,1), col = "grey40", lty = "dashed")
abline(lm(y_temp ~ x_temp), col = "red", lty = "dotted")


# Time Course -------------------------------------------------------------
# plot Ct values with day of illness
# cases only
p_cases <- subset(p, p$Group == "case")
library(ggplot2)

# ORF1/a for the FLOQSwab
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS1_FLOQ_A)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "ORF1/a Ct value", pch=20, col = alpha("black", 0.4),
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     main = "FLOQSwab")
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")
# include negative cases
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS1_FLOQ_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "ORF1/a Ct value",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(4,16)[as.numeric(as.factor(p_cases$COBAS1_FLOQ_A_cat))], col = alpha("black", 0.4),
     main = "FLOQSwab")
legend("bottomright", pch = c(16,4), legend = c("positive", "negative"), cex=1)
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")

# ORF1/a for the Python swab
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS1_POLY_A)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "ORF1/a Ct value", pch=20, col = alpha("black", 0.4),
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     main = "Python swab")
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")
# include negative cases
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS1_POLY_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "ORF1/a Ct value",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(4,16)[as.numeric(as.factor(p_cases$COBAS1_POLY_A_cat))], col = alpha("black", 0.4),
     main = "Python swab")
legend("bottomright", pch = c(16,4), legend = c("positive", "negative"), cex=1)
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")


# E-gene for the FLOQSwab
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS2_FLOQ_A)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "E-gene Ct value", pch=20, col = alpha("black", 0.4),
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     main = "FLOQSwab")
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")
# include negative cases
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS2_FLOQ_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "E-gene Ct value",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(4,16)[as.numeric(as.factor(p_cases$COBAS2_FLOQ_A_cat))], col = alpha("black", 0.4),
     main = "FLOQSwab")
legend("bottomright", pch = c(16,4), legend = c("positive", "negative"), cex=1)
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")

# E-gene for the Python swab
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS2_POLY_A)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "E-gene Ct value", pch=20, col = alpha("black", 0.4),
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     main = "Python swab)")
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")
# include negative cases
x_temp <- p_cases$Day_of_illness
y_temp <- as.numeric(p_cases$COBAS2_POLY_A_num)
cor_temp <- cor.test(x_temp,y_temp, method="pearson")
plot(x_temp,y_temp, xlab = "Day of Illness", ylab = "E-gene Ct value",
     sub = paste("r =", signif(cor_temp$estimate,3), "p =", signif(cor_temp$p.value,3)), font.lab = 2,
     pch = c(4,16)[as.numeric(as.factor(p_cases$COBAS2_POLY_A_cat))], col = alpha("black", 0.4),
     main = "Python swab")
legend("bottomright", pch = c(16,4), legend = c("positive", "negative"), cex=1)
abline(lm(y_temp ~ x_temp), col = "red", lty="dotted")


# Descriptive statistics ------------------------------------------------------
# compare mean Ct values
p_ORF1a <- p[p$COBAS1_FLOQ_A_cat == "Pos" & p$COBAS1_POLY_A_cat == "Pos",]
t.test(as.numeric(p_ORF1a$COBAS1_FLOQ_A_num), as.numeric(p_ORF1a$COBAS1_POLY_A_num), paired=TRUE)
p_egene <- p[p$COBAS2_FLOQ_A_cat == "Pos" & p$COBAS2_POLY_A_cat == "Pos",]
t.test(as.numeric(p_egene$COBAS2_FLOQ_A_num), as.numeric(p_egene$COBAS2_POLY_A_num), paired=TRUE)

# Day of illness
range(as.numeric(p$Day_of_illness), na.rm=T)
median(as.numeric(p$Day_of_illness), na.rm=T)



