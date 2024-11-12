pdf("Biomarkers and Pain Analysis Report.pdf")
setwd("/Users/DI/R practice/")
library(tidyverse)
library(readxl)
library(ggplot2)
library(randomForest)
library(e1071)
library(modelr)
library(ggcorrplot)
library(fitdistrplus)
library(stargazer)
rm(list = ls())
# Clear the current workspace

load_data <- function(){
    biomarkers <- read_excel("biomarkers.xlsx") 
    convariates <- read_excel("covariates.xlsx")
    convariates <- set_names(convariates, c("PatientID", "Age", "Sex", "Smoker", "Vas-at-inclusion", "Vas-12months"))
    biomarkers <- separate(biomarkers, "Biomarker", c("PatientID", "timepoint"), sep = "-", convert = TRUE)
    data <- left_join(biomarkers, convariates, by = "PatientID")
    data <- drop_na(data)
    names(data) <- make.names(names(data))
    data$Smoker <- as.factor(data$Smoker)
    return(data)
}

#Data grouping: According to certain conditions, the data is divided into two groups, a and b, and the timepoint is 0weeks(that is, at inclusion).
division_smoker_nonsmoker <- function(data){
    smoker <- filter(data, Smoker == 1, timepoint == '0weeks')
    nonsmoker <- filter(data, Smoker == 2, timepoint == '0weeks')
    output <- list(a=smoker, b=nonsmoker)
    return(output)
}

## Plot the histogram of data distribution
plot_data_distribute <- function(data, title=""){
    img <- ggplot(data.frame(x = data), aes(x)) +
            geom_histogram(bins = 20, fill = "lightblue", color = "black") +
            labs(title = title, x = "Value", y = "Frequency")
    print(img)
    return(img)
}

## Plot a scatter plot of two variables
plot_scatter <- function(data, x, y){
    img <- ggplot(data, aes_string(x = x, y = y)) +
            geom_point() + geom_smooth(method = "loess", color = "black", fill = "lightgray") +
            labs(title = sprintf("%s vs %s", x, y), x = x, y = y)
    print(img)
}

#Mapping heat map
plot_heat_map <- function(data,biomarkers){
    correlation <- round(cor(data[, c(biomarkers, "Vas.12months")]),3)
    p <- cor_pmat(data[, c(biomarkers, "Vas.12months")])
    img <- ggcorrplot(correlation, p.mat = p)
    print(img)
}

# Check normality
check_noml <- function(x, y, p = 0.05) {
  z = c(x, y)
  normality_test = shapiro.test(z)
  if(normality_test$p.value > p) {
    return(TRUE)  
  } else {
    return(FALSE) 
  }
}


#Check whether the two samples are from the same distribution
ks_test <- function(data1, data2, distri){
    fit1 <- fitdist(data1, distri)
    fit2 <- fitdist(data2, distri)
    if(fit1$loglik < fit2$loglik){
        fit <- fit1
        data <- data2
    }else {
        fit <- fit2
        data <- data1
    }
    eval(parse(text = sprintf("result <- ks_test_%s(fit, data)", distri)))
    return(result)
}

ks_test_gamma <- function(fit1, data2){
    result = ks.test(data2, 'pgamma', shape = fit1$estimate["shape"], rate = fit1$estimate["rate"])
    return(result)
}

ks_test_lnorm <- function(fit1, data2){
    result = ks.test(data2, 'plnorm', meanlog = fit1$estimate["meanlog"], sdlog = fit1$estimate["sdlog"])
    return(result)
}

ks_test_weibull <- function(fit1, data2){
    result = ks.test(data2, 'pweibull', shape = fit1$estimate["shape"], scale = fit1$estimate["scale"])
    return(result)
}

ks_test_exp <- function(fit1, data2){
    result = ks.test(data2, 'pexp', rate = fit1$estimate["rate"])
    return(result)
}

ks_test_logis <- function(fit1, data2){
    result = ks.test(data2, 'plogis', location = fit1$estimate["location"], scale = fit1$estimate["scale"])
    return(result)
}

ks_test_cauchy <- function(fit1, data2){
    result = ks.test(data2, 'pcauchy', location = fit1$estimate["location"], scale = fit1$estimate["scale"])
    return(result)
}

ks_test_unif <- function(fit1, data2){
    result = ks.test(data2, 'punif', min = fit1$estimate["min"], max = fit1$estimate["max"])
    return(result)
}

#Check all distributions, return the most likely distribution, and plot a fit for all distributions
test_all_distribution <- function(input1, input2, name){
    most_likelihood <- -10000
    # log0 = -inf
    most_likelihood_distri <- ""
    prop_distri <- c("gamma", "lnorm", "weibull", "exp", "logis","cauchy","unif")
    input <- c(input1, input2)
    for(distri in prop_distri){
        result <- fitdist(input, distri)
        log_likelihood <- result$loglik
        if(log_likelihood > most_likelihood){
            most_likelihood <- log_likelihood
            most_likelihood_distri <- distri
        }
        plot(result, hist = TRUE)
        title(main=sprintf("%s fit for %s", distri, name),col.main="red")
    }
    return(most_likelihood_distri)
}

#Check for significant differences in individual features
hypothesis_testing <- function(a, b, name, p = 0.05, use_para = TRUE){
    print(sprintf("check Biomarker: %s", name))
    a <- a[[name]]
    b <- b[[name]]
    plot_data_distribute(a, name)
    plot_data_distribute(b, name)
    if(check_noml(a, b, p)){
        print("two sample are from normal distribution")
        t_result <- t.test(a, b)
        f_result <- var.test(a, b)
        t_test_p <- t_result$p.value
        f_test_p <- f_result$p.value
        print(sprintf("t-test p-value is %f, f-test p-value is %f", t_test_p, f_test_p))
        p_value = t_test_p
    }else{
        most_likelihood_distri <- test_all_distribution(a, b, name)
        hypothesis_testing_by_para <- ks_test(a, b, most_likelihood_distri)
        print(sprintf("most likely distribution is %s, from the same distribution p-value is %f",
                                     most_likelihood_distri, hypothesis_testing_by_para$p.value))
        non_para_test <- wilcox.test(a, b)
        print(sprintf("non-parametric test p-value is %f", non_para_test$p.value))
        if(use_para){
            p_value = hypothesis_testing_by_para$p.value
        }else{
            p_value = non_para_test$p.value
        }
    }
    return(p_value)
}


#Test for significant differences in each feature (9) in a and b
hypothesis_testing_all <- function(a, b, p = 0.05){
    biomarkers <- c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1")
    for(biomarker in biomarkers){
        p_value <- hypothesis_testing(a, b, biomarker, p)
        if(p_value < p){
            # error_1 <- error_1 * (1 - p_value)
            print(sprintf("biomarker %s is significant different ", biomarker))
        }else{
            print(sprintf("biomarker %s is not significant different", biomarker))
        }
        cat("\n")
    }
}

#Analyze all the data
analysis <- function(data, divider){
    output <- divider(data)
    a <- output$a
    b <- output$b
    hypothesis_testing_all(a, b, 0.05/9)
}

# Fitting model
# （1）linear model
linear_model <- function(train, ...){
    model <- lm(Vas.12months ~ 0 + OPG + IL.6 + CSF.1 + Vas.at.inclusion, data = train)
    alpha <- 0.05
    print(confint(model, level = 1 - alpha))
    stargazer(model, type = "text")
    plot(model$fitted.values, resid(model), 
         xlab = "Fitted Values", 
         ylab = "Residuals")
    abline(h = 0, col = "red", lwd = 2)
    qqnorm(resid(model))
    qqline(resid(model), col = "red", lwd = 2)
    return(model)
}
# （2）random Forest model
randomForest_model <- function(train, ...){
    model <- randomForest(Vas.12months ~ IL.8 + TGF.beta.1  + VEGF.A + IL.6 +  CSF.1 + Vas.at.inclusion,
                             data = train, ntree = 300, ,important = TRUE, proximity = TRUE)
    return(model)
}
# （3）SVM model
svm_model <- function(train, cost = 10){
    model <- svm(Vas.12months ~ IL.6 +  OPG + CSF.1 + Vas.at.inclusion,
                 data = train, kernel = 'linear', cost = cost)
    return(model)
}


#Model training and evaluation (Train and test the model)
train_test_model <- function(dataset, model, ...){
    sample <- sample(c(TRUE, FALSE), nrow(dataset), replace = TRUE, prob = c(0.8, 0.2))
    train <- dataset[sample,]
    test <- dataset[!sample,]
    model <- model(train, ...)
    prediction <- predict(model, newdata = test)
    r_squared <- cor(prediction, test$Vas.12months)^2
    rmse <- sqrt(mean((prediction - test$Vas.12months)^2))
    return(rmse)
}

data <- load_data()

# q1
analysis(data, divider = division_smoker_nonsmoker) 

#q2
biomarkers_convar <- c("IL.8", "VEGF.A", "OPG", "TGF.beta.1", "IL.6", "CXCL9", "CXCL1", "IL.18", "CSF.1", "Vas.at.inclusion")

for(biomarker in biomarkers_convar){
    plot_scatter(data, biomarker, "Vas.12months")
}

dataset <- filter(data, timepoint == '0weeks')

plot_heat_map(dataset, biomarkers_convar)

average_rmse <- 0
#The training test was repeated 1000 times and the mean root mean square error was calculated
# rounds <- 10000
rounds <- 1
for(i in c(1:rounds)){
    #linear_model, randomForest_model, svm_model
    rmse <- train_test_model(dataset, linear_model)
    average_rmse <- average_rmse + rmse
}
average_rmse <- average_rmse / rounds
print(sprintf("average rmse: %f", average_rmse))

dev.off()

