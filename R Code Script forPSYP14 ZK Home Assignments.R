#Statistical Report Home Assignments PSYP14-HT20  Zoltan Kekecs
#R Code Script for the statstical analysis
#Author: Marcel Borgsen

# Packages  -----------------------------------------------------

#Here are listed all used Packages for all assignments

library(gridExtra) #for the error correction visualization
library(psych) # for the describe	function
library(tidyverse) # for the tidy coding
library(lm.beta) # for lm.beta
library(car) # for residualPlots, vif, pairs.panels, ncvTest
library(lmtest) # bptest
library(sandwich) # for coeftest vcovHC estimator
library(boot) # for bootstrapping
library(lmboot) # for wild bootsrapping
library(cAIC4) # for cAIC
library(r2glmm) # for r2beta
library(lme4) # for lmer
library(lmerTest) # to get singificance test in lmer
library(MuMIn) # for r.squaredGLMM
library(optimx) # for optimx optimizer
library(apaTables) # used for making APA Style tables in word

# Costum Functions --------------------------------------------------------

#Using costum function for coef_table
coef_table = function(model) {
  require(lm.beta)
  mod_sum = summary(model)
  mod_sum_p_values = as.character(round(mod_sum$coefficients[,
                                                             4], 3))
  mod_sum_p_values[mod_sum_p_values != "0" & mod_sum_p_values !=
                     "1"] = substr(mod_sum_p_values[mod_sum_p_values != "0" &
                                                      mod_sum_p_values != "1"], 2, nchar(mod_sum_p_values[mod_sum_p_values !=
                                                                                                            "0" & mod_sum_p_values != "1"]))
  mod_sum_p_values[mod_sum_p_values == "0"] = "<.001"
  mod_sum_table = cbind(as.data.frame(round(cbind(coef(model),
                                                  confint(model), c(0, lm.beta(model)$standardized.coefficients[c(2:length(model$coefficients))])),
                                            2)), mod_sum_p_values)
  names(mod_sum_table) = c("b", "95%CI lb", "95%CI ub", "Std.Beta",
                           "p-value")
  mod_sum_table["(Intercept)", "Std.Beta"] = "0"
  return(mod_sum_table)
}


#Using costum function for stdCoef.merMod

stdCoef.merMod <- function(object) {
  sdy <- sd(getME(object,"y"))
  sdx <- apply(getME(object,"X"), 2, sd)
  sc <- fixef(object)*sdx/sdy
  se.fixef <- coef(summary(object))[,"Std. Error"]
  se <- se.fixef*sdx/sdy
  return(data.frame(stdcoef=sc, stdse=se))}

# Research Question 1 Start -----------------------------------------------------

data_sample_1 = read.csv("https://tinyurl.com/ha-dataset1")


# RQ1 Data Diagnostics and Vizualisation ----------------------------------------------------

View(data_sample_1)
str(data_sample_1)
summary(data_sample_1)
describe(data_sample_1)

#After using the Summary function, following coding errors are visible: 
#There is an coding error in the age variable with an age of 444, which is at the moment biologically impossible for humans.
#There is an coding error in the STAI_trait variable with a value of 3.9, since this variable start at minimum with a value of 20.

#Visualization of the data distribution of age, sex, STAI, pain catastrophizing, mindfulness, cortisol measures and pain.

data_sample_1 %>%
  ggplot() +
  aes(x = age) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = sex) +
  geom_bar()

data_sample_1 %>%
  ggplot() +
  aes(x = STAI_trait) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = pain_cat) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = mindfulness) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = cortisol_serum) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = cortisol_saliva) +
  geom_histogram()

data_sample_1 %>%
  ggplot() +
  aes(x = pain) +
  geom_histogram()

#Visualization of the linear correlations between the predictors and the outcome 

#For a better overview, here I used the geom bar function and filled pain with sex 
data_sample_1 %>%
  ggplot() + 
  aes(x = pain, fill = sex) + 
  geom_bar()


#The coding error influences the correlation btw. age and pain
data_sample_1 %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

#The coding error influences the correlation btw. STAI_trait and pain
data_sample_1 %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_1 %>%
  ggplot() + 
  aes(x = pain_cat, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_1 %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_1 %>%
  ggplot() + 
  aes(x = cortisol_serum, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_1 %>%
  ggplot() + 
  aes(x = cortisol_saliva, y = pain) + 
  geom_point() +
  geom_smooth()

# RQ 1 Coding Error Correction -----------------------------------------------

#Coding Error Correction, in order to keep the original data set, a new data set is build up
data_corrected <- data_sample_1 %>%
  mutate(age = replace(age, age=="444", NA)) %>%
  mutate(STAI_trait = replace(STAI_trait, STAI_trait == "3.9", NA))


#Since STRAIT_trait is still numeric after correction, here the transformation to integer vector type
data_corrected$STAI_trait <-as.integer(data_corrected$STAI_trait)
str(data_corrected$STAI_trait)

#For further analysis the rows with NA values will be removed (exclude cases with the coding errors)
data_corrected <- drop_na(data_corrected) 
view(data_corrected)
summary(data_corrected)
str(data_corrected)

#Plot for comparing the differences for and after error correction in age
old_plot_age_smooth <-
  data_sample_1 %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

new_plot_age_smooth <-
  data_corrected %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

grid.arrange(old_plot_age_smooth, new_plot_age_smooth, ncol=2)

#Plot for comparing the differences for and after error correction in STAI_trait
old_plot_STAI_trait_smooth <-
  data_sample_1 %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

new_plot_STAI_trait_smooth <-
  data_corrected %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

grid.arrange(old_plot_STAI_trait_smooth, new_plot_STAI_trait_smooth, ncol=2)

# RQ1 Regressionmodels ---------------------------------------------------

#Regressionmodel 1
mod_1 = lm(pain ~ age + sex, data = data_corrected)
summary(mod_1)
confint(mod_1)
lm.beta(mod_1)
coef_table(mod_1)

#Regressionmodel 2
mod_2 = lm(pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum + cortisol_saliva + mindfulness, data = data_corrected)
summary(mod_2)
confint(mod_2)
lm.beta(mod_2)
coef_table(mod_2)

# RQ1 Model Diagnostics ---------------------------------------------------

#Assumption of normality
residuals_mod_2 = enframe(residuals(mod_2))
residuals_mod_2 %>% 
ggplot() + aes(x =value) + geom_histogram()
describe(residuals(mod_2))
#Skew and Kurtosis between -1 and +1 

# Assumption of linearity
mod_2 %>% residualPlots()


#Assumption of homoscedasticity
mod_2 %>% plot(which = 3)
mod_2 %>% ncvTest()
mod_2 %>% bptest()

#Assumption of no excess in multicollinearity
mod_2 %>% vif()
data_corrected %>% 
select(pain, age, sex, STAI_trait,pain_cat,cortisol_serum,cortisol_saliva, mindfulness) %>%
 pairs.panels(col = "red", lm = T)


#Checking for influental outliner
mod_2 %>% plot(which = 5)
mod_2 %>% plot(which = 4)
mod_2 %>% plot(which = 2)


# RQ1 Final Model ---------------------------------------------------------

#cortisol serum and saliva are highly correlated, therefore saliva will be excluded
mod_final = lm(pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum + mindfulness, data = data_corrected)
summary(mod_final)
confint(mod_final)
lm.beta(mod_final)
coef_table(mod_final)

#Re-Check Model diagnostics

#Assumption of normality
residuals_mod_final = enframe(residuals(mod_final))
residuals_mod_final %>% 
  ggplot() + aes(x =value) + geom_histogram()
describe(residuals(mod_final))
#Skew and Kurtosis between -1 and +1 

#Assumption of linearity
mod_final %>% residualPlots()

#Assumption of homoscedasticity
mod_final %>% plot(which = 3)
mod_final %>% ncvTest()
mod_final %>% bptest()

#Assumption of no excess in multicollinearity
mod_final %>% vif()

#Checking for influental outliner
mod_final %>% plot(which = 5)
mod_final %>% plot(which = 4)
mod_final %>% plot(which = 2)

#Visualization outliners of Age and pain

data_corrected %>% 
  mutate(rownum = row.names(data_corrected)) %>% 
  ggplot() +
  aes(x = age, y = pain, label = rownum) + geom_label()

data_corrected %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth(method = "lm")

#Check for differences after excluding cortisol saliva
summary(mod_2)$adj.r.squared
summary(mod_final)$adj.r.squared
AIC(mod_2)
AIC(mod_final)
summary(mod_2)
summary(mod_final)
# Slightly but no significant difference found, therefore the decision of excluding cortisol saliva was right. 


# RQ1 Model Comparison ----------------------------------------------------

summary(mod_1)
summary(mod_final)
AIC(mod_1)
AIC(mod_final)
coef_table(mod_1)
coef_table(mod_final)

anova(mod_1,mod_final)

#APA style tables for the report
apa.reg.table(mod_1, filename = "Table1_APA.doc", table.number = 1)
apa.reg.table(mod_final, filename = "Table2_APA.doc", table.number = 2)


# Research Question 2 Start -----------------------------------------------------

#Since in Research question 1 two cases were excluded, here the corrected data set will be used again. 
#Data Diagnostics
View(data_corrected)
str(data_corrected)
summary(data_corrected)
describe(data_corrected)

#The variables age, sex, STAI_trait, pain catastrophizing, mindfulness, serum cortisol were therefore already checked

#Visualization of the new variables weight, IQ, household income.

data_corrected %>%
  ggplot() +
  aes(x = weight) +
  geom_histogram()

data_corrected %>%
  ggplot() +
  aes(x = IQ) +
  geom_histogram()

data_corrected %>%
  ggplot() +
  aes(x = household_income) +
  geom_histogram()

#There is a coding error in the household income variable

# RQ 2 Coding Error Correction -----------------------------------------------

#Coding Error Correction, in order to keep the original corrected data set, a new corrected data set is build up

data_corrected2 <- data_corrected %>%
  mutate(household_income = replace(household_income, household_income=="-3732", NA)) 


#For further analysis the rows with NA values will be removed (exclude cases with the coding errors)

data_corrected2 <- drop_na(data_corrected2) 

view(data_corrected2)
summary(data_corrected2)
str(data_corrected2)

#Plot for comparing the differences for and after error correction in household_income
old_plot_household_income_smooth <-
  data_corrected %>%
  ggplot() + 
  aes(x = household_income, y = pain) + 
  geom_point() +
  geom_smooth()

new_plot_household_income_smooth <-
  data_corrected2 %>%
  ggplot() + 
  aes(x = household_income, y = pain) + 
  geom_point() +
  geom_smooth()

grid.arrange(old_plot_household_income_smooth, new_plot_household_income_smooth, ncol=2)


# RQ 2 Backward regression model ------------------------------------------

mod_backward_initial = lm(pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum + mindfulness + weight + IQ + household_income, data = data_corrected2)
summary(mod_backward_initial)



# RQ 2 Model Diagnostics --------------------------------------------------

mod_backward_initial%>% plot(which = 2)

#Assumption of normality
residuals_mod_backward_initial = enframe(residuals(mod_backward_initial))
residuals_mod_backward_initial %>% 
  ggplot() + aes(x =value) + geom_histogram()
# checking skew and kurtosis
describe(residuals(mod_backward_initial))


# #Assumption of linearity
mod_backward_initial %>% residualPlots()


##Assumption of homodescasity
mod_backward_initial %>% plot(which = 3)
mod_backward_initial %>% ncvTest()
mod_backward_initial %>% bptest()


#Assumption of no excess in mulitcorrelatility

mod_backward_initial %>% vif()

#Running the backward regression 

backward_model = step(mod_backward_initial, direction = "backward")

#Since there is a new case excluded, the final model of Research Question 1 refeitted on the corrected data set

theory_based_model = lm(pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum + cortisol_saliva + mindfulness, data = data_corrected2)

# RQ2 Model Comparison ----------------------------------------------

summary(mod_backward_initial)
AIC(mod_backward_initial)
summary(backward_model)
AIC(backward_model)
summary(theory_based_model)
AIC(theory_based_model)

anova(mod_backward_initial,backward_model)
anova(backward_model,theory_based_model)

#APA Table for the report
coef_table(backward_model)
apa.reg.table(backward_model, filename = "Table3_APA.doc", table.number = 3)

# RQ2 New Dataset Diagnostics and Error Correction---------------------------------------------------------

data_sample_2 = read.csv("https://tinyurl.com/ha-dataset2")

#Data Diagnostics
View(data_sample_2)
str(data_sample_2)
summary(data_sample_2)
describe(data_sample_2)

#Visualization of the new variables weight, IQ, household income.

data_sample_2 %>%
  ggplot() +
  aes(x = weight) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = IQ) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = household_income) +
  geom_histogram()

#Visualization of the data distribution of age, sex, STAI, pain catastrophizing, mindfulness, cortisol measures and pain.

data_sample_2 %>%
  ggplot() +
  aes(x = age) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = sex) +
  geom_bar()

data_sample_2 %>%
  ggplot() +
  aes(x = STAI_trait) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = pain_cat) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = mindfulness) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = cortisol_serum) +
  geom_histogram()

data_sample_2 %>%
  ggplot() +
  aes(x = pain) +
  geom_histogram()

#Visualization of the linear correlations between the predictors and the outcome 

#For a better overview, here I used the geom bar function and filled pain with sex 
data_sample_2 %>%
  ggplot() + 
  aes(x = pain, fill = sex) + 
  geom_bar()

data_sample_2 %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = pain_cat, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = cortisol_serum, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = IQ, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = household_income, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_2 %>%
  ggplot() + 
  aes(x = weight, y = pain) + 
  geom_point() +
  geom_smooth()

# RQ2 Error Correction Dataset 2 ------------------------------------------

#Coding Error Correction, in order to keep the original corrected data set, a new corrected data set is build up

data2_corrected <- data_sample_2 %>%
  mutate(mindfulness = replace(mindfulness, mindfulness=="7.17", NA)) 

#For further analysis the rows with NA values will be removed (exclude cases with the coding errors)

data2_corrected <- drop_na(data2_corrected) 

view(data2_corrected)
summary(data2_corrected)
str(data2_corrected)

#Plot for comparing the differences for and after error correction in mindfulness
old_plot_mindfulness_smooth <-
  data2_corrected %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

new_plot_mindfulness_smooth <-
  data2_corrected %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

grid.arrange(old_plot_mindfulness_smooth,new_plot_mindfulness_smooth, ncol=2)




# RQ2 Prediction Comparison -----------------------------------------------

# calculate predicted values
pred_theory <- predict(theory_based_model, data2_corrected)
pred_back <- predict(backward_model, data2_corrected)
                        
# now we calculate the sum of squared residuals
RSS_theory  = sum((data2_corrected[, "pain"] - pred_theory)^2)
RSS_back = sum((data2_corrected[, "pain"] - pred_back)^2)

RSS_theory 
RSS_back 

# Research Question 3 Start -----------------------------------------------

#Loading data set

data_sample_3 = read.csv("https://tinyurl.com/ha-dataset3")
data_sample_4 = read.csv("https://tinyurl.com/ha-dataset4")

# RQ 3 Data Correction  ---------------------------------------------------

#Data Diagnostics Data Set 3 
View(data_sample_3)
str(data_sample_3)
summary(data_sample_3)
describe(data_sample_3)

#Visualization of the data distribution of age, sex, STAI, pain catastrophizing, mindfulness, cortisol measures and pain.

data_sample_3 %>%
  ggplot() +
  aes(x = age) +
  geom_histogram()

data_sample_3 %>%
  ggplot() +
  aes(x = sex) +
  geom_bar()

data_sample_3 %>%
  ggplot() +
  aes(x = STAI_trait) +
  geom_histogram()

data_sample_3 %>%
  ggplot() +
  aes(x = pain_cat) +
  geom_histogram()

data_sample_3 %>%
  ggplot() +
  aes(x = mindfulness) +
  geom_histogram()

data_sample_3 %>%
  ggplot() +
  aes(x = cortisol_serum) +
  geom_histogram()

data_sample_3 %>%
  ggplot() +
  aes(x = pain) +
  geom_histogram()

#Visualization of the linear correlations between the predictors and the outcome 

#For a better overview, here I used the geom bar function and filled pain with sex 
data_sample_3 %>%
  ggplot() + 
  aes(x = pain, fill = sex) + 
  geom_bar()

data_sample_3 %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_3 %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_3 %>%
  ggplot() + 
  aes(x = pain_cat, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_3 %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_3 %>%
  ggplot() + 
  aes(x = cortisol_serum, y = pain) + 
  geom_point() +
  geom_smooth()

#Coding Error Correction, in order to keep the original corrected data set, a new corrected data set is build up

data_sample_3_corrected <- data_sample_3%>%
mutate(sex = replace(sex, sex=="femlae", "female")) 

view(data_sample_3_corrected)
summary(data_sample_3_corrected)
str(data_sample_3_corrected)

#Data Diagnostics Data Set 4 
View(data_sample_4)
str(data_sample_4)
summary(data_sample_4)
describe(data_sample_4)


#Visualization of the data distribution of age, sex, STAI, pain catastrophizing, mindfulness, cortisol measures and pain.

data_sample_4 %>%
  ggplot() +
  aes(x = age) +
  geom_histogram()

data_sample_4 %>%
  ggplot() +
  aes(x = sex) +
  geom_bar()

data_sample_4 %>%
  ggplot() +
  aes(x = STAI_trait) +
  geom_histogram()

data_sample_4 %>%
  ggplot() +
  aes(x = pain_cat) +
  geom_histogram()

data_sample_4 %>%
  ggplot() +
  aes(x = mindfulness) +
  geom_histogram()

data_sample_4 %>%
  ggplot() +
  aes(x = cortisol_serum) +
  geom_histogram()

data_sample_4 %>%
  ggplot() +
  aes(x = pain) +
  geom_histogram()

#Visualization of the linear correlations between the predictors and the outcome 

#For a better overview, here I used the geom bar function and filled pain with sex 
data_sample_4 %>%
  ggplot() + 
  aes(x = pain, fill = sex) + 
  geom_bar()

data_sample_4 %>%
  ggplot() + 
  aes(x = age, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_4 %>%
  ggplot() + 
  aes(x = STAI_trait, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_4 %>%
  ggplot() + 
  aes(x = pain_cat, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_4 %>%
  ggplot() + 
  aes(x = mindfulness, y = pain) + 
  geom_point() +
  geom_smooth()

data_sample_4 %>%
  ggplot() + 
  aes(x = cortisol_serum, y = pain) + 
  geom_point() +
  geom_smooth()


#Coding Error Correction, in order to keep the original corrected data set, a new corrected data set is build up

data_sample_4_corrected <- data_sample_4 %>%
  mutate(mindfulness = replace(mindfulness, mindfulness =="6.05" , NA))

#For further analysis the rows with NA values will be removed (exclude cases with the coding errors)

data_sample_4_corrected <- drop_na(data_sample_4_corrected) 

view(data_sample_4_corrected)
summary(data_sample_4_corrected)
str(data_sample_4_corrected)



# RQ 3 First linear mixed model -------------------------------------------------

# Grouping the hospitals
data_sample_3_corrected %>%
  mutate(hospital = factor(hospital))

#Creating the first linear mixed model with the fixed predictors
mod_random = lmer(pain ~ age + sex + STAI_trait + pain_cat + cortisol_serum + mindfulness + (1 | hospital), data = data_sample_3_corrected)
summary(mod_random)
confint(mod_random)
stdCoef.merMod(mod_random)
# Calculating the marginal R squared with confidence intervals
r2beta(mod_random, method = "nsj", data = data_sample_3_corrected)
# marginal and conditional R squared values
r.squaredGLMM(mod_random)


# RQ2 Predicting  ---------------------------------------------------------

#calculating the predicted values
pred_data4 <- predict(mod_random,allow.new.levels = TRUE, data_sample_4_corrected)
pred_data4_mean<- lm(pain ~ 1, data = data_sample_4_corrected)

# calculating the sum of residual squared residuals
RSS_data4 = sum((data_sample_4_corrected[, "pain"] - pred_data4)^2)
list(RSS_data4)
TSS_data4 = sum((data_sample_4_corrected$pain - predict(pred_data4_mean))^2)
list(TSS_data4)
#Variance explained
1-(RSS_data4/TSS_data4)

# RQ 3 Second linear mixed model -------------------------------------------------

#Using cortisol serum as the most 
mod_random2_int = lmer(pain ~ cortisol_serum +(1| hospital), control = lmerControl(optimizer = "Nelder_Mead"), data = data_sample_3_corrected)
mod_random2_slope= lmer(pain ~ cortisol_serum +(cortisol_serum| hospital), control = lmerControl(optimizer = "Nelder_Mead"), data = data_sample_3_corrected)
 
#Building a prediction variable
data_sample_3_corrected_mod = data_sample_3_corrected %>%
  mutate(pred_slope = predict(mod_random2), pred_int =predict(mod_random2_int))
  view(data_sample_3_corrected_mod)
  
#Visualization of the separate hospitals 
data_sample_3_corrected_mod %>%
  ggplot() +
  aes(y = pain, x = cortisol_serum, group = hospital)+
  geom_point(aes(color = hospital), size = 4) +
  geom_line(color='red', aes(y=pred_slope, x=pain))+
  facet_wrap( ~ hospital, ncol = 2)

data_sample_3_corrected_mod %>%
  ggplot() +
  aes(y = pain, x = cortisol_serum, group = hospital)+
  geom_point(aes(color = hospital), size = 4) +
  geom_line(color='red', aes(y=pred_int, x=pain))+
  facet_wrap( ~ hospital, ncol = 2)
