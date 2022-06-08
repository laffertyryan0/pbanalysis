devtools::load_all(path="C:/Users/laffertyrm/Documents/work/PB/pbanalysis")
pb.example("logistic")

#test code for multinomial


setwd(r"(C:\Users\laffertyrm\Documents\work\PB)")
data = read.csv("bmi_cat.csv")
data$race = array("other",c(nrow(data)))
data$race[data$deltaR0==1] = "white"
data$race[data$deltaR1==1] = "black"

out = pb.fit(bmi_cat ~ age + age_square + pir + insurance + phy.act + alc.consump + smoke2 + smoke3,
             data = data,
             weights = data$sample_weight,
             disparity.group = "race",
             majority.group = "white",
             minority.group = c("black","other"),
             prop.odds.fail = NULL, #c("phy.act","alc.consump","smoke2","smoke3"),
             family = "ordinal")

print(out)


#test code for linear


setwd(r"(C:\Users\laffertyrm\Documents\work\PB)")
data = read.csv("bmi.csv")
data$race = array("other",c(nrow(data)))
data$race[data$deltaR0==1] = "white"
data$race[data$deltaR1==1] = "black"

out = pb.fit(bmi ~ age + age_square + pir + insurance + phy.act + alc.consump + smoke2 + smoke3,
             data = data,
             weights = data$sample_weight,
             disparity.group = "race",
             majority.group = "white",
             minority.group = c("black","other"),
             prop.odds.fail = NULL, #c("phy.act","alc.consump","smoke2","smoke3"),
             family = "gaussian")

print(out)


#test code for binary logistic


setwd(r"(C:\Users\laffertyrm\Documents\work\PB)")
data = read.csv("bmi_cat.csv")
data$race = array("other",c(nrow(data)))
data$race[data$deltaR0==1] = "white"
data$race[data$deltaR1==1] = "black"
data$bmi_cat[data$bmi_cat == 3] = 2

out = pb.fit(bmi_cat ~ age + age_square + pir + insurance + phy.act + alc.consump + smoke2 + smoke3,
             data = data,
             weights = data$sample_weight,
             disparity.group = "race",
             majority.group = "white",
             minority.group = c("black","other"),
             prop.odds.fail = NULL, #c("phy.act","alc.consump","smoke2","smoke3"),
             family = "multinomial")

print(out)
