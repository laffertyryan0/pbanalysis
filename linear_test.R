library(devtools)

load_all('C:\Users\laffertyrm\Documents\work\PB\pbanalysis') #directory where package is located

setwd(r"(C:\Users\laffertyrm\Documents\work\PB)")            #directory where data is located
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
             family = "multinomial")

print(out)
