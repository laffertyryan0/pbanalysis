#test code for ordinal


data = data.frame(bmi_cat = NHANES::NHANES$BMI_WHO,
             race = NHANES::NHANES$Race1,
             age = NHANES::NHANES$Age,
             poverty = NHANES::NHANES$Poverty
             )
data = na.omit(data)

out = pb.fit(bmi_cat ~ age + poverty,
             data = data,
             disparity.group = "race",
             majority.group = "White",
             minority.group = c("Black"),
             prop.odds.fail = c("age"),
             family = "ordinal")

print(out)

