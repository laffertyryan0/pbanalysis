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


#test code for binary outcome


data = data.frame(bmi_bin = as.integer(NHANES::NHANES$BMI > 20) ,
                  race = NHANES::NHANES$Race1,
                  age = NHANES::NHANES$Age,
                  poverty = NHANES::NHANES$Poverty
)
data = na.omit(data)

out = pb.fit(bmi_bin ~ age + poverty,
             data = data,
             disparity.group = "race",
             majority.group = "White",
             minority.group = c("Black"),
             family = "multinomial")

print(out)

#test code for continuous outcome


data = data.frame(bmi_con = NHANES::NHANES$BMI ,
                  race = NHANES::NHANES$Race1,
                  age = NHANES::NHANES$Age,
                  poverty = NHANES::NHANES$Poverty
)
data = na.omit(data)

out = pb.fit(bmi_con ~ age + poverty,
             data = data,
             disparity.group = "race",
             majority.group = "White",
             minority.group = c("Black"),
             family = "gaussian")

print(out)

