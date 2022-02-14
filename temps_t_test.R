library(dplyr)

temp.liz_F <- na.omit(diagnosed.forest$TEMP.LIZ)
temp.sub_F <- na.omit(diagnosed.forest$TEMP.SUB)
temp.liz_U <- na.omit(diagnosed.urban$TEMP.LIZ)
temp.sub_U <- na.omit(diagnosed.urban$TEMP.SUB)

t.test(temp.liz_F, temp.liz_U)
t.test(temp.sub_F, temp.sub_U)
