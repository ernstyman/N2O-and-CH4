#Loading Packages
library(ggplot2)
library(plyr)
library(dplyr)
library(GGally)
library(broom)
library(ggpubr)
library(gridExtra)

library(srvyr)

#Loading Datasets
Data0 <- read.csv("Physicalchem.csv")

Data0 <- within(Data0, {Soiltype <- factor(Soiltype,
                                           levels = c("Acidic",
                                                      "Neutral",
                                                      "Alkaline"))
Treat <- factor(Treat,
                levels = c("Control",
                           "Biochar",
                           "Straw"))})

str(Data0)

ggplot(Data0, aes(x= PH)) +
  geom_histogram(binwidth = 1.2,
                 color = "black",
                 fill = "white")

ggplot(Data0, aes(x= WFPS)) +
  geom_histogram(binwidth = 5,
                 color = "black",
                 fill = "white")

ggplot(Data0, aes(x= NH4)) +
  geom_histogram(binwidth = 1.5,
                 color = "black",
                 fill = "white")

ggplot(Data0, aes(x= NO3)) +
  geom_histogram(binwidth = 5,
                 color = "black",
                 fill = "white")

ggplot(Data0, aes(x= DOC)) +
  geom_histogram(binwidth = 10,
                 color = "black",
                 fill = "white")

Aggreg0 <- Data0 %>%
  group_by(Soiltype, Treat) %>%
  summarize (meanPH = mean(PH),
             SDPH = sd(PH),
             n = n(),
             SEPH = SDPH / sqrt(n),
             meanWFPS = mean(WFPS),
             SDWFPS = sd(WFPS),
             SEWFPS = SDWFPS / sqrt(n),
             meanNH4 = mean(NH4),
             SDNH4 = sd(NH4),
             SENH4 = SDNH4 / sqrt(n),
             meanNO3 = mean(NO3),
             SDNO3 = sd(NO3),
             SENO3 = SDNO3 / sqrt(n),
             meanDOC = mean(DOC),
             SDDOC = sd(DOC),
             SEDOC = SDDOC / sqrt(n))

Aggreg0 %>% print(n = Inf)

print(select(Aggreg0, c(Treat, meanDOC, SDDOC, SEDOC)), n = Inf)

#Analysis of Variance
SoilpH <- aov(PH ~ Soiltype * Treat, Data0)
summary(SoilpH)
plot(SoilpH, 1)
TukeyHSD(SoilpH)



SoilWFPS <- aov(WFPS ~ Soiltype * Treat, Data0)
summary(SoilWFPS)
plot(SoilWFPS, 1)
TukeyHSD(SoilWFPS)

SoilNH4 <- aov(NH4 ~ Soiltype * Treat, Data0)
summary(SoilNH4)
plot(SoilNH4, 1)
TukeyHSD(SoilNH4)

SoilNO3 <- aov(NO3 ~ Soiltype * Treat, Data0)
summary(SoilNO3)
plot(SoilNO3, 1)
TukeyHSD(SoilNO3)

SoilDOC <- aov(DOC ~ Soiltype * Treat, Data0)
summary(SoilDOC)
plot(SoilDOC, 1)
TukeyHSD(SoilDOC)
ggplot(augment(SoilNO3), aes(.fitted, .resid)) +
  geom_point() +
  geom_smooth(se = F)

Data1 <- read.csv("cumufs.csv")


#Data wrangling
Data1 <- within(Data1, {Soiltype <- factor(Soiltype,
                                           levels = c("Acidic",
                                                      "Neutral",
                                                      "Alkaline"))
                        Treat <- factor(Treat,
                                        levels = c("Control",
                                                   "Biochar",
                                                   "Straw"))})
#Datasets structure
str(Data1)

#Exploratory analysis
#Density plot of N02 by Soil type
ggplot(Data1, aes(x = Nitros, fill = Soiltype)) +
  geom_density(alpha = 0.5)

#Density plot of NO2 by Soil type and Treatment
ggplot(Data1, aes(x = Nitros, fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Soiltype)


#Density plot of Methane by Soil type
ggplot(Data1, aes(x = Methan, fill = Soiltype)) +
  geom_density(alpha = 0.5)

#Density plot of Methane by Soil type and Treatment
ggplot(Data1, aes(x = Methan, fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Soiltype)

#Aggregating data
Aggreg <- Data1 %>%
  group_by(Soiltype, Treat) %>%
  summarize (meanNitros = mean(Nitros),
             SDNitros = sd(Nitros),
             n = n(),
             SENitros = SDNitros / sqrt(n),
             meanMethan = mean(Methan),
             SDMethan = sd(Methan),
             SEMethan = SDMethan / sqrt(n))
             
Aggreg
print(select(Aggreg, -c(meanNitros, SDNitros, SENitros)), n = Inf)

#Anova Nitrous oxide
anova1 <- aov(Nitros ~ Soiltype * Treat, Data1)

#Model anova1 summary
summary(anova1)

#Posthoc of anova1 model by Tukey test
TukeyHSD(anova1)

#Model anova1 diagnostic
plot(anova1)

#Anova Methane
anova3 <- aov(Methan ~ Soiltype * Treat, Data1)

#Model anov3 summary
summary(anova3)

#Model anova3 diagnostic
plot(anova3)

#Plotting Marginal means for Nitrous oxide
p1 <- ggplot(Aggreg, aes(x = Soiltype, y = meanNitros, group = Treat)) +
  geom_point(aes(shape = Treat, color = Treat), size = 3, 
             position = position_dodge(width = 0.3)) +
  stat_summary(geom = "bar", position = "dodge") +
  geom_text(aes(label = c("a", "b", "ab", "c", "c", "c", "c", "c", "c")), position = position_dodge(width = 0.3)) +
  geom_errorbar(aes(ymin = meanNitros - SENitros,
                    ymax = meanNitros + SENitros, color = Treat),
                width = 0.1, 
               position = position_dodge(width = 0.3)) +
  ylim(c(0, 15)) +
  labs(x = "Soil type", 
       y = "Cumulative emissions of Nitrous oxide (Kg N m-2 ha-1)") #+
  
  #scale_color_discrete(labels = c("With biochar", "Without biochar")) +
  #theme_classic() +
  #guides(color = guide_legend(title = NULL)) +
  #theme(legend.position = c(0.9, 0.9)) +
  #theme(legend.background = element_blank()) +
  #theme(legend.key = element_blank())
 p1  

#Plotting Marginal means for Methane
p2 <- ggplot(Aggreg, aes(x = Soiltype, y = meanMethan, 
                   color = Treat)) +
  geom_point(aes(shape = Treat), size = 3, 
             position = position_dodge(width = 0.3)) +
  stat_summary(geom = "bar", position = "dodge") +
  geom_errorbar(aes(ymin = meanMethan - SEMethan, 
                    ymax = meanMethan + SEMethan), width = 0.15,
                position = position_dodge(width = 0.3)) +
  ylim(c(0, 0.15)) +
  labs(x = "Soil type", 
       y = "Cumulative emissions of Nitrous oxide (ug N m-2 h-1)") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(legend.title = element_blank(),
        legend.position = "bottom") +
  annotate("segment", x = 1, xend = 3, y = 0.12, yend = 0.12) + #, arrow = arrow(ends = "both", angle = 30, length = unit(.2, "cm"))) +
  annotate("segment", x = 1, xend = 1, y = 0.116, yend = 0.12) +
  annotate("segment", x = 3, xend = 3, y = 0.116, yend = 0.12) +
  annotate("text", x = 2, y = 0.125, label = "p < 0.05")
  
p2
#scale_color_discrete(labels = c("With biochar", "Without biochar")) +
  #theme_classic() +
  #guides(color = guide_legend(title = NULL)) +
  #theme(legend.position = c(0.9, 0.9)) +
  #theme(legend.background = element_blank()) +
  #theme(legend.key = element_blank())
  
figuremeans <- ggarrange(p1, p2, common.legend = T,
                         nrow = 1, legend = "none", align = "hv", labels = c("a", "b"), font.label = list(size = 12), heights = 1:1)

Marginplot <- annotate_figure(figuremeans,
                              bottom = text_grob("Typical soil pH status", size = 10, vjust = -4, hjust = -0.05))

Marginplot %>%
  gridExtra::grid.arrange(get_legend(p2), 
                          heights = unit(c(80, 5), "mm"))


#Data_frame2
Data2 <- read.csv("Gas_data.csv")
Data2 <- within(Data2, {Soiltype <- factor(Soiltype,
                                           levels = c("Acidic", 
                                                      "Neutral", 
                                                      "Alkaline"))
                        Treat <- factor(Treat,
                                        levels = c("Control",
                                                   "Biochar",
                                                   "Straw"))
                        Time <- as.Date.factor(Time, format = "%m/%d/%y")})
str(Data2)

Aggreg2 <- Data2 %>%
  group_by(Time, Soiltype, Treat)%>%
  summarize(meanNitros = mean(Nitros, na.rm = T),
            sdNitros = sd(Nitros, na.rm = T),
            n = n(),
            SENitros = sdNitros/sqrt(n),
            meanMethan = mean(Methan, na.rm = T),
            sdMethan = sd(Methan, na.rm = T),
            SEMethan = sdMethan/sqrt(n))

AcidsoilN<- Aggreg2 %>% filter(Soiltype == "Acidic")
NeutralsoilN <- Aggreg2 %>% filter(Soiltype == "Neutral")
AlkalinesoilN <- Aggreg2 %>% filter(Soiltype == "Alkaline")


#Plotting Nitrous Oxide Dynamics
ggplot(Aggreg2, aes(Time, meanNitros, 
                     group = interaction(Soiltype, Treat))) +
  geom_errorbar(aes(ymin = meanNitros - SENitros, 
                    ymax = meanNitros + SENitros, color = Soiltype, linetype = Treat), width = 0.2) +
  geom_line(aes(color = Soiltype, linetype = Treat), size = 0.6) +
  geom_point(aes(shape = Treat, color = Soiltype)) +
  theme(axis.text.x = element_text(angle = 45)) #+
  #scale_x_date()

N <- Acidsoil %>% ggplot (aes(Time, meanNitros,
                         color = Treat, group = Treat)) +
  geom_point() +
  geom_line(aes(color = Treat)) +
  geom_errorbar(aes(ymin = meanNitros - SENitros,
                    ymax = meanNitros + SENitros, color = Treat), width = 0.2) +
  theme(aspect.ratio = 1/1) +
  ylab(parse(text = "CH[4]* ' '*emissions*' '*from*' '*acidic*' '*soil* ' '*(mg*' '*C*' '*m^-2* ' '*hr^-1)")) +
  theme_classic()
N

NitrosDynamics <- 
  ggplot(Aggreg2, aes(Time, meanNitros, color = Treat, group = Treat)) +
  geom_point() +
  geom_line(aes(color = Treat)) +
  geom_errorbar(aes(ymin = meanNitros - SENitros, 
                    ymax = meanNitros + SENitros, color = Treat), width = 0.2) +
  facet_wrap(~Soiltype, nrow = 1, scales = "free_y") +
  theme_classic() +
  scale_x_date() +
  ylab(parse(text = "N[2]*O* ' '*emissions* ' '*(Kg*' '*N*' '*m^-2* ' '*ha^-1)")) +
  theme(aspect.ratio = 0.3/0.5) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.title = element_text(size = rel(0.8)))
  #theme(plot.margin = unit(c(1, 1, -0.5, 1),"lines"))
  #theme(axis.text.x = element_text(angle = 45))
NitrosDynamics

M <- Acidsoil %>% ggplot (aes(Time, meanMethan,
                         color = Treat, group = Treat)) +
  geom_point() +
  geom_line(aes(color = Treat)) +
  geom_errorbar(aes(ymin = meanMethan - SEMethan,
                    ymax = meanMethan + SEMethan, color = Treat), width = 0.2) +
  theme(aspect.ratio = 1/1) +
  theme_classic()

ggarrange(N, M, common.legend = T,
          nrow = 1, legend = "bottom", align = "hv",
          heights = 1:1)

#Plotting Methane Dynamics 
MethanDynamics <- 
  ggplot(Aggreg2, aes(Time, meanMethan, color = Treat, group = Treat)) +
  geom_point() +
  geom_line(aes(color = Treat)) +
  geom_errorbar(aes(ymin = meanMethan - SEMethan, 
                    ymax = meanMethan + SEMethan, color = Treat), width = 0.2) +
  facet_wrap(~Soiltype, nrow = 1, scales = "free_y") +
  theme_classic() +
  scale_x_date() +
  ylab(parse(text = "CH[4]* ' '*emissions* ' '*(Kg*' '*C*' '*m^-2* ' '*ha^-1)")) +
  theme(aspect.ratio = 0.3/0.5) +
  theme(strip.background = element_blank(),
        strip.text = element_blank()) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank()
        ) +
  theme(legend.title = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.title = element_text(size = rel(0.8)))
  #theme(plot.margin = unit(c(-0.5, 1, 1, 1),"lines"))
#theme(axis.text.x = element_text(angle = 45))

MethanDynamics

figuredy <- ggarrange(MethanDynamics, NitrosDynamics, common.legend = T,
          nrow = 2, legend = "bottom", align = "hv", labels = c("a", "b"), hjust = c(-4, -3.5), vjust = -0.01, font.label = list(size = 13))

annotate_figure(figuredy,
                top = text_grob(paste("Acidic Soil                                   ", 
                                      "Neutral Soil                                   ", 
                                      "Alkaline Soil",
                                       
                                      sep = "                    "), hjust = 0.47))


#Exploratory data analysis
#Data frame 3

Data3 <- read.csv("Soildata_ModelingNitrous.csv")
Data3 <- within(Data3, {Soiltype <- factor(Soiltype,
                                           labels = c("Acidic", "Neutral", "Alkaline"))

Treat <- factor (Treat,
                 levels = c("Control", "Biochar", "Straw"))})
head(Data3)
str(Data3)



#Histogram for Nitrous oxide
ggplot(Data3, aes(Nitros)) +
  geom_histogram(binwidth = 500, 
                 color = "black", 
                 fill = "white")
#Nitrous oxide is right skewed

#Log transform for Nitrous oxide
ggplot(Data3, aes(log(Nitros))) +
  geom_histogram(binwidth = 1, 
                 color = "black", 
                 fill = "white")
#Log transforming Nitrous oxide can make remove the skewness
#and rather produce a more symmetrical shape

#Density plot for log(N2O) by Treatments
ggplot(Data3, aes(log(Nitros), fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Soiltype, scales = "free")
#It turns out that Nitrous oxide emission did not differ did not differ 
#by Biochar application

#Density plot for log(N2O) by Soil types
ggplot(Data3, aes(log(Nitros_2), fill = Soiltype)) +
  geom_density(alpha = 0.5)
#It turns out that acidic soil contribute more to Nitrous oxide emission

#Correlation among soil factors
ggpairs(Data3, columns = 6:12)

#Histogram pH

ggplot(Data3, aes(SoilpH)) +
  geom_histogram(binwidth = 0.5, 
                 color = "black", 
                 fill = "white")
#Soil pH is explained by a binormal distribution

#Density plot for Soil pH by treatment
ggplot(Data3, aes(SoilpH, fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Soiltype,  scales = "free")
#Treatments can not explain bimodality

#Density plot for soil pH by Soil types
ggplot(Data3, aes(SoilpH, fill = Soiltype)) +
  geom_density(alpha = 0.5)
#It turns out that the bimodal shape is explained by different types of Soil 


Regression_labelspH <- function(dat) {
  model <- lm(log(Nitros) ~ SoilpH, data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(dat$SoilpH, log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
      notation <- print("**")
  }else if (p_value <= 0.05) {
        notation <- print("*")
  } else {
        notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == %.2f", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}
 
labelspH <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelspH)
labelspH


#Scatter plot for log(N2O) against soil pH
PlotpHNitros <- 
  ggplot(Data3, aes(x = SoilpH, y = log(Nitros), color = Treat)) +
  geom_point(aes(shape = Treat), size = 1) +
  geom_text(x = c(6, 7.5, 8), y = c(8.2, 7.8, 10.6), aes(label = c(formula)), data = filter(labelspH, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8), y = c(7.8, 7.1, 9.9), aes(label = c(formula)), data = filter(labelspH, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8), y = c(7.4, 6.4, 9.2), aes(label = c(formula)), data = filter(labelspH, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8), y = c(8, 7.5, 10.3), aes(label = c(r2)), data = filter(labelspH, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8), y = c(7.6, 6.8, 9.6), aes(label = c(r2)), data = filter(labelspH, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8), y = c(7.2, 6.1, 8.9), aes(label = c(r2)), data = filter(labelspH, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(6, 7.5, 8.13), y = c(7.6, 6.8, 9.6), aes(label = c(p)), data = filter(labelspH, Treat == "Biochar"), size = 2, show.legend = F) +
  geom_text(x = c(6, 7.5, 8.13), y = c(7.2, 6.1, 8.9), aes(label = c(p)), data = filter(labelspH, Treat == "Straw"), size = 2, show.legend = F) +
  geom_smooth(method = "lm", se = F, alpha = 0.25) +
  theme_classic() +
  scale_color_manual(values = c("darkred", "darkgreen", "darkblue")) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  theme(aspect.ratio = 0.01/0.01) +
  xlab("Soil pH") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8))) +
  theme(legend.title = element_blank()
        )
PlotpHNitros


#Explaining Nitrous oxide emission by soil pH
MpH <- lm(log(Nitros) ~ SoilpH * Soiltype * Treat, Data3)
summary(MpH)  
plot(MpH)

library(kableExtra)
 library(data.table)
  

tM <- tidy(MpH, conf.int = T) #%>%
  kable(col.names = c("term",  "estimate", "std.error", "statistic", "p.value", "conf.low", "conf.high"))
dt <- setDT(tM)
dt %>% kable()
augpH <- augment(MpH)
head(augpH)
residupH <- resid(MpH)
head(residupH)
dfpH <- data.frame(augpH, residupH)
head(dfpH)

ggplot(dfpH, aes(.fitted, residupH)) +
  geom_point() +
  #geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F)

#Histogram for WFPS 
ggplot(Data3, aes(WFPS)) +
  geom_histogram(binwidth = 3, 
                 color = "black", 
                 fill = "white") #+
#facet_wrap(Soil_type ~ Treatment)

#Density plot for WFPS by treatment
ggplot(Data3, aes(WFPS, fill = Treat)) +
  geom_density(alpha = 0.5)
#Biochar application was likely to increase water content

#Density plot for WFPS by Soil types
ggplot(Data3, aes(WFPS, fill = Soiltype)) +
  geom_density(alpha = 0.5)
#Acidic soil retained water more than Alkaline and Neutral did

Regression_labelsWFPS <- function(dat) {
  model <- lm(log(Nitros) ~ WFPS, data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(dat$WFPS, log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
    notation <- print("**")
  }else if (p_value <= 0.05) {
    notation <- print("*")
  } else {
    notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == %.2f", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}

labelsWFPS <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelsWFPS)
labelsWFPS

#Scatterplot for log(N2O) against WFPS
plotWFPSNitros <- 
  ggplot(Data3, aes(x = WFPS, y = log(Nitros), color = Treat)) +
  geom_point(aes(shape = Treat), size = 1) +
  geom_text(x = c(65, 42, 55), y = c(5.8, 3.5, 9.5), aes(label = c(formula)), data = filter(labelsWFPS, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 55), y = c(5.4, 2.8, 8.8), aes(label = c(formula)), data = filter(labelsWFPS, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 55), y = c(5.0, 2.1, 8.1), aes(label = c(formula)), data = filter(labelsWFPS, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 55), y = c(5.6, 3.2, 9.2), aes(label = c(r2)), data = filter(labelsWFPS, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 55), y = c(5.2, 2.5, 8.5), aes(label = c(r2)), data = filter(labelsWFPS, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 55), y = c(4.8, 1.8, 7.8), aes(label = c(r2)), data = filter(labelsWFPS, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(65, 42, 56), y = c(5.2, 2.5, 8.5), aes(label = c(p)), data = filter(labelsWFPS, Treat == "Biochar"), size = 2, show.legend = F) +
  geom_text(x = c(65, 42, 56), y = c(4.8, 1.8, 7.8), aes(label = c(p)), data = filter(labelsWFPS, Treat == "Straw"), size = 2, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme(aspect.ratio = 1/1) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("WFPS") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))


MWFPS <- lm(log(Nitros) ~ WFPS * Soiltype * Treat, Data3)
summary(MWFPS)
plot(MWFPS)

#aug <- 
  augment(MWFPS)
head(aug)
residuwt <- resid(MWFPS)
head(residuwt)
df <- data.frame(aug, residuwt)
head(df)

ggplot(df, aes(.fitted, residuwt)) +
  geom_point() +
  geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F)


#Histogram for Ammonium
ggplot(Data3, aes(x = Ammo)) +
  geom_histogram(binwidth = 10,
                 color = "black",
                 fill = "white")
#The distribution of Ammonium is right skewed

#Log transformation of Ammonium
ggplot(Data3, aes(x = log(Ammo))) +
  geom_histogram(binwidth = 0.5,
                 color = "black",
                 fill = "white")
#The shape is binormal distributed

#Density plot for Ammonium by Treat
ggplot(Data3, aes(x = log(Ammo), fill = Treat)) +
  geom_density(alpha = 0.5)

#Density plot for Ammonium by Soil type
ggplot(Data3, aes(x = log(Ammo), fill = Soiltype)) +
  geom_density(alpha = 0.5)
#Acidity soil has higher levels of Ammonium


Regression_labelsAmmo <- function(dat) {
  model <- lm(log(Nitros) ~ log(Ammo), data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(log(dat$Ammo), log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
    notation <- print("**")
  }else if (p_value <= 0.05) {
    notation <- print("*")
  } else {
    notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == %.2f", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}

labelsAmmo <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelsAmmo)
labelsAmmo

#Scatterplot for Nitrous oxide against Ammonium
PlotAmmoNitros <- 
  ggplot(Data3, aes(x = log(Ammo), y = log(Nitros),
                                    color = Treat)) +
  geom_point(aes(shape = Treat), size = 1) +
  geom_text(x = c(0, 3.0, 3.0), y = c(8.0, 2.8, 2.8), aes(label = c(formula)), data = filter(labelsAmmo, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.0), y = c(7.6, 2.1, 2.1), aes(label = c(formula)), data = filter(labelsAmmo, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.0), y = c(7.2, 1.4, 1.4), aes(label = c(formula)), data = filter(labelsAmmo, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.0), y = c(7.8, 2.5, 2.5), aes(label = c(r2)), data = filter(labelsAmmo, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.0), y = c(7.4, 1.8, 1.8), aes(label = c(r2)), data = filter(labelsAmmo, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.0), y = c(7.0, 1.1, 1.1), aes(label = c(r2)), data = filter(labelsAmmo, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.13), y = c(7.8, 2.5, 2.5), aes(label = c(p)), data = filter(labelsAmmo, Treat == "Control"), size = 2, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.13), y = c(7.4, 1.8, 1.8), aes(label = c(p)), data = filter(labelsAmmo, Treat == "Biochar"), size = 2, show.legend = F) +
  geom_text(x = c(0, 3.0, 3.13), y = c(7.0, 1.1, 1.1), aes(label = c(p)), data = filter(labelsAmmo, Treat == "Straw"), size = 2, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme(aspect.ratio = 0.01/0.01) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Ammonium") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))


MAmm <- lm(log(Nitros) ~ log(Ammo) * Soiltype * Treat, Data3)
summary(MAmm)
plot(MAmm)

augment(MAmm)
#Histogram for Nitrate
ggplot(Data3, aes(x = Nitrate)) +
  geom_histogram(binwidth = 4,
                 color = "black",
                 fill = "white")

#Log transformation of Nitrate
ggplot(Data3, aes(x = log(Nitrate))) +
  geom_histogram(binwidth = 0.5,
                 color = "black",
                 fill = "white")
#Density plot for Nitrate by Treatment
ggplot(Data3, aes(x = Nitrate, fill = Treat)) +
  geom_density(alpha = 0.5)

#Density plot for Nitrate by Soiltype
ggplot(Data3, aes(x = log(Nitrate), fill = Soiltype)) +
  geom_density(alpha = 0.5)

Regression_labelsNitrate <- function(dat) {
  model <- lm(log(Nitros) ~ log(Nitrate), data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(log(dat$Nitrate), log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
    notation <- print("**")
  }else if (p_value <= 0.05) {
    notation <- print("*")
  } else {
    notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == %.2f", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}

labelsNitrate <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelsNitrate)
labelsNitrate

#Scatterplot for Nitrous oxide against Nitrate
PlotNitrateNitros <- 
  ggplot(Data3, aes(x = log(Nitrate), y = log(Nitros),
                                       color = Treat)) +
  geom_point(aes(shape = Treat), size = 1) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(5.5, 2.8, 2.8), aes(label = c(formula)), data = filter(labelsNitrate, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(5.1, 2.1, 2.1), aes(label = c(formula)), data = filter(labelsNitrate, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(4.7, 1.4, 1.4), aes(label = c(formula)), data = filter(labelsNitrate, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(5.3, 2.5, 2.5), aes(label = c(r2)), data = filter(labelsNitrate, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(4.9, 1.8, 1.8), aes(label = c(r2)), data = filter(labelsNitrate, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.0), y = c(4.5, 1.1, 1.1), aes(label = c(r2)), data = filter(labelsNitrate, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.13), y = c(5.3, 2.5, 2.5), aes(label = c(p)), data = filter(labelsNitrate, Treat == "Control"), size = 2, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.13), y = c(4.9, 1.8, 1.8), aes(label = c(p)), data = filter(labelsNitrate, Treat == "Biochar"), size = 2, show.legend = F) +
  geom_text(x = c(-2, 3.0, 3.13), y = c(4.7, 1.1, 1.1), aes(label = c(p)), data = filter(labelsNitrate, Treat == "Straw"), size = 2, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme(aspect.ratio = 0.01/0.01) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Nitrate") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))

MNitrate <- lm(log(Nitros) ~ log(Nitrate) * Soiltype * Treat, Data3)
summary(MNitrate)
plot(MNitrate)

augment(MNitrate)
#Histogram or Nitrite
ggplot(Data3, aes(x = DOC)) +
  geom_histogram(binwidth = 4,
                 color = "black",
                 fill = "white")

#log transformation for Nitrite
ggplot(Data3, aes(x = log(Nitrite))) +
  geom_histogram(binwidth = 0.5,
                 color = "black",
                 fill = "white")

#Density plot for Nitrite by Treatment
ggplot(Data3, aes(x = DOC, fill = Treat)) +
  geom_density(alpha = 0.5)

#Density plot for Nitrite by Soil type
ggplot(Data3, aes(x = DOC, fill = Soiltype)) +
  geom_density(alpha = 0.5)

Regression_labelsDOC <- function(dat) {
  model <- lm(log(Nitros) ~ DOC, data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(dat$DOC, log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
    notation <- print("**")
  }else if (p_value <= 0.05) {
    notation <- print("*")
  } else {
    notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == '%.2f'", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}

labelsDOC <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelsDOC)
labelsDOC

#Scatterplot for Nitrous oxide against Nitrite
#PlotDOCNitros <- 
  ggplot(Data3, aes(x = DOC, y = log(Nitros),
                                       color = Treat)) +
  geom_point(aes(shape = Treat)) +
  #geom_text(aes(label = Number), vjust = -2, size = 2) +
  geom_text(x = c(-5, 0, -0.5), y = c(5.5, 2.8, 8.5), aes(label = c(formula)), data = filter(labelsDOC, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(5.1, 2.1, 7.8), aes(label = c(formula)), data = filter(labelsDOC, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(4.7, 1.4, 7.1), aes(label = c(formula)), data = filter(labelsDOC, Treat == "Straw"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(5.3, 2.5, 8.2), aes(label = c(r2)), data = filter(labelsDOC, Treat == "Control"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(4.8, 1.8, 7.5), aes(label = c(r2)), data = filter(labelsDOC, Treat == "Biochar"), parse = T, size = 1, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(4.5, 1.1, 6.8), aes(label = c(r2)), data = filter(labelsDOC, Treat == "Straw"), parse = T, size = 2, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(5.3, 2.5, 8.2), aes(label = c(p)), data = filter(labelsDOC, Treat == "Control"), size = 2, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(4.8, 1.8, 7.5), aes(label = c(p)), data = filter(labelsDOC, Treat == "Biochar"), size = 2, show.legend = F) +
  geom_text(x = c(-5, 0, -0.5), y = c(4.5, 1.1, 6.8), aes(label = c(p)), data = filter(labelsDOC, Treat == "Straw"), size = 2, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme(aspect.ratio = 1/1) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("DOC") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))

PlotDOCNitros

MDOC <- lm(log(Nitros) ~ DOC * Soiltype * Treat, Data3)
summary(MDOC)
plot(MDOC)

aug <- augment(MDOC)
head(aug)
residu <- resid(MDOC)
head(residu)
df <- data.frame(aug, residu)
head(df)

ggplot(df, aes(.fitted, residu, color = Treat)) +
  geom_point(aes(shape = Treat)) +
  geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F) #+
  facet_wrap(.~Soiltype, scales = "free", nrow = 2)

#Histogram for Temperature
ggplot(Data3, aes(x = Temp)) +
  geom_histogram(binwidth = 1.5,
                 color = "black",
                 fill = "white")

#Density plot for Temperature by Treatment
ggplot(Data3, aes(x = Temp, fill = Treat)) +
  geom_density(alpha = 0.5)

#Density plot for Temperature by Soil type
ggplot(Data3, aes(x = Temp, fill = Soiltype)) +
  geom_density(alpha = 0.5)

Regression_labelsTemp <- function(dat) {
  model <- lm(log(Nitros) ~ Temp, data = dat)
  
  formula <- sprintf("italic(y) == %.2f %+.2f * italic(x)",
                     round(coef(model)[1], 2), round(coef(model)[2], 2))
  
  cor <- cor.test(log(dat$Temp), log(dat$Nitros))
  r <- as.numeric(cor[4])
  p_value <- as.numeric(cor[3])
  if(p_value <= 0.001) {
    
    notation <- print("***")
  } else if (p_value <= 0.01) {
    notation <- print("**")
  }else if (p_value <= 0.05) {
    notation <- print("*")
  } else {
    notation <- print("")
  }
  r2 <- sprintf("italic(R^2) == %.2f", r^2)
  p <- sprintf("%.3s", notation)
  data.frame(formula = formula, r2 = r2, p = p, stringsAsFactors = FALSE)
}

labelsTemp <- ddply(Data3, c("Soiltype", "Treat"), Regression_labelsTemp)
labelsTemp

#Scatterplot for Nitrous oxide against temperature
PlotTempNitros <- 
  ggplot(Data3, aes(x = Temp, y = log(Nitros),
                                    color = Treat)) +
  geom_point(aes(shape = Treat), alpha = 1) +
    geom_text(x = c(28, 32, 32), y = c(5.3, 2.8, 2.8), aes(label = c(formula)), data = filter(labelsTemp, Treat == "Control"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(28, 32, 32), y = c(4.9, 2.1, 2.1), aes(label = c(formula)), data = filter(labelsTemp, Treat == "Biochar"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(28, 32, 32), y = c(4.5, 1.4, 1.4), aes(label = c(formula)), data = filter(labelsTemp, Treat == "Straw"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(28, 32, 32), y = c(5.1, 2.5, 2.5), aes(label = c(r2)), data = filter(labelsTemp, Treat == "Control"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(28, 32, 32), y = c(4.7, 1.8, 1.8), aes(label = c(r2)), data = filter(labelsTemp, Treat == "Biochar"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(28, 32, 32), y = c(4.3, 1.1, 1.1), aes(label = c(r2)), data = filter(labelsTemp, Treat == "Straw"), parse = T, size = 2, show.legend = F) +
    geom_text(x = c(25, 32, 32.13), y = c(5.1, 2.5, 2.5), aes(label = c(p)), data = filter(labelsTemp, Treat == "Control"), size = 3, show.legend = F) +
    geom_text(x = c(25, 32, 32.13), y = c(4.7, 1.8, 1.8), aes(label = c(p)), data = filter(labelsTemp, Treat == "Biochar"), size = 3, show.legend = F) +
    geom_text(x = c(25, 32, 32.13), y = c(4.3, 1.1, 1.1), aes(label = c(p)), data = filter(labelsTemp, Treat == "Straw"), size = 3, show.legend = F) +
  geom_smooth(method = "lm", se = F) +
  theme_classic() +
  facet_wrap(~ Soiltype, nrow = 2, scales = "free") +
  #theme(aspect.ratio = 1/1) +
  theme(strip.background = element_blank(),
        legend.position = c(0.6, 0.35))+
  theme(strip.text = element_text()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  #theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Temperature") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))
PlotTempNitros

figure <- ggarrange(plotWFPSNitros, PlotTempNitros,  
                 ncol = 1, common.legend = T, legend = "bottom",
                    labels = c("a", "b"), font.label = list(size = 13), heights = 1:1)

figure

figure <- ggarrange(PlotpHNitros, PlotAmmoNitros,
                    PlotNitrateNitros, PlotNitriteNitros,
                    ncol = 2, nrow = 2, common.legend = T, legend = "bottom",
                    labels = c("a", "b", "c", "d"), font.label = list(size = 13), heights = 1:1)
annotate_figure(figure, left = text_grob(parse(text = "N[2]*O* ' '*emissions* ' '*(Kg*' '*N*' '*m^-2* ' '*ha^-1)"), rot = 90, size = 10, vjust = 2.5),
                top = text_grob(paste("Acidic Soil           ", 
                                      "Neutral Soil           ", 
                                      "Alkaline Soil                                 ", 
 sep = " "), hjust = 0.47))
annotate_figure(figure, left = text_grob(parse(text = "N[2]*O* ' '*emissions* ' '*(Kg*' '*N*' '*m^-2* ' '*ha^-1)"), rot = 90, size = 10, vjust = 2.5),
                top = text_grob(paste("Acidic Soil           ", 
                                      "Neutral Soil           ", 
                                      "Alkaline Soil                                 ", 
                                      "Acidic Soil          ", 
                                      "Neutral Soil          ", 
                                      "Alkaline Soil", sep = " "), hjust = 0.47))


MTemp <- lm(log(Nitros) ~ Temp * Soiltype * Treat, Data3)
summary(MTemp)
plot(MTemp)

ff <- Data3 %>% filter(Soiltype == "Acidic", Treat == "Control")
m1 <- lm(log(Nitros) ~ Temp, ff)
summary(m1)

ee <- Data3 %>% filter(Soiltype == "Alkaline", Treat == "Control")
dim(ee)
m2 <- lm(log(Nitros) ~ Temp, ee)
summary(m2)
anova(m1, m2)

#Exploratory data analysis
#Data_frame4

Data4 <- read.csv("Soildata_ModelingMethane.csv")
Data4 <- within(Data4, {Soiltype <- factor(Soiltype,
                                           labels = c("Acidic", "Neutral", "Alkaline"))
                        
                        Treat <- factor (Treat,
                                         levels = c("Control", "Biochar", "Straw"))})
head(Data4)
str(Data4)



#Histogram for Methane
histmeth <- ggplot(Data4, aes(Methan)) +
  geom_histogram(binwidth = 5, 
                 color = "black", 
                 fill = "white") +
  theme_classic()
  #theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank()) #+
  
histmeth
#Nitrous oxide is right skewed
#Inspecting Outliers
boxmeth <- ggplot(Data4, aes(x = Methan, y= "")) +
  geom_boxplot() +
  theme_classic() +
  theme(axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        axis.text = element_blank())
boxmeth
ggarrange(boxmeth,histmeth, heights = 1:2, nrow = 2, align = "v")

#Which outliers Methan > 25

Data4out <- Data4 %>% filter(Methan > 25)
Data4out

Data4 <- Data4 %>% filter(Number > 9)

#It looks like the first day  


#Density plot for Methane by Treatments & Soil type
ggplot(Data4, aes(Methan, fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Soiltype, scales = "free")
#It turns out that Nitrous oxide emission did not differ did not differ 
#by Biochar application

#Density plot for log(N2O) by Soil types
ggplot(Data3, aes(log(Nitros_2), fill = Soiltype)) +
  geom_density(alpha = 0.5)
#It turns out that acidic soil contribute more to Nitrous oxide emission


#Histogram pH

ggplot(Data3, aes(SoilpH)) +
  geom_histogram(binwidth = 0.5, 
                 color = "black", 
                 fill = "white")
#Soil pH is explained by a binormal distribution

#Density plot for Soil pH by treatment
ggplot(Data3, aes(SoilpH, fill = Treat)) +
  geom_density(alpha = 0.5) +
  facet_wrap(~Soiltype,  scales = "free")
#Treatments can not explain bimodality

#Density plot for soil pH by Soil types
ggplot(Data3, aes(SoilpH, fill = Soiltype)) +
  geom_density(alpha = 0.5)
#It turns out that the bimodal shape is explained by different types of Soil 

#Scatter plot for Methane against soil pH
PlotpHMethan <- 
  ggplot(Data4, aes(x = SoilpH, y = Methan, color = Treat)) +
  geom_point(aes(shape = Treat), size = 2) +
  #geom_text(aes(label = Number), vjust = - 2, size = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") + 
    theme_classic() +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  theme(aspect.ratio = 0.02/0.02) +
  xlab("Soil pH") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8))) +
  theme(legend.title = element_blank())
PlotpHMethan


#Explaining Mewthane emission by soil pH
MpH2 <- lm(Methan ~ SoilpH * Soiltype * Treat, Data4)
summary(MpH2)  
plot(MpH2)


augMpH2 <- augment(MpH2)
head(aug)
residuMpH2 <- resid(MpH2)
head(residuwt)
dfMpH2diag <- data.frame(augMpH2, residuMpH2)
head(df)

ggplot(dfMpH2diag, aes(.fitted, residuMpH2)) +
  geom_point() +
  geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F)

#Scatterplot for Methane against WFPS
plotWFPSMethan <- 
  ggplot(Data4, aes(x = WFPS, y = Methan, color = Treat)) +
  geom_point(aes(shape = Treat), size = 2) +
  #geom_text(aes(label = Number), size = 2, vjust = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme_classic() +
  theme(aspect.ratio = 0.05/0.05) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("WFPS") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))


MWFPS2 <- lm(Methan ~ WFPS * Soiltype * Treat, Data4)
summary(MWFPS2)
plot(MWFPS2)

aug <- augment(MWFPS)
head(aug)
residuwt <- resid(MWFPS)
head(residuwt)
df <- data.frame(aug, residuwt)
head(df)

ggplot(df, aes(.fitted, residuwt)) +
  geom_point() +
  geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F)

#Scatterplot for Methane against Ammonium
PlotAmmoMethan <- 
  ggplot(Data4, aes(x = log(Ammo), y = Methan,
                  color = Treat)) +
  geom_point(aes(shape = Treat), size = 2) +
  #geom_text(aes(label = Number), size = 2, vjust = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme_classic() +
  theme(aspect.ratio = 0.01/0.01) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Ammonium") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))


MAmm2 <- lm(Methan ~ log(Ammo) * Soiltype * Treat, Data4)
summary(MAmm2)
plot(MAmm2)

#Scatterplot for Methane against Nitrate
PlotNitrateMethan <- 
  ggplot(Data4, aes(x = log(Nitrate), y = Methan,
                  color = Treat)) +
  geom_point(aes(shape = Treat), size = 2) +
  #geom_text(aes(label = Number), size = 2, vjust = -2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme_classic() +
  theme(aspect.ratio = 0.01/0.01) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Nitrate") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))
  
MNitrate2 <- lm(Methan ~ log(Nitrate) * Soiltype * Treat, Data4)
summary(MNitrate2)
plot(MNitrate2)

#Scatterplot for Methane against Nitrite
PlotDOCMethan <- 
  ggplot(Data4, aes(x = DOC, y = Methan,
                  color = Treat)) +
  geom_point(aes(shape = Treat)) +
  geom_text(aes(label = Number), size = 2, vjust = -2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme_classic() +
  theme(aspect.ratio = 1/1) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("DOC") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))
  
PlotDOCMethan
MDOC2 <- lm(Methan ~ DOC * Soiltype * Treat, Data4)
summary(MDOC2)
plot(MDOC2)

aug <- augment(MDOC2)
head(aug)
residu <- resid(MDOC2)
head(residu)
df <- data.frame(aug, residu)
head(df)

ggplot(df, aes(.fitted, residu, color = Treat)) +
  geom_point(aes(shape = Treat)) +
  geom_text(aes(label = .rownames), vjust = -2, size = 2) +
  geom_smooth(se = F) +
  #facet_wrap(~ Soiltype, scales = "free") +
  theme(aspect.ratio = 1/1)

#Scatterplot for Methan against temperature
PlotTempMethan <- 
  ggplot(Data4, aes(x = Temp, y = Methan,
                  color = Treat)) +
  geom_point(aes(shape = Treat), size = 2) +
  #geom_text(aes(label = Number), size = 2, vjust = 2) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~ Soiltype, nrow = 1, scales = "free") +
  theme_classic() +
  theme(aspect.ratio = 0.01/0.01) +
  theme(strip.background = element_blank())+
  theme(strip.text = element_blank()) +
  theme(axis.text = element_text(size = rel(0.7))) +
  theme(axis.ticks = element_line(size = rel(0.01))) +
  theme(axis.line = element_line(size = rel(0.2))) +
  xlab("Temperature") +
  ylab("") +
  theme(axis.title = element_text(size = rel(0.8)))
  
  MTemp2 <- lm(Methan ~ Temp * Soiltype * Treat, Data4)
  summary(MTemp2)
  plot(MTemp2)
  
  augMTemp2 <- augment(MTemp2)
  head(augMTemp2)
  residuMTemp2 <- resid(MTemp2)
  head(residu)
  dfTempdiag <- data.frame(augMTemp2, residuMTemp2)
  head(df)
  
  ggplot(dfTempdiag, aes(.fitted, residuMTemp2)) +
    geom_point(aes(shape = Treat)) +
    geom_text(aes(label = .rownames), vjust = -2, size = 2) +
    geom_smooth(se = F)

figure <- ggarrange(PlotpHMethan, plotWFPSMethan, PlotAmmoMethan,
          PlotNitrateMethan, PlotTempMethan,PlotNitriteMethan,  
          nrow = 3, ncol = 2, common.legend = T, legend = "bottom",
          labels = c("a", "b", "c", "d", "e", "f"), vjust = 1, hjust = c(-4.9, -4.7, -4.9, -4.7, -4.9, -10), font.label = list(size = 13), align = "hv")
  
figure2 <- ggarrange(PlotpHMethan, plotWFPSMethan,
                    nrow = 2, ncol = 1, common.legend = T, legend = "bottom",
                    labels = c("a", "b"), vjust = 1, hjust = c(-4.9, -4.7), font.label = list(size = 13), align = "hv")

save <- annotate_figure(figure2, left = text_grob(parse(text = "CH[4]* ' '*emissions* ' '*(Kg*' '*C*' '*m^-2* ' '*ha^-1)"), rot = 90, size = 10, vjust = 2.5),
                top = text_grob(paste("Acidic Soil           ", 
                                      "Neutral Soil           ", 
                                      "Alkaline Soil                                 ", 
                                      "Acidic Soil          ", 
                                      "Neutral Soil          ", 
                                      "Alkaline Soil", sep = " "), hjust = 0.47))
save
ggsave("save.png", width = 30, height = 30, units = "cm", dpi = 1600)


Nit <- Data3 %>% filter(Soiltype == "Neutral", Treat == "Control")
summary(lm(log(Nitros) ~ SoilpH, Nit))

datann <- read.csv("annual fluxes.csv")
colnames(datann) <- c("Date", "Maximum temperature", 
                      "Minimum temperature", "Precipitation")

datann$Date <- as.Date(datann$Date, format = "%m/%d/%Y")
str(datann)
summary(datann)


library(lubridate)
summary(datann)
g <- ggplot(datann) +
  geom_bar(aes(x = Date, y = Precipitation / 1.5 - 6), stat = "identity") +
  geom_line(aes(x = Date, y = `Maximum temperature`,linetype = "Maximum temperature"), color = "red") +
  geom_line(aes(x = Date, y = `Minimum temperature`, linetype = "Minimum temperature"), color = "red", na.rm = T) +
  scale_x_date(name = "Month", breaks = seq.Date(as.Date("2020-01-01"), as.Date("2021-01-01"), by = "1 month"),
               labels = function(date){return(month(date, label = T))}) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(~ (. + 6) * 1.5, name = "Precipitation"), breaks = seq(-5, 40, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.line.y.left = element_line(colour = "red"),
        axis.line.y.right = element_line(colour = "blue"),
        axis.ticks.y.left = element_line(colour = "red"),
        axis.ticks.y.right = element_line(colour = "blue"),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9))
  
g$layers[[1]] <- geom_segment(aes(x = Date, y = Precipitation / 1.5 - 6, xend = Date, yend = -6), color = "blue", lineend = "butt")
g  


datastudy <- read.csv("Precipitation_temp.csv")
colnames(datastudy) <- c("Date", "Maximum temperature", 
                      "Minimum temperature", "Precipitation")

datastudy$Date <- as.Date(datastudy$Date, format = "%m/%d/%Y")
str(datastudy)
summary(datastudy)


P <- ggplot(datastudy) +
  geom_bar(aes(x = Date, y = Precipitation * 40/70), fill = "blue", stat = "identity") +
  geom_line(aes(x = Date, y = `Maximum temperature`,linetype = "Maximum temperature"), color = "red") +
  geom_line(aes(x = Date, y = `Minimum temperature`, linetype = "Minimum temperature"), color = "red", na.rm = T) +
  scale_x_date(name = "Time", breaks = seq.Date(as.Date("2020-07-20"), as.Date("2020-09-21"), by = "1 week"), #+
               labels = date_format("%b %d")) +
  scale_y_continuous(name = "Temperature",
                     sec.axis = sec_axis(~ . * 70/40, name = "Precipitation"), breaks = seq(0, 40, 10)) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  theme(axis.line.y.left = element_line(colour = "red"),
        axis.line.y.right = element_line(colour = "blue"),
        axis.ticks.y.left = element_line(colour = "red"),
        axis.ticks.y.right = element_line(colour = "blue"),
        legend.title = element_blank(),
        legend.position = c(0.9, 0.9))

#g$layers[[1]] <- geom_segment(aes(x = Date, y = Precipitation / 1.5 - 6, xend = Date, yend = -6), color = "blue", lineend = "butt")
P  

