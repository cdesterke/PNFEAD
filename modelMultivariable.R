##model


data<-read.table("PNF_EAD.csv",sep=",",h=T,na.string="NA")
library(dplyr)


data%>%select(outcome_PNF_EAD,BILI,PLQ,TEMPS_INTERVENTION,REPERFUSION_STEATOSE_MICRO)->small



library(multirocauc)
small$outcome_PNF_EAD<-as.numeric(small$outcome_PNF_EAD)

geneset<-colnames(small[,-1])
roc.list<-roclist(small,geneset,outcome="outcome_PNF_EAD")
roc.list


rocplot(roc.list,line=1,title="PNF+EAD",police=14)

paste(geneset,collapse="+")



library(corrplot)
library(Hmisc)

df<-small[,-1]

res <- rcorr(as.matrix(df))  # attention, rcorr nécessite une matrice numérique
cor_matrix <- res$r
p_matrix <- res$P




# corrplot
corrplot(cor_matrix,
         method = "number",      
         type = "upper",
         tl.col = "black",
         tl.cex = 0.7,
         number.cex = 2,       
         p.mat = p_matrix,
         sig.level = 0.05,
         insig = "blank")      










library(Epi)

      
ROC(form = outcome_PNF_EAD~BILI+PLQ+TEMPS_INTERVENTION+REPERFUSION_STEATOSE_MICRO, plot="ROC", data=small)


      



model<-glm(outcome_PNF_EAD~BILI+PLQ+TEMPS_INTERVENTION+REPERFUSION_STEATOSE_MICRO,data=small,family="binomial")
summary(model)

library(broom)
library(broom.helpers)
library(GGally)

ggcoef_model(model,exponentiate=T,colour_guide=TRUE)+
scale_color_brewer(palette="Set2")+theme(text=element_text(size=18))+theme(legend.position="none")

res <- tidy(model, conf.int = TRUE, conf.level = 0.95,exponentiate=T)
write.csv(res,file="multivariable4variables.csv",row.names=F)

small%>%select(outcome_PNF_EAD,BILI,PLQ,TEMPS_INTERVENTION,REPERFUSION_STEATOSE_MICRO)->selection

selection%>%mutate(score_EAD_PNF=((BILI*0.0027964)+(PLQ*-0.0045604)+(TEMPS_INTERVENTION*0.0034242)+(REPERFUSION_STEATOSE_MICRO*0.0225287)))->selection
selection$outcome_PNF_EAD<-as.factor(selection$outcome_PNF_EAD)
library(ggplot2)
library(ggbeeswarm)
ggplot(selection,aes(outcome_PNF_EAD,score_EAD_PNF))+geom_boxplot(outlier.shape=NA) + scale_fill_brewer(palette="Pastel1")+
  geom_point(aes(fill=factor(outcome_PNF_EAD)),size=2,shape = 21, alpha = .8, position = position_dodge2(width = .5))+
  theme_classic(base_size=18) +
  theme(legend.position = "none")+xlab("")+ggtitle("")




library(ggpubr)

ggplot(selection, aes(outcome_PNF_EAD, score_EAD_PNF)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(fill = factor(outcome_PNF_EAD)), size = 2.5, shape = 21, alpha = 1, position = position_dodge2(width = .5)) +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  xlab("outcome_EAD+PNF") +ylab("SCORE") +
  ggtitle("") +
  stat_compare_means(method = "t.test", label = "p.format")+ geom_hline(yintercept = 2.19860, 
             linetype = "dashed", 
             color = "red", 
             size = 1.2) 



data$outcome_PNF_EAD<-as.factor(data$outcome_PNF_EAD)

data$score_EAD_PNF<-selection$score_EAD_PNF

ggplot(data, aes(outcome_PNF_EAD, score_EAD_PNF)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(fill = factor(SEXE_DONNEUR)), size = 3, shape = 21, alpha = 1, position = position_dodge2(width = .5)) +
  theme_classic(base_size = 18) +
  theme(legend.position = "bottom") +
  xlab("outcome_EAD+PNF") +ylab("SCORE") +
  ggtitle("") +
  stat_compare_means(method = "t.test", label = "p.format") 


ggplot(selection, aes(outcome_PNF_EAD, score_EAD_PNF)) +
  geom_boxplot(outlier.shape = NA) +
  scale_fill_brewer(palette = "Pastel1") +
  geom_point(aes(fill = factor(outcome_PNF_EAD)), 
             size = 2.5, shape = 21, alpha = 1, 
             position = position_dodge2(width = .5)) +
  geom_hline(yintercept = 2.19860, 
             linetype = "dashed", 
             color = "red", 
             size = 1.2) +
  annotate("text", x = 1.5, y = 2.19860 + 2.7, 
           label = "threshold = 2.19860", 
           color = "red", size = 5, fontface = "italic") +
  theme_classic(base_size = 16) +
  theme(legend.position = "none") +
  xlab("outcome_EAD+PNF") +
  ylab("SCORE") +
  ggtitle("") +
  stat_compare_means(method = "t.test", label = "p.format")






library(rms)

m <- lrm(outcome_PNF_EAD ~ BILI + PLQ +TEMPS_INTERVENTION + REPERFUSION_STEATOSE_MICRO, 
             data = small, x = TRUE, y = TRUE)


cal <- calibrate(m, method = "boot", B = 100)
plot(cal)

library(regplot)
regplot(m, clickable=F, points=T, droplines=T,rank="sd")


save(selection,file="selection.rda")
save(small,file="small.rda")
save(data,file="data.rda")








library(cutpointr)

sel<-selection[complete.cases(selection$score_EAD_PNF),]
cp <- cutpointr(sel, score_EAD_PNF,outcome_PNF_EAD,method = maximize_metric, metric = sum_sens_spec)

plot(cp)
cp





sel%>%mutate(class_EAD_PNF=ifelse(score_EAD_PNF>=2.19860,1,0))->sel

df<-data[complete.cases(data$class_EAD_PNF),]


mini<-glm(outcome_PNF_EAD~class_EAD_PNF,data=df,family="binomial")
exp(coef(mini))

mini<-glm(outcome_PNF_EAD~score_EAD_PNF,data=df,family="binomial")
exp(coef(mini))

ROC(form = outcome_PNF_EAD~score_EAD_PNF, plot="ROC", data=df)



draw_confusion_matrix <- function(cm) {

  total <- sum(cm$table)
  res <- as.numeric(cm$table)

  # Generate color gradients. Palettes come from RColorBrewer.
  greenPalette <- c("#F7FCF5","#E5F5E0","#C7E9C0","#A1D99B","#74C476","#41AB5D","#238B45","#006D2C","#00441B")
  redPalette <- c("#FFF5F0","#FEE0D2","#FCBBA1","#FC9272","#FB6A4A","#EF3B2C","#CB181D","#A50F15","#67000D")
  getColor <- function (greenOrRed = "green", amount = 0) {
    if (amount == 0)
      return("#FFFFFF")
    palette <- greenPalette
    if (greenOrRed == "red")
      palette <- redPalette
    colorRampPalette(palette)(100)[10 + ceiling(90 * amount / total)]
  }

  # set the basic layout
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CONFUSION MATRIX', cex.main=2)

  # create the matrix 
  classes = colnames(cm$table)
  rect(150, 430, 240, 370, col=getColor("green", res[1]))
  text(195, 435, classes[1], cex=1.2)
  rect(250, 430, 340, 370, col=getColor("red", res[3]))
  text(295, 435, classes[2], cex=1.2)
  text(125, 370, 'Predicted', cex=1.3, srt=90, font=2)
  text(245, 450, 'Actual', cex=1.3, font=2)
  rect(150, 305, 240, 365, col=getColor("red", res[2]))
  rect(250, 305, 340, 365, col=getColor("green", res[4]))
  text(140, 400, classes[1], cex=1.2, srt=90)
  text(140, 335, classes[2], cex=1.2, srt=90)

  # add in the cm results
  text(195, 400, res[1], cex=1.6, font=2, col='white')
  text(195, 335, res[2], cex=1.6, font=2, col='white')
  text(295, 400, res[3], cex=1.6, font=2, col='white')
  text(295, 335, res[4], cex=1.6, font=2, col='white')

  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = "DETAILS", xaxt='n', yaxt='n')
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=2)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=2)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 85, names(cm$byClass[5]), cex=1.2, font=2)
  text(50, 70, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 85, names(cm$byClass[6]), cex=1.2, font=2)
  text(70, 70, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 85, names(cm$byClass[7]), cex=1.2, font=2)
  text(90, 70, round(as.numeric(cm$byClass[7]), 3), cex=1.2)

  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=2)
  text(30, 20, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=2)
  text(70, 20, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}


library(caret)
sel$class_EAD_PNF<-as.factor(sel$class_EAD_PNF)
cm <- confusionMatrix(data = sel$class_EAD_PNF, reference = sel$outcome_PNF_EAD)



draw_confusion_matrix(cm)




library(vcd)

struct <- structable(~ class_EAD_PNF+outcome_PNF_EAD, data = sel)
mosaic(struct, , direction = "v", pop = FALSE,colorize = T, shade = TRUE)
labeling_cells(text = as.table(struct), margin = 0)(as.table(struct))
chisq.test(struct)


save(sel,file="selection_model_data_computed.rda")

save(model,file="model_EADPNF.rda")
